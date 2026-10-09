#####
##### TRACER-MIP prescribed two-mode aerosol profiles (Tier 1 control), ported from the
##### reference notebook `Pyplt.TRACER_Aerosol_Profiles_MIP.ipynb` (ARM-Synergy/tracer-mip @ 416423ac).
#####
##### Notebook algorithm, reproduced exactly:
#####   1. surface numbers N_k [cm⁻³] (Table 3) → mixing ratio n_k = N_k / ρ_sfc [mg⁻¹], ρ_sfc = 1.159 kg m⁻³
#####      (the notebook's `density[1]`, the base-state density at the first level above ground);
#####   2. the lidar-derived shape f(z) = x0 e^{-γz} - xerf erf((z - rm)/s) + xc (z in km),
#####      normalized by its maximum (at z = 0) → f̂(z) ∈ (0, 1];
#####   3. n_k(z) = n_k f̂(z), floored at r_k · 50 mg⁻¹ with r_k = n_k / Σ n (so the TOTAL never
#####      drops below 50 mg⁻¹ and the modal ratio is preserved); above 6 km AGL the floor is imposed.
##### The notebook's `profile_aug7` branch for z > 6 km only affects the shape there, which the
##### floor overrides, so it has no effect on the control profiles and is not reproduced.
#####

using SpecialFunctions: erf

struct TracerMIPAerosolMode{FT}
    number_cm3 :: FT          # ambient surface concentration [cm⁻³]
    median_diameter :: FT     # [m]
    sigma :: FT               # geometric standard deviation [-]
end

struct TracerMIPAerosolProfile{FT, M}
    modes :: M                              # Tuple of TracerMIPAerosolMode
    reference_density :: FT                 # [kg m⁻³], cm⁻³ → kg⁻¹ conversion
    minimum_total_number :: FT              # [kg⁻¹] (50 mg⁻¹ = 5e7 kg⁻¹ in the protocol)
    floor_height :: FT                      # [m] above which only the floor applies (6 km)
    shape :: NTuple{6, FT}                  # (x0, γ, xerf, rm, s, xc), z in km
    multiplier :: FT                        # sensitivity factor on surface numbers (CTRL = 1)
end

"""
    tracer_mip_aerosol_profile(case; protocol = tracer_mip_protocol(), FT = Float64, multiplier = 1, scale_floor = false)

The prescribed two-mode aerosol profile of `case` (`:aug07` or `:jun17`). `multiplier` scales the
surface numbers of every mode (3 for the HIGH variant; the LOW variant's 0.3 versus 1/3 is
unresolved in the roadmap — pass it explicitly). The 50 mg⁻¹ total floor is applied as in the
notebook (unscaled) unless `scale_floor = true` (the roadmap text says the floor "scales with the
surface aerosol number concentrations"; the notebook does not scale it).
"""
function tracer_mip_aerosol_profile(case; protocol = tracer_mip_protocol(), FT = Float64,
                                    multiplier = 1, scale_floor = false)
    a = protocol["aerosol"]
    modes_dict = a[String(case)]
    modes = Tuple(TracerMIPAerosolMode{FT}(modes_dict[k]["number_cm3"], 1e-9 * modes_dict[k]["median_diameter_nm"], modes_dict[k]["sigma"])
                  for k in ("mode1", "mode2"))
    fit = a["shape_fit"]
    shape = FT.((fit["x0"], fit["gam"], fit["xerf"], fit["rm"], fit["s"], fit["xc"]))
    floor = 1e6 * a["minimum_total_number_mg"] * (scale_floor ? multiplier : 1)
    return TracerMIPAerosolProfile{FT, typeof(modes)}(modes, a["reference_density_kg_m3"], floor, 6000, shape, multiplier)
end

# Unnormalized lidar shape fit; z in meters.
@inline function aerosol_shape_fit(shape, z)
    x0, γ, xerf, rm, s, xc = shape
    zkm = z / 1000
    return x0 * exp(-γ * zkm) - xerf * erf((zkm - rm) / s) + xc
end

"""
    aerosol_shape(profile, z)

Normalized vertical shape f̂(z) ∈ (0, 1] at height `z` [m AGL] (1 at the surface).
"""
@inline aerosol_shape(p::TracerMIPAerosolProfile, z) = aerosol_shape_fit(p.shape, z) / aerosol_shape_fit(p.shape, zero(z))

"""
    surface_number_mixing_ratios(profile)

Surface number mixing ratio of each mode [kg⁻¹] (`multiplier · N_k / ρ_sfc`, with 1 cm⁻³ = 10⁶ m⁻³).
"""
surface_number_mixing_ratios(p::TracerMIPAerosolProfile) =
    map(m -> p.multiplier * 1e6 * m.number_cm3 / p.reference_density, p.modes)

"""
    mode_number_mixing_ratios(profile, z)

Number mixing ratio of each mode [kg⁻¹] at height `z` [m AGL]: the surface value times the
shape, floored at the mode's share of the 50 mg⁻¹ total floor, and the floor alone above 6 km.
"""
@inline function mode_number_mixing_ratios(p::TracerMIPAerosolProfile, z)
    n₀ = surface_number_mixing_ratios(p)
    total = sum(n₀)
    f = aerosol_shape(p, z)
    above = z > p.floor_height
    return map(n₀) do nₖ
        floorₖ = p.minimum_total_number * nₖ / total
        ifelse(above, floorₖ, max(nₖ * f, floorₖ))
    end
end

"""
    total_number_mixing_ratio(profile, z)

Total (all modes) aerosol number mixing ratio [kg⁻¹] at `z` [m AGL].
"""
@inline total_number_mixing_ratio(p::TracerMIPAerosolProfile, z) = sum(mode_number_mixing_ratios(p, z))

"""
    aerosol_profile_table(profile, heights)

Rows `(z, n₁, n₂, total)` in m and mg⁻¹, matching the notebook's printed table.
"""
aerosol_profile_table(p::TracerMIPAerosolProfile, heights) =
    [(z = z, modes = 1e-6 .* mode_number_mixing_ratios(p, z), total = 1e-6 * total_number_mixing_ratio(p, z)) for z in heights]

"""
    tracer_mip_p3_aerosol_modes(profile; kappa = 0.26, aerosol_density = 1770, molecular_weight_aerosol = 0.132)

Breeze P3 `AerosolMode`s for the two TRACER-MIP modes at their *surface* number mixing ratios,
geometric mean radius = median diameter / 2 and the protocol σ. The protocol prescribes a bulk
hygroscopicity κ = 0.26 instead of a composition; P3's lognormal activation uses the solute
parameter β = ν φ εₘ M_w ρₐ / (Mₐ ρ_w), which equals κ for a κ-Köhler particle (Petters and
Kreidenweis 2007: κ = ν Φ ρₛ M_w / (ρ_w Mₛ)), so the van 't Hoff factor is chosen to make
β = κ given the (otherwise inconsequential) density and molecular weight: ν = κ Mₐ ρ_w / (M_w ρₐ).
Precision follows `Oceananigans.defaults.FloatType`.
"""
function tracer_mip_p3_aerosol_modes(p::TracerMIPAerosolProfile; kappa = 0.26,
                                     aerosol_density = 1770, molecular_weight_aerosol = 0.132,
                                     thermodynamic_constants = ThermodynamicConstants())
    ρ_w = thermodynamic_constants.liquid.density
    M_w = thermodynamic_constants.vapor.molar_mass
    vant_hoff_factor = kappa * molecular_weight_aerosol * ρ_w / (M_w * aerosol_density)
    n₀ = surface_number_mixing_ratios(p)
    modes = map(p.modes, n₀) do m, n
        AerosolMode(; number_mixing_ratio = n, mean_radius = m.median_diameter / 2, geometric_std = m.sigma,
                    vant_hoff_factor, osmotic_potential = 1, mass_fraction_soluble = 1,
                    aerosol_density, molecular_weight_aerosol, thermodynamic_constants)
    end
    return modes
end
