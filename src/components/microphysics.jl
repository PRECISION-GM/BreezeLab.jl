#####
##### Microphysics helpers. Case constructors take a Breeze microphysics object directly
##### (e.g. `P3Microphysics(; cloud, aerosol)`); the helpers here build the aerosol parts the
##### ENA protocols prescribe and return Breeze types. Precision follows
##### `Oceananigans.defaults.FloatType`.
#####

"""
    lasso_aerosol(; setting=:aer2, reference_density, maximum_supersaturation=nothing,
                    aerosol_density=1790, molecular_weight_aerosol=0.115, vant_hoff_factor=3,
                    mass_fraction_soluble=1)

The two lognormal aerosol modes of the LASSO-ENA spectral-bin configuration as a Breeze
`AerosolActivation` with a prognostic (depleting) reservoir, as the SBM's. LASSO quotes
number *concentrations* per cm³:

| setting | mode 1 (r = 0.018 μm, σ = 1.53) | mode 2 (r = 0.066 μm, σ = 1.78) |
|---------|---------------------------------|---------------------------------|
| aer1    | 138 cm⁻³                        | 140.5 cm⁻³                      |
| aer2    | 276 cm⁻³                        | 281 cm⁻³                        |
| aer3    | 552 cm⁻³                        | 562 cm⁻³                        |

P3 takes the number *per unit mass* of air, `nᵃ = 10⁶ N / ρ` [kg⁻¹]. The LASSO HUJI-SBM
initializes its CCN with a constant *mixing ratio* with height
(`MICRO_HUJISBM/microphysics.f90`: `FCCN0 = FCCNR_mp * rhocgs(k)/rhocgs(1)`), so the
conversion uses the first-level reference density (see
[`first_level_reference_density`](@ref)): `number_mixing_ratio = 1e6 N / ρ₁`, uniform in z.
The chemistry is the HUJI-SBM's (aerosol density 1790 kg m⁻³, molecular weight
0.115 kg mol⁻¹, van 't Hoff factor 3, fully soluble: `micro_prm.f90` has no insoluble
fraction, so Breeze's default `mass_fraction_soluble = 0.9` is overridden).

`maximum_supersaturation` emulates the SBM's activation cap (`ss_max = 0.003` in
`MICRO_HUJISBM/microphysics.f90`: no particle activates above 0.3 %): each mode's number is
reduced to the fraction of its particles with critical supersaturation below the cap
([`activated_fraction`](@ref)). Breeze's Morrison–Grabowski activation itself has no cap
(see `docs/cases/ena_aerosol_audit.md`). The cap is an approximation that keeps the
lognormal shape and scales the number; it is not the SBM.
"""
function lasso_aerosol(; setting = :aer2, reference_density, maximum_supersaturation = nothing,
                         aerosol_density = 1790, molecular_weight_aerosol = 0.115, vant_hoff_factor = 3,
                         mass_fraction_soluble = 1)
    modes = lasso_aerosol_modes(; setting, reference_density, maximum_supersaturation,
                                aerosol_density, molecular_weight_aerosol, vant_hoff_factor, mass_fraction_soluble)
    return AerosolActivation(modes...; prognostic = true)
end

function lasso_aerosol_modes(; setting, reference_density, maximum_supersaturation, chemistry...)
    N₁, N₂ = setting === :aer1 ? (138.0, 140.5) :
             setting === :aer2 ? (276.0, 281.0) :
             setting === :aer3 ? (552.0, 562.0) :
             throw(ArgumentError("unknown LASSO aerosol setting $setting (aer1, aer2, aer3)"))
    geometry = ((; mean_radius = 0.018e-6, geometric_std = 1.53), (; mean_radius = 0.066e-6, geometric_std = 1.78))
    fractions = isnothing(maximum_supersaturation) ? (1.0, 1.0) :
        Tuple(activated_fraction(g.mean_radius, g.geometric_std, maximum_supersaturation; chemistry...) for g in geometry)
    n₁ = fractions[1] * 1e6 * N₁ / reference_density
    n₂ = fractions[2] * 1e6 * N₂ / reference_density
    return (AerosolMode(; number_mixing_ratio = n₁, geometry[1]..., chemistry...),
            AerosolMode(; number_mixing_ratio = n₂, geometry[2]..., chemistry...))
end

"""
    covert_aerosol(; reference_density, droplet_number=75e6, maximum_supersaturation=0.003)

A Covert-consistent prognostic aerosol for the public Covert et al. (2022) ENA case, which
prescribes the *observed* droplet number `Nc = 75 cm⁻³` and publishes no aerosol spectrum.
The modes keep the LASSO-ENA aer2 shapes and SBM chemistry (see [`lasso_aerosol`](@ref)) and
are scaled by one factor so that the number activatable below `maximum_supersaturation`
equals `droplet_number` [m⁻³] at the first-level reference density. With the diagnostic-CCN
projection this bounds the in-cloud droplet number at the observed value. It is a labelled
configuration of this package, not an ARM or Covert prescription.
"""
function covert_aerosol(; reference_density, droplet_number = 75e6, maximum_supersaturation = 0.003)
    aer2 = lasso_aerosol_modes(; setting = :aer2, reference_density, maximum_supersaturation,
                               aerosol_density = 1790, molecular_weight_aerosol = 0.115, vant_hoff_factor = 3,
                               mass_fraction_soluble = 1)
    activatable = aer2[1].number_mixing_ratio + aer2[2].number_mixing_ratio          # kg⁻¹, already capped
    factor = droplet_number / reference_density / activatable
    modes = Tuple(AerosolMode(; number_mixing_ratio = factor * m.number_mixing_ratio, mean_radius = m.mean_radius,
                              geometric_std = m.geometric_std, vant_hoff_factor = m.vant_hoff_factor,
                              osmotic_potential = m.osmotic_potential, mass_fraction_soluble = m.mass_fraction_soluble,
                              aerosol_density = m.aerosol_density, molecular_weight_aerosol = m.molecular_weight_aerosol) for m in aer2)
    return AerosolActivation(modes...; prognostic = true)
end

"""
    activated_fraction(mean_radius, geometric_std, S; T=285, aerosol_density=1790,
                       molecular_weight_aerosol=0.115, vant_hoff_factor=3, mass_fraction_soluble=1,
                       osmotic_potential=1)

Fraction of a lognormal aerosol mode activated at supersaturation `S` (fraction, not %)
and temperature `T` in the Morrison & Grabowski (2007) Köhler closure Breeze's P3 uses,
evaluated with Breeze's own `activated_number` on a unit-number mode so that the two can
never disagree.
"""
function activated_fraction(mean_radius, geometric_std, S; T=285.0, aerosol_density=1790, molecular_weight_aerosol=0.115,
                            vant_hoff_factor=3, mass_fraction_soluble=1, osmotic_potential=1)
    mode = AerosolMode(Float64; number_mixing_ratio=1.0, mean_radius, geometric_std, aerosol_density,
                       molecular_weight_aerosol, vant_hoff_factor, mass_fraction_soluble, osmotic_potential)
    return activated_number(mode, AerosolActivation(mode), Float64(T), Float64(S))
end

"""
    first_level_reference_density(data_dir, z_faces; moisture_basis=:mixing_ratio)

The anelastic reference density [kg m⁻³] at the first cell center of a case whose SAM
`snd` (and `prm`, for `day0`) are in `data_dir`, on the vertical grid `z_faces`, at the
precision `Oceananigans.defaults.FloatType`. It is the value the ENA constructors compute
from the same inputs, so aerosol numbers converted with it (see [`lasso_aerosol`](@ref),
[`covert_aerosol`](@ref)) match the protocol defaults exactly.
"""
function first_level_reference_density(data_dir, z_faces; moisture_basis = :mixing_ratio)
    soundings = read_sam_sounding(joinpath(data_dir, "snd"))
    prm = joinpath(data_dir, "prm")
    day0 = isfile(prm) ? Float64(get(read_sam_namelist(prm), "day0", soundings[1].day)) : soundings[1].day
    profiles = SoundingProfiles(initial_sounding(soundings, day0); moisture_basis)
    grid = RectilinearGrid(CPU(); size = length(z_faces) - 1, z = z_faces, halo = 5, topology = (Flat, Flat, Bounded))
    reference_state = sounding_reference_state(grid, profiles, ThermodynamicConstants())
    return Array(interior(reference_state.density, 1, 1, :))[1]
end

"""
    sounding_reference_state(grid, profiles, constants)

Breeze `ReferenceState` of a SAM sounding: base pressure at the sea surface (z = 0) and the
sounding's potential temperature and total water.
"""
sounding_reference_state(grid, profiles::SoundingProfiles, constants) =
    ReferenceState(grid, constants;
                   base_pressure = profiles.surface_pressure,
                   potential_temperature = z -> profiles(:θ, z),
                   vapor_mass_fraction = z -> profiles(:qᵗ, z))

is_p3(microphysics) = microphysics isa PredictedParticlePropertiesMicrophysics

"""
    has_aerosol_reservoir(microphysics)

`true` for P3 with a prognostic (depleting) aerosol reservoir `ρnᵃ`.
"""
has_aerosol_reservoir(microphysics) =
    is_p3(microphysics) && !isnothing(microphysics.aerosol) && has_prognostic_aerosol(microphysics.aerosol)

"""
    microphysics_record(microphysics; reference_density)

Provenance of a microphysics object: its type, the P3 cloud droplet number and, for each
aerosol mode `i`, the number mixing ratio `nᵢ` [kg⁻¹], the equivalent first-level
concentration `Nᵢ` [cm⁻³], and the mode's shape and chemistry parameters.
"""
function microphysics_record(microphysics; reference_density)
    record = (; microphysics = string(nameof(typeof(microphysics))))
    is_p3(microphysics) || return record
    record = merge(record, (; cloud_droplet_number = Float64(microphysics.cloud.number_concentration),
                              aerosol = isnothing(microphysics.aerosol) ? "none" : summary(microphysics.aerosol)))
    (isnothing(microphysics.aerosol) || !hasproperty(microphysics.aerosol, :modes)) && return record
    for (i, m) in enumerate(microphysics.aerosol.modes)
        n = Float64(m.number_mixing_ratio)
        names = [Symbol("n", subscript(i)), Symbol("N", subscript(i))]
        values = Any[n, n * reference_density / 1e6]
        for property in (:mean_radius, :geometric_std, :kappa, :aerosol_density, :molecular_weight_aerosol,
                         :vant_hoff_factor, :mass_fraction_soluble)
            hasproperty(m, property) || continue
            push!(names, Symbol(property, "_", i))
            push!(values, Float64(getproperty(m, property)))
        end
        record = merge(record, NamedTuple{Tuple(names)}(Tuple(values)))
    end
    return merge(record, (; reference_density = Float64(reference_density)))
end

subscript(i::Integer) = join(Char(0x2080 + d) for d in reverse(digits(i)))

"""
    specific_prognostic_names(microphysics)

The specific names (`qᶜˡ`, `nᶜˡ`, ...) of the scheme's prognostic density fields (`ρqᶜˡ`, ...),
the keys under which Breeze accepts their forcings.
"""
specific_prognostic_names(microphysics) =
    Tuple(Symbol(string(ρname)[nextind(string(ρname), 1):end])
          for ρname in Breeze.AtmosphereModels.prognostic_field_names(microphysics))

# The floating-point type a Breeze microphysics object was built with: its first floating-point
# field, searched breadth-first through nested parameter structs (`nothing` when it has none).
function microphysics_precision(microphysics; depth = 4)
    level = Any[microphysics]
    for _ in 1:depth
        next = Any[]
        for x in level
            x isa AbstractFloat && return typeof(x)
            x isa Union{AbstractArray, Function, Module, Type, Symbol, AbstractString, Nothing} && continue
            isstructtype(typeof(x)) || continue
            for name in fieldnames(typeof(x))
                isdefined(x, name) && push!(next, getfield(x, name))
            end
        end
        isempty(next) && return nothing
        level = next
    end
    return nothing
end

"""
    check_precision(microphysics, FT)

Throw when `microphysics` was built at a precision other than the grid's `FT` (for example
when `Oceananigans.defaults.FloatType` changed between building the two).
"""
function check_precision(microphysics, FT)
    P = microphysics_precision(microphysics)
    isnothing(P) || P === FT ||
        throw(ArgumentError("the microphysics was built with $P but the grid uses $FT; build both after setting Oceananigans.defaults.FloatType"))
    return nothing
end
