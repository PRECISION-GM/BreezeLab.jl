# SEA STARR — driver inventory and Breeze mapping

SEA STARR (Southeast Atlantic Stratocumulus Transitions with Aerosol-Rain-Radiation
interactions) is an LES/SCM intercomparison of a 66-hour Lagrangian stratocumulus-to-cumulus
transition under a biomass-burning smoke layer, driven by a composite of GEOS-5/SEVIRI
trajectories (Diamond et al., in prep.). Official inputs: Zenodo record
[22241697](https://zenodo.org/records/22241697) (DOI 10.5281/zenodo.22241697), pinned in
`cases/seastarr/inputs.toml`. This page records, per plan
step 1, what the setup PDF and the DEPHY driver NetCDF actually contain and how each quantity
maps onto Breeze. It is an inventory of the inputs and of the implementation choices; it is
**not** evidence that the Breeze CTRL reproduces any other model.

## Protocol summary (setup PDF, `SEA_STARR_Model_Setup_Specifications.pdf`)

| Item | Specification |
| --- | --- |
| Horizontal grid | Δx = Δy = 50 m; 384² (19.2 km) or, if prohibitive, 192² (9.6 km) |
| Vertical grid | 10 m spacing from the surface to 2500 m, then each level 10 % thicker until 6500 m is reached |
| Time step | 1 s |
| Period | 2017-08-15 21:00 UTC → 2017-08-18 15:00 UTC (66 h; the first 3 h are spin-up) |
| Nudging | only above the *domain maximum* inversion height (maximum ∂θₗ/∂z) + 100 m, with a Blossey et al. (2013) ramp over the next 200 m; T, qₜ, Nₐ with τ = 1800 s; u, v with τ = 10800 s |
| Advection | no horizontal advective tendencies for free-tropospheric variables (the nudging accounts for them); vertical advection from the prescribed vertical velocity (`wa` or `wap`) |
| Geostrophic forcing | on, with the Coriolis parameter from the (single) latitude |
| Radiation | on, interactive aerosol with single-scatter albedo 0.85 at 550 nm |
| Aerosol | single lognormal accumulation mode: geometric mean diameter 185 nm, σ = 1.5, κ = 0.2; surface source 70 cm⁻² s⁻¹ (Yamaguchi et al. 2017) |
| Surface | ocean; time-varying skin temperature (`ts`); no prescribed moisture/wind forcing (bulk fluxes); z₀ = 10⁻⁴ m |

## Driver files (`SEA_STARR_{CTRL,N100,N030}_SCM_driver.nc`, DEPHY SCM format v1)

Dimensions: `time` = 67 hourly records from 2017-08-15 21:00 to 2017-08-18 15:00 UTC,
`t0` = 1, `zh` = 336 mid-levels stored **top-down** (76 255 m → 5 m; 10 m spacing for the 250
levels below 2500 m, stretched above), `lat` = `lon` = 1 (−11.649° N, −8.336° E), `orog` = 0.
Global attributes: `forc_wap = forc_wa = forc_geo = 1`, `adv_ta = adv_thetal = 0`,
`radiation = "on"`, `surface_forcing_temp = "ts"`, `surface_forcing_moisture = "none"`,
`surface_forcing_wind = "none"`, `surface_source_aerosol = "70 cm-2 s-1"`,
`nudging_* = 1800 s / 10800 s`, `zh_nudge_gap = "200 m"`.

The members differ **only** in `na` and `na_nud` (verified by element-wise comparison of every
numeric variable and attribute): CTRL has 100 mg⁻¹ in the boundary layer and a smoke plume of
≈1000–1100 mg⁻¹ in the free troposphere (nudging target up to 1481 mg⁻¹ later); N100 is 100 mg⁻¹
everywhere and at all times; N030 is 30 mg⁻¹ everywhere and at all times.

### Variable table

Column "role": IC = initial condition (dimension `t0`), F = time-varying forcing (dimension
`time`), N = nudging target. "Breeze target" names the quantity in the BreezeLab
implementation (`src/sea_starr.jl`).

| Variable | Dims | Units | Role | Interpolation | Breeze target / treatment |
| --- | --- | --- | --- | --- | --- |
| `zh` | zh | m | coordinate (mid-levels, top-down) | — | reversed to ascending by the reader; the LES grid is built from the protocol (10 m / 10 % stretch), **not** from `zh` |
| `lat`, `lon` | 1 | deg | geometry | — | Coriolis `FPlane(latitude = -11.649)`; solar position (see departures) |
| `ps` | t0 | Pa | IC | — | `ReferenceState(base_pressure = ps)` at z = 0 |
| `pa` | zh, t0 | Pa | IC (diagnostic) | — | not used for the reference state (Breeze integrates hydrostatically from `ps`); recorded. Note `pa(5 m) = 101 659 Pa` versus `ps = 101 984 Pa`, a 325 Pa drop over 5 m, so the driver column is not hydrostatic with `ps` near the surface |
| `thetal` | zh, t0 | K | IC | linear in z | liquid-water potential temperature → `ReferenceState(potential_temperature)` and `set!(model; θ)` (`:LiquidIcePotentialTemperature` formulation) |
| `ta` | zh, t0 | K | IC (diagnostic) | — | "liquid water temperature"; not used (θₗ is the prognostic basis) |
| `qt` | zh, t0 | 1 (mass fraction) | IC | linear in z | total water → `ReferenceState(vapor_mass_fraction)` and `set!(model; qᵗ)`; **no mixing-ratio conversion** (`moisture_basis = :mass_fraction`) |
| `na` | zh, t0 | mg⁻¹ (per mg of air) | IC | linear in z | aerosol number mixing ratio `nᵃ = 10⁶ na` [kg⁻¹]; set as a height-varying initial field of P3's prognostic reservoir `ρnᵃ` |
| `ua`, `va` | zh, t0 | m s⁻¹ | IC | linear in z | `set!(model; u, v)` |
| `ts`, `z0` | t0 | K, m | IC | — | initial SST; momentum roughness length for the bulk formulae |
| `wa` | zh, time | m s⁻¹ | F | linear in z, linear in time | large-scale vertical velocity → `LargeScaleVerticalAdvection` (full-field upwind −w ∂ₓϕ) on u, v, θ, qᵛ and every microphysical species, **exactly once per field**. `wa` is identical at all 67 times (a single profile); ≈ −4.5 mm s⁻¹ at 2500 m |
| `wap` | zh, time | Pa s⁻¹ | F (alternative) | — | not used (`wa` is used) |
| `pa_force` | zh, time | Pa | F (diagnostic) | — | not used by the anelastic model (fixed reference pressure); recorded |
| `ug`, `vg` | zh, time | m s⁻¹ | F | linear in z and time | `TimeVaryingGeostrophicForcing` (`-f vg`, `+f ug`). `ug ≡ ua_nud`, `vg ≡ va_nud` in all three files |
| `ta_nud` | zh, time | K | N | — | not used (θₗ is nudged) |
| `thetal_nud` | zh, time | K | N | linear in z and time | target of the inversion-following mean-profile nudging of θ (τ = 1/`nudging_constant_thetal` = 1800 s) |
| `qt_nud` | zh, time | 1 | N | linear in z and time | target of the nudging of the moisture variable (P3: `qᵛ`; above the inversion there is no condensate so qₜ = qᵛ) |
| `na_nud` | zh, time | mg⁻¹ | N | linear in z and time | target of the nudging of `nᵃ` (× 10⁶) |
| `ua_nud`, `va_nud` | zh, time | m s⁻¹ | N | linear in z and time | targets of the wind nudging (τ = 10800 s) |
| `nudging_constant_*` | zh, time | s⁻¹ | N | — | uniform in height and time: 1/1800 s⁻¹ (ta, thetal, qt, na), 1/10800 s⁻¹ (ua, va). The inversion mask is therefore **not** in the file and must be diagnosed online |
| `ts_force` | time | K | F | linear in time | SST → shared surface-temperature field (bulk fluxes and radiation); 293.33 → 298.06 K |
| `ps_force` | time | Pa | F (diagnostic) | — | not used; 101 984 → 101 523 Pa |
| `o3` | zh, time | mol mol⁻¹ | F (radiation) | — | ozone profile for the radiation column; see departures |

Does any thermodynamic tendency already include vertical transport or radiation? There are
**no** horizontal advective tendencies in the file (`adv_ta = adv_thetal = 0`; no `tnta_*`,
`tnqt_*` variables). The only large-scale thermodynamic forcings are (i) vertical advection by
`wa`, which the LES must apply itself, and (ii) the free-tropospheric nudging, which the PDF
states "already accounts for" advection of free-tropospheric variables. Radiation is fully
interactive in the LES; none of the driver tendencies contains it.

### Initial state

At t₀ the column is a well-mixed layer (θₗ = 290.0 K, qₜ = 9.0 g kg⁻¹, `na` = 100 mg⁻¹) up to the
inversion at 1160 m (maximum ∂θₗ/∂z = 0.064 K m⁻¹, jump of ≈ +10 K), capped by the smoke layer
(θₗ increasing ≈ 10 K km⁻¹, qₜ = 2.36 g kg⁻¹, `na` ≈ 995–1100 mg⁻¹ up to ≈ 3 km, decaying to
≈ 130 mg⁻¹ by 5 km). The boundary layer is saturated above ≈ 850 m (qₛₐₜ(1000 m) ≈ 7.7 g kg⁻¹ <
9 g kg⁻¹), i.e. the initial state contains a stratocumulus deck. The nudging targets'
"inversion" (maximum gradient of `thetal_nud`) jumps between 600 m and 1790 m during the
period, which is why the protocol restricts nudging to above the LES's own inversion.

## Raw trajectories (`SEA_STARR_Raw_Trajectories.nc`)

40 GEOS-5 trajectories (`ind`), 67 hourly times, 72 GEOS levels, with SEVIRI cloud retrievals
(Nd, CF, LWP, rₑ, cloud-top T/p/z, SZA) and GEOS state/aerosol species. The composite's mean
path moves from (−16.4° N, −1.2° E) at t = 0 to (−7.6° N, −14.9° E) at t = 66 h; the driver's
single (`lat`, `lon`) = (−11.65, −8.34) is the mid-period position (t ≈ 30 h). The mean `TS`
of the trajectories equals `ts_force`. The SEVIRI retrievals are the observational reference
for evaluation; no comparison has been made yet.

## Standard variables (`SEA_STARR_Standard_Variables.xlsx`)

Three sheets: statistical output (15 min; profiles and scalars, tier 1/2), 2D output (15 min)
and 3D output (1 h). Time axis "seconds since 2017-08-16 00:00:00". The tier-1 statistical
set is: `pa, ta, qv, qt, hur, ua, va, wa, theta, thetal, tke, hfss, hfls, ps, ts, zi, w2, ql,
qr, cl, clt, pr, cwp, lwp, rwp, nt, na, nc, nr, rlh, rsh, rsdt, rsut, rlut, rsutcs, rlutcs`;
2D tier 1: `prw, zi, pr, cwp, lwp, rwp, optc, optr, zct, zcb, rsdt, rsut, rlut, rsutcs, rlutcs`.
Number concentrations are per kg (`kg-1`). Cloud fraction `cl` uses a 0.01 g kg⁻¹ liquid
threshold; `clt` counts columns with cloud optical thickness > 2. `nt` is total aerosol
(dry + in droplets + in rain); `na` is dry (unactivated) aerosol.

## Implementation (`src/dephy_driver.jl`, `src/sea_starr_aerosol.jl`, `src/sea_starr_forcings.jl`, `src/sea_starr.jl`)

`sea_starr(; member, arch, FT, Nx, Ny, ...)` reads the driver, builds the protocol grid
(`sea_starr_vertical_faces`: 10 m × 250 to 2500 m, then 10 % growth, last interface clipped
to 6500 m, Nz = 288), the anelastic reference state (hydrostatic from `ps`, θₗ and qₜ of the
driver), P3 microphysics with one `KappaAerosolMode` and a prognostic reservoir, full-field
upwind subsidence (`LargeScaleVerticalAdvection` with the driver `wa`) on u, v, θ, qᵛ and
every microphysical prognostic exactly once, the time-varying geostrophic wind, bulk surface
fluxes from `ts_force` with the driver `z0`, RRTMGP LW+SW with the driver ozone profile,
inversion-following mean-profile nudging of θ, qᵛ, nᵃ (τ = 1800 s) and u, v (τ = 10800 s),
the surface aerosol source, evaporation regeneration, an equilibrium initial cloud, JLD2
output (15-min statistics, 1-min time series with the aerosol number budget columns, 15-min 2D
fields, hourly 3D fields) and a 3-hourly checkpoint. `cases/sea_starr.jl` constructs, writes
provenance, runs and plots; `execution/sea_starr.sbatch` is the UUID-pinned launcher;
`test/sea_starr.jl` holds the reference-value tests.

Aerosol processes (plan step 4), and what the pinned Breeze P3 provides:

| Process | Pinned Breeze P3 | BreezeLab |
| --- | --- | --- |
| Activation | M&G2007 erf form with the mode's *constant* number, only capped by the local reservoir | `aerosol_activation_rate` specialized for `KappaAerosolMode`: target `f(S) (nᶜˡ + nᵃ)` from the κ-Köhler critical supersaturation of the *local* reservoir; one aerosol consumed per droplet (Breeze) |
| Surface source | none (`FluxBoundaryCondition` on `ρnᵃ` is ignored by P3's field materialization) | `SurfaceAerosolSource`: 7×10⁵ m⁻² s⁻¹ as a bottom-cell tendency |
| Regeneration on droplet evaporation | none (`nᶜˡ` kept until `qᶜˡ` is exhausted, then discarded) | `EvaporationRegeneration`: nᶜˡ removed ∝ P3's own diagnosed cloud-evaporation rate and returned to nᵃ (one-step lag) |
| Collision scavenging | droplet/rain number lost in autoconversion, accretion, self-collection; nothing returned | unchanged (this *is* the coalescence sink of aerosol) |
| Rain removal | rain number sediments through the surface | unchanged |
| Regeneration on rain evaporation | P3 reduces `nʳ`; nothing returned | not implemented (nʳ ≪ nᵃ; labelled) |
| Interstitial/impaction scavenging | none | not implemented (labelled) |
| Transport | — | advection, subsidence, nudging on `nᵃ` like every other scalar |

The column budget `∫(ρnᵃ + ρnᶜˡ + ρnʳ) dz` is written every minute (`n_total_column`); it
changes only through the surface source, nudging, subsidence, coalescence and rain removal.

## Implementation choices and departures

The following are either not specified by the PDF/driver or not available in the pinned
Breeze revision. Each is listed in `case.config.departures` and written to the provenance file.
A run with them is an **exploratory CTRL**, not a protocol-complete one.

1. **Nudging ramp.** Blossey et al. (2013) could not be retrieved from this cluster (paywall).
   The implemented ramp is linear in height: weight 0 at z_b = max(zᵢ) + 100 m, 1 at z_b + 200 m
   (the driver's `zh_nudge_gap`). The inversion height is the per-column height of maximum
   ∂θₗ/∂z (searched below 4 km), and `max(zᵢ)` is the maximum over all columns ("domain maximum
   inversion height"). The nudging acts on the horizontal mean so that resolved
   free-tropospheric eddies are not damped.
2. **Nudged thermodynamic variable.** θₗ (the prognostic) is nudged to `thetal_nud`; `ta_nud`
   is not used. The P3 moisture variable is the vapor mass fraction `qᵛ`, nudged to `qt_nud`.
3. **Aerosol chemistry.** κ enters the activation directly (κ-Köhler critical supersaturation;
   identical to Morrison & Grabowski's β_act form with β = κ).
4. **Aerosol processes.** See the table above (proportional evaporation regeneration; no
   regeneration from rain evaporation; no interstitial scavenging). ENA's
   `DiagnosticCCNProjection` reset is **not** used.
5. **Radiation.** Breeze's RRTMGP has no aerosol optics (single-scatter albedo 0.85 not
   representable, so the smoke-layer absorption is absent), its column ends at the LES top
   (6.5 km; no atmosphere above, zero downwelling longwave at the top), and `ApparentSolarPosition`
   takes a fixed coordinate (the driver's single lat/lon rather than the moving composite
   trajectory, which spans ≈ 9° of latitude and 14° of longitude over the period).
6. **Vertical grid top.** 38 stretched levels reach 6504 m; the last interface is clipped to
   exactly 6500 m (Nz = 288).
7. **Upper sponge.** `SAMSponge` over the top 15 % of the domain (z > 5525 m), not in the protocol.
8. **Initial cloud.** Warm-phase saturation adjustment of the driver (θₗ, qₜ) gives the initial
   cloud; in cloudy cells the whole aerosol reservoir is initially activated (what the activation
   scheme would do at the large initial supersaturation), elsewhere `nᵃ` is the driver profile.
9. **Surface fluxes.** Breeze bulk formulae (Large & Yeager polynomials, Li et al. stability) with
   z₀ = 10⁻⁴ m and 0.1 m s⁻¹ gustiness; the protocol only prescribes `ts` and `z0`.
10. **Time step.** Fixed 1 s as specified; the CFL wizard may only shorten it.
