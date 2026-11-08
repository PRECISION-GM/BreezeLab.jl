# SAM `oceflx` surface fluxes in BreezeLab (LASSO-ENA `flxsst`)

LASSO-ENA runs SAM with `SFC_FLX_FXD = .false.`, `OCEAN = .true.` and
`UNIFORM_SFC_FLX = .false.`, so `surface.f90` calls `oceflx.f90` once per column and
staggered point. BreezeLab's `surface_flux_law = :sam_oceflx` (the `ena_lasso` default) now
evaluates that routine itself: `sam_ocean_fluxes` in
`src/components/sam_ocean_surface_fluxes.jl` is a line-by-line transcription, and
`SAMOceanSurfaceFlux` feeds it per column. The previous approximation (SAM's neutral laws
fitted to Breeze's polynomial form) is kept as `:sam_oceflx_fitted`.

Source: `~/lasso_sam_sbm/SRC/oceflx.f90`, `surface.f90`, `sat.f90`, `params.f90` at public
HEAD `12d0244` (identical in these files to the LASSO model hash `f83adf58`).

## `surface.f90`: what is passed (per column `i, j`)

| `oceflx` argument | scalar call (cell centre) | x-stress call (u point) | y-stress call (v point) |
|---|---|---|---|
| `pmidm1` | `pres(1)` [mb] | same | same |
| `ubot` | `½(u(i+1)+u(i)) + ug` | `u(i) + ug` | 4-point average of `u` + `ug` |
| `vbot` | `½(v(j+1)+v(j)) + vg` | 4-point average of `v` + `vg` | `v(j) + vg` |
| `thbot` | `t(i,j,1)` | `½(t(i-1)+t(i))` | `½(t(j-1)+t(j))` |
| `tbot` | `t - gamaz(1)` | averaged `t - gamaz(1)` | averaged `t - gamaz(1)` |
| `qbot` | `qv(i,j,1)` | `½(qv(i-1)+qv(i))` | `½(qv(j-1)+qv(j))` |
| `zbot` | `z(1)` | same | same |
| `ts` | `sstxy + t00` | averaged | averaged |

`t` is SAM's prognostic liquid-ice static energy over `cp`,
`t = T + gz/cp − (lcond/cp)(qcl+qpl) − (lsub/cp)(qci+qpi)` with `cp = 1004`, `ggr = 9.81`,
`lcond = 2.5104e6`, `lsub = 2.8440e6`. So SAM's "potential temperature" `thbot` is `t`, and its
"temperature" `tbot` is `T − (L/cp) q_condensate`. `(u + ug, v + vg)` is the ground-relative
wind: the model wind plus the namelist translation velocity. No gustiness is added.

Outputs are used as `fluxbt = shf`, `fluxbq = lhf` (kinematic, positive upward), and
`fluxbu = taux/rhow(1)`, `fluxbv = tauy/rhow(1)`. The diffusion routines multiply every
bottom flux by `rhow(1)`, so the applied mass fluxes are `rhow(1)·shf` (of `t`),
`rhow(1)·lhf` (of `qv`) and `taux`, `tauy` (of momentum). `rhow(1) = 2 rhow(2) − rhow(3)` with
the hydrostatic interface densities `rhow(k) = (pres(k−1) − pres(k)) / (z(k) − z(k−1)) / ggr`.
(`SFC_TAU_FXD` is false for `flxsst`, so the stress is not overridden.)

## `oceflx.f90`

Constants: `κ = 0.4`, `zref = 10 m`, `umin = 1 m/s`, `g = 9.81`; neutral 10-m laws
`cdn(U) = 0.0027/U + 0.000142 + 0.0000764 U`, Stanton root `rhn = 0.0327` (unstable) or
`0.018` (stable), Dalton root `ren = 0.0346` (the coefficients are `rhn·√cdn`, `ren·√cdn`).

Setup:

- `ρbot = 100 p / (287 tbot)`, `U = max(umin, √(u² + v²))`
- `θv = thbot (1 + 0.606 qbot)`, `Δθ = thbot − ts`, `Δq = qbot − ssq`
- `ssq = qsatw(ts, p)` if `ts > 271 K`, else `qsati(ts, p)` (Flatau et al. 1992 polynomials in
  `sat.f90`; `q = 0.622 e / max(e, p − e)`)
- `ℓ = ln(zbot/zref)`

Initial guess (z/L = 0, U₁₀ₙ = U): `stable = ½ + sign(½, Δθ)`, `rdn = √cdn(U)`,
`rhn(stable)`, `ren`; `u★ = rdn U`, `θ★ = rhn Δθ`, `q★ = ren Δq`.

Two identical iterations, each:

1. `ζ = κ g zbot (θ★/θv + q★/(1/0.606 + qbot)) / u★²`, clipped to `|ζ| ≤ 10`;
   `stable = ½ + sign(½, ζ)`.
2. `x = max(√|1 − 16ζ|, 1)^½`;
   `ψm = −5ζ` (stable) or `ln((1 + x(2 + x))(1 + x²)/8) − 2 atan x + 1.571` (unstable);
   `ψh = −5ζ` (stable) or `2 ln((1 + x²)/2)` (unstable).
3. 10-m neutral wind with the previous `rdn`: `rd = rdn / (1 + rdn/κ (ℓ − ψm))`,
   `U₁₀ₙ = U rd/rdn`.
4. New neutral coefficients `rdn = √cdn(U₁₀ₙ)`, `rhn(stable)`, `ren = 0.0346`, shifted to `zbot`:
   `rd = rdn / (1 + rdn/κ (ℓ − ψm))`, `rh = rhn / (1 + rhn/κ (ℓ − ψh))`,
   `re = ren / (1 + ren/κ (ℓ − ψh))`.
5. `u★ = rd U`, `θ★ = rh Δθ`, `q★ = re Δq`.

Outputs: `τ = ρbot u★²`, `taux = −τ u/U`, `tauy = −τ v/U`, `shf = −u★θ★` [K m/s],
`lhf = −u★q★` [m/s], `ssq`. (The routine also interpolates a 2-m reference temperature that
SAM discards.)

## BreezeLab mapping (`SAMOceanSurfaceFlux`)

- `p = pᵣ(z₁)/100`: Breeze's anelastic reference pressure at the first centre (SAM's `pres(1)`).
- `zbot = z₁`, the first cell centre.
- `t = T + (9.81/1004) z₁ − (2.5104e6/1004) qˡ − (2.8440e6/1004) qⁱ` from Breeze's temperature and
  condensate (liquid includes rain, ice all ice), with SAM's constants; `thbot = t`,
  `tbot = t − (9.81/1004) z₁`; `qbot = qᵛ`.
- Wind: Breeze's ground-relative `(u, v)` (the LASSO adapter does not use a translating frame,
  so this is SAM's `u + ug`), averaged exactly as in `surface.f90` (scalar fluxes at centres,
  each stress component at its own face with averaged scalars and SST, the other component
  averaged over four points).
- `ts`: the shared SST field (its halos are filled so the face averages are valid).
- Applied fluxes: `ρu`, `ρv` bottom flux `taux`, `tauy`; moisture `ρʷ lhf`; energy key `ρE`
  `cᵖᵐ ρʷ shf`, which Breeze converts to `Jᶿ = 𝒬/(cᵖᵐ Π) = ρʷ shf / Π`, i.e. the first cell's
  `T` changes by `ρʷ shf`, as SAM's `t` does at fixed condensate. `ρʷ` is SAM's `rhow(1)`
  recipe applied to Breeze's reference pressure (`sam_surface_interface_density`).

## Verification

`~/ena-les/lasso/oceflx_parity/`: the unchanged `oceflx.f90` and the `sat.f90` saturation
functions are compiled with gfortran behind `driver.f90` and evaluated on 8736 input
combinations (U = 0.5–25 m/s, air–sea ΔT = −5 to +3 K, RH = 60–100 %, z₁ = 10–40 m,
SST = 285, 292 and 269.5 K, p = 1015 mb), in double (`-fdefault-real-8`) and default single
precision; `parity.jl` evaluates `sam_ocean_fluxes` on the same inputs. Results are recorded in
`plans/status/ENA-LASSO.md`. `test/ena_lasso.jl` ("SAM oceflx surface law: parity with
oceflx.f90") holds six reference rows from the double-precision run and checks the model-level
boundary conditions against the function.

## Remaining differences from SAM

- **Inputs, not the algorithm.** Breeze's anelastic `ρθ` dynamics, advection and microphysics
  produce the lowest-level state; given the same state the fluxes agree to round-off.
- **Pressure.** SAM refreshes `pres` from the domain-mean state (`pressz.f90`); Breeze's
  reference pressure is fixed. A 1 hPa difference changes `ssq` by 0.1 % and `ρbot` by 0.1 %.
- **Vapor flux at fixed θ.** Breeze adds vapor at fixed `θ`, which changes `T` through the
  moist Exner function by a few per cent of `(cᵖᵛ − cᵖᵈ) T δqᵛ / cᵖᵐ`; SAM adds `qv` at fixed `t`.
- **Precision.** SAM is compiled in single precision (`FFLAGS = -O3`, no `-r8`); BreezeLab
  runs at `Oceananigans.defaults.FloatType`. `sam_ocean_fluxes` is generic: in `Float32` it
  rounds its constants as SAM does and agrees with the single-precision Fortran to Float32
  round-off (see the status note); in `Float64` it agrees with the `-fdefault-real-8` build to
  double round-off.
- **Model version.** The LASSO bundle README cites SAM hash `f83adf58`; this checkout's HEAD is
  `12d0244`. `f83adf58` is in the checkout's history (five commits back), and
  `git diff f83adf58 12d0244` leaves `oceflx.f90`, `surface.f90`, the diffusion routines,
  `setdata.f90`, `pressz.f90` and the `params.f90` constants unchanged; `sat.f90` only gains
  `get_inversion_height` and `saturation_adjustment_liquid`, and the other changes are to
  nudging and output. The version difference therefore does not touch the surface fluxes.
