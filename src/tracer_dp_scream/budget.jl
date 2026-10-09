#####
##### Column water and moist-enthalpy budgets of a TRACER–DP-SCREAM run from its saved output.
#####
##### Water: W = PW + LWP + RWP + IWP (domain means, kg m⁻²), dW/dt = E − P + ∫ρ qls dz.
#####
##### Energy: the column static energy S = ∫ρ s dz with s = cᵖᵐ T + g z − ℒˡᵣ qˡ − ℒⁱᵣ qⁱ is *not* conserved
##### by phase changes or by precipitation leaving the column (removing condensate removes its
##### −ℒᵣ q contribution), so the earlier S-budget (which omitted those terms) could not close in a
##### precipitating case. The budget here uses the moist enthalpy of Breeze's thermodynamics,
#####
#####     h = cᵖᵐ (T − Tᵣ) + ℒˡᵣ qᵛ + (ℒˡᵣ − ℒⁱᵣ) qⁱ,
#####
##### which phase changes conserve (ℒˡ(T) = ℒˡᵣ + (cᵖᵛ − cˡ)(T − Tᵣ), likewise for ice), so that
#####
#####     H = ∫ρ (h + g z) dz = S + ℒˡᵣ W − Tᵣ [(cᵖᵛ − cᵖᵈ) PW + (cˡ − cᵖᵈ)(LWP + RWP) + (cⁱ − cᵖᵈ) IWP] + const
#####
##### is computed exactly from the 60-s time series. Its sources (anelastic, fixed total density:
##### departing or arriving water is exchanged with dry air at the local temperature) are
#####
#####   surface   F_E + E [ℒˡᵣ + (cᵖᵛ − cᵖᵈ)(T₁ − Tᵣ)]          (F_E: the heating the ρE flux applies)
#####   radiation ∫ −∂F/∂z dz
#####   forcing   ∫ρ [F_tls + (ℒˡᵣ + (cᵖᵛ − cᵖᵈ)(T − Tᵣ)) R] dz   (F_tls: the heating of the energy forcing,
#####                                                             R: its vapor source)
#####   precip.   −(cˡ − cᵖᵈ)(T₁ − Tᵣ) Pʳ − [(cⁱ − cᵖᵈ)(T₁ − Tᵣ) − (ℒˡᵣ − ℒⁱᵣ)] Pⁱ
#####
##### with T₁ the lowest-level temperature. In the potential-temperature formulation an energy
##### source F heats as cᵖᵐ dT = F and a vapor source at fixed θ leaves T unchanged, so a
##### `(cᵖᵛ − cᵖᵈ) T ×` vapor-rate "temperature-neutral" term in F (designed for the static-energy
##### formulation) is a net heating there; `applied_vapor_cross_term` states whether the run's
##### forcing and surface flux contained it, and the budget reports its integral.
##### Profiles (30-min means) supply T and the condensate for cᵖᵐ and the cross terms.
#####

using JLD2: jldopen
using Breeze: ThermodynamicConstants

function saved_series(file, names; interior_range = nothing)
    jldopen(file) do f
        ts = f["timeseries"]
        iterations = sort(parse.(Int, collect(keys(ts["t"]))))
        times = Float64[ts["t/$i"] for i in iterations]
        values = map(names) do name
            haskey(ts, string(name)) || return nothing
            if isnothing(interior_range)
                Float64[ts["$name/$i"][1, 1, 1] for i in iterations]
            else
                reduce(hcat, [Float64.(vec(ts["$name/$i"])[interior_range]) for i in iterations])
            end
        end
        return times, NamedTuple{Tuple(names)}(Tuple(values))
    end
end

held(times, values, t) = values[clamp(searchsortedlast(times, t), 1, length(values))]

"""
    tracer_dp_scream_budget(case, timeseries_file, profiles_file; applied_vapor_cross_term = true,
                            interval = 86400)

Water and moist-enthalpy budgets (see the file header) of a TRACER–DP-SCREAM run, totals and per
`interval` (default daily). `case` is a `tracer_dp_scream(...)` bundle with the run's vertical
grid and window (it supplies ρᵣ, the forcing profiles and the surface series; any horizontal size).
`applied_vapor_cross_term` states whether the run's energy forcing and surface energy flux included
the `(cᵖᵛ − cᵖᵈ) T ×` vapor-rate term; by default it is read from `case` (pass `true` explicitly when
analysing a run made before the formulation-aware fix with a `case` built from newer code). Also returns
the legacy static-energy residual (precipitation and phase-change terms excluded) for comparison.
"""
# Whether the case's energy forcing applies the vapor heat-capacity cross term (it did in every
# θ-formulation model before the forcing became formulation-aware; the surface flux was fixed together).
function applies_vapor_cross_term(case)
    found = Ref{Any}(nothing)
    visit(f) = f isa LargeScaleEnergyForcing ? (found[] = f) :
               f isa Breeze.Forcings.SpecificForcing ? visit(f.forcing) :
               f isa Oceananigans.Forcings.MultipleForcings ? foreach(visit, f.forcings) : nothing
    foreach(visit, values(case.model.forcing))
    f = found[]
    isnothing(f) && return false
    return hasfield(typeof(f), :static_energy) ? f.static_energy : true
end

function tracer_dp_scream_budget(case, timeseries_file, profiles_file; applied_vapor_cross_term = applies_vapor_cross_term(case),
                                 interval = 86400)
    constants = ThermodynamicConstants(Float64)
    ℒˡ = constants.liquid.reference_latent_heat
    ℒⁱ = constants.ice.reference_latent_heat
    cᵖᵈ = constants.dry_air.heat_capacity
    cᵖᵛ = constants.vapor.heat_capacity
    cˡ = constants.liquid.heat_capacity
    cⁱ = constants.ice.heat_capacity
    Tᵣ = constants.energy_reference_temperature

    t, s = saved_series(timeseries_file, (:lwp, :rwp, :iwp, :precipitable_water, :column_static_energy,
                                          :column_radiative_heating, :rain_flux, :ice_flux))
    zero_series = zeros(length(t))
    IWP = something(s.iwp, zero_series); Pⁱ = something(s.ice_flux, zero_series)
    R = something(s.column_radiative_heating, zero_series)
    PW, LWP, RWP, S, Pʳ = s.precipitable_water, s.lwp, s.rwp, s.column_static_energy, s.rain_flux
    W = PW .+ LWP .+ RWP .+ IWP
    H = S .+ ℒˡ .* W .- Tᵣ .* ((cᵖᵛ - cᵖᵈ) .* PW .+ (cˡ - cᵖᵈ) .* (LWP .+ RWP) .+ (cⁱ - cᵖᵈ) .* IWP)

    grid = case.grid
    Nz = size(grid, 3)
    Hz = grid.Hz
    ρᵣ = Float64.(Array(interior(case.model.dynamics.reference_state.density, 1, 1, :)))
    Δz = diff(Float64.(Array(znodes(grid, Face()))))
    tp, p = saved_series(profiles_file, (:T, :qᵛ, :qᶜˡ, :qʳ, :qⁱ); interior_range = (Hz + 1):(Hz + Nz))
    qⁱp = something(p.qⁱ, zeros(Nz, length(tp)))
    profile_at(tt) = clamp(searchsortedfirst(tp, tt), 1, length(tp))      # window (tpₘ₋₁, tpₘ] mean

    tls, qls = case.forcing_profiles.tls, case.forcing_profiles.qls
    forcing_column(fts, n) = Float64.(Array(interior(fts[n], 1, 1, :)))
    tls_columns = [forcing_column(tls, n) for n in eachindex(tls.times)]
    qls_columns = [forcing_column(qls, n) for n in eachindex(qls.times)]
    sfc = case.sfc
    sfc_times = day_to_seconds.(sfc.day, case.config.day0)
    cross = applied_vapor_cross_term ? 1.0 : 0.0

    N = length(t)
    Sq = zeros(N); Fls = zeros(N); Hls = zeros(N); Xls = zeros(N)
    Fsfc = zeros(N); Hsfc = zeros(N); Xsfc = zeros(N); E = zeros(N); Hpr = zeros(N); Fls_cpd = zeros(N)
    for (n, tt) in enumerate(t)
        nf = clamp(searchsortedlast(tls.times, tt), 1, length(tls.times))
        m = profile_at(tt)
        T = p.T[:, m]; qᵛ = p.qᵛ[:, m]; qˡ = p.qᶜˡ[:, m] .+ p.qʳ[:, m]; qⁱ = qⁱp[:, m]
        cᵖᵐ = cᵖᵈ .+ (cᵖᵛ - cᵖᵈ) .* qᵛ .+ (cˡ - cᵖᵈ) .* qˡ .+ (cⁱ - cᵖᵈ) .* qⁱ
        τ = tls_columns[nf]; r = qls_columns[nf]           # mass-fraction basis: the vapor rate is qls
        w = ρᵣ .* Δz
        Sq[n] = sum(w .* r)
        Xls[n] = sum(w .* (cᵖᵛ - cᵖᵈ) .* T .* r)
        Fls[n] = sum(w .* cᵖᵐ .* τ) + cross * Xls[n]
        Fls_cpd[n] = cᵖᵈ * sum(w .* τ)
        Hls[n] = Fls[n] + sum(w .* (ℒˡ .+ (cᵖᵛ - cᵖᵈ) .* (T .- Tᵣ)) .* r)
        T₁ = T[1]
        Hs = held(sfc_times, sfc.sensible_heat_flux, tt); LE = held(sfc_times, sfc.latent_heat_flux, tt)
        Tg = held(sfc_times, sfc.sst, tt)
        E[n] = LE / ℒˡ
        Xsfc[n] = (cᵖᵛ - cᵖᵈ) * Tg * E[n]
        Fsfc[n] = Hs + cross * Xsfc[n]
        Hsfc[n] = Fsfc[n] + E[n] * (ℒˡ + (cᵖᵛ - cᵖᵈ) * (T₁ - Tᵣ))
        Hpr[n] = -(cˡ - cᵖᵈ) * (T₁ - Tᵣ) * Pʳ[n] - ((cⁱ - cᵖᵈ) * (T₁ - Tᵣ) - (ℒˡ - ℒⁱ)) * Pⁱ[n]
    end
    P = Pʳ .+ Pⁱ

    integral(f, a, b) = sum((f[k] + f[k + 1]) / 2 * (t[k + 1] - t[k]) for k in a:b-1; init = 0.0)
    function window(a, b)
        sources_h = (; surface = integral(Hsfc, a, b), radiation = integral(R, a, b),
                       large_scale = integral(Hls, a, b), precipitation = integral(Hpr, a, b))
        ΔH = H[b] - H[a]
        residual_h = ΔH - sum(values(sources_h))
        scale_h = integral(abs.(Hsfc), a, b) + integral(abs.(R), a, b) + integral(abs.(Hls), a, b) + integral(abs.(Hpr), a, b)
        ΔW = W[b] - W[a]
        water = (; ΔW, evaporation = integral(E, a, b), precipitation = integral(P, a, b), large_scale = integral(Sq, a, b))
        residual_w = ΔW - (water.evaporation - water.precipitation + water.large_scale)
        legacy_sources = integral(Fsfc, a, b) + integral(R, a, b) + integral(Fls_cpd, a, b)
        legacy_residual = (S[b] - S[a]) - legacy_sources
        return (; t_start = t[a], t_stop = t[b],
                  energy = merge((; ΔH, residual = residual_h, relative = abs(residual_h) / max(scale_h, eps())), sources_h),
                  water = merge(water, (; residual = residual_w,
                                          relative = abs(residual_w) / max(water.evaporation + water.precipitation + abs(water.large_scale), eps()))),
                  legacy_static_energy_residual = legacy_residual,
                  latent_heat_of_precipitation = ℒˡ * water.precipitation,
                  vapor_cross_term_surface = integral(Xsfc, a, b), vapor_cross_term_large_scale = integral(Xls, a, b),
                  mean_surface_heating_W_m2 = integral(Fsfc, a, b) / (t[b] - t[a]))
    end
    edges = [findfirst(≥(e), t) for e in 0:interval:t[end]]
    edges = unique(filter(!isnothing, edges))
    last(edges) == N || push!(edges, N)
    periods = [window(edges[k], edges[k + 1]) for k in 1:length(edges) - 1]
    return (; total = window(1, N), periods, applied_vapor_cross_term,
              constants = (; ℒˡ, ℒⁱ, cᵖᵈ, cᵖᵛ, cˡ, cⁱ, Tᵣ))
end
