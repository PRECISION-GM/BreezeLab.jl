#####
##### Checkpoint behaviour of BreezeLab callbacks
#####
#
# `Oceananigans.Checkpointer` saves `prognostic_state(simulation)`, which descends into every
# callback's `func`. The generic fallback `prognostic_state(obj) = obj` would serialize these
# callback structs (with their GPU fields) and then has no matching `restore_prognostic_state!`
# method on pickup. The three BreezeLab callbacks hold no evolving state of their own: they
# recompute surface stress / SST / the diagnostic CCN projection from the model clock and the
# prognostic fields on their next call. Declaring them stateless makes checkpoints skip them and
# makes `run!(simulation; pickup=true)` restore cleanly. After a pickup the first call of each
# callback refreshes its field (see docs/cases/ena.md, restart notes).

Oceananigans.prognostic_state(::PrescribedStressUpdater) = nothing
Oceananigans.prognostic_state(::SeaSurfaceTemperatureUpdater) = nothing
Oceananigans.prognostic_state(::DiagnosticCCNProjection) = nothing
Oceananigans.prognostic_state(::ProgressMessenger) = nothing   # wall-clock reference only
