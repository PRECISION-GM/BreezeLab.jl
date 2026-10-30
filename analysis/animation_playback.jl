# Frame rate of an animation from the spacing of its saved frames and a target playback speed.
#
# Movie standard (Greg, 2026-10-13): cloud-scale movies play at no more than 15 simulated minutes
# per second of video and are encoded at 10–15 fps, which needs 2-D slices saved every 1 min.
#
#   ANIMATION_SPEED    target playback speed in simulated minutes per video second (default 15);
#                      an upper bound: the frame rate is rounded down so the movie is never faster
#   ANIMATION_FPS_MAX  highest frame rate used (default 15); ANIMATION_FPS_MIN lowest acceptable (10)
#   ANIMATION_FPS      an explicit frame rate overriding the rule (the resulting speed is reported)
#
# frames saved every Δt seconds at `fps` frames per second play at fps × Δt / 60 simulated minutes
# per second, so the rule is fps = min(fps_max, floor(60 speed / Δt)). Frames sparser than
# 60 speed / fps_min seconds cannot meet both bounds; the speed bound wins and a warning names the
# slice interval needed.
using Statistics: median
using Printf: @sprintf

function playback_framerate(times; speed = parse(Float64, get(ENV, "ANIMATION_SPEED", "15")),
                                   fps_max = parse(Int, get(ENV, "ANIMATION_FPS_MAX", "15")),
                                   fps_min = parse(Int, get(ENV, "ANIMATION_FPS_MIN", "10")),
                                   explicit = get(ENV, "ANIMATION_FPS", nothing))
    length(times) ≥ 2 || return (; framerate = fps_min, frame_spacing = NaN, playback_speed = NaN)
    Δt = median(diff(collect(Float64, times)))                     # s between saved frames
    framerate = if isnothing(explicit)
        max(1, min(fps_max, floor(Int, 60speed / Δt + 1e-9)))
    else
        parse(Int, explicit)
    end
    playback_speed = framerate * Δt / 60                            # simulated minutes per video second
    if isnothing(explicit) && framerate < fps_min
        @warn "frames every $(Δt) s are too sparse for ≥ $fps_min fps at ≤ $speed sim-min/s; " *
              "encoding at $framerate fps ($(round(playback_speed, digits=1)) sim-min/s). " *
              "Save slices every ≤ $(round(60speed / fps_min, digits=1)) s to meet the standard."
    elseif playback_speed > speed + 1e-9
        @warn "ANIMATION_FPS = $framerate plays at $(round(playback_speed, digits=1)) sim-min/s (> $speed)"
    end
    return (; framerate, frame_spacing = Δt, playback_speed)
end

playback_note(p) = @sprintf("%d fps, frames every %.0f s → %.1f simulated min per second of video",
                            p.framerate, p.frame_spacing, p.playback_speed)
