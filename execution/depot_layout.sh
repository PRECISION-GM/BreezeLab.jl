# Shared package sources across per-CPU compiled-cache depots. Source this after CPU_TAG is set.
# Julia writes compiled caches AND installs package sources into DEPOT_PATH[1]; caches must be per-CPU
# (login/AMD/Intel nodes) but sources, artifacts, registries and dataset scratch caches must be shared,
# so those subdirectories of the per-CPU depot are symlinks into $W/depot/shared.
W=${W:-/shared/home/greg/breezelab-work/tracer-mip}
SHARED="$W/depot/shared"
DEPOT="$W/depot/$CPU_TAG"
mkdir -p "$SHARED" "$DEPOT"
for d in packages artifacts clones registries scratchspaces; do
    mkdir -p "$SHARED/$d"
    if [ -d "$DEPOT/$d" ] && [ ! -L "$DEPOT/$d" ]; then
        cp -an "$DEPOT/$d/." "$SHARED/$d/" 2>/dev/null || true   # merge what this depot already holds
        mv "$DEPOT/$d" "$DEPOT/$d.merged.$(date +%s)"
    fi
    [ -L "$DEPOT/$d" ] || ln -s "$SHARED/$d" "$DEPOT/$d"
done
