#!/bin/bash
# Palace launcher used by tile_lib.run_palace (Linux / WSL).  usage: bash palace.sh [-np N] config.json
# Palace 0.18.1 was installed with Spack; adapt the two lines below to your install (or set PALACE_SH to your own launcher).
export HWLOC_COMPONENTS=-gl      # WSL2: stops hwloc/OpenMPI hanging on X11 display probing
source ~/spack/share/spack/setup-env.sh && spack load palace
exec palace "$@"
