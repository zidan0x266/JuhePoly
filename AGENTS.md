# AGENTS.md

JuhePoly is a C++23 simulation toolkit for polymerization (FRP and ATRP chemistry, end-linked networks). The repo is scientific code, so the umbrella rules in `~/Work/Science/AGENTS.md` apply: never fabricate results, and record command, config, seed and commit for any experiment-backed claim.

## Build and test

CMake >= 3.20 and a C++23 compiler (GCC 14+ or Clang 18+). Everything builds and tests from the repo root.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build                          # all unit + smoke tests
ctest --test-dir build -R '^dll\.'              # one module
ctest --test-dir build -R 'statistical.smoke_frp' --output-on-failure   # one test
build/src/3dbfm/bfm_react --print-config        # list every config key with defaults
```

Useful options: `-DPOLY_WARNINGS_AS_ERRORS=ON`, `-DPOLY_CXX_STANDARD=26`, `-DPOLY_NATIVE_ARCH=ON`, `-DPOLY_LTO=OFF`, `-DPOLY_BUILD_TESTS=OFF`, `-DPOLY_PGO=GENERATE|USE`. Warnings are strict (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion`); there is no linter or formatter config. Executables land in `build/src/<module>/`.

## Architecture

One directory per model under `src/`, each with `include/<ns>/`, `lib/`, `apps/` (the executable) and a README documenting its config keys. Tests live in `tests/<module>/test_main.cpp`.

- `src/common/`: code used by at least two modules and ignorant of any single one (`poly::` namespaces). Code moves here only when a second module needs it. Static libs: `poly::chemistry`, `poly::chains`, `poly::endlink`. The rest is header-only (`poly::common`): RNG (xoshiro256++), INI reader, FCC lattice geometry, union-find molecules, bond tables, PDB/PSF/LAMMPS-data writers.
- `3dbfm` (`bfm_react`) and `dll` (`dll_react`): reactive lattice Monte Carlo (bond fluctuation model, dynamic lattice liquid). `statistical` (`stat_react`): coordinate-free Gillespie kinetics. All three share the reaction schemes in `poly/chem/`.
- `randomwalk` (`rw_pack`) and `gusev` (`gusev_build`): build end-linking starting configurations for LAMMPS. They share `poly/endlink/` and the chain-length distributions in `poly/chains/`.
- `src/lammps/`: `fix poly/react` for LAMMPS, **GPL-2.0-only** (the rest is PolyForm Noncommercial), and it must not use code from `src/common/`. `polyreact_core` (reaction-file parser and rule engine) is tested without LAMMPS. The plugin `polyreactplugin` is only built when `-DPOLYREACT_LAMMPS_SOURCE=<lammps/src>` is set. It must match that LAMMPS's compiler, `LAMMPS_SIZES` and C++ standard (`POLYREACT_LAMMPS_CXX_STANDARD`). The `lammps.plugin` ctest also needs `-DPOLYREACT_LAMMPS_EXECUTABLE=...`.

Adding things: a new simulation method is a new module under `src/` (register it in the root `CMakeLists.txt` and `tests/CMakeLists.txt`). A new chemistry is a new reaction scheme in `poly/chem/` plus its rate expressions in the statistical model. A new polymer system is a new example, not new code.

## Configuration and examples

Every executable is driven by an INI file with `key=value` command-line overrides (`bfm_react examples/example_1/atrp_paper.ini seed=12345`). Each run writes `run_config.ini` with the full resolved config including the seed, so runs are exactly repeatable. Inputs live in `examples/example_<n>/` (1 = 3dbfm, 2 = randomwalk, 3 = gusev, 4 = dll, 5 = statistical, 6 = LAMMPS), each with a README and a `smoke_small.ini` used by the smoke tests.

## Tests

Unit tests are plain executables (`<module>_tests`) that return non-zero on failure, with no test framework. Smoke tests in `tests/CMakeLists.txt` run the real binaries on the small `examples/*/smoke_*.ini` systems and write into `build/smoke/`. A new module needs both a unit test executable and smoke tests registered there.
