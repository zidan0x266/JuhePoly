# common — code shared by the modules

Everything here is used by at least two modules and knows nothing about any one of them. Code
moves here when a second module needs it; until then it stays in its module. The LAMMPS fix in
[`src/lammps`](../lammps) is licensed separately and uses nothing from this directory.

| Folder | Namespace | Contents | Used by |
|:--|:--|:--|:--|
| `poly/` | `poly`, `poly::ini` | random numbers (xoshiro256++), configuration-file reader, index set | 3dbfm, randomwalk, gusev, dll, statistical |
| `poly/geometry/` | `poly::geometry` | integer and real vectors, periodic face-centered cubic (FCC) lattice | 3dbfm, randomwalk, gusev, dll |
| `poly/topology/` | `poly::topology` | molecules as disjoint sets of units (union-find), bond lists of fixed capacity | 3dbfm, dll, statistical |
| `poly/chem/` | `poly::chem` | unit states, reaction schemes (FRP, ATRP), reaction bookkeeping, per-sample output tables | 3dbfm, dll, statistical |
| `poly/chains/` | `poly::chains` | chain-length distributions (Gaussian, Schulz–Zimm, log-normal, Flory, Poisson, Weibull–Nycander–Gold) of given number average and dispersity, their configuration keys and the table of lengths drawn | randomwalk, gusev |
| `poly/endlink/` | `poly::endlink` | the end-linking system that the packing modules build (chains, cross-linkers, links with their periodic images, lattice or continuous positions), its statistics and check, its LAMMPS data, XYZ and link writers, the shared keys and the steps of a run | randomwalk, gusev |
| `poly/io/` | `poly::io` | output files that report failure, PDB and PSF writers, LAMMPS data header | 3dbfm, randomwalk, gusev, dll, statistical |

The chemistry, the chain-length distributions and the end-linking systems are the static libraries
`poly::chemistry`, `poly::chains` and `poly::endlink` (`lib/`); the rest is header-only and comes
with `poly::common`.

The three directions the toolkit grows in map onto this layout: a new simulation method is a new
module under `src/`, a new chemistry is a new reaction scheme in `poly/chem/` (with its rate
expressions in the statistical model), and a new polymer system (linear, branched, network,
copolymer) is a new example, not new code. Unit tests are in [`tests/common`](../../tests/common).
