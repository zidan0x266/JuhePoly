# randomwalk — chain packing for end-linked networks

Precursor chains are grown as self-avoiding random walks on a face-centered cubic (FCC) lattice
and single-bead cross-linkers are put on random empty sites. The chains all have one length, or
their lengths follow a distribution (Gaussian, Schulz–Zimm, log-normal, Flory, Poisson or
Weibull–Nycander–Gold) with a given number-average length and, for the two-parameter shapes, a
given dispersity. The result is written as a LAMMPS data file, the starting point of
coarse-grained end-linking simulations with `fix bond/create`. The original program, the
*Amorphous Polymer Generator*, was obtained from Mark A. Tschopp and Don K. Ward, Center for
Advanced Vehicular Systems, Mississippi State University, who used it in [1, 2]; the end-linking
protocol follows Riggleman and co-workers [3–5]. This module built the monodisperse, end-linked
networks of the 2026 fracture paper listed in the [repository README](../../README.md). The
[gusev](../gusev/README.md) module builds such networks already linked, by Gusev's Monte Carlo
procedure, with the same keys for the chains and the same atom types.

## Run

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini chain_length=8 n_chains=12416 n_crosslinkers=6208 directory=endlink_No10
build/src/randomwalk/rw_pack examples/example_2/endlink_polydisperse.ini
build/src/randomwalk/rw_pack --print-config        # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly. The example inputs and the LAMMPS scripts that relax, compress and
end-link the packed system are described in [`examples/example_2`](../../examples/example_2/README.md).

## Configuration

Keys can be grouped under optional `[system]`, `[walk]`, `[run]` and `[output]` headers.

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `chain_length` | `18` | beads per precursor chain (at least 2); with a length distribution, the number-average chain length L<sub>n</sub>, which need not be a whole number. Both end beads of every chain are reactive |
| | `length_distribution` | `monodisperse` | `monodisperse`, `gaussian`, `schulz_zimm`, `lognormal`, `flory`, `poisson` or `weibull_nycander_gold`; see *Chain lengths* |
| | `dispersity` | `1` | dispersity Đ = L<sub>w</sub>/L<sub>n</sub> of the chain lengths; `1` gives monodisperse chains with a two-parameter distribution. Left out for `flory` and `poisson`, whose dispersity follows from `chain_length` |
| | `min_chain_length` | `2` | shortest chain of a distribution (at least 2) |
| | `max_chain_length` | `0` | longest chain of a distribution (`0` or at least `min_chain_length`); `0` limits it only by the number of chain beads, `n_chains` × `chain_length` |
| | `n_chains` | `5680` | number of precursor chains |
| | `n_crosslinkers` | `2840` | number of single-bead cross-linkers |
| | `functionality` | `4` | sites per cross-linker (1 to 255); the end-linking input of example_2 assumes 4 |
| | `bond_length` | `1.0` | bond length b, the nearest-neighbor distance of the lattice |
| | `mass` | `1.0` | mass of every atom type |
| | `occupancy` | `0.05` | sets the box: the smallest number of cells in which the chain beads fill less than this fraction (between 0 and 1) of the lattice sites |
| | `cells` | `0` | cubic cells per box edge (2 to 1000); `0` derives it from `occupancy` |
| **Walk** | `bond_angles` | `60, 90, 120` | bond angles a chain may form, from 60, 90, 120 and 180 degrees |
| | `site_weighting` | `legacy` | `legacy`, `free_neighbors` or `uniform`; see *Model* |
| | `max_attempts` | `1000` | growth attempts per chain before the run stops with an error |
| **Run** | `seed` | `0` | random seed; `0` picks a fresh one |
| **Output** | `directory` | `.` | output directory |
| | `data_file` | `config.lmp` | name of the LAMMPS data file |
| | `write_xyz` | `true` | also write an XYZ file for visualization |

## Output

- `config.lmp`: LAMMPS data file, atom style `molecular`, with molecule IDs and image flags.
  The chains come first, one after another: chain *c* (counted from 0) is molecule *c* + 1 and,
  with equal lengths *L*, holds atoms *c*·*L* + 1 to (*c* + 1)·*L*. Then come the cross-linkers,
  each a molecule of its own. Atom types carry the labels `bead` (1, chain interior), `end` (2,
  reactive chain end), `xlink` (3, cross-linker) and `xlink_full` (4, cross-linker with all its
  sites reacted, assigned during end-linking). Bonds have type 1.
- `chain_lengths.dat` (with a length distribution only): for every chain length, the number of
  chains drawn and the number expected from the distribution, after a header with L<sub>n</sub>,
  L<sub>w</sub> and Đ of both.
- `config.xyz`: extended XYZ file with the atom type and the wrapped position of every atom.
- `run_config.ini`: the configuration of the run.

The log reports the chain lengths (with a distribution, its solved parameters and the lengths
drawn), the lattice with its box edge, occupancy and number density, and the number of trapped
chains that were grown again.

## Model

Sites of the FCC lattice are one bond length apart, so every bond has length b and no two beads
are closer than b. The box edge is √2·b times the number of cells. Each chain starts on a random
empty site; every further bead goes to an empty neighbor of the last one that forms an allowed
bond angle with the previous bond, chosen according to `site_weighting`. A chain whose last bead
has no such neighbor is removed and grown again from a new random site. The cross-linkers take
random empty sites after all chains are in place.

The defaults reproduce the program used in the publications, including two of its quirks:

- `bond_angles`: the original program meant to allow 90 and 120 degrees only, but its
  floating-point angle test also lets 60-degree folds through, so the default is `60, 90, 120`.
- `site_weighting`: the original program prefers sites in less crowded regions. It weights each
  candidate by 12 minus the number of neighboring beads from which a chain has already grown,
  but it takes the weights, in its scan order, from the list of all empty neighbors and applies
  them to the shorter list of those with an allowed angle. `legacy` (default) reproduces this;
  `free_neighbors` weights each candidate by its own number of empty neighbors, as intended;
  `uniform` gives every candidate the same probability.

With the defaults, every candidate list and weight matches the original program step by step
(checked on traced runs), and bond-angle fractions, end-to-end distances and bead contacts agree
within statistical error. Three of its faults are not reproduced:

- Image flags are now written, so LAMMPS reads the packed system without warning about
  inconsistent image flags (after end-linking, bonds across the box bring the harmless warning
  back).
- A trapped chain is grown again elsewhere. The original program backed up one bead but kept
  growing next to the removed bead, which can leave a bond of twice the bond length, and stopped
  when this happened twice. In 40 test runs with the paper's settings, 3 contained such a bond
  and 1 stopped.
- A cross-linker can no longer land on an occupied site.

### Chain lengths

With a `length_distribution` other than `monodisperse`, the chain length *N* (in beads) follows a
discrete distribution *P*(*N*) on `min_chain_length` ≤ *N* ≤ `max_chain_length`, of one of six
shapes:

| `length_distribution` | *P*(*N*) ∝ | Parameters | Unbounded shape | Typical origin |
|:--|:--|:--|:--|:--|
| `gaussian` | exp[−(*N* − μ)²/(2σ²)] | μ, σ | Đ = 1 + σ²/L<sub>n</sub>² | narrow, symmetric distributions |
| `schulz_zimm` | *N*<sup>*k*−1</sup> exp(−*N*/θ) | *k*, θ | L<sub>n</sub> = *k*θ, Đ = 1 + 1/*k* | Schulz [6], Zimm [7] |
| `lognormal` | *N*<sup>−1</sup> exp[−(ln *N* − μ)²/(2σ²)] | μ, σ | L<sub>n</sub> = exp(μ + σ²/2), Đ = exp(σ²) | broad distributions (Wesslau [8]) |
| `flory` | *a*<sup>*N*</sup> | *a* | L<sub>n</sub> = 1/(1 − *a*), Đ = 1 + *a* | the most probable distribution of linear step growth or of radical chain growth with termination by disproportionation (Flory [9]) |
| `poisson` | *a*<sup>*N*−1</sup>/(*N* − 1)! | *a* | L<sub>n</sub> = *a* + 1, Đ = 1 + *a*/(*a* + 1)² | living polymerization with fast initiation (Flory [10]) |
| `weibull_nycander_gold` | (1/*r*) ∫<sub>0</sub><sup>τ</sup> e<sup>−(τ−*w*)/*r*</sup> Poisson(*N* − 1; *w*) d*w* | *v*, *r* | Đ from 1 up to 4/3 | living polymerization with slow initiation, *r* = *k*<sub>p</sub>/*k*<sub>i</sub> (Weibull and Nycander [11], Gold [12]) |

The Flory and Poisson distributions have one parameter, so `chain_length` alone sets them and
their dispersity follows: leave `dispersity` out. The other four have two, solved so that this
discrete, bounded distribution has exactly the number-average length `chain_length` and the
dispersity `dispersity`. The relations in the fourth column hold for the unbounded shapes only
(lengths from 1 on for Flory, Poisson and Weibull–Nycander–Gold, continuous for the others), and
they drift for short
chains: for L<sub>n</sub> = 18 and Đ = 2, *k* = 1 and θ = 18 give L<sub>n</sub> = 19.5 and
Đ = 1.85 once lengths below 2 are left out, while the solved Schulz–Zimm distribution has
*k* = 0.76; a Flory distribution of L<sub>n</sub> = 10 from 2 on has Đ = 1.72, not 1.9. The solved
parameters are written to the log and to `chain_lengths.dat`.

The Weibull–Nycander–Gold distribution describes a living polymerization with initiation and
propagation only. With τ = *k*<sub>p</sub>∫[M]d*t*, an initiator molecule that started a time *w*
before the end carries 1 + Poisson(*w*) monomer units, and initiation slows down as the initiator
is used up, which gives the integral in the table; *v* is the number-average length of all
initiator molecules, those not initiated counting as 0, and fixes τ through
*v* = τ − (*r* − 1)(1 − e<sup>−τ/*r*</sup>). Fast initiation (*r* → 0) gives the
Poisson distribution, slow initiation (*r* → ∞) broadens it up to Đ = 4/3 for long chains (Gold
[12]). *P*(*N*) is computed from exact recurrences between neighboring lengths, each run in its
stable direction; it reproduces Gold's table and agrees with numerical integration to 10<sup>−8</sup>
for *r* from 0.05 to 3.

At a given L<sub>n</sub>, each two-parameter shape reaches only a range of dispersities: with
*N* ≥ 2 and L<sub>n</sub> = 10, for instance, a Gaussian reaches up to Đ = 1.72, a
Weibull–Nycander–Gold distribution from 1.09 to 1.30 and a Schulz–Zimm distribution up to 2.52,
while a log-normal distribution goes further. A `chain_length` that is not a whole number also
sets a lowest dispersity. Beyond the range, `rw_pack` stops with a message that names the limit.

The lengths are drawn by stratified sampling: chain *i* of *n* takes the length at the quantile
(*i* + *u*)/*n* of the distribution, *u* uniform in [0, 1), and the lengths are then shuffled and
the chains grown in that order. The number of chains of every length then differs from its
expected number *n P*(*N*) by less than 2, so the chains drawn have nearly the requested
L<sub>n</sub> and Đ even when there are only a few thousand of them. The tail beyond the
1 − 1/*n* quantile gets a single chain, so with few chains of a broad distribution Đ usually comes
out slightly low and occasionally well above the target. Monodisperse lengths, including
`dispersity = 1`, take no random numbers, so they leave the random sequence of the packing
unchanged.

The box is sized from the chain beads drawn. The number of cross-linkers stays as given: for
end-linking with tetrafunctional cross-linkers it is half the number of chains, whatever their
lengths.

## Performance

Wall time and peak memory of a complete run, including the output files, on one core of an
Apple M4 with Apple Clang 21; the original program was compiled with the same optimization
flags and its screen output discarded.

| System | Original program | This code |
|:--|--:|--:|
| fracture paper, N<sub>o</sub> = 20 (105,080 beads) | 0.15 s, 100 MB | 0.024 s, 9 MB |
| ten times larger (1,050,800 beads) | 1.46 s, 470 MB | 0.26 s, 78 MB |

The original program is limited to 300 cells per edge and stops when a chain is trapped a
second time, which happened in every test at occupancies of 0.2 and above or with chains of 1000
beads. This code handles these cases and boxes of up to 1000 cells per edge.

## Code structure

```
include/rw/
  model/   the FCC lattice with its bond angles
  sim/     packing: chain growth and cross-linker placement, with the check of the result; the driver
  io/      configuration
lib/       implementation
apps/      rw_pack
```

randomwalk and [gusev](../gusev/README.md) are both packing modules and share their structure: the
end-linking system they build, with its statistics and consistency check, its LAMMPS data, XYZ and
link writers, the keys of the chains, cross-linkers, seed and output, and the steps of a run are
in [`src/common`](../common/README.md) (`poly/endlink/`), so each module holds only its method;
randomwalk keeps its positions as lattice points. Unit tests are in
[`tests/randomwalk`](../../tests/randomwalk). The chain-length distributions are shared with
gusev, the FCC geometry with dll (which takes its random-walk starting configurations from this
module), and the random number generator and the configuration reader with all modules.

## References

1. D. Hossain, M. A. Tschopp, D. K. Ward, J. L. Bouvard, P. Wang, M. F. Horstemeyer. Molecular
   dynamics simulations of deformation mechanisms of amorphous polyethylene. *Polymer* **2010**,
   51, 6071–6083. [10.1016/j.polymer.2010.10.009](https://doi.org/10.1016/j.polymer.2010.10.009)
2. M. A. Tschopp, J. L. Bouvard, D. K. Ward, M. F. Horstemeyer. Atomic scale deformation
   mechanisms of amorphous polyethylene under tensile loading. In *Supplemental Proceedings:
   Materials Fabrication, Properties, Characterization, and Modeling* (TMS 2011); Wiley, **2011**;
   pp 789–794. [10.1002/9781118062142.ch95](https://doi.org/10.1002/9781118062142.ch95)
3. Z. Ye, R. A. Riggleman. Molecular View of Cavitation in Model-Solvated Polymer Networks.
   *Macromolecules* **2020**, 53, 7825–7834.
   [10.1021/acs.macromol.0c00994](https://doi.org/10.1021/acs.macromol.0c00994)
4. C. W. Barney, Z. Ye, I. Sacligil, K. R. McLeod, H. Zhang, G. N. Tew, R. A. Riggleman,
   A. J. Crosby. Fracture of model end-linked networks. *Proc. Natl. Acad. Sci. U.S.A.* **2022**,
   119, e2112389119. [10.1073/pnas.2112389119](https://doi.org/10.1073/pnas.2112389119)
5. H. Zhang, R. A. Riggleman. Predicting failure locations in model end-linked polymer networks.
   *Phys. Rev. Mater.* **2024**, 8, 035604.
   [10.1103/PhysRevMaterials.8.035604](https://doi.org/10.1103/PhysRevMaterials.8.035604)
6. G. V. Schulz. Über die Kinetik der Kettenpolymerisationen. V. *Z. Phys. Chem.* **1939**, 43B,
   25–46. [10.1515/zpch-1939-4304](https://doi.org/10.1515/zpch-1939-4304)
7. B. H. Zimm. Apparatus and Methods for Measurement and Interpretation of the Angular Variation of
   Light Scattering; Preliminary Results on Polystyrene Solutions. *J. Chem. Phys.* **1948**, 16,
   1099–1116. [10.1063/1.1746740](https://doi.org/10.1063/1.1746740)
8. H. Wesslau. Die Molekulargewichtsverteilung einiger Niederdruckpolyäthylene. *Makromol. Chem.*
   **1956**, 20, 111–142. [10.1002/macp.1956.020200109](https://doi.org/10.1002/macp.1956.020200109)
9. P. J. Flory. Molecular Size Distribution in Linear Condensation Polymers. *J. Am. Chem. Soc.*
   **1936**, 58, 1877–1885. [10.1021/ja01301a016](https://doi.org/10.1021/ja01301a016)
10. P. J. Flory. Molecular Size Distribution in Ethylene Oxide Polymers. *J. Am. Chem. Soc.*
    **1940**, 62, 1561–1565. [10.1021/ja01863a066](https://doi.org/10.1021/ja01863a066)
11. B. Weibull, B. Nycander. The Distribution of Compounds Formed in the Reaction between Ethylene
    Oxide and Water, Ethanol, Ethylene Glycol, or Ethylene Glycol Monoethyl Ether. *Acta Chem.
    Scand.* **1954**, 8, 847–858. [10.3891/acta.chem.scand.08-0847](https://doi.org/10.3891/acta.chem.scand.08-0847)
12. L. Gold. Statistics of Polymer Molecular Size Distribution for an Invariant Number of Propagating
    Chains. *J. Chem. Phys.* **1958**, 28, 91–99. [10.1063/1.1744088](https://doi.org/10.1063/1.1744088)
