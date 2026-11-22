# gusev — end-linked networks by Gusev's Monte Carlo procedure

This module builds an end-linked polymer network directly, as a starting structure for molecular
dynamics: f-functional cross-linkers are placed in a periodic box, linear strands with two reactive
ends are linked to their sites by the Monte Carlo procedure of Gusev [1, 2], and the beads of every
connected strand are then placed as a Brownian bridge between its two cross-linkers [3, 4], those
of the other strands as random walks. It serves the same purpose as the
[randomwalk](../randomwalk/README.md) module, which packs the unlinked precursors for end-linking
in the simulation itself: both take the same keys for the chains and write LAMMPS data files with
the same atom types, so the two routes can be compared directly. The code
implements the published descriptions and is checked against the networks that Argun and Statt
built with the same procedure [5]; Sen and Olsen compared Gusev networks with those of other
network generators [6].

## Run

```bash
build/src/gusev/gusev_build examples/example_3/endlink_melt.ini
build/src/gusev/gusev_build examples/example_3/endlink_melt.ini procedure=strands conversion=0.99 directory=gusev_strands
build/src/gusev/gusev_build --print-config        # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly. The LAMMPS input that relaxes and equilibrates the network is described
in [`examples/example_3`](../../examples/example_3/README.md).

## Configuration

Keys can be grouped under optional `[system]`, `[network]`, `[run]` and `[output]` headers.

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `chain_length` | `18` | beads per strand (at least 2); with a length distribution, the number-average length |
| | `length_distribution` | `monodisperse` | `monodisperse`, `gaussian`, `schulz_zimm`, `lognormal`, `flory`, `poisson` or `weibull_nycander_gold`, as in randomwalk (see its [*Chain lengths*](../randomwalk/README.md#chain-lengths)) |
| | `dispersity` | `1` | dispersity Đ = L<sub>w</sub>/L<sub>n</sub> of the strand lengths; left out for `flory` and `poisson` |
| | `min_chain_length` | `2` | shortest strand of a distribution (at least 2) |
| | `max_chain_length` | `0` | longest strand of a distribution (`0` or at least `min_chain_length`); `0` limits it only by the number of strand beads |
| | `n_chains` | `5680` | number of strands |
| | `n_crosslinkers` | `2840` | number of cross-linkers (at least 1); the stoichiometric imbalance is r = f `n_crosslinkers` / (2 `n_chains`) |
| | `functionality` | `4` | sites per cross-linker, f (1 to 255) |
| | `bond_length` | `1.0` | b: a strand of N segments has the mean-square end-to-end distance N b<sup>2</sup> and reaches at most N b; the steps of the Brownian bridges have mean-square length b<sup>2</sup> |
| | `mass` | `1.0` | mass of every atom type |
| | `density` | `0.85` | number density of all beads, strands and cross-linkers; sets the box unless `box` is given |
| | `box` | `0` | box edge; `0` derives it from `density` |
| **Network** | `procedure` | `ends` | `ends` or `strands`; see *Linking* |
| | `conversion` | `0.95` | `ends`: the target conversion of the minority groups; `strands`: the probability ξ that an end reacts; above 0 and at most 1 |
| **Run** | `seed` | `0` | random seed; `0` picks a fresh one |
| **Output** | `directory` | `.` | output directory |
| | `data_file` | `config.lmp` | name of the LAMMPS data file |
| | `write_xyz` | `true` | also write an XYZ file for visualization |

## Output

- `config.lmp`: LAMMPS data file, atom style `molecular`, with molecule IDs and image flags.
  The strands come first, one after another from end 0 to end 1: strand *s* (counted from 0) is
  molecule *s* + 1 and, with equal lengths *L*, holds atoms *s*·*L* + 1 to (*s* + 1)·*L*. Then come
  the cross-linkers, each a molecule of its own. Atom types are those of example_2 after
  end-linking: `bead` (1, strand beads, including linked ends), `end` (2, unlinked strand ends),
  `xlink` (3, cross-linker) and `xlink_full` (4, cross-linker with all f sites linked). Bonds have
  type 1: along the strands and from every linked end to its cross-linker. LAMMPS warns about
  inconsistent image flags: an end linked to a periodic image of its cross-linker is bonded across
  the box, and a network that spans the box cannot be unwrapped consistently. Bonds use the
  minimum image, so the warning is harmless.
- `strands.dat`: for every strand, its number (*s* + 1), its number of beads, the atom IDs of the
  cross-linkers of its two ends (0 for an unlinked end) and the vector from the first cross-linker
  to the second through the periodic image it was linked to (0 0 0 unless both ends are linked).
- `chain_lengths.dat` (with a length distribution only): the strand lengths drawn against the
  distribution, as in [randomwalk](../randomwalk/README.md#output).
- `config.xyz`: extended XYZ file with the atom type and the wrapped position of every atom.
- `run_config.ini`: the configuration of the run.

The log reports the conversion reached and the ends that found no active site within reach; the
numbers of connected, dangling and free strands, of primary loops and loops through a periodic
image, and of fully linked cross-linkers; and the mean R<sup>2</sup>/(N b<sup>2</sup>) of the
connected strands, primary loops (R = 0) included.

## Model

### Linking

A strand of n beads has N = n + 1 segments between its two cross-linkers. The cross-linkers take
random positions in the box and stay there, and excluded volume is ignored, as in Gusev's
procedure [1, 2]:

- An end whose strand has no linked end yet reacts with an active site chosen at random from all
  active sites in the box, so a cross-linker is chosen in proportion to its active sites.
- An end whose other end is linked at X reacts with an active site within the contour length N b
  of X, chosen with probability proportional to exp(−3d<sup>2</sup>/(2N b<sup>2</sup>)), the Gaussian
  distribution of the end-to-end vector. The candidates include the remaining sites of X itself
  (primary loops, d = 0) and all periodic images within N b, which matter once a strand can reach
  across half the box. Images whose weight is below e<sup>−37</sup> (10<sup>−16</sup>) of the
  closest one's are left out; even thousands of them together change the probabilities by less
  than 10<sup>−11</sup>.
- An end that finds no active site within reach stays unlinked; the log counts such ends.

Every strand is linked with its own N, so with a length distribution the topology follows the
lengths: short strands find their own cross-linker more often and long ones reach farther. In
mini networks of 2000 Schulz–Zimm strands (L<sub>n</sub> = 10, Đ = 1.5) with 1000 tetrafunctional
cross-linkers at conversion 0.95 (ten seeds), the fraction of strands that form a primary loop
falls from 0.31 for strands of 2 to 5 beads to 0.16 for 6 to 9, 0.11 for 10 to 14, 0.07 for 15
to 24 and 0.03 for longer ones, while R<sup>2</sup>/(N b<sup>2</sup>) of the other connected
strands falls from 1.31 to 1.15, 1.09, 1.07 and 1.01. With the primary loops included, the mean
in the log is 0.91 to 0.97.

`procedure` selects the order in which the ends react:

- `ends` (Gusev and Schwarz [2]): every step reacts an end chosen at random from all remaining
  active ends, until round(`conversion` × M) ends are linked, M = min(2 `n_chains`, f
  `n_crosslinkers`) being the reactive groups of the minority component, or until no end can react.
- `strands` (Gusev [1], as used by Argun and Statt [5]): the strands are taken one after another.
  A random active site is chosen and end 0 reacts with it with probability ξ = `conversion`; the
  site for end 1 is then chosen around that cross-linker, whether end 0 reacted or not, and end 1
  reacts with probability ξ. With active sites in reach, the free, dangling and connected strands
  are (1 − ξ)<sup>2</sup>, 2ξ(1 − ξ) and ξ<sup>2</sup> of all strands.

Sen and Olsen's implementation [6] differs in three points that this code does not follow: it
weights cross-linkers rather than their active sites, it does not try end 1 when end 0 has not
reacted, and its network file leaves out loops, dangling and free strands.

### Beads

The beads of a strand follow a Gaussian random walk whose steps have mean-square length
b<sup>2</sup> (each component normal with variance b<sup>2</sup>/3). A connected strand is a
Brownian bridge from its first cross-linker to the image of the second [3, 4], sampled exactly:
with m steps left from x to the target y, the next bead is normal with mean x + (y − x)/m and
variance (b<sup>2</sup>/3)(m − 1)/m per component. A dangling strand is a free walk from its
cross-linker, a free strand one from a random point. The bond lengths therefore scatter around b
and beads overlap, so the network needs a push-off before a Lennard-Jones model can run; the
LAMMPS input of example_3 ramps up a soft repulsion under capped displacements first.

### Checks

The unit tests compare the image search with brute force, the choice of the second site with its
exact probabilities in a small configuration, the conversion of the `ends` procedure with its
target, the strand fractions of the `strands` procedure with (1 − ξ)<sup>2</sup>, 2ξ(1 − ξ) and
ξ<sup>2</sup>, and the Brownian bridges with their exact means and variances.

Argun and Statt [5] built their homogeneous gels with the same procedure (`strands`, ξ = 0.99,
the second cross-linker drawn from the Gaussian distribution of distances, Brownian-bridge
strands), setting the density by the box edge; their primary-loop fractions (strands of 100
beads, f = 4) are reproduced with an eighth of their system at the same densities (250 strands,
125 cross-linkers, six seeds each):

| Box (theirs) | Box (here) | Primary-loop fraction, ref. 5 | Here |
|--:|--:|:-:|:-:|
| 80 | 40 | ≈ 0.11 | 0.113 ± 0.015 |
| 100 | 50 | | 0.195 ± 0.023 |
| 120 | 60 | ≈ 0.30 | 0.305 ± 0.026 |

## Code structure

```
include/gusev/
  sim/     packing: cross-linker placement and Brownian-bridge beads; linking: Gusev's procedure;
           the driver
  io/      configuration
lib/       implementation
apps/      gusev_build
```

gusev and [randomwalk](../randomwalk/README.md) are both packing modules and share their
structure: the end-linking system they build, with its statistics and consistency check, its
LAMMPS data, XYZ and link writers (`strands.dat`), the keys of the chains, cross-linkers, seed and
output, and the steps of a run are in [`src/common`](../common/README.md) (`poly/endlink/`), so
each module holds only its method; gusev keeps its positions as continuous coordinates. Unit
tests are in [`tests/gusev`](../../tests/gusev). The chain-length distributions are shared with
randomwalk, and the random number generator and the configuration reader with all modules.

## References

1. A. A. Gusev. Numerical Estimates of the Topological Effects in the Elasticity of Gaussian
   Polymer Networks and Their Exact Theoretical Description. *Macromolecules* **2019**, 52,
   3244–3251. [10.1021/acs.macromol.9b00262](https://doi.org/10.1021/acs.macromol.9b00262)
2. A. A. Gusev, F. Schwarz. Molecular Dynamics Validation and Applications of the Maximum Entropy
   Homogenization Procedure for Predicting the Elastic Properties of Gaussian Polymer Networks.
   *Macromolecules* **2019**, 52, 9445–9455.
   [10.1021/acs.macromol.9b01766](https://doi.org/10.1021/acs.macromol.9b01766)
3. S. Wang, D. Ramkrishna, V. Narsimhan. Exact sampling of polymer conformations using Brownian
   bridges. *J. Chem. Phys.* **2020**, 153, 034901. [10.1063/5.0010368](https://doi.org/10.1063/5.0010368)
4. T. Bernhard, F. Schwarz, A. A. Gusev. pylimer-tools: A Python Package for Generating and
   Analyzing Bead-Spring Polymer Networks. arXiv:2508.11509 (**2025**).
   [10.48550/arXiv.2508.11509](https://doi.org/10.48550/arXiv.2508.11509)
5. B. R. Argun, A. Statt. Interplay of Spatial and Topological Defects in Polymer Networks. *ACS
   Eng. Au* **2024**, 4, 351–358. [10.1021/acsengineeringau.3c00072](https://doi.org/10.1021/acsengineeringau.3c00072)
6. D. Sen, B. D. Olsen. Quantification of cyclic topology in polymer networks using 3D nets.
   *Phys. Rev. Mater.* **2025**, 9, 015602.
   [10.1103/PhysRevMaterials.9.015602](https://doi.org/10.1103/PhysRevMaterials.9.015602)
