# Example 2 — randomwalk

| File | Purpose |
|:--|:--|
| `endlink_paper.ini` | precursors of the end-linked network with N<sub>o</sub> = 20 (Zhang et al., *Macromolecules* 2026, Table 1) |
| `endlink_polydisperse.ini` | the same composition with Schulz–Zimm chain lengths, L<sub>n</sub> = 18 and Đ = 1.5 |
| `smoke_small.ini` | small system that packs in milliseconds |
| `in.prepare` | LAMMPS: relaxation, equilibration and compression to number density 0.85 |
| `in.react` | LAMMPS: end-linking with `fix bond/create` |

From the repository root:

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini
cd endlink_No20
lmp -in ../examples/example_2/in.prepare -var data config.lmp
lmp -in ../examples/example_2/in.react
```

N<sub>o</sub> is the paper's label of a network, close to its mean strand length. The other
systems of Table 1 differ in three keys, given on the command line:

| N<sub>o</sub> | `chain_length` | `n_chains` | `n_crosslinkers` |
|--:|--:|--:|--:|
| 20 | 18 | 5680 | 2840 |
| 10 | 8 | 12416 | 6208 |
| 7 | 6 | 19520 | 9760 |
| 6 | 5 | 27154 | 13577 |
| 5 | 4 | 35020 | 17510 |

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini chain_length=8 n_chains=12416 n_crosslinkers=6208 directory=endlink_No10
```

Polydisperse precursors keep the number of chains and cross-linkers; `chain_length` becomes the
number-average length, and `length_distribution` and `dispersity` set the shape and width of the
distribution (see [`src/randomwalk/README.md`](../../src/randomwalk/README.md), *Chain lengths*):

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_polydisperse.ini
build/src/randomwalk/rw_pack examples/example_2/endlink_polydisperse.ini length_distribution=gaussian dispersity=1.1 directory=endlink_G1.1
```

The run writes `chain_lengths.dat` with the lengths drawn and those expected. The LAMMPS scripts
work for any chain lengths, as they find the reactive chain ends by atom type.

The LAMMPS scripts need LAMMPS 30 Sep 2026 or newer, built with the MOLECULE and MC packages.
They use the step counts of the paper (610,000 steps of preparation and 120 million steps of
end-linking, after which more than 96 % of the cross-linker sites have reacted); shorter
runs are set with `-var relax_steps …`, `-var pre_steps …`, `-var compress_steps …` or
`-var react_steps …`, and the seed with `-var seed …`. All keys of `rw_pack` are described in
[`src/randomwalk/README.md`](../../src/randomwalk/README.md).
