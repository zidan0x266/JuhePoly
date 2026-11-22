# Example 3 — gusev

| File | Purpose |
|:--|:--|
| `endlink_melt.ini` | network with the composition of N<sub>o</sub> = 20 (Zhang et al., *Macromolecules* 2026, Table 1) at melt density, linked to conversion 0.96 |
| `endlink_polydisperse.ini` | the same network with Schulz–Zimm strand lengths, L<sub>n</sub> = 18 and Đ = 1.5 |
| `smoke_small.ini` | small network that builds in milliseconds |
| `in.relax` | LAMMPS: push-off, relaxation and equilibration of a built network |

From the repository root:

```bash
build/src/gusev/gusev_build examples/example_3/endlink_melt.ini
cd gusev_No20
lmp -in ../examples/example_3/in.relax -var data config.lmp
```

`gusev_build` writes the network already linked, so the end-linking run of example_2 is not
needed: `in.relax` removes the overlaps of the Brownian-bridge strands with a soft repulsion that
is ramped up under capped displacements, switches to the bead-spring model of example_2 (shifted
Lennard-Jones beads, cutoff 2.5, harmonic bonds K = 1000, r0 = 1), relaxes it and equilibrates it
at T = 1. Shorter runs are set with `-var push_steps …`, `-var relax_steps …` or
`-var equil_steps …`, and the seed with `-var seed …`. The input needs LAMMPS 30 Sep 2026 or newer
with the MOLECULE package. On the small network of `smoke_small.ini`, with
`-var push_steps 5000 -var relax_steps 2000 -var equil_steps 5000` (12,000 steps), it brings the
longest bond from 2.5 to 1.08 and the pair energy to −4.20 per atom.

The other systems of the fracture paper differ in `chain_length`, `n_chains` and
`n_crosslinkers` as listed in [example_2](../example_2/README.md), and
`procedure=strands conversion=0.99` selects Gusev's original procedure. Polydisperse strands
(`endlink_polydisperse.ini`) use the same keys as in example_2; every strand is linked with its
own number of segments, and the run also writes `chain_lengths.dat`.
All keys of `gusev_build` are described in [`src/gusev/README.md`](../../src/gusev/README.md).
