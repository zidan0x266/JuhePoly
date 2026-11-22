# Example 1 — 3dbfm

| File | Purpose |
|:--|:--|
| `atrp_paper.ini` | ATRP network, 4 % cross-linker (Zhang et al., *Macromolecules* 2025) |
| `frp_paper.ini` | FRP network, 4 % cross-linker (same paper) |
| `smoke_small.ini` | small 32³ system; the smoke tests (`tests/CMakeLists.txt`) run it to conversion 0.3 with faster rates in under a second, and as given it takes about three minutes |

From the repository root:

```bash
build/src/3dbfm/bfm_react examples/example_1/atrp_paper.ini
```

Other cross-linker contents are set with `monovinyl_fraction`, e.g.
`monovinyl_fraction=0.92 directory=atrp_X08` for 8 %; every run needs its own `directory`, since
the tables are appended to. All keys are described in
[`src/3dbfm/README.md`](../../src/3dbfm/README.md).
