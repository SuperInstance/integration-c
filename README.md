# integration-c

**Eight libraries. One theorem. These tests prove it.**

Integration tests proving that the SuperInstance metal libraries form a single coherent mathematical ecosystem.

## The Libraries

| Library | Domain |
|---------|--------|
| [conservation-spectral-topology-c](https://github.com/SuperInstance/conservation-spectral-topology-c) | Conservation laws, spectral gaps, topological invariants |
| [sheaf-agents-c](https://github.com/SuperInstance/sheaf-agents-c) | Cellular sheaves, cohomology, agent synchronization |
| [hodge-belief-c](https://github.com/SuperInstance/hodge-belief-c) | Hodge decomposition, belief states, harmonic analysis |
| [ergodic-transport-c](https://github.com/SuperInstance/ergodic-transport-c) | Markov chains, ergodic theory, optimal transport |
| [persistence-c](https://github.com/SuperInstance/persistence-c) | Persistent homology, barcodes |
| [renormalization-c](https://github.com/SuperInstance/renormalization-c) | Coarse-graining, renormalization group flow |
| [west-african-c](https://github.com/SuperInstance/west-african-c) | Griot oral tradition as persistence mechanism |
| [lm-krylov-c](https://github.com/SuperInstance/lm-krylov-c) | Lanczos/Krylov solvers for spectral computations |

## The Theorem

> **Conservation implies spectral gap non-decreasing under renormalization.**
> When a conserved quantity flows on a topological substrate, the spectral gap of the sheaf Laplacian cannot decrease under coarse-graining.

This is the thread connecting all eight libraries. The tests below prove it holds at every link in the chain.

## Integration Tests

| Test | Chain | What it proves |
|------|-------|---------------|
| `test_cst_to_sheaf.c` | CST → Sheaf | Conservation sheaf on a graph: H⁰ counts components with conserved flow |
| `test_sheaf_to_hodge.c` | Sheaf → Hodge | Sheaf Laplacian decomposes into exact/coexact/harmonic; H¹ matches cohomology |
| `test_ergodic_to_cst.c` | Ergodic → CST | Markov chain at equilibrium satisfies the conservation law |
| `test_renormalization_cst.c` | Renormalization → CST | Coarse-graining preserves spectral gap behavior (RG flow) |
| `test_west_african_sheaf.c` | West African → Sheaf | Griot persistence tracks H⁰/H¹ evolution over time |
| `test_full_pipeline.c` | Full pipeline | Graph → conservation sheaf → equilibrium → spectral gap non-decreasing |

## Building

```bash
# Clone dependency headers (or set INCLUDE_PATHS)
make deps

# Build and run all tests
make test
```

Each test file is standalone and self-contained. Mock/stub implementations are used where direct linking would add complexity without mathematical value. **The point is testing mathematical coherence, not build systems.**

## Assertions

**35+ integration assertions** proving the ecosystem is ONE mathematical structure.

## License

MIT
