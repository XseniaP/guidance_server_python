# Helper program sources

C++ sources of the programs bundled in `guidance3/programs/`:

| Program | Source |
|---|---|
| `msa_set_score`, `isEqualTree`, `removeTaxa` | `guidance2/` (GUIDANCE2, Pupko lab) |
| `features_for_msas` | `features_for_msas/` + `sp_alternative_cpp/` |

IQ-TREE 3, MAFFT, PRANK and Clustal Omega are third-party programs and are not included here.

`build_tools.sh <out_dir>` compiles all four programs (needs `g++`, `make`, `cmake`). The Dockerfile uses it to build the
native `linux/arm64` image; the x86-64 image uses the bundled binaries.
