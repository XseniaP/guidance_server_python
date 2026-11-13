### GUIDANCE3

> **Supported platforms: macOS (arm64 / Intel) and Linux (x86-64) only.**
> Windows is not supported — no Windows wheel is published on PyPI and the bundled compiled binaries are not available for Windows.

---

#### Installation — CLI

```bash
pip install guidance3
```

All Python dependencies and bundled binaries (iqtree3, msa_set_score, clustalo, etc.) are installed automatically. The Command Line Interface (CLI) detects your platform and uses the correct executables.

**MAFFT and PRANK are not on PyPI** and must be installed separately, e.g.:

```bash
conda install -c bioconda mafft prank      # or: brew install mafft / apt install mafft
```

> **Installation time:** a few minutes; TensorFlow (~500 MB) is the largest download. The install is a one-time cost.

The development version can be installed directly from GitHub:
`pip install git+https://github.com/XseniaP/guidance_server_python.git`

---

#### Quick start

```bash
guidance3 --seqFile <path_to_fasta> --msaProgram MAFFT --seqType aa \
          --outDir <path_to_output_dir> --program GUIDANCE3 \
          --bootstraps 100 --proc_num 8
```

Use `--seqType nuc` for nucleotides or `--seqType codon` for codon-aware analysis. Run `guidance3 -h` for the full option list.

---

#### Container (Docker / Podman)

```bash
docker build --platform linux/amd64 -t guidance3 .
docker run --rm -v "$PWD":/data guidance3 --seqFile /data/seqs.fasta --msaProgram MAFFT --seqType aa \
           --outDir /data/guidance3_out --proc_num 4
```

The image includes MAFFT and PRANK; nothing else needs to be installed.

---

#### Installation — Web server / full conda environment

Required if you want to run the web server, or need GPU/CUDA support for the DL model.

```bash
git clone https://github.com/XseniaP/guidance_server_python.git
cd guidance_server_python
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda env create -f environment.yml
conda activate guidance3
export PYTHONPATH=$(pwd)
python app/__init__.py
```

The SECRET_KEY, RECAPTCHA_SITE_KEY, and RECAPTCHA_SECRET_KEY should be saved in a **.env** file in the **./app** folder:

```
SECRET_KEY = ''
RECAPTCHA_SITE_KEY = ''
RECAPTCHA_SECRET_KEY = ''
```

These keys are loaded by **dotenv** from `app/__init__.py`.

After you finish working with the conda environment, deactivate it:

```bash
conda deactivate
```

---

#### Development / Contributing

To work on the code itself:

```bash
git clone https://github.com/XseniaP/guidance_server_python.git
cd guidance_server_python
pip install -e .
```

The `-e` (editable) flag means changes to the source are reflected immediately without reinstalling. You can also run the pipeline directly without installing:

```bash
python3 guidance3/pipeline/main.py --seqFile <fasta> --msaProgram MAFFT --seqType aa \
        --outDir <outDir> --program GUIDANCE3 --bootstraps 100 --proc_num 8
```

---

#### Platform notes

Pre-compiled binaries for the following programs are bundled under `guidance3/programs/mac/` and `guidance3/programs/linux/`:

- msa_set_score, removeTaxa, isEqualTree, features_for_msas, iqtree3 (IQ-TREE 3.1.4), clustalo (Clustal Omega)

The CLI detects your platform at runtime and uses the correct executables automatically. If you are on an unsupported platform, you will need to obtain the binaries from the respective upstream projects and update the paths in **SharedConsts.py** and **guidance3/constants.py**.

---

#### Citations

If you use GUIDANCE3 in your research, please cite:

- **GUIDANCE2** — Sela I, Ashkenazy H, Katoh K, Pupko T (2015) GUIDANCE2: accurate detection of unreliable alignment regions accounting for the uncertainty of multiple parameters. *Nucleic Acids Res* 43:W7–W14. https://doi.org/10.1093/nar/gkv318
- **GUIDANCE** — Penn O, Privman E, Ashkenazy H, Landan G, Graur D, Pupko T (2010) GUIDANCE: a web server for assessing alignment confidence scores. *Nucleic Acids Res* 38:W23–W28. https://doi.org/10.1093/nar/gkq443
- **HoT** — Landan G, Graur D (2007) Heads or tails: a simple reliability measure for multiple sequence alignments. *Mol Biol Evol* 24:1380–1391. https://doi.org/10.1093/molbev/msm060 
- **HoT** — Landan G, Graur D (2008) Local reliability measures from sets of co-optimal multiple sequence alignments. *Pac Symp Biocomput* 13:15-24. PMID: 18229673. https://pubmed.ncbi.nlm.nih.gov/18229673/ 

---

#### Third-party tools

GUIDANCE3 relies on the following external tools. Please cite them alongside GUIDANCE3 when using this software in your research.

**Bundled binaries** (installed automatically with `pip install guidance3`):

- **IQ-TREE 3** — phylogenetic tree inference (bundled v3.1.4; set `GUIDANCE_IQTREE` to use another binary, otherwise `iqtree3` on the PATH is used if the bundled one cannot run)  
  Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Mol Biol Evol* 43. https://doi.org/10.1093/molbev/msag117

- **Clustal Omega** — multiple sequence alignment  
  Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011) Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Mol Syst Biol* 7:539. https://doi.org/10.1038/msb.2011.75

**External tools** (must be installed separately via conda or your system package manager):

- **MAFFT** — multiple sequence alignment  
  Katoh K, Standley DM (2013) MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Mol Biol Evol* 30:772–780. https://doi.org/10.1093/molbev/mst010

- **PRANK** — phylogeny-aware multiple sequence alignment  
  Löytynoja A (2014) Phylogeny-aware alignment with PRANK. *Methods Mol Biol* 1079:155–170. https://doi.org/10.1007/978-1-62703-646-7_10

**Key Python libraries:**

- **Biopython** — Cock PJ et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422–1423. https://doi.org/10.1093/bioinformatics/btp163

---

#### Command-line reference

**`guidance3`** — main pipeline

#### GUIDANCE3 Options

**Required:**

  --seqFile USRSEQ_FILE
                        Specify the sequence file.

  --msaProgram {MAFFT,PRANK}
                        Specify the MSA program. Default=""

  --seqType {aa,nuc,codon}
                        Specify the sequence type: aa, nuc, or codon.

  --outDir OUTDIR       Specify the full path to the output directory.

**Algorithm:**

  --program {HoT,GUIDANCE3}
                        Specify the algorithm to run. Default is GUIDANCE3.

  --inputType {seq,re_align,msa}
                        Specify the type of input provided: seq (unaligned sequences), re_align (re-align an existing MSA), or msa (evaluate a
                        user-provided MSA without re-aligning). Default is seq.

  --bootstraps BOOTSTRAPS
                        Number of bootstrap iterations. Default is 100.

  --disableConvergence
                        Disable the convergence stopping criterion. When set, GUIDANCE3 always runs the exact number of bootstraps specified by
                        --bootstraps rather than stopping early when scores stabilise. Useful for reproducibility or debugging.

  --proc_num PROC_NUM   Number of processors to use. Default=2.

**Phylogeny across alternative MSAs:**

  --noConsensusTree     Do not build the GUIDANCE3 tree. By default GUIDANCE3 infers a maximum-likelihood tree
                        (IQ-TREE 3, -fast; LG+G for proteins, GTR+G for nucleotides) from each of --nAltTrees random
                        alternative MSAs, builds their extended majority-rule consensus and resolves it by a
                        constrained search on the base MSA. Each branch gets an alignment support value: the % of
                        the alternative-MSA trees that contain it.

  --nAltTrees N_ALT_TREES
                        Number of alternative MSAs used for the trees. Default is 20.

  --userTree USER_TREE  A Newick tree with the input sequence names; GUIDANCE3 reports the alignment support of
                        each of its branches.

**Set of MSAs as input:**

  --msaSet MSA_SET      A folder, or .zip / .tar.gz archive, of aligned FASTA files of the same sequences (same
                        names; 2 to 2000 MSAs and at least 4 sequences). The deep-learning model ranks them by predicted
                        distance from the true MSA, and GUIDANCE3 evaluates the best one without re-alignment
                        (as --inputType msa). Use instead of --seqFile. The ranking is written to MSA_set_ranking.csv.

**Filtering thresholds:**

  --seqCutoff SP_SEQ_CUTOFF
                        Confidence cutoff for sequence removal (0–1). Sequences scored below this value are removed. Default is 0.6.

  --colCutoff SP_COL_CUTOFF
                        Confidence cutoff for column removal (0–1). Columns scored below this value are removed. Default is 0.93.

  --Z_Seq_Cutoff Z_SEQ_CUTOFF
                        Z-score threshold as an additional criterion to filter sequences. EXPERIMENTAL. Default is NA (not active).

  --Z_Col_Cutoff Z_COL_CUTOFF
                        Z-score threshold as an additional criterion to filter columns. EXPERIMENTAL. Default is NA (not active).

**Output:**

  --outOrder {aligned,as_input}
                        Order of sequences in output: aligned (default) or as_input.

  --dataset DATASET     Prefix used for all output file names. Default=MSA.

  --msaFile USERMSA_FILE
                        Provide an existing MSA file instead of computing one (use with --inputType msa). Default=None.

**Codon options:**

  --genCode {1,15,6,10,2,5,3,13,9,14,4}
                        Genetic code table. Default is 1 (Nuclear Standard).
                        1  = Nuclear Standard
                        15 = Nuclear Blepharisma
                        6  = Nuclear Ciliate
                        10 = Nuclear Euplotid
                        2  = Mitochondria Vertebrate
                        5  = Mitochondria Invertebrate
                        3  = Mitochondria Yeast
                        13 = Mitochondria Ascidian
                        9  = Mitochondria Echinoderm
                        14 = Mitochondria Flatworm
                        4  = Mitochondria Protozoan

**MSA program parameters:**

  --MSA_Param ALIGN_PARAM
                        Parameters to pass to the MSA program. To include a flag with a leading `-`, prefix it with `\`, e.g. `\-F` for PRANK.

  --outOrder {aligned,as_input}
                        Output sequence order. Default is aligned.

**Advanced / alignment program paths:**

  --mafft MAFFT_PROG    Path to mafft executable. Default=mafft.

  --prank PRANK_PROG    Path to prank executable. Default=prank.

  --clustalo CLUSTALW_PROG
                        Path to clustalo executable. Default=clustalo.

  --muscle MUSCLE_PROG  Path to muscle executable. Default=muscle.

  --pagan PAGAN_PROG    Path to pagan executable. Default=pagan.

  --ruby RUBY_PROG      Path to ruby executable (required by PAGAN). Default=ruby.

  --RootingType {BioPerl,MidPoint}
                        Guide-tree rooting method. Default=BioPerl.

  --BBL {YES,NO}        Perform branch-length optimisation on guide trees. Default=NO.

  --GapPenDist {UNIF,EMP}
                        Distribution from which gap penalties are sampled: uniform (UNIF) or empirical (EMP). Default=UNIF.
                        Relevant only for GUIDANCE3.

---

**`guidance3-rank`** — rank a set of MSAs with the bundled DL model and return the best one

```bash
guidance3-rank --msas-dir <folder | set.zip | set.tar.gz> --out-dir <output_dir>
```

The DL model predicts each MSA's distance from the (unknown) true MSA; the MSA with the lowest predicted score ranks first. The MSAs must be alignments of the same sequences (checked); accepted extensions: `.fasta .fas .fa .faa .fna .aln .afa .mfa`. Outputs in `<output_dir>`:

- `ranking.csv` — `rank`, `code` (MSA file name) and `predicted_score`, best first
- `best_msa.fasta` — the top-ranked MSA, unchanged
- `features.csv` — the alignment features the model used

The sequence type (amino acids or nucleotides) is detected automatically; use `--seq-type codon` for codon alignments, which are ranked on their amino-acid translation. Add `--verbose` for detailed progress. The MSAs can come from any source, e.g. different aligners or settings, or the alternative MSAs of a `guidance3` run.

---

**`guidance3-tree`** — the GUIDANCE3 tree and alignment-aware branch support from any set of alternative MSAs

The tree step of the `guidance3` pipeline as a stand-alone command, for MSAs made elsewhere (other aligners, other settings, or the alternative MSAs of an earlier run):

```bash
# GUIDANCE3 tree (+ supports of your own tree) from a folder or .zip / .tar.gz of MSAs of the same sequences
guidance3-tree --msas-dir <alternative_msas> --base-msa <default_msa.fasta> \
               --user-tree <my_tree.nwk> --out-dir <output_dir> --proc-num 4

# re-use trees already inferred (e.g. <prefix>_AltMSA_trees.nwk of a guidance3 run): takes seconds
guidance3-tree --alt-trees <prefix>_AltMSA_trees.nwk --base-msa <default_msa.fasta> \
               --user-tree <my_tree.nwk> --out-dir <output_dir>
```

- An ML tree (IQ-TREE 3, `-fast`; LG+G for proteins, GTR+G for nucleotides) is inferred from each of `--n-trees` (default 20) randomly chosen MSAs of the set.
- **GUIDANCE3 tree**: the extended majority-rule consensus of these trees, resolved by an ML search on the base MSA constrained to it.
- **Branch supports** (one table row per branch): alignment support (% of the trees containing the branch), aBayes and SH-aLRT (1000 replicates) on the base MSA with the topology fixed, and combined = aBayes × alignment support / 100.
- `--base-msa` defaults to the MSA of the set that the DL model ranks best (ranking written to `<prefix>_DL_ranking.csv`). With `--alt-trees` it is required for the GUIDANCE3 tree; with `--no-consensus --user-tree` and no base MSA, only alignment support is reported (aBayes / SH-aLRT columns are `NA`).
- Outputs (`--prefix`, default `MSA`): `<prefix>_GUIDANCE3_tree.nwk` / `_branch_support.tsv` / `.html` (interactive viewer), `<prefix>_UserTree_alignment_support.nwk` / `.tsv` / `.html`, and `<prefix>_AltMSA_trees.nwk`.

---

**`guidance3-concat-msa`** — concatenate N randomly sampled MSAs from a folder into a super-MSA

```bash
guidance3-concat-msa --msas-dir <path/to/alternative_msas> --n 20 --out-dir <output_dir>
```

Writes the result to `<output_dir>/concatenated_msa.fasta`.

Picks `--n` MSAs at random from the folder and concatenates them column-wise. All MSAs must contain the same set of sequences.

---

#### Licence

Copyright (C) 2026 Tal Pupko. This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version. This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
The full text is in [LICENSE](LICENSE).
