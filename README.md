### GUIDANCE3

> **Supported platforms: macOS (arm64 / Intel) and Linux (x86-64) only.**
> Windows is not supported — the bundled compiled binaries are not available for Windows.

#### Installation — CLI (recommended)

The quickest way to get the `guidance3` command available system-wide:

```bash
git clone https://github.com/XseniaP/guidance_server_python.git
cd guidance_server_python
pip install -e .
```

All Python dependencies are installed automatically. Executables for MacOS and Linux are bundled inside the package — the CLI detects your platform and uses the correct ones.

*Simple example* — amino-acid sequences, GUIDANCE3 algorithm, 100 bootstraps, 8 CPUs:

```bash
guidance3 --seqFile <path_to_fasta> --msaProgram MAFFT --seqType aa \
          --outDir <path_to_output_dir> --program GUIDANCE3 \
          --bootstraps 100 --proc_num 8
```

Use `-h` at any time to print the full option list.

The `--seqType` should be changed to `nuc` for nucleotides or `codon` for codon-aware analysis.

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

#### Platform notes

Executables for the following programs are bundled for both MacOS-arm64 (M1) and Ubuntu Linux:

- msa_set_score, semphy, removeTaxa, isEqualTree, iqtree, clustalo (ClustalOmega)

Each program's source and makefile are in `./script/programs/<program>/`. If you are on a different platform, build replacements using the makefiles and update the paths in **SharedConsts.py**.

---

#### Command-line usage

All commands below become available after `pip install -e .`. Without pip install you can run the main pipeline directly from the project root:

```bash
python3 guidance3/pipeline/main.py --seqFile <fasta> --msaProgram MAFFT --seqType aa \
        --outDir <outDir> --program GUIDANCE3 --bootstraps 100 --proc_num 8
```

---

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

**`guidance3-predict`** — select the best MSA from a folder of alternatives using the bundled DL model

```bash
guidance3-predict \
  --msas-dir  <path/to/folder_with_alternative_msas> \
  --seq-type  aa \
  --out-dir   <output_dir>
```

Use `--seq-type nuc` for nucleotide sequences. The command handles everything internally: feature extraction, model loading (bundled pretrained weights), prediction, and output. The best MSA is written to `<out-dir>/best_msa.fasta`. Add `--verbose` for detailed progress.

The alternative MSAs folder should contain `.fasta` files — typically the `BP/GUIDANCE3_MSA/` sub-directory produced by a `guidance3` run, or any set of MSAs you want to rank.

---

**`guidance3-concat-msa`** — concatenate N randomly sampled MSAs from a folder into a super-MSA

```bash
guidance3-concat-msa --msas-dir <path/to/alternative_msas> --n 20 --out-dir <output_dir>
```

Writes the result to `<output_dir>/concatenated_msa.fasta`.

Picks `--n` MSAs at random from the folder and concatenates them column-wise. All MSAs must contain the same set of sequences.
