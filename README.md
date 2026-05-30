### GUIDANCE3

#### Local run on MacOS-arm64 or Ubuntu Linux

**Prerequisites:**

* Download or clone the project directory from GitHub
* Install all prerequisites by creating a conda environment using the **environment.yml** file provided in the project directory

**Using environment.yml**

`git clone https://github.com/XseniaP/guidance_server_python.git`

`cd guidance_server_python`

`conda config --add channels conda-forge`

`conda config --add channels bioconda`

`conda config --add channels defaults`

`conda env create -f environment.yml`

`conda activate Guidance_py_env`

`export PYTHONPATH=$(pwd)`

`python app/__init__.py`   OR    `python3 app/__init__.py`


Details behind the .yml environment: the environment installs all dependencies, including MAFFT v7.525, prank v.170427, and pip dependencies.
Other required programs' executables (for both Linux and MacOS) are located in the **./script/programs** folder. No changes are required if either Linux or MacOS is used to run the web server.

The SECRET_KEY, RECAPTCHA_SITE_KEY, and RECAPTCHA_SECRET_KEY should be saved in **.env** file in the **./app** folder in the following format:

SECRET_KEY = '' <br>
RECAPTCHA_SITE_KEY = '' <br>
RECAPTCHA_SECRET_KEY = ''

These keys are loaded by **dotenv** from the __init__.py file.  
________________________________________________

* The **./script/programs** folder has multiple subfolders with the .cpp programs' source code which require platform-specific builds to be performed. Builds for MacOS and Linux are located in this folder for the following programs:
  - semphy  
  - removeTaxa  
  - msa_set_score  
  - isEqualTree  
  - iqtree
  - clustalo (ClustalOmega)

Each program's makefile is located in its subfolder. The existing executables are built for MacOS-arm64 (M1) and Ubuntu Linux. If you are on a different platform, delete them, build replacements using the makefiles, and update the paths in **SharedConsts.py**.


#### Command-line usage

The canonical entry point is `guidance3/pipeline/main.py`. The legacy wrapper `script/guidance_main.py` delegates to it and is used by the web server internally.

*Simple example* — amino-acid sequences, GUIDANCE3 algorithm, 100 bootstraps, 8 CPUs:

`cd <base_directory_of_the_project>`

`python3 guidance3/pipeline/main.py --seqFile <path_to_fasta> --msaProgram MAFFT --seqType aa --outDir <path_to_output_dir> --program GUIDANCE3 --bootstraps 100 --proc_num 8`

After you finish working with the conda environment on Linux/Ubuntu, deactivate it:

`conda deactivate`

The `--seqType` should be changed to `nuc` for nucleotides or `codon` for codon-aware analysis.

Use `-h` at any time to print the full option list.

#### GUIDANCE3 Options

**Required:**

  --seqFile USRSEQ_FILE
                        Specify the sequence file.

  --msaProgram {MAFFT,MAFFT_LINSI,PRANK,CLUSTALO,MUSCLE,PAGAN}
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
