### Guidance2.0.3_Beta_Version

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

OPI

* The **./script/programs** folder has multiple subfolders with the .cpp programs' source code which require platform-specific builds to be performed. Build for MacOS and Linux are located in this folder for the following programs:
  - semphy  
  - removeTaxa  
  - msa_set_score  
  - isEqualTree  
  - iqtree
  - clustalo (ClustalOmega)

Each program makefile is located in this program subfolder accordingly. The existing executables in the folders are built for MacOS-arm64 (M1) and Ubuntu Linux accordingly, if this is not the platform you are working on, they should be deleted and replaced with the executables which you build on your platform using the makefiles, paths should be updated in the SharedConst.py file.


#### To Run Guidance2.0.3 Beta Version:

*Simple example* of running the program from the command line:

`cd <base_directory_of_the_project>`

`**python3** **script/guidance_main.py** --seqFile <path_to_the_fasta_file_with_sequences> --msaProgram MAFFT --seqType aa --outDir <path_to_the_output_directory> --program GUIDANCE2 --bootstraps 100 --proc_num 8`

After you finished working with the conda virtual environment on Linux/Ubuntu please **deactivate** the environment by running the following command line:

`conda deactivate`

In this sample run, it is assumed that the input is amino-acids (aa) sequences, 100 bootstrap trees are created, and 8 CPUs are used   
The '--seqType' should be changed to 'nuc' in case of nucleotides.

-h option can be used at any time to see the full list of program options (listed below)

#### Guidance2.1 Beta Version Options:

 -h, --help            show this help message and exit
 
  --seqFile USRSEQ_FILE
                        Specify the sequence file (required).
                        
  --msaProgram {MAFFT,PRANK,CLUSTALO,MUSCLE,PAGAN}
                        Specify the MSA program (Required). <MAFFT|PRANK|CLUSTALO|MUSCLE|PAGAN>. Default=""
                        
  --seqType {aa,nuc,codon}
                        Specify the sequence type: aa, nuc, or codon (Required
                        
  --outDir OUTDIR       Specify the full path to the output directory (required).
  
  --program {GUIDANCE,HoT,GUIDANCE2}
                        Specify the program to run (optional): GUIDANCE, HoT or GUIDANCE2. Default is GUIDANCE2.
                        
  --inputType {seq,re_align,msa}
                        Specify the type of input provided (optional): seq, re_align or msa. Default is seq.
                        
  --bootstraps BOOTSTRAPS
                        Specify the number of bootstrap iterations. Default is 100.
                        
  --genCode {1,15,6,10,2,5,3,13,9,14,4}
                        Specify the codon table. Default is 1 (Nuclear Standard). <option value=1> Nuclear Standard, <option value=15> Nuclear
                        Blepharisma, <option value=6> Nuclear Ciliate, <option value=10> Nuclear Euplotid, <option value=2> Mitochondria
                        Vertebrate, <option value=5> Mitochondria Invertebrate, <option value=3> Mitochondria Yeast, <option value=13> Mitochondria
                        Ascidian <option value=9> Mitochondria Echinoderm <option value=14> Mitochondria Flatworm <option value=4> Mitochondria
                        Protozoan
                        
  --outOrder {aligned,as_input}
                        Specify the output order (optional). Default is aligned.
                        
  --msaFile USERMSA_FILE
                        Specify the MSA file (optional). Not recommended, see documentation online guidance.tau.ac.il. Default=None
                        
  --seqCutoff SP_SEQ_CUTOFF
                        Specify confidence cutoff between 0 to 1. Default is 0.6.
                        
  --colCutoff SP_COL_CUTOFF
                        Specify confidence cutoff between 0 to 1. Default is 0.93.
                        
  --Z_Seq_Cutoff Z_SEQ_CUTOFF
                        Specify Z score as additional criteria to filter sequences. EXPERIMENTAL. Default is NA (not active).
                        
  --Z_Col_Cutoff Z_COL_CUTOFF
                        Specify Z score as additional criteria to filter position. EXPERIMENTAL. Default is NA (not active).
                        
  --mafft MAFFT_PROG    Specify path to mafft executable. Default=mafft.
  
  --prank PRANK_PROG    Specify path to prank executable. Default=prank.
  
  --clustalo CLUSTALW_PROG
                        Specify path to clustalo executable. Default=clustalo.
                        
  --muscle MUSCLE_PROG  Specify path to muscle executable. Default=muscle.
  
  --pagan PAGAN_PROG    Specify path to pagan executable. Default=pagan.
  
  --ruby RUBY_PROG      Specify path to ruby executable. Default=ruby.
  
  --dataset DATASET     Specify a unique name for the Dataset - will be used as prefix to outputs. Default=MSA.
  
  --MSA_Param ALIGN_PARAM
                        Specify the parameters for the alignment program. To pass parameter containing - in it, add \ before each - e.g. \-F for
                        PRANK
                        
  --proc_num PROC_NUM   Specify num of processors to use. Default=1.
  
  --RootingType {BioPerl,MidPoint}
                        Specify Rooting Type: BioPerl or MidPoint. Default=BioPerl
                        
  --BBL {YES,NO}        Specify if to do branch length optimization (BBL): YES or NO. Default=NO
  
  --GapPenDist {UNIF,EMP}
                        Specify if to sample gap penalties from uniform (UNIF) or empirical (EMP) distribution. Default = UNIF => RELEVANT ONLY FOR
                        GUIDANCE 2
