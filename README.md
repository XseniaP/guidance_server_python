### Guidance2.1_Beta_Version

#### Running with Docker:

You need to have a docker engine installed https://www.docker.com/products/docker-desktop/. 
The following commands can be run in Terminal (Bash) to build and run Docker image: 

`git clone https://github.com/XseniaP/Guidance_mid.git`

`cd Guidance_mid`

`docker build -t guidance .`

`docker run -v <path_to_folder_with_seq_input_file>:/input -v <path_to_out_dir>:/output guidance --seqFile "/input/<name_of_seq_input_file>" --msaProgram <MAFFT or PRANK> --seqType aa --outDir "/output/<out_dir_name>/" --program GUIDANCE2 --bootstraps <integer number of bootstrap trees> --proc_num <number of CPUs>`

*Sample run* for the following file structure, s.t. fasta sequence file is located in data folder and the results of this run are expected to be saved to ABD_results folder:

|-- user <br />
| &nbsp;  &nbsp; |-- data <br />
| &nbsp;  &nbsp; | &nbsp;  &nbsp; |-- ABD.fasta <br />
| &nbsp;  &nbsp; |-- ABD_results <br />


`docker run -v /user/data:/input -v /user/ABD_results:/output guidance --seqFile "/input/ABD.fasta" --msaProgram MAFFT --seqType aa --outDir "/output/ABD_results/" --program GUIDANCE2 --bootstraps 100 --proc_num 8`

*Sample run* for the following file structure, s.t. fasta sequence file has an absolute path /user/Downloads/data/AGMAT.fas  and the results of this run will be saved in /user/Downloads/results/AGMAT_results/ folder:

|-- user <br />
| &nbsp;  &nbsp; |-- Downloads <br />
| &nbsp;  &nbsp; | &nbsp;  &nbsp; |-- data <br />
| &nbsp;  &nbsp; | &nbsp;  &nbsp; | &nbsp;  &nbsp; |-- AGMAT.fas <br />
| &nbsp;  &nbsp; | &nbsp;  &nbsp; |-- results <br />
| &nbsp;  &nbsp; | &nbsp;  &nbsp; | &nbsp;  &nbsp; |-- AGMAT_results <br />


`docker run -v /Users/user/Downloads/data:/input -v /Users/user/Downloads/results:/output guidance --seqFile "/input/AGMAT.fas" --msaProgram MAFFT --seqType aa --outDir "/output/AGMAT_results/" --program GUIDANCE2 --bootstraps 100 --proc_num 8`


#### Local run on MacOS-arm64 or Ubuntu Linux

**Prerequisites:**

* Create python project and save all the files and folders either from **guidance_MacOS-arm64**  OR **guidance_Linux** folder  into your project folder accordingly
* Install all prerequisites listed in requirements.txt

**On MacOS (via pip or pip3):**

`git clone https://github.com/XseniaP/Guidance_mid.git`

`cd Guidance_mid`

`pip install -r ./guidance_MacOS-arm64/requirements.txt`

**On Ubuntu Linux:**

`sudo apt install python3.12-venv`

`sudo apt install python3-pip`

`sudo apt install git` 

`git clone https://github.com/XseniaP/Guidance_mid.git`

`cd Guidance_mid`

`python3 -m venv .venv`

`source .venv/bin/activate`

`python3 -m pip install -r ./guidance_Linux/requirements.txt`

* The **./script/programs** folder has multiple subfolders with the .cpp programs' source code which require platform-specific builds to be performed. 
  - semphy - **./script/programs/semphy/semphy** path is assumed
  - removeTaxa - **script/programs/removeTaxa/removeTaxa** path is assumed
  - msa_set_score - **script/programs/msa_set_score/msa_set_score** path is assumed 
  - isEqualTree - **script/programs/isEqualTree/isEqualTree** path is assumed
  - iqtree - **./script/programs/iqtree/bin/iqtree2** path is assumed // executable for your platform (in case it's not Ubuntu Linux or MacOS-arm64/M1) can be downloaded from [http://www.iqtree.org/doc/Quickstart](http://www.iqtree.org/doc/Quickstart) and named `iqtree2` and located in the following path: 
  ./script/programs/iqtree/bin/ 

Each program makefile is located in this program subfolder accordingly. The existing executables in the folders are built for MacOS-arm64 (M1) and Ubuntu Linux accordingly, if this is not the platform you are working on, they should be deleted and replaced with the executables which you build on your platform using the makefiles.

* Other prerequisites to be installed:

  - MAFFT v7.525 should be installed from [https://mafft.cbrc.jp/alignment/software/source.html](https://mafft.cbrc.jp/alignment/software/source.html) and globally callable with `mafft` command line; if local executable is used (not recommended) then path to it should be updated in **./script/config.py** under MAFFT_GUIDANCE variable and **./script/guidance_args_library.py**
     
     - on Linux/Ubuntu can be installed with  `sudo apt install mafft` command line
    
  - prank v.170427 should be installed from [http://wasabiapp.org/software/prank/prank_installation/](http://wasabiapp.org/software/prank/prank_installation/) and globally callable with `prank` command line; if local executable is used (not recommended) then path to it should be updated in **./script/config.py** under PRANK_LECS and PRANK variables and **./script/guidance_args_library.py**
  
    - on Linux/Ubuntu can be installed with  `sudo apt install prank`  command line

  NOT YET SUPPORTED:
  
  - ClustalOmega executable should be located in the script/programs/ under the name `clustalo`. The final path assumed is **script/programs/clustalo**
  - muscle 5.1 should be installed and globally callable with `muscle` command line; if local executable is used then path to it should be updated in **./script/config.py** under MUSCLE variable and **./script/guidance_args_library.py**

#### To Run Guidance2.1 Beta Version:

*Simple example* of running the program from the command line:

`cd <base_directory_of_the_project>`  or  `cd guidance_Linux` or `cd guidance_MacOS-arm64` in case you used git pull and are currently in the Guidance_mid folder

`python3 script/guidance_main.py --seqFile <path_to_the_fasta_file_with_sequences> --msaProgram MAFFT --seqType aa --outDir <path_to_the_output_directory> --program GUIDANCE2 --bootstraps 100 --proc_num 8`

After you finished working with the virtual environment on Linux/Ubuntu please deactivate the environment by running the following command line:

`deactivate`

In this sample run it is assumed that the input is amino-acids (aa) sequences, 100 bootstrap trees are created and 8 CPUs are used   
The '--seqType' should be changed to 'nuc' in case of nucleotides

-h option can be used at any time to see full list of program options (listed below)

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
