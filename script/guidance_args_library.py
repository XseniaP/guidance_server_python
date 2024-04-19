import json
import pickle
import argparse
from shutil import copy

import guidance_CONSTANTS
from guidance_sequence_functions import *
from guidance_tree_functions import *
Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)

class Library:
    def __init__(self):
        self.isServer = 0
        self.input_type = "seq"
        # self.proc_num = None
        self.outDir = ""
        self.stored_data_file = ""
        self.stored_form_data = ""
        # self.semphy_prog = os.path.join(os.getcwd(), 'programs/semphy/semphy')
        # self.semphy_prog = '/Users/kpolonsky/PycharmProjects/guidance_server_python/script/programs/semphy/semphy'
        self.semphy_prog = os.path.join(BIN_DIR, 'script/programs/semphy/semphy')
        # self.iqtree_prog = os.path.join(os.getcwd(), 'programs/iqtree-2.2.2.6-MacOSX/bin/iqtree2')
        # self.iqtree_prog = os.path.join('/Users/kpolonsky/PycharmProjects/guidance_server_python/script/programs/iqtree-2.2.2.6-MacOSX/bin/iqtree2')
        self.iqtree_prog = os.path.join(BIN_DIR, 'script/programs/iqtree-2.2.2.6-MacOSX/bin/iqtree2')
        self.mafft_prog = "mafft"
        self.prank_prog = "prank"
        # self.clustalw_prog = "clustalw"
        self.clustalw_prog = os.path.join(BIN_DIR, 'script/programs/clustalo')
        self.ruby_prog = "ruby"
        self.muscle_prog = "muscle"
        # self.msa_set_score_prog = os.path.join(os.getcwd(), 'programs/msa_set_score/msa_set_score')
        self.msa_set_score_prog = os.path.join(BIN_DIR, 'script/programs/msa_set_score/msa_set_score')
        self.pagan_prog = "pagan"
        # self.remove_taxa_prog = os.path.join(os.getcwd(), 'programs/removeTaxa/removeTaxa')
        self.remove_taxa_prog = os.path.join(BIN_DIR, 'script/programs/removeTaxa/removeTaxa')
        self.orig_argv = sys.argv[1:]  # Skip the script name itself
        self.overview_URL = ""
        self.gallery_URL = ""
        self.home_URL = ""
        self.status_file = ""  # will follow the status of alternative MSA creation on the web-server
        self.DNA_AA = {}
        self.align_param = ""
        self.SP_SEQ_CUTOFF = "0.6"
        self.SP_COL_CUTOFF = "0.93"
        self.Z_Col_Cutoff = "NA"
        self.Z_Seq_Cutoff = "NA"
        self.dataset = "MSA"
        self.rooting_type = "BioPerl"
        self.BBL = "NO"
        self.GapPenDist = "UNIF"
        self.proc_num = 1
        self.userMSA_File = ""
        self.PROGRAM = "GUIDANCE2"
        self.Bootstraps = 100
        self.CodonTable = 1  # Nuclear standard code
        self.Seq_Type = ""
        self.Align_Order = "aligned"
        self.usrSeq_File = ""
        self.MSA_Program = ""
        self.server_output = ""
        # self.usage_ = "FILL IN USAGE LATER"
        self.usage_ = "USAGE: --seqFile <seqFile> --msaProgram <MAFFT|PRANK|CLUSTALO|MUSCLE|PAGAN> --seqType <aa|nuc|codon> --outDir <full path outDir> \
    Optional parameters: \
  --program <GUIDANCE|HoT|GUIDANCE2> default=GUIDANCE2 \
  --bootstraps <number of bootstrap iterations> default=100 \
  --genCode <option value> default=1 \
                <option value=1>  Nuclear Standard \
                <option value=15> Nuclear Blepharisma \
           	    <option value=6>  Nuclear Ciliate \
                <option value=10> Nuclear Euplotid \
                <option value=2>  Mitochondria Vertebrate \
                <option value=5>  Mitochondria Invertebrate \
                <option value=3>  Mitochondria Yeast \
                <option value=13> Mitochondria Ascidian \
                <option value=9>  Mitochondria Echinoderm \
                <option value=14> Mitochondria Flatworm \
                <option value=4>  Mitochondria Protozoan \
  --outOrder <aligned|as_input> default=aligned \
  --msaFile <msaFile> - not recommended, see documentation online guidance.tau.ac.il \
  --seqCutoff <confidence cutoff between 0 to 1> default=0.6 \
  --colCutoff <confidence cutoff between 0 to 1> default=0.93 \
  --Z_Seq_Cutoff <Z score as additional criteria to filter sequences> EXPERIMENTAL, default=NA (not active) \
  --Z_Col_Cutoff <Z score as additional criteria to filter position> EXPERIMENTAL, default=NA (not active) \
  --mafft <path to mafft executable> default=mafft \
  --prank <path to prank executable> default=prank \
  --clustalw <path to clustalw executable> default=clustalw \
  --muscle <path to muscle executable> default=muscle \
  --pagan <path to pagan executable> default=pagan \
  --ruby <path to ruby executable> default=ruby \
  --dataset Unique name for the Dataset - will be used as prefix to outputs (default=MSA) \
  --MSA_Param passing parameters for the alignment program. To pass parameter containning '-' in it, add \\ before each '-' e.g. \\-F for PRANK \
  --proc_num <num of processors to use> default=1"


    def check_and_set_input_and_output_variables(self, arguments):
        """get variables from input arguments and fill out the Variable Class properties"""
        if not len(sys.argv) >= 2:
            print("At least one --seqFile is required.\n" + self.usage_)
            sys.exit()

        if len(sys.argv) > 1 and sys.argv[1].startswith('-'):
            self.parse_arguments()
            self.check_arguments_for_errors()
            self.check_align_parameters()
            self.check_output_files()
            self.validate_sequences()
        else:
            self.server_input_mode_arguments()
        self.set_rest_of_the_variables()
        self.user_provided_MSA()
        self.convert_fasta_names_to_codes()
        self.if_codons()


    def parse_arguments(self):
        """parse arguments and fill out the relevant Variable Class properties"""
        parser = argparse.ArgumentParser(prog="GUIDANCE.v2.02", description="Program description: NA")

        # Required arguments
        parser.add_argument('--seqFile', dest='usrSeq_File', type=str, required=True,
                            help='Specify the sequence file (required).')
        parser.add_argument('--msaProgram', dest='MSA_Program', type=str,
                            choices=['MAFFT', 'PRANK', 'CLUSTALO', 'MUSCLE', 'PAGAN'], required=True, default="",
                            help='Specify the MSA program (Required). <MAFFT|PRANK|CLUSTALO|MUSCLE|PAGAN>. Default=""')
        parser.add_argument('--seqType', dest='Seq_Type', type=str, choices=['aa', 'nuc', 'codon'],
                            required=True, default="", help='Specify the sequence type: aa, nuc, or codon (Required')
        parser.add_argument('--outDir', required=True, dest='outDir', type=str,
                            help='Specify the full path to the output directory (required).')

        # Optional parameters
        parser.add_argument('--program', dest='PROGRAM', type=str, choices=['GUIDANCE', 'HoT', 'GUIDANCE2'],
                            default='GUIDANCE2',
                            help='Specify the program to run (optional): GUIDANCE, HoT or GUIDANCE2. Default is GUIDANCE2.')
        parser.add_argument('--inputType', dest='input_type', type=str, choices=['seq', 're_align', 'msa'],
                            default='seq',
                            help='Specify the type of input provided (optional): seq, re_align or msa. Default is seq.')
        parser.add_argument('--bootstraps', type=int, dest='Bootstraps', default=100,
                            help='Specify the number of bootstrap iterations. Default is 100.')
        parser.add_argument('--genCode', type=int, dest='CodonTable', default=1,
                            choices=[1, 15, 6, 10, 2, 5, 3, 13, 9, 14, 4],
                            help='Specify the codon table. Default is 1 (Nuclear Standard).\n'
                                 '<option value=1>  Nuclear Standard, \n'
                                 '<option value=15> Nuclear Blepharisma, \n'
                                 '<option value=6>  Nuclear Ciliate, \n'
                                 '<option value=10> Nuclear Euplotid, \n'
                                 '<option value=2>  Mitochondria Vertebrate, \n'
                                 '<option value=5>  Mitochondria Invertebrate, \n'
                                 '<option value=3>  Mitochondria Yeast, \n'
                                 '<option value=13> Mitochondria Ascidian \n'
                                 '<option value=9>  Mitochondria Echinoderm \n'
                                 '<option value=14> Mitochondria Flatworm \n'
                                 '<option value=4>  Mitochondria Protozoan \n')
        parser.add_argument('--outOrder', dest='Align_Order', type=str, choices=['aligned', 'as_input'],
                            default='aligned', help='Specify the output order (optional). Default is aligned.')
        parser.add_argument('--msaFile', dest='userMSA_File', type=str, default="",
                            help='Specify the MSA file (optional). Not recommended, see documentation online guidance.tau.ac.il. Default=None')
        parser.add_argument('--seqCutoff', dest='SP_SEQ_CUTOFF', type=str, default="0.6",
                            help='Specify confidence cutoff between 0 to 1. Default is 0.6.')
        parser.add_argument('--colCutoff', dest='SP_COL_CUTOFF', type=str, default="0.93",
                            help='Specify confidence cutoff between 0 to 1. Default is 0.93.')
        parser.add_argument('--Z_Seq_Cutoff', default='NA', dest='Z_Seq_Cutoff', type=str,
                            help='Specify Z score as additional criteria to filter sequences. EXPERIMENTAL. Default is NA (not active).')
        parser.add_argument('--Z_Col_Cutoff', default='NA', dest='Z_Col_Cutoff', type=str,
                            help='Specify Z score as additional criteria to filter position. EXPERIMENTAL. Default is NA (not active).')
        parser.add_argument('--mafft', default='mafft', dest='mafft_prog', type=str,
                            help='Specify path to mafft executable. Default=mafft.')
        parser.add_argument('--prank', default='prank', dest='prank_prog', type=str,
                            help='Specify path to prank executable. Default=prank.')
        parser.add_argument('--clustalo', default='clustalo', dest='clustalw_prog', type=str,
                            help='Specify path to clustalo executable. Default=clustalo.')
        parser.add_argument('--muscle', default='muscle', dest='muscle_prog', type=str,
                            help='Specify path to muscle executable. Default=muscle.')
        parser.add_argument('--pagan', default='pagan', dest='pagan_prog', type=str,
                            help='Specify path to pagan executable. Default=pagan.')
        parser.add_argument('--ruby', default='ruby', dest='ruby_prog', type=str,
                            help='Specify path to ruby executable. Default=ruby.')
        parser.add_argument('--dataset', default='MSA', dest='dataset', type=str,
                            help='Specify a unique name for the Dataset - will be used as prefix to outputs. Default=MSA.')
        parser.add_argument('--MSA_Param', dest='align_param', type=str, default="",
                            help='Specify the parameters for the alignment program. To pass parameter containing - in it, add \\ before each - e.g. \\-F for PRANK')
        parser.add_argument('--proc_num', dest='proc_num', type=int, default=1,
                            help='Specify num of processors to use. Default=1.')

        ### EXPERIMENTAL FEATURES.... MOST ACTIVE ONLY LOCAL
        parser.add_argument('--RootingType', dest='rooting_type', choices=['BioPerl', 'MidPoint'],
                            default='BioPerl', type=str,
                            help='Specify Rooting Type: BioPerl or MidPoint. Default=BioPerl')
        parser.add_argument('--BBL', dest='BBL', choices=['YES', 'NO'], default='NO', type=str,
                            help='Specify if to do branch length optimization (BBL): YES or NO. Default=NO')
        parser.add_argument('--GapPenDist', dest='GapPenDist', choices=['UNIF', 'EMP'], default='UNIF',
                            type=str,
                            help='Specify if to sample gap penalties from uniform (UNIF) or empirical (EMP) distribution. Default = UNIF => RELEVANT ONLY FOR GUIDANCE 2')

        # Parse the command-line arguments
        args = parser.parse_args()
        for arg_name, arg_value in vars(args).items():
            setattr(self, arg_name, arg_value)
        return args

    # def unalign(self):
    #     with open(self.userMSA_File, "r") as file_in, open(self.usrSeq_File, 'w') as file_out:
    #         lines = file_in.readlines()
    #         new_line = ""
    #         for line in lines:
    #             if line.startswith(">"):
    #                 if new_line != "":
    #                     file_out.write(new_line.strip() + "\n")
    #                 file_out.write(line.strip() + "\n")
    #                 new_line = ""
    #             else:
    #                 new_line += line.replace("-", "").strip()
    #         file_out.write(new_line)

    def check_arguments_for_errors(self):

        if self.BBL.upper() == "YES":
            self.semphy_prog = "/groups/pupko/haim/Programs/semphy_test_clean_log_BBL/programs/semphy/semphy"

        if not isinstance(self.Bootstraps, int) or not str(self.Bootstraps).isdigit():
            raise ValueError("ERROR: Bootstraps parameter must be a number")

        if self.outDir == None or self.outDir == "":
            raise ValueError("ERROR: No path for output\n")

        if self.usrSeq_File == None or self.usrSeq_File == "":
            raise ValueError("ERROR: seqFile is required\n")

        if self.proc_num < 1:
            raise ValueError("ERROR: Number of processors must be >= 1\n")

        if self.MSA_Program == "PAGAN":
            print("WARNING: PAGAN requires that MAFFT is also installed on your system. Otherwise GUIDANCE will fail\n")

        self.MSA_Program = self.MSA_Program.upper()

        if self.MSA_Program not in ["MAFFT", "PRANK", "CLUSTALO", "MUSCLE", "PAGAN"]:
            raise ValueError("ERROR: msaProgram should be MAFFT or PRANK or CLUSTALO or MUSCLE or PAGAN (case sensitive)\n")


        # if self.input_type == "msa" or self.input_type == "re_align":
            # self.userMSA_File = self.usrSeq_File
            # self.userMSA_File = self.Alignment_File
            # os.remove(self.usrSeq_File)
            # self.usrSeq_File = f"{self.usrSeq_File}_seq"
            # self.unalign()
            # if self.input_type == "re_align":
            #     self.userMSA_File = ""

        if not self.outDir.endswith("/"):
            self.outDir += "/"

        print(f"outDir: {self.outDir}\n")
        self.WorkingDir = self.outDir

        if not os.path.exists(self.outDir):
            os.system(f"mkdir {self.outDir}")

        if self.GapPenDist == "" and self.PROGRAM == "GUIDANCE3":
            self.GapPenDist = "UNIF"
        elif self.GapPenDist == "" and self.PROGRAM == "GUIDANCE3_HOT":
            self.GapPenDist = "UNIF"
        elif self.GapPenDist == "" and self.PROGRAM == "GUIDANCE2":
            self.GapPenDist = "UNIF"

        if self.PROGRAM == "GUIDANCE2":
            self.overview_URL = "http://guidance.tau.ac.il/ver2/overview.php"
            self.gallery_URL = "http://guidance.tau.ac.il/ver2/Gallery.php"
            self.home_URL = "http://guidance.tau.ac.il/ver2/"
            self.credits_URL = "http://guidance.tau.ac.il/ver2/credits.php"
        else:
            self.overview_URL = "http://guidance.tau.ac.il/overview.html"
            self.gallery_URL = "http://guidance.tau.ac.il/Gallery.htm"
            self.home_URL = "http://guidance.tau.ac.il/"
            self.credits_URL = "http://guidance.tau.ac.il/credits.html"

        # self.OutLogFile = f"{self.WorkingDir}/log"
        self.OutLogFile = f"{self.outDir}/log"
        self.Output = self.OutLogFile
        try:
            with open(self.OutLogFile, "a") as log_file:
                log_file.write(
                    "\n\n========================================= NEW GUIDANCE RUN STARTED ===========================================\n")
                log_file.write(f"GUIDANCE COMMAND: python {sys.argv[0]} {' '.join(sys.argv[1:])}\n")
        except IOError as e:
            # raise FileNotFoundError(f"ERROR: Can't open Log File: {VARS['OutLogFile']}")
            exit_on_error('sys_error', f"Can't open Log File: {self.OutLogFile}: {e}", self)

        # INPUT FILES
        self.SeqsFile = "Seqs.Orig.fas"  # generic fixed name for the sequence file. the user file will be copied to this file.
        self.SeqsFile_Codons = "Seqs.Orig_DNA.fas"  # generic fixed name for the DNA CODONS sequence file. the user file will be copied to this file.
        if self.usrSeq_File != "" and self.usrSeq_File != None:
            shutil.copy(self.usrSeq_File, os.path.join(f"{self.WorkingDir}", f"{self.SeqsFile}"))
            if self.isServer == 1:
                convertNewline(self.SeqsFile, self.WorkingDir)
            removeEndLineExtraChars(self.SeqsFile, self.WorkingDir)
        if self.Seq_Type == 'aa':
            self.Seq_Type = 'AminoAcids'
        elif self.Seq_Type == 'nuc':
            self.Seq_Type = 'Nucleotides'
        elif self.Seq_Type == 'codon':
            self.Seq_Type = 'Codons'
            shutil.move(os.path.join(self.WorkingDir, self.SeqsFile), os.path.join(self.WorkingDir, self.SeqsFile_Codons))
        else:
            raise ValueError("ERROR: --seqType must be: aa or nuc or codon\n")

    def check_align_parameters(self):

        if 'addfragments' in self.align_param and self.MSA_Program == 'MAFFT':     # support --addfragments option of MAFFT, seqFile is the seq for the coreMSA
            user_fragments_file = re.search(r'addfragments (\S+)', self.align_param).group(1)
            self.fragments_file_name = "add_fragments_file"
            self.fragments_file_name_seqName_coded = f"{self.fragments_file_name}.SeqNameCoded"
            try:
                copy(user_fragments_file, os.path.join(self.WorkingDir, self.fragments_file_name))
            except Exception as e:
                # raise RuntimeError(f"Can't copy user fragments file: {user_fragments_file} to {self.VARS['WorkingDir']}{self.VARS['fragments_file_name']}. Error: {e}\n")
                exit_on_error("sys_error",f"Can't copy user fragments file: {user_fragments_file} to {self.WorkingDir}{self.fragments_file_name}\n", self)
            with open(f'{self.OutLogFile}', "a") as log_file:
                validation_message = f"Validating fragments: Guidance::validate_Seqs({self.WorkingDir}, {self.fragments_file_name}, {self.Seq_Type}, No): \n"
                print(validation_message, end="")
                log_file.write(validation_message)

            ans = validate_seqs(self.WorkingDir, self.fragments_file_name, self.Seq_Type, "No", "")

            if "sys_error" in ans:
                exit_on_error('sys_error', " ".join(ans), self)
            elif ans[0] != "OK":
                exit_on_error('user_error', "".join(ans), self)

            if ans[0] == "OK" and ans[1] != "":
                with open(f'{self.OutLogFile}', "a") as log_file:
                    log_file.write(f"Warning: {ans[1]}; Nevertheless, calculation is continued\n")
                print(f"Warning: {ans[1]}; Nevertheless, calculation is continued\n")

            with open(f'{self.OutLogFile}', "a") as log_file:
                log_file.write("return: " + ans + "\n")

            self.fragments_file_name = ans[2]
            self.NumOfFragments = ans[3]

            tmp = re.split(re.escape("--"), self.align_param)
            # tmp = re.split(r'\-\-', self.align_param)
            tmp_size = len(tmp)
            if tmp_size >= 1:  # command line type
                for i in range(len(tmp)):
                    if 'addfragments' in tmp[i]:
                        tmp[i] = f"addfragments {os.path.join(self.WorkingDir, self.fragments_file_name_seqName_coded)}"
                self.align_param = "--".join(tmp)
                # self.align_param = "\-\-".join(tmp)
            else:  # server type
                tmp = re.split(re.escape("--"), self.align_param)
                # tmp = re.split(r'--', self.align_param)
                for i in range(len(tmp)):
                    if 'addfragments' in tmp[i]:
                        tmp[i] = f"addfragments {os.path.join(self.WorkingDir, self.fragments_file_name_seqName_coded)}"
                self.align_param = "--".join(tmp)
                # self.align_param = "\-\-".join(tmp)

            # Check if need to remove reorder with fragments
            if '--reorder' in self.align_param or '\-\-reorder' in self.align_param:
                self.align_param = self.align_param.replace('--reorder', '')        # if seed is provided reorder must be removed so the seeds will be first
                print(
                    "WARNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored, and the output order will be the same as input (with seeds first)\n")
                with open(f'{self.OutLogFile}', "a") as log_file:
                    log_file.write("WARNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored, and the output order will be the same as input (with seeds first)\n")

        if "--reorder" in self.align_param and "seed" in self.align_param:
            self.align_param = re.sub(r"--reorder", "", self.align_param)
            print(
                "WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be the same as input (with seeds first)\n")
            with open(f'{self.OutLogFile}', "a") as log_file:
                log_file.write(
                    "WARNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored, and the output order will be the same as input (with seeds first)\n")

        retree_match = re.search(r"--retree ([0-9]+)", self.align_param)
        if retree_match:
            retree_value = int(retree_match.group(1))
            if retree_value > 1:
                print(f"WARNNING: --retree {retree_value} is not supported in GUIDANCE, therefore this argument is ignored.\n")
                with open(f'{self.OutLogFile}', "a") as log_file:
                    log_file.write(
                        f"WARNNING: --retree {retree_value} is not supported in GUIDANCE, therefore this argument is ignored.\n")
                self.align_param = re.sub(r"--retree ([0-9]+)", "", self.align_param)

    def check_output_files(self):

        if self.userMSA_File != "" and self.userMSA_File != None:       # The user gave the base MSA as input
            self.Alignment_File = "UserMSA"
            if not os.path.exists(os.path.join(self.WorkingDir, self.Alignment_File)):
                shutil.copy(self.userMSA_File, os.path.join(self.WorkingDir, self.Alignment_File))
                if self.isServer == 1:
                    convertNewline(self.Alignment_File, self.WorkingDir)
                removeEndLineExtraChars(self.Alignment_File, self.WorkingDir)
        else:
            self.Alignment_File = f"{self.dataset}.{self.MSA_Program}.aln"

        if "addfragments" in self.align_param:
            self.Core_Alignment_File = f"{self.dataset}.{self.MSA_Program}.CORE.aln"

        self.output_page = "log"

    def validate_sequences(self):

        if self.userMSA_File == "":  # Only seq file provided, no MSA
            with open(self.OutLogFile, 'a') as log_file:
                ans = []
                if self.Seq_Type != "Codons":
                    log_file.write(f"Guidance::validate_Seqs({self.WorkingDir}, {self.SeqsFile}, {self.Seq_Type}, No): \n")
                    ans = validate_seqs(self.WorkingDir, self.SeqsFile, self.Seq_Type, "No", "")
                else:
                    log_file.write(f"Guidance::validate_Seqs({self.WorkingDir}{self.SeqsFile_Codons}, {self.Seq_Type}, No): \n")
                    ans = validate_seqs(self.WorkingDir, self.SeqsFile_Codons, self.Seq_Type, "No","")

                if ans[0] == "sys_error":
                    exit_on_error('sys_error', ans[1], self)
                elif ans[0] != "OK":
                    exit_on_error('user_error', ''.join(ans), self)

                if ans[0] == "OK" and ans[1] != "":
                    log_file.write(f"Warning: {ans[1]}; Nevertheless calculation is continued\n")
                    print(f"Warning: {ans[1]}; Nevertheless calculation is continued\n")

                ans_joined = " ".join(ans)
                log_file.write(f"return: {ans_joined}\n")

            self.SeqsFile = ans[2]
            self.NumOfSeq = int(ans[3])


        elif os.path.exists(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}")) and os.path.getsize(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}")) > 0:   # Alignment provided
            try:
                with open(f"{self.OutLogFile}", 'a') as log_file:
                    ans = []
                    if self.Seq_Type != "Codons":
                        log_file.write(
                            f"Guidance::validate_Seqs({self.WorkingDir}, {self.Alignment_File}, {self.Seq_Type}, Yes): \n")
                        ans = validate_seqs(self.WorkingDir, self.Alignment_File, self.Seq_Type,
                                                     "Yes", "")
                    else:
                        log_file.write(
                            f"Guidance::validate_Seqs({self.WorkingDir},{self.Alignment_File},{self.Seq_Type},Yes,{self.CodonTable}) \n")
                        ans = validate_seqs(self.WorkingDir, self.Alignment_File, self.Seq_Type,
                                                     "Yes", self.CodonTable)

                    if ans[0] == "sys_error":
                        exit_on_error('sys_error', ans[1], self)
                    elif ans[0] != "OK":
                        exit_on_error('user_error', ''.join(ans), self)
                    if ans[0] == "OK" and ans[1] != "":
                        print(f"Warning: {ans[1]}; Nevertheless calculation is continued")
                        log_file.write(f"Warning: {ans[1]}; Nevertheless calculation is continued\n")
                        # self.Alignment_File = ans[2]
                        # self.NumOfSeq = int(ans[3])

                    log_file.write(f"return: {' '.join(ans)}\n")

            except IOError as e:
                exit_on_error('sys_error', f"Can't open Log File: {self.OutLogFile}: {e}", self)

            self.Alignment_File = ans[2]
            self.NumOfSeq = int(ans[3])

        if self.NumOfSeq < 4 and self.PROGRAM == "GUIDANCE":
            exit_on_error('user_error',
                          f"Only {self.NumOfSeq} sequences were provided, however at least 4 sequences are required for GUIDANCE<br>You can run HoT algorithm instead.\n", self)

        if self.NumOfSeq < 4 and self.PROGRAM == "GUIDANCE2":
            exit_on_error('user_error',
                          f"Only {self.NumOfSeq} sequences were provided, however at least 4 sequences are required for GUIDANCE2<br>You can run HoT algorithm instead.\n", self)


    def server_input_mode_arguments(self):
        # Server input mode
        self.isServer = 1

        self.stored_data_file = sys.argv[1]
        self.stored_form_data = sys.argv[2]

        # Data from files
        with open(self.stored_data_file, "r") as file:
            json_string = file.read()
            vars_data = json.loads(json_string)

        with open(self.stored_form_data, "r") as file:
            json_string = file.read()
            form_data = json.loads(json_string)

        for key, value in vars_data.items():
            setattr(self, key, value)

        for key, value in form_data.items():
            setattr(self, key, value)

        self.dataset = "MSA"

        # GUIDANCE2 server defaults
        self.rooting_type = "BioPerl"
        self.BBL = "NO"
        self.GapPenDist = "UNIF"

        if self.BBL.upper() == "NO":
            self.semphy_prog = guidance_CONSTANTS.SEMPHY
        else:
            self.semphy_prog = guidance_CONSTANTS.SEMPHY_BBL  # TO DO: Change its location to a more stable one

        self.mafft_prog = guidance_CONSTANTS.MAFFT_GUIDANCE
        self.prank_prog = guidance_CONSTANTS.PRANK_LECS
        # self.clustalw_prog = guidance_CONSTANTS.CLUSTALW_LECS
        self.clustalw_prog = guidance_CONSTANTS.CLUSTAL_OMEGA
        self.muscle_prog = guidance_CONSTANTS.MUSCLE
        self.pagan_prog = guidance_CONSTANTS.PAGAN_LECS
        self.ruby_prog = guidance_CONSTANTS.RUBY
        self.msa_set_score_prog = MSA_SET_SCORE

        # Defaults (still not supported by the web server implementation, experimental feature)
        self.Z_Col_Cutoff ='NA'
        self.Z_Seq_Cutoff ='NA'

        if self.CALLING_SERVER == "GUIDANCE2":
            self.overview_URL = "http://guidance.tau.ac.il/ver2/overview.php"
            self.gallery_URL = "http://guidance.tau.ac.il/ver2/Gallery.php"
            self.home_URL = "http://guidance.tau.ac.il/"
            self.credits_URL = "http://guidance.tau.ac.il/credits.php"
        else:
            self.overview_URL = "http://guidance.tau.ac.il/overview.html"
            self.gallery_URL = "http://guidance.tau.ac.il/Gallery.htm"
            self.home_URL = "http://guidance.tau.ac.il/"
            self.credits_URL = "http://guidance.tau.ac.il/credits.html"

        output_file_path = os.path.join(self.WorkingDir, self.output_page)
        self.server_output = output_file_path

    def set_rest_of_the_variables(self):
        # codons handling
        self.TranslateErrors = "xCodons.html"
        # Output Files
        self.code_fileName = "Seqs.Codes"
        self.codded_seq_fileName = "Seqs.numberd.fas"
        self.codded_seq_fileName_Codons = "Seqs_DNA.numberd.fas"
        self.Alignment_File_PROT = f"{self.dataset}.{self.MSA_Program}.PROT.aln"
        self.BootStrap_Dir = f"{self.WorkingDir}BP/"
        self.BootStrap_MSA_Dir = f"{self.BootStrap_Dir}BP_MSA/"
        self.HoT_MSAs_Dir = "COS_MSA"
        self.GUIDANCE2_MSAs_Dir = f"{self.BootStrap_Dir}GUIDANCE2_MSA/"
        self.Tree_File = ""
        self.Semphy_OutFile = ""
        self.Semphy_LogFile = ""
        self.Semphy_StdFile = ""
        self.COL_SCORES_FIGURE = "Col_Scores_Graph.png"
        self.Scoring_Alignments_Dir = ""  # The dir with the alignment used to create the score
        self.send_email_dir = guidance_CONSTANTS.SEND_EMAIL_DIR_IBIS
        self.DNA_AA = {}

        # if self.isServer == 1:
        #     output_file_path = os.path.join(self.WorkingDir, self.output_page)
        #     self.Output = output_file_path
        #     # with open(output_file_path, "a") as OUTPUT:
            #     pass


    def user_provided_MSA(self):
        log_file_path = f"{self.OutLogFile}"
        try:
            with open(log_file_path, 'a') as log_file:
                if os.path.exists(os.path.join(self.WorkingDir, self.Alignment_File)) and os.path.getsize(
                        os.path.join(self.WorkingDir, self.Alignment_File)) > 0:
                    self.code_fileName_aln = self.code_fileName + ".fromALN"
                    if self.Seq_Type != "Codons":
                        log_file.write(f"==== Alignment_File:{self.Alignment_File}\tSeq_Type: {self.Seq_Type}\n")
                        log_file.write(
                            f"extract_seq_from_MSA({self.WorkingDir}{self.Alignment_File},{self.WorkingDir}{self.SeqsFile})\n")
                        extract_seq_from_MSA(f"{self.WorkingDir}{self.Alignment_File}",
                                             f"{self.WorkingDir}{self.SeqsFile}", self)
                    else:
                        log_file.write(f"==== Alignment_File:{self.Alignment_File}\tSeq_Type: {self.Seq_Type}\n")
                        log_file.write(
                            f"extract_seq_from_MSA({self.WorkingDir}{self.Alignment_File},{self.WorkingDir}{self.SeqsFile_Codons})\n")
                        extract_seq_from_MSA(f"{self.WorkingDir}{self.Alignment_File}",
                                             f"{self.WorkingDir}{self.SeqsFile_Codons}", self)

                    if self.PROGRAM == "GUIDANCE2" or self.PROGRAM == "GUIDANCE3_HOT":
                        log_file.write(
                            f"Guidance::name2codeFastaFrom1({self.WorkingDir}{self.Alignment_File}, {self.WorkingDir}{self.code_fileName_aln}, {self.WorkingDir}{self.Alignment_File}.WithCodesName,0,seqNum);\n")
                        ans = name2code_fasta_from1(f"{self.WorkingDir}{self.Alignment_File}",
                                                            f"{self.WorkingDir}{self.code_fileName_aln}",
                                                            f"{self.WorkingDir}{self.Alignment_File}.WithCodesName", 0,
                                                            "seqNum")
                        if ans[0] != "ok":
                            exit_on_error("sys_error", f"Guidance::name2codeFastaFrom1: {' '.join(ans)}\n", self)
                    else:
                        log_file.write(
                            f"Guidance::name2codeFastaFrom1({self.WorkingDir}{self.Alignment_File}, {self.WorkingDir}{self.code_fileName_aln}, {self.WorkingDir}{self.Alignment_File}.WithCodesName);\n")
                        ans = name2code_fasta_from1(f"{self.WorkingDir}{self.Alignment_File}",
                                                            f"{self.WorkingDir}{self.code_fileName_aln}",
                                                            f"{self.WorkingDir}{self.Alignment_File}.WithCodesName")
                        if ans[0] != "ok":
                            exit_on_error("sys_error", f"Guidance::name2codeFastaFrom1: {' '.join(ans)}\n", self)

                    if self.Seq_Type == "Nucleotides":
                        if self.MSA_Program == "MAFFT":
                            convert_fs_to_lower_case(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}.WithCodesName"))
                        elif self.MSA_Program == "PRANK":
                            convert_fs_to_upper_case(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}.WithCodesName"))
                        elif self.MSA_Program == "CLUSTALO":
                            convert_fs_to_upper_case(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}.WithCodesName"))
                        elif self.MSA_Program == "MUSCLE":
                            convert_fs_to_upper_case(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}.WithCodesName"))
                        elif self.MSA_Program == "PAGAN":
                            convert_fs_to_upper_case(os.path.join(f"{self.WorkingDir}",f"{self.Alignment_File}.WithCodesName"))

                    if os.path.getsize(
                            f"{self.WorkingDir}{self.Alignment_File}.WithCodesName") == 0 or not os.path.exists(
                            f"{self.WorkingDir}{self.Alignment_File}.WithCodesName"):
                        exit_on_error("user_error",
                                      f"Sequences were not found on the <A htef=\"{self.Alignment_File}\">uploaded file</A><br>Make sure it is a Plain text FASTA Format<br>", self)

                    if self.PROGRAM == "HoT":
                        names_according_cos(f"{self.WorkingDir}{self.Alignment_File}.WithCodesName")

                    os.system(f"cp {self.WorkingDir}{self.Alignment_File} {self.WorkingDir}{self.Alignment_File}.ORIG")
                    shutil.move(f"{self.WorkingDir}{self.Alignment_File}.WithCodesName",
                                f"{self.WorkingDir}{self.Alignment_File}")
        except FileNotFoundError as e:
            print(f"Couldn't open Log file: {e}")
            sys.exit()


    def convert_fasta_names_to_codes(self):
        log_file_path = f"{self.OutLogFile}"
        try:
            with open(log_file_path, 'a') as log_file:
                seq_counter = 0
                if '--seed' in self.align_param:
                    seeds = re.findall(r'--seed ({})'.format(self.WorkingDir + 'seedFile[0-9]+'), self.align_param)
                    for seed_file in seeds:
                        log_file.write("Guidance::name2codeFasta_without_codded_out({}, {}, {})\n".format(seed_file,
                                                                                               self.WorkingDir + self.code_fileName,
                                                                                               seq_counter))
                        ans = name2code_fasta_without_codded_out(seed_file, self.WorkingDir + self.code_fileName,
                                                                         seq_counter)
                        if ans[0] != "ok":
                            exit_on_error("sys_error", "Guidance::name2codeFasta_without_codded_out: " + " ".join(ans), self)
                        else:
                            seq_counter = ans[1]

                if self.Seq_Type != "Codons":
                    # ans = []
                    if self.PROGRAM in ["GUIDANCE2", "GUIDANCE3_HOT"]:
                        log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {}, seqNum;\n)".format(
                            self.WorkingDir + self.SeqsFile, self.WorkingDir + self.code_fileName,
                            self.WorkingDir + self.codded_seq_fileName, seq_counter))
                        ans = name2code_fasta_from1(self.WorkingDir + self.SeqsFile,
                                                           self.WorkingDir + self.code_fileName,
                                                           self.WorkingDir + self.codded_seq_fileName, seq_counter,
                                                           "seqNum")
                    else:
                        log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {});\n".format(self.WorkingDir + self.SeqsFile,
                                                                                     self.WorkingDir + self.code_fileName,
                                                                                     self.WorkingDir + self.codded_seq_fileName, seq_counter))
                        ans = name2code_fasta_from1(self.WorkingDir + self.SeqsFile,
                                                           self.WorkingDir + self.code_fileName,
                                                           self.WorkingDir + self.codded_seq_fileName, seq_counter)
                    if ans[0] != "ok":
                        exit_on_error("sys_error", "Guidance::name2codeFastaFrom1: " + " ".join(ans) + "\n", self)

                    if os.path.getsize(self.WorkingDir + self.codded_seq_fileName) == 0 or not os.path.exists(
                            self.WorkingDir + self.codded_seq_fileName):
                        exit_on_error("user_error",
                                      "Sequences were not found on the <A htef=\"{}\">uploaded file</A><br>Make sure it is a Plain text FASTA Format<br>".format(
                                          self.SeqsFile), self)

                    if '--addfragments' in self.align_param and self.MSA_Program == "MAFFT":
                        offset = self.NumOfSeq + 1
                        if self.PROGRAM in ["GUIDANCE2", "GUIDANCE3_HOT"]:
                            log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {}, seqNum)\n".format(
                                self.WorkingDir + self.fragments_file_name, self.WorkingDir + self.code_fileName,
                                self.WorkingDir + self.fragments_file_name_seqName_coded, offset))
                            ans = name2code_fasta_from1(self.WorkingDir + self.fragments_file_name,
                                                               self.WorkingDir + self.code_fileName,
                                                               self.WorkingDir + self.fragments_file_name_seqName_coded,
                                                               offset, "seqNum")
                        else:
                            log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {})\n".format(
                                self.WorkingDir + self.fragments_file_name, self.WorkingDir + self.code_fileName,
                                self.WorkingDir + self.fragments_file_name_seqName_coded, offset))
                            ans = name2code_fasta_from1(self.WorkingDir + self.fragments_file_name,
                                                               self.WorkingDir + self.code_fileName,
                                                               self.WorkingDir + self.fragments_file_name_seqName_coded,
                                                               offset)
                        if ans[0] != "ok":
                            exit_on_error("sys_error", "Guidance::name2codeFastaFrom1: " + " ".join(ans) + "\n", self)
                        if os.path.getsize(
                                self.WorkingDir + self.fragments_file_name_seqName_coded) == 0 or not os.path.exists(
                                self.WorkingDir + self.fragments_file_name_seqName_coded):
                            exit_on_error("user_error",
                                          "Sequences were not found on the <A htef=\"{}\">uploaded file</A><br>Make sure it is a Plain text FASTA Format<br>".format(
                                              self.fragments_file_name), self)

                        # create file with the list of fragments name
                        self.fragments_codes = "add_fragments.seqCodes"
                        try:
                            with open(self.WorkingDir + self.fragments_codes, 'w') as fragments_names:
                                for i in range(offset, offset + self.NumOfFragments):
                                    fragments_names.write(str(i) + "\n")
                        except IOError as e:
                            exit_on_error("sys_error", f"Can't open fragments codes file name: {self.WorkingDir}{self.fragments_codes} {e}", self)

                else:
                    if self.PROGRAM in ["GUIDANCE2", "GUIDANCE3", "GUIDANCE3_HOT"]:
                        log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {}, seqNum)\n".format(
                            self.WorkingDir + self.SeqsFile_Codons, self.WorkingDir + self.code_fileName,
                            self.WorkingDir + self.codded_seq_fileName, seq_counter))
                        ans = name2code_fasta_from1(self.WorkingDir + self.SeqsFile_Codons,
                                                           self.WorkingDir + self.code_fileName,
                                                           self.WorkingDir + self.codded_seq_fileName_Codons,
                                                           seq_counter, "seqNum")
                    else:
                        log_file.write("Guidance::name2codeFastaFrom1({}, {}, {}, {})\n".format(
                            self.WorkingDir + self.SeqsFile_Codons, self.WorkingDir + self.code_fileName,
                            self.WorkingDir + self.codded_seq_fileName, seq_counter))
                        ans = name2code_fasta_from1(self.WorkingDir + self.SeqsFile_Codons,
                                                           self.WorkingDir + self.code_fileName,
                                                           self.WorkingDir + self.codded_seq_fileName_Codons,
                                                           seq_counter)
                    if ans[0] != "ok":
                        exit_on_error("sys_error", "Guidance::name2codeFastaFrom1: " + " ".join(ans) + "\n", self)
                    if os.path.getsize(
                            self.WorkingDir + self.codded_seq_fileName_Codons) == 0 or not os.path.exists(
                            self.WorkingDir + self.codded_seq_fileName_Codons):
                        exit_on_error("user_error",
                                      "Sequences were not found on the <A htef=\"{}\">uploaded file</A><br>Make sure it is a Plain text in FASTA Format<br>".format(
                                          self.SeqsFile_Codons), self)

        except FileNotFoundError as e:
            print(f"Couldn't open Log file: {e}")
            sys.exit()

    def if_codons(self):
        # IF CODONS TRANSLATE TO AA
        if self.Seq_Type == "Codons":
            ans = []
            if self.PROGRAM == "GUIDANCE2" or self.PROGRAM == "GUIDANCE3" or self.PROGRAM == "GUIDANCE3_HOT":
                with open(self.OutLogFile, "a") as logFile:
                    logFile.write(
                        f"Translating DNA to AA: codonAlign::translate_DNA_to_AA('{self.WorkingDir}{self.codded_seq_fileName_Codons}', '{self.WorkingDir}{self.codded_seq_fileName}', '{self.CodonTable}', '{self.WorkingDir}{self.TranslateErrors}', '{self.WorkingDir}{self.output_page}', {self.DNA_AA}, '', '', seqNum);\n")
                ans = translate_DNA_to_AA(f"{self.WorkingDir}{self.codded_seq_fileName_Codons}",
                                                      f"{self.WorkingDir}{self.codded_seq_fileName}",
                                                      self.CodonTable,
                                                      f"{self.WorkingDir}{self.TranslateErrors}",
                                                      f"{self.WorkingDir}{self.output_page}", self.DNA_AA, '', '',
                                                      "seqNum")
            else:
                with open(self.OutLogFile, "a") as logFile:
                    logFile.write(
                    f"Translating DNA to AA: codonAlign::translate_DNA_to_AA('{self.WorkingDir}{self.codded_seq_fileName_Codons}', '{self.WorkingDir}{self.codded_seq_fileName}', '{self.CodonTable}', '{self.WorkingDir}{self.TranslateErrors}', '{self.WorkingDir}{self.output_page}', {self.DNA_AA});\n")
                ans = translate_DNA_to_AA(f"{self.WorkingDir}{self.codded_seq_fileName_Codons}",
                                                      f"{self.WorkingDir}{self.codded_seq_fileName}",
                                                      self.CodonTable,
                                                      f"{self.WorkingDir}{self.TranslateErrors}",
                                                      f"{self.WorkingDir}{self.output_page}", self.DNA_AA)

            if ans[0] != "ok":
                if ans[0] == "user":
                    exit_on_error("user_error", " ".join(ans[1]), self)
                elif ans[0] == "sys":
                    exit_on_error("sys_error", " ".join(ans), self)
