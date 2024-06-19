from guidance_sequence_functions import *
from guidance_tree_functions import *
from guidance_common_functions import *
from guidance_scoring_and_visualization import calculate_sp_scores, check_convergence, \
    calculate_sp_scores_convergence, add_scores_to_dict
from multiprocessing import Process, Manager
import uuid
import sys
import multiprocessing as mp
from time_decorator import timeit
from multiprocessing.sharedctypes import Value, Array
from multiprocessing import Process, Manager, Lock

#@timeit
def run_hot_internal(args_library, op_vals_arr_ref, countTrees, tree_good_BranchLength, Branch):
    HOT_COS_GUIDANCE2_cmd = f"cd {args_library.WorkingDir}; python3 {HOT_GUIDANCE2_PROGRAM} {args_library.dataset}_{countTrees} {args_library.HoT_MSA_Program}"
    print(HOT_COS_GUIDANCE2_cmd)

    if args_library.Seq_Type in ["AminoAcids", "Codons"]:
        HOT_COS_GUIDANCE2_cmd += " aa"
    elif args_library.Seq_Type == "Nucleotides":
        HOT_COS_GUIDANCE2_cmd += " nt"

    HOT_COS_GUIDANCE2_cmd += f" {args_library.codded_seq_fileName} . \"\" 0 {args_library.HoT_MSA_Program_path} {tree_good_BranchLength} {Branch}"

    if args_library.MSA_Program == "MAFFT":
        if args_library.PROGRAM == "GUIDANCE2":
            HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} --op {op_vals_arr_ref[countTrees]}"
        elif args_library.PROGRAM == "GUIDANCE3_HOT":
            HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} --op {op_vals_arr_ref[countTrees]} --ep {ep_vals_arr_ref[countTrees]}"
    elif args_library.MSA_Program == "PRANK":
        HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} -gaprate={op_vals_arr_ref[countTrees]}"
    elif args_library.MSA_Program == "CLUSTALO":
        # HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} -GAPOPEN={op_vals_arr_ref[countTrees]}"
        HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param}"

    HOT_COS_GUIDANCE2_cmd += " >> COS.std"

    # log_file.write(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
    # print(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
    os.system(HOT_COS_GUIDANCE2_cmd)
    return HOT_COS_GUIDANCE2_cmd

#@timeit
def run_hot_process_on_tree(args_library, epsilon, proc, RandomBranches,op_vals_arr_ref, ep_vals_arr_ref, Num_of_Aln_from_HoT_per_Run, lock):
    try:
        log_file = open(args_library.OutLogFile, "a")
    except OSError:
        print("run_guidance() could not open log file\n")
        sys.exit()

    bp_per_proc = (args_library.Bootstraps // args_library.proc_num) + 1

    for tree_num in range(bp_per_proc):

        convergence = 0

        # countTrees = proc * bp_per_proc + tree_num
        countTrees = proc + args_library.proc_num * tree_num
        print(f"proc num {proc}\ttree num {tree_num} --> global tree index {countTrees}\n")

        if countTrees >= args_library.Bootstraps:
            break
        # alt_msas = len(os.listdir(args_library.Scoring_Alignments_Dir))
        # if alt_msas >= args_library.convergence:
        #     break

        Branch = RandomBranches[countTrees]

        if args_library.MSA_Program == "MAFFT":
            if 'addfragments' in args_library.align_param:
                tree = f"{args_library.prune_BootStrap_Dir}tree_{countTrees}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{countTrees}CORE.rooted"
            else:
                tree = f"{args_library.BootStrap_Dir}tree_{countTrees}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{countTrees}.rooted"

        elif args_library.MSA_Program == "CLUSTALO":
            tree = f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{countTrees}.rooted"

        else:
            tree = f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{countTrees}"

        if args_library.MSA_Program == "PRANK":
            if 'iterate' in args_library.align_param:
                args_library.align_param = args_library.align_param.replace(r'\-iterate=\d+', '')
                log_file.write("[WARNING] -iterate argument is ignored for the perturbed alignments stage\n")
                print("[WARNING] -iterate argument is ignored when reconstructing the perturbed alignments\n")
            if 'once' not in args_library.align_param:
                args_library.align_param += " -once"

        tree_good_BranchLength = f"{tree}.GoodBranchLength"
        reformat_trees_branch_length(tree, tree_good_BranchLength)

        # HOT_COS_GUIDANCE2_cmd = f"cd {args_library.WorkingDir}; python3 {HOT_GUIDANCE2_PROGRAM} {args_library.dataset}_{countTrees} {args_library.HoT_MSA_Program}"
        # print(HOT_COS_GUIDANCE2_cmd)
        #
        # if args_library.Seq_Type in ["AminoAcids", "Codons"]:
        #     HOT_COS_GUIDANCE2_cmd += " aa"
        # elif args_library.Seq_Type == "Nucleotides":
        #     HOT_COS_GUIDANCE2_cmd += " nt"
        #
        # HOT_COS_GUIDANCE2_cmd += f" {args_library.codded_seq_fileName} . \"\" 0 {args_library.HoT_MSA_Program_path} {tree_good_BranchLength} {Branch}"
        #
        # if args_library.MSA_Program == "MAFFT":
        #     if args_library.PROGRAM == "GUIDANCE2":
        #         HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} --op {op_vals_arr_ref[countTrees]}"
        #     elif args_library.PROGRAM == "GUIDANCE3_HOT":
        #         HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} --op {op_vals_arr_ref[countTrees]} --ep {ep_vals_arr_ref[countTrees]}"
        # elif args_library.MSA_Program == "PRANK":
        #     HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} -gaprate={op_vals_arr_ref[countTrees]}"
        # elif args_library.MSA_Program == "CLUSTALW":
        #     HOT_COS_GUIDANCE2_cmd += f" --- {args_library.align_param} -GAPOPEN={op_vals_arr_ref[countTrees]}"
        #
        # HOT_COS_GUIDANCE2_cmd += " >> COS.std"
        #
        # log_file.write(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
        # print(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
        # os.system(HOT_COS_GUIDANCE2_cmd)

        HOT_COS_GUIDANCE2_cmd = run_hot_internal(args_library, op_vals_arr_ref, countTrees, tree_good_BranchLength, Branch)
        log_file.write(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
        print(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")

        # pertubed_aln = []

        if args_library.NumOfSeq > 2:
            local_dataset = f"{args_library.dataset}_{countTrees}_cos_"
            pathname = os.path.join(
                f"{args_library.WorkingDir}", f"{local_dataset}{args_library.HoT_MSA_Program}/b[01]*.fasta")
            pertubed_aln = glob.glob(pathname)

        else:
            local_dataset = f"{args_library.dataset}_{countTrees}_cos_"
            pathname = os.path.join(
                f"{args_library.WorkingDir}", f"{local_dataset}{args_library.HoT_MSA_Program}/b[01]*.fasta")
            pertubed_aln = glob.glob(pathname)

        shuffled = random.sample(pertubed_aln, len(pertubed_aln))

        for j in range(Num_of_Aln_from_HoT_per_Run):
            aln = shuffled[j]
            base = os.path.basename(aln).split(".fasta")[0]

            if args_library.PROGRAM == "GUIDANCE2":
                cp_from = f"{args_library.WorkingDir}{args_library.dataset}_{countTrees}_cos_{args_library.HoT_MSA_Program}/{base}.fasta"
                cp_to = f"{args_library.Scoring_Alignments_Dir}{base}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_Split_{Branch}.fasta"
                # cmd = f"cp {args_library.WorkingDir}{args_library.dataset}_{countTrees}_cos_{args_library.HoT_MSA_Program}/{base}.fasta {args_library.Scoring_Alignments_Dir}{base}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_Split_{Branch}.fasta"
                cmd = f"cp {cp_from} {cp_to}"
            elif args_library.PROGRAM == "GUIDANCE3_HOT":
                cmd = f"cp {args_library.WorkingDir}{args_library.dataset}_{countTrees}_cos_{args_library.HoT_MSA_Program}/{base}.fasta {args_library.Scoring_Alignments_Dir}{base}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_EP_{ep_vals_arr_ref[countTrees]}_Split_{Branch}.fasta"

            os.system(cmd)

        cmd = ""

        if args_library.PROGRAM == "GUIDANCE2":
            cmd = f"mv {args_library.WorkingDir}{args_library.dataset}_{countTrees}_cos_{args_library.HoT_MSA_Program} {args_library.BootStrap_Dir}{args_library.dataset}_cos_{args_library.HoT_MSA_Program}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_Split_{Branch}"
        elif args_library.PROGRAM == "GUIDANCE3_HOT":
            cmd = f"mv {args_library.WorkingDir}{args_library.dataset}_{countTrees}_cos_{args_library.HoT_MSA_Program} {args_library.BootStrap_Dir}{args_library.dataset}_cos_{args_library.HoT_MSA_Program}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}EP_{ep_vals_arr_ref[countTrees]}_Split_{Branch}"

        os.system(cmd)


        # check convergence starting from the 20th tree (starting from 80 MSAs)
        # if countTrees != 0:
        if args_library.input_type != "msa":
            if countTrees >= 20:
                # check the convergence only for every nth tree
                if args_library.proc_num >=2 or (args_library.proc_num == 1 and countTrees % 3 == 0):
                # if countTrees % 1 == 0:
                    try:
                        alt_msas = calculate_sp_scores_convergence(args_library, countTrees)
                        add_scores_to_dict(args_library, epsilon, countTrees, lock)
                        # print(args_library.mean_res_pair_score)
                        # print(args_library.mean_col_score)
                        os.system(
                            f'rm {os.path.join(args_library.WorkingDir, args_library.Output_Prefix + f"_tree_{countTrees}_*.scr")}')
                        convergence = check_convergence(args_library, epsilon)
                        print(
                                f"convergence of proc num {proc}\ttree num {tree_num} --> global tree index {countTrees} is {convergence} \n")
                    except Exception as e:
                        log_file.write(f"failed to calculate scores for convergence of proc num {proc}\ttree num {tree_num} \t# of alternative MSAs {alt_msas} error {e}\n")
                        print(f"failed to calculate scores for convergence of proc num {proc}\ttree num {tree_num} \t# of alternative MSAs {alt_msas} \n")
            if convergence == 1:
                # print(f"run_HOT_COS_GUIDANCE2 converged at tree #{alt_msas}\n")
                # if alt_msas < args_library.convergence:
                #     with lock:
                #         args_library.convergence = alt_msas
                #         args_library.count_convergence.value += 1
                with lock:
                    # args_library.convergence = alt_msas
                    args_library.count_convergence.value += 1
                # print(f'.done {proc}, generated {alt_msas}', flush=True)
                # print(args_library.count_convergence.value)
                if args_library.proc_num > 2:
                    print(f'.done {proc}', flush=True)
                    break
                else:
                    if args_library.count_convergence.value >= 2:
                        print(f'.done {proc}', flush=True)
                        break

    # end of child // end of process
    sys.exit(0)


# else:
#     e = Exception
#     raise Exception(f"ERROR: fork failed: {e}\n")
#@timeit
def run_guidance(args_library):
    # Align
    ##############
    if args_library.form_user_MSA_File == None or args_library.form_user_MSA_File == "":  # align if user did not supplied alignment
        align(args_library)  # ADD ALIGN

    # handle the adjustdirection
    if "adjustdirection" in args_library.align_param or "adjustdirectionaccurately" in args_library.align_param:  # support the --adjustdirectionaccurately and --adjustdirection option in MAFFT
        if args_library.form_MSA_Program == "MAFFT" and args_library.form_Seq_Type == "Nucleotides" and not args_library.form_userMSA_File:
            # back up the original seq provided
            shutil.copy(f"{args_library.WorkingDir}{args_library.codded_seq_fileName}",
                        f"{args_library.WorkingDir}{args_library.codded_seq_fileName}.OrigDirection")
            shutil.copy(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                        f"{args_library.WorkingDir}{args_library.Alignment_File}.OrigMAFFT")

            try:
                log_file = open(args_library.OutLogFile, "a")
            except OSError:
                print("run_guidance() could not open log file\n")
                sys.exit()
            # After the base alignment, take all the sequences out of the alignment and continue with them
            log_file.write("adjustdirection type option is in use; going to:\n")
            log_file.write(
                f"\t\t(1) copy {args_library.WorkingDir}{args_library.codded_seq_fileName} to {args_library.WorkingDir}{args_library.codded_seq_fileName}.OrigDirection\n")
            log_file.write(
                f"\t\t(2) copy {args_library.WorkingDir}{args_library.Alignment_File} to {args_library.WorkingDir}{args_library.Alignment_File}.OrigMAFFT\n")
            log_file.write(
                f"\t\t(3) Fix {args_library.WorkingDir}{args_library.codded_seq_fileName} to contain the sequences in the correct direction\n")
            log_file.write(
                f"\t\t(4) Fix {args_library.WorkingDir}{args_library.Alignment_File} to remove the '_R_' prefix added by MAFFT to the flipped seq names\n")

            try:
                MSA = open(f"{args_library.WorkingDir}{args_library.Alignment_File}.OrigMAFFT")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {args_library.WorkingDir}{args_library.Alignment_File}' for reading {e}\n", args_library)
            try:
                SEQ_FILE = open(f"{args_library.WorkingDir}{args_library.codded_seq_fileName}", "w")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {args_library.WorkingDir}{args_library.codded_seq_fileName}' for writing {e}\n", args_library)
            try:
                MSA_NEW = open(f"{args_library.WorkingDir}{args_library.Alignment_File}", "w")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {args_library.WorkingDir}{args_library.Alignment_File} for writing {e}", args_library)

            FRAGMENT_FILE = None
            NEW_CORE_ALN = None
            if "addfragments" in args_library.align_param:  # alignment with addfragments
                log_file.write(
                    f"\t\t(5) Copy original core alignment {args_library.WorkingDir}{args_library.Core_Alignment_File} to {args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT")
                shutil.copy(f"{args_library.WorkingDir}{args_library.Core_Alignment_File}",
                            f"{args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT")
                log_file.write(
                    f"\t\t(6) Copy {args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded} to {args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}.OrigDirection")
                shutil.copy(f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}",
                            f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}.OrigDirection")
                log_file.write(
                    f"\t\t(7) Fix {args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded} to contain the fragments in the correct direction")
                log_file.write(
                    f"\t\t(8) Fix {args_library.WorkingDir}{args_library.Core_Alignment_File} to remove the '_R_' prefix added by MAFFT to the flipped seq names")

                try:
                    FRAGMENT_FILE = open(
                        f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}", "w")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded} for writing {e}", args_library)
                try:
                    CORE_ALN = open(f"{args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT",
                                    "r")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT for reading {e}", args_library)
                try:
                    NEW_CORE_ALN = open(f"{args_library.WorkingDir}{args_library.Core_Alignment_File}", "w")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {args_library.WorkingDir}{args_library.Core_Alignment_File} for writing {e}", args_library)

                log_file.write(f"TOTAL SEQ IN CORE ALN: {args_library.NumOfSeq}")

                for line in CORE_ALN:
                    if match := re.match(r'^>_R_(.*)', line):
                        NEW_CORE_ALN.write(f">{match.group(1)}")
                        log_file.write(f"[NOTICE] MAFFT REVERSE SEQUENCE {match.group(1)} in the CORE ALIGNMENT\n")
                    elif match := re.match(r'^>', line):
                        NEW_CORE_ALN.write(line + "\n")
                    else:
                        NEW_CORE_ALN.write(line + "\n")

                CORE_ALN.close()
                NEW_CORE_ALN.close()

            IsReversedSeq = 0
            i = 0
            for line in MSA:
                line = line.strip()
                if match := re.match(r'^>_R_(.*)', line):
                    reverse_seq_name = match.group(1)
                    i += 1
                    IsReversedSeq = 1
                    if i <= args_library.NumOfSeq:
                        SEQ_FILE.write(f">{reverse_seq_name}\n")
                        log_file.write(f"[NOTICE] MAFFT REVERSE SEQUENCE {reverse_seq_name}\n")
                    elif "addfragments" in args_library.align_param:
                        FRAGMENT_FILE.write(f">{reverse_seq_name}\n")
                        log_file.write(f"[NOTICE] MAFFT REVERSE FRAGMENT SEQUENCE {reverse_seq_name}\n")
                    MSA_NEW.write(f">{reverse_seq_name}\n")
                elif match := re.match(r'^>', line):
                    i += 1
                    if i <= args_library.NumOfSeq:
                        SEQ_FILE.write(line + "\n")
                    elif "addfragments" in args_library.align_param:  # fragment
                        FRAGMENT_FILE.write(line + "\n")
                    MSA_NEW.write(line + "\n")
                else:
                    MSA_NEW.write(line + "\n")
                    line = line.replace("-", "")
                    if line != "" and i <= args_library.NumOfSeq:
                        SEQ_FILE.write(line + "\n")
                    elif line != "" and "addfragments" in args_library.align_param:  # fragment
                        FRAGMENT_FILE.write(line + "\n")

            MSA.close()
            SEQ_FILE.close()
            MSA_NEW.close()
            if "addfragments" in args_library.align_param:
                FRAGMENT_FILE.close()

            if IsReversedSeq == 0:
                shutil.copy(f"{args_library.WorkingDir}{args_library.Alignment_File}.OrigMAFFT",
                            f"{args_library.WorkingDir}{args_library.Alignment_File}")
                shutil.copy(f"{args_library.WorkingDir}{args_library.codded_seq_fileName}.OrigDirection",
                            f"{args_library.WorkingDir}{args_library.codded_seq_fileName}")
                os.unlink(f"{args_library.WorkingDir}{args_library.codded_seq_fileName}.OrigDirection")
                os.unlink(f"{args_library.WorkingDir}{args_library.Alignment_File}.OrigMAFFT")

                if '--addfragments' in args_library.align_param:
                    shutil.copy(f"{args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT",
                                f"{args_library.WorkingDir}{args_library.Core_Alignment_File}")
                    shutil.copy(
                        f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}.OrigDirection",
                        f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}")
                    os.unlink(
                        f"{args_library.WorkingDir}{args_library.fragments_file_name_seqName_coded}.OrigDirection")
                    os.unlink(f"{args_library.WorkingDir}{args_library.Core_Alignment_File}.OrigMAFFT")

                log_file.write("No flipped sequences were found... NOTHING OF THE ABOVE MENTIONED WAS NEEDED!!\n")

            log_file.close()

        ans = convert_names_of_align_with_seed(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                               f"{args_library.WorkingDir}{args_library.Alignment_File}_new")
        if ans[0] != "ok":
            exit_on_error("sys_error",
                          f"ConvertNamesOfAlignWithSeed({args_library.WorkingDir}{args_library.Alignment_File}, {args_library.WorkingDir}{args_library.Alignment_File}_new): {' '.join(ans)}\n", args_library)
        else:
            os.rename(f"{args_library.WorkingDir}{args_library.Alignment_File}_new",
                      f"{args_library.WorkingDir}{args_library.Alignment_File}")

#@timeit
def run_guidance2(args_library):
    args_library.Scoring_Alignments_Dir = args_library.GUIDANCE2_MSAs_Dir

    # INIT
    if args_library.MSA_Program in ["MAFFT", "MAFFT_LINSI"]:
        args_library.HoT_MSA_Program = "MFT"
        args_library.HoT_MSA_Program_path = args_library.mafft_prog
        check_mafft_profile = subprocess.getoutput(f"which {args_library.HoT_MSA_Program_path}-profile")
        if "Command not found" in check_mafft_profile:
            raise Exception(
                f"It seems that {args_library.HoT_MSA_Program_path}-profile is not properly installed or found in PATH args_library. "
                "Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n")

    if args_library.MSA_Program == "MUSCLE":
        exit_on_error('user_error', "GUIDANCE2 currently does not support MUSCLE, please run GUIDANCE<br>")
    if args_library.MSA_Program == "PAGAN":
        exit_on_error('user_error', "GUIDANCE2 currently does not support PAGAN, please run GUIDANCE<br>")
    elif args_library.MSA_Program == "CLUSTALO":
        args_library.HoT_MSA_Program = "CLO"
        args_library.HoT_MSA_Program_path = args_library.clustalw_prog
    elif args_library.MSA_Program == "PRANK":
        args_library.HoT_MSA_Program = "PRK"
        args_library.HoT_MSA_Program_path = args_library.prank_prog

    try:
        log_file = open(f'{args_library.OutLogFile}', "a")
    except Exception as e:
        print(f"Could not open Log file: {e}\n")
        sys.exit()
    # number of alignment to sample out of HoT (GUIDANCE2) perturbed alignment (per split,tree and op value, 4 or 8 HoT outputs are created)
    Num_of_Aln_from_HoT_per_Run = 4
    # HoT assumes that all the sequences are upper case
    log_file.write(f"convert_fs_to_upper_case({args_library.WorkingDir}{args_library.codded_seq_fileName})\n")
    convert_fs_to_upper_case(f"{args_library.WorkingDir}{args_library.codded_seq_fileName}")

    # Align
    ##############
    if args_library.userMSA_File == "" or args_library.userMSA_File is None:  # align if user did not supply alignment
        align(args_library)
        # HOT ASSUME THAT THE SEQUENCES ARE ALL UPPER CASE, SO WE CONVERT THE ALN TO UPPER CASE
        convert_fs_to_upper_case(
            f"{args_library.WorkingDir}{args_library.Alignment_File}")

    # TO DO: handle the adjustdirection

    # BootStrap Trees
    ##################
    Bootstrap_Trees(args_library)

    # pull out the trees
    ######################
    numUniqueTrees = ""
    numRepeats4UniqueTree = []

    if str(args_library.BBL).upper() == "YES":
        log_file.write(
            f"Guidance::pullOutBPtrees_BBL({args_library.WorkingDir}, {args_library.dataset}, {args_library.Bootstraps}, {args_library.MSA_Program});\n")
        ans = pull_out_bp_trees_bbl(args_library.WorkingDir, args_library.dataset, args_library.Bootstraps,
                                    args_library.MSA_Program)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::pullOutBPtrees_BBL: {' '.join(ans)}\n")
        if args_library.MSA_Program != "MAFFT":
            numUniqueTrees = ans[1]
            numRepeats4UniqueTree = ans[2]
    else:
        log_file.write(
            f"Guidance::pullOutBPtrees({args_library.WorkingDir}, {args_library.dataset}, {args_library.Bootstraps}, {args_library.MSA_Program});\n")
        ans = pull_out_bp_trees(args_library.WorkingDir, args_library.dataset, args_library.Bootstraps,
                                args_library.MSA_Program, args_library)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::pullOutBPtrees: {' '.join(ans)}\n", args_library)
        if args_library.MSA_Program != "MAFFT":
            numUniqueTrees = ans[1]
            numRepeats4UniqueTree = ans[2]

    # Convert trees to MAFFT format
    #################################
    if args_library.MSA_Program == "MAFFT" or args_library.MSA_Program == "MAFFT_LINSI":  # FOR MAFFT BUILDED ALIGNMENT FOR ALL TREES
        # if addfragments is used first prune the trees for the 'core' alignment
        # and afterward will use the created 'BP core' and the 'full tree' to create the full alignment
        if args_library.align_param.find('addfragments') != -1:
            args_library.prune_BootStrap_Dir = args_library.BootStrap_Dir + "PRUNE_BP_FOR_CORE_ALN/"
            # prune trees
            if not os.path.exists(args_library.prune_BootStrap_Dir):
                os.system(f"mkdir {args_library.prune_BootStrap_Dir}")
            for full_tree_dir in glob.glob(f"{args_library.BootStrap_Dir}tree*"):
                if os.path.isdir(full_tree_dir):
                    match = re.search(r'tree_([0-9]+)', full_tree_dir)
                    if match:
                        tree_num = match.group(1)
                        pruned_tree_dir = f"{args_library.prune_BootStrap_Dir}{tree_num}"
                        if not os.path.exists(pruned_tree_dir):
                            os.system(f"mkdir {pruned_tree_dir}")
                        full_tree = f"{full_tree_dir}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.{tree_num}"
                        pruned_tree = f"{pruned_tree_dir}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.{tree_num}CORE"
                        os.system(
                            f"{args_library.remove_taxa_prog} {full_tree} {args_library.WorkingDir}{args_library.fragments_codes} {pruned_tree}")
            # root pruned trees
            log_file.write(
                f"Guidance::root_BP_trees({args_library.prune_BootStrap_Dir},{args_library.dataset}, {args_library.MSA_Program}, {args_library.Bootstraps},'',{args_library.rooting_type});\n")
            ans = root_BP_trees(args_library.prune_BootStrap_Dir, args_library.dataset,
                                args_library.MSA_Program,
                                args_library.Bootstraps, "", args_library.rooting_type)
            if ans[0] != "ok":
                exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", args_library)
            if not os.path.exists(
                    f"{args_library.prune_BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}CORE.rooted") or os.path.getsize(
                f"{args_library.prune_BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}CORE.rooted") == 0:
                exit_on_error("sys_error",
                              f"{args_library.prune_BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}CORE.rooted does not exist/empty\n", args_library)  # TO DO: Consider to be numUniqueTrees instead Bootstraps
        # prepare the trees
        log_file.write(
            f"Guidance::root_BP_trees({args_library.BootStrap_Dir},{args_library.dataset}, {args_library.MSA_Program}, {args_library.Bootstraps},'',{args_library.rooting_type});\n")
        ans = root_BP_trees(args_library.BootStrap_Dir, args_library.dataset, args_library.MSA_Program,
                            args_library.Bootstraps,
                            "", args_library.rooting_type)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", args_library)
        if not os.path.exists(
                f"{args_library.BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted") or os.path.getsize(
            f"{args_library.BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted") == 0:
            exit_on_error("sys_error",
                          f"{args_library.BootStrap_Dir}tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted does not exist/empty\n", args_library)

    if args_library.MSA_Program == "CLUSTALO":
        # prepare the trees
        log_file.write(
            f"Guidance::root_BP_trees({args_library.BootStrap_Dir},{args_library.dataset}, {args_library.MSA_Program}, {args_library.Bootstraps},'',{args_library.rooting_type});\n")
        ans = root_BP_trees(f"{args_library.BootStrap_Dir}nonUniqueTrees/", args_library.dataset, args_library.MSA_Program,
                            args_library.Bootstraps,
                            "", args_library.rooting_type)

        # tree = f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{countTrees}.rooted"

        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", args_library)
        if not os.path.exists(
                f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted") or os.path.getsize(
            f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted") == 0:
            exit_on_error("sys_error",
                          f"{args_library.BootStrap_Dir}nonUniqueTrees/tree_{args_library.Bootstraps - 1}/{args_library.dataset}.{args_library.MSA_Program}.iqtree.tree_{args_library.Bootstraps - 1}.rooted does not exist/empty\n",
                          args_library)

    args_library.MSA_Depth = calculate_msa_depth(args_library.WorkingDir + args_library.Alignment_File, args_library)
    # Sample OP
    if (args_library.MSA_Program == "MAFFT" and args_library.PROGRAM == "GUIDANCE2"):
        OP_DistFile = MAFFT_OP_DIST
        # OP_DistFile = MAFFT_OP_DIST_0_25
    elif (args_library.MSA_Program == "MAFFT" and args_library.PROGRAM == "GUIDANCE3_HOT"):
        OP_DistFile = MAFFT_OP_DIST_0_25
        EP_DistFile = MAFFT_EP_DIST_0_25

    OutEP = os.path.join(args_library.WorkingDir, "SampledEPVals.log")
    OutOP = os.path.join(args_library.WorkingDir, "SampledOPVals.log")
    op_vals_arr_ref = ""
    ep_vals_arr_ref = ""

    if args_library.PROGRAM == "GUIDANCE2":
        if args_library.GapPenDist.upper() == "EMP":
            log_file.write(
                f"Sample op according to empiric distribution: Guidance::SampelFromEmpiricDistribution({OP_DistFile},{OutOP},{args_library.Bootstraps})\n")
            op_vals_arr_ref = sample_from_empirical_distribution(OP_DistFile, OutOP, args_library.FORM['Bootstraps'])
        elif args_library.GapPenDist.upper() == "UNIF":
            if args_library.MSA_Program == "MAFFT":
                log_file.write(
                    f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(1,3,{OutOP},{args_library.Bootstraps})\n")
                # according to mafft web-site defaults: http://mafft.cbrc.jp/alignment/server/index.html
                op_vals_arr_ref = sample_from_uniform_dist(1, 3, OutOP, args_library.Bootstraps)
            elif args_library.MSA_Program == "PRANK":
                log_file.write(
                    f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(0,0.5,{OutOP},{args_library.Bootstraps})\n")
                op_vals_arr_ref = sample_from_uniform_dist(0, 0.5, OutOP, args_library.Bootstraps)  # for prank v.140110 the defaults are: dna 0.025 / prot 0.005
            elif args_library.MSA_Program == "CLUSTALO":
                log_file.write(
                    f"Sample gap opening panelty according to uniform distribution: Guidance::SampleFromUniformDist(4,16,{OutOP},{args_library.Bootstraps})\n")
                op_vals_arr_ref = sample_from_uniform_dist(4, 16, OutOP, args_library.Bootstraps)

    if args_library.PROGRAM == "GUIDANCE3_HOT":
        if args_library.GapPenDist.upper() == "EMP":
            log_file.write(
                f"Sample op according to empiric distribution: Guidance::SampleFromEmpiricDistribution({OP_DistFile},{OutOP},{args_library.Bootstraps})\n")
            op_vals_arr_ref = sample_from_empirical_distribution(OP_DistFile, OutOP, args_library.Bootstraps)
            log_file.write(
                f"Sample ep according to empiric distribution: Guidance::SampleFromEmpiricDistribution({EP_DistFile},{OutEP},{args_library.Bootstraps})\n")
            ep_vals_arr_ref = sample_from_empirical_distribution(EP_DistFile, OutEP, args_library.FORM['Bootstraps'])
        if args_library.GapPenDist.upper() == "UNIF":
            log_file.write(
                f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(0,6,{OutOP},{args_library.Bootstraps})\n")
            op_vals_arr_ref = sample_from_uniform_dist(0, 6, OutOP, args_library.Bootstraps)
            log_file.write(
                f"Sample ep according to uniform distribution: Guidance::SampleFromUnifomDist(0,4,{OutEP},{args_library.Bootstraps})\n")
            ep_vals_arr_ref = sample_from_uniform_dist(0, 4, OutEP, args_library.Bootstraps)

    if args_library.isServer == 1:
        args_library.status_file = args_library.WorkingDir + "MSA_STATUS.txt"
        # with open(args_library.status_file, "w") as STATUS:
        #     STATUS.write("<ul><li><p><font face=Verdana size=2>Start creating alternative alignments<br></li></ul>\n")
        update_progress(f"{args_library.WorkingDir}{args_library.progress_report}", "Started generating alternative alignments")

        with open(f"{args_library.server_output}", "a") as OUTPUT:
        # with open(f"{args_library.WorkingDir}{args_library.server_output}", "a") as OUTPUT:
            OUTPUT.write(
                "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t$fil = fopen('MSA_STATUS.txt', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

    # Get random branches for HoT
    NumOfBranches = (2 * args_library.MSA_Depth) - 3
    RandomBranches = [random.randint(0, NumOfBranches - 1) for _ in range(args_library.Bootstraps)]

    # CREATE THE PERTURBED ALN DIR
    os.mkdir(args_library.Scoring_Alignments_Dir)

    countTrees = 0
    epsilon = 0.0006
    manager = Manager()
    lock = manager.Lock()
    args_library.mean_res_pair_score = manager.list()
    args_library.mean_col_score = manager.list()
    args_library.count_convergence = manager.Value('i', 0)
    # this args_library.convergence value might be unnecessary at the end and can be removed
    args_library.convergence = args_library.Bootstraps * Num_of_Aln_from_HoT_per_Run

    # Running parallel processes using multiprocessing, each will run an equal share of the BP alignments (?)
    processes = [Process(target=run_hot_process_on_tree, args=(args_library, epsilon, proc, RandomBranches, op_vals_arr_ref, ep_vals_arr_ref, Num_of_Aln_from_HoT_per_Run, lock)) for proc in range(args_library.proc_num)]
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    # this is the number of alternative MSAs produced at the end
    alt_msas = len(os.listdir(args_library.Scoring_Alignments_Dir))
    args_library.convergence = alt_msas
    log_file.write(f"run_HOT_COS_GUIDANCE2 converged at tree #{alt_msas}\n")
    log_file.close()


    # Ksenia removed this part
    if args_library.isServer == 1:
        update_progress(f"{args_library.WorkingDir}{args_library.progress_report}",
                        f"Finished generating {alt_msas} alternative alignments")
    #     with open(args_library.status_file, "w") as PROGRESS:
    #         PROGRESS.write(
    #             f"\n<ul><li>{alt_msas} out of {args_library.Bootstraps * 4} alternative alignments were created</li></ul>\n")
    #         PROGRESS.write(
    #             f"\n<ul><li>{alt_msas} alternative alignments were created</li></ul>\n")

    # To validate all alns were created
    # aln_count = len(os.listdir(args_library.Scoring_Alignments_Dir))
    # expected_count = Num_of_Aln_from_HoT_per_Run * args_library.Bootstraps
    # expected_count = (args_library.convergence) * Num_of_Aln_from_HoT_per_Run
    print(f"the convergence final number is {args_library.convergence}")
    # if aln_count < expected_count:
    #     exit_on_error("sys_error",
    #                   f"run_Guidance2: Only {aln_count} alignments were created on {args_library.Scoring_Alignments_Dir} while expecting {expected_count}\n", args_library)
    # else:
    #     print("\nSUCCESS!\n")
    # print(args_library.mean_res_pair_score)
    # print(args_library.mean_col_score)
    print("\nSUCCESS!\n")

#@timeit
def run_hot(args_library):
    #	python3 ../hot_cos_main.py caseID msa_method seq_type input_fasta_file . output_dir >& COS.std
    #msa_method: MA0 = mafft ; CW2 = clustalW.
    #seq_type: aa = amino-acid ; nt = nucleotides
    #input_fasta_file = the input sequences file

    #The base MSA is: output_dir_cos_msa_method /hot_H.fasta
    #So you can copy it:
    #cp output_dir_cos_msa_method/hot_H.fasta ./caseID_mafft.fasta

    #The MSA sets should be copied to a directory:
    #mkdir COS_MSA
    #cp ./output_dir_cos_msa_method/b0#*.fasta ./COS_MSA/
    try:
        with open(args_library.OutLogFile, "a") as log_file:
            if args_library.MSA_Program in ["MAFFT", "MAFFT_LINSI"]:
                args_library.HoT_MSA_Program = "MFT"
                args_library.HoT_MSA_Program_path = args_library.mafft_prog
                check_mafft_profile = subprocess.getoutput(f"which {args_library.HoT_MSA_Program_path} -profile")
                if "Command not found" in check_mafft_profile:
                    raise ValueError(
                        f"It seems that {args_library.HoT_MSA_Program_path}-profile is not properly installed or found in PATH args_library. Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n")

            elif args_library.MSA_Program == "MUSCLE":
                exit_on_error('user_error', "HoT currently does not support MUSCLE, please run GUIDANCE<br>", args_library)

            elif args_library.MSA_Program == "PAGAN":
                exit_on_error('user_error', "HoT currently does not support PAGAN, please run GUIDANCE<br>",args_library)

            elif args_library.MSA_Program == "CLUSTALO":
                args_library.HoT_MSA_Program = "CLO"
                args_library.HoT_MSA_Program_path = args_library.clustalw_prog

            elif args_library.MSA_Program == "PRANK":
                args_library.HoT_MSA_Program = "PRK"
                args_library.HoT_MSA_Program_path = args_library.prank_prog

            cmd = f"cd {args_library.WorkingDir}; python3 {HOT_PROGRAM} {args_library.dataset} {args_library.HoT_MSA_Program}"

            if args_library.Seq_Type in ["AminoAcids", "Codons"]:
                cmd += " aa"
            elif args_library.Seq_Type == "Nucleotides":
                cmd += " nt"

            log_file.write(f"convert_fs_to_upper_case({args_library.WorkingDir}{args_library.codded_seq_fileName})")
            convert_fs_to_upper_case(
                f"{args_library.WorkingDir}{args_library.codded_seq_fileName}")  # HoT assumes that all the sequences are upper case

            if args_library.align_param == "":
                cmd += f" {args_library.codded_seq_fileName} . {args_library.WorkingDir}MSA_STATUS.txt 0 {args_library.HoT_MSA_Program_path} \"\" all > COS.std"
            else:
                cmd += f" {args_library.codded_seq_fileName} . {args_library.WorkingDir}MSA_STATUS.txt 0 {args_library.HoT_MSA_Program_path} \"\" all --- {args_library.align_param} > COS.std"

            log_file.write(f"run_HoT: {cmd}\n")

            if args_library.isServer == 1:
                status_file = f"{args_library.WorkingDir}MSA_STATUS.txt"
                with open(status_file, "w") as STATUS:
                    STATUS.write("\n")

                with open(f"{status_file}.0", "w") as STATUS0:
                    STATUS0.write("\n")

                with open(f"{args_library.server_output}", "a") as OUTPUT:
                    OUTPUT.write(
                    "<?php\n\tif (file_exists('MSA_STATUS.txt.0'))\n\t{\n\t\t$fil =fopen('MSA_STATUS.txt.0', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt.0'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

                    OUTPUT.write(
                    "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t$fil = fopen('MSA_STATUS.txt', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

            os.system(cmd)

            if not os.path.exists(f"{args_library.WorkingDir}{args_library.Alignment_File}"):
                cmd = f"cp {args_library.WorkingDir}{args_library.dataset}_cos_{args_library.HoT_MSA_Program}/hot_H.fasta {args_library.WorkingDir}{args_library.Alignment_File}"
                log_file.write(f"run_HoT: {cmd}\n")  # Copy Alignment
                os.system(cmd)

            if args_library.Align_Order == "as_input":
                print(f"MSA_parser::sort_alignment({args_library.WorkingDir}{args_library.Alignment_File},fasta);\n")
                ans = sort_alignment(f"{args_library.WorkingDir}{args_library.Alignment_File}", "fasta")
                print("".join(ans))
                args_library.Alignment_File_NOT_SORTED = args_library.Alignment_File
                args_library.Alignment_File = args_library.Alignment_File + ".Sorted"

            os.system(f"mkdir {args_library.WorkingDir}{args_library.HoT_MSAs_Dir}")

            if args_library.NumOfSeq > 2:
                cmd = f"cp {args_library.WorkingDir}{args_library.dataset}_cos_{args_library.HoT_MSA_Program}/b[01]*.fasta {args_library.WorkingDir}{args_library.HoT_MSAs_Dir}"
            else:
                cmd = f"cp {args_library.WorkingDir}{args_library.dataset}_cos_{args_library.HoT_MSA_Program}/hot*.fasta {args_library.WorkingDir}{args_library.HoT_MSAs_Dir}"

            log_file.write(f"run_HoT: {cmd}\n")
            os.system(cmd)
            # HoT assumes that all the sequences are upper case
            if os.path.exists(f"{args_library.WorkingDir}{args_library.Alignment_File}.WithCodesName"):
                convert_fs_to_upper_case(f"{args_library.WorkingDir}{args_library.Alignment_File}.WithCodesName")
            args_library.Scoring_Alignments_Dir = f"{args_library.WorkingDir}{args_library.HoT_MSAs_Dir}"
    except Exception as e:
        sys.exit("run_hot() Error: " + str(e) + "\n")


