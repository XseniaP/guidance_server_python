from guidance3.sequences.fasta import *
from guidance3.trees.phylo import *
from guidance3.utils.common import *
from guidance3.pipeline.scoring import calculate_sp_scores, check_convergence, \
    calculate_sp_scores_convergence, add_scores_to_dict
import uuid
import sys
import re
import subprocess
import multiprocessing as mp
from guidance3.utils.timing import timeit
from multiprocessing.sharedctypes import Value, Array
from multiprocessing import Process, Manager, Lock
# import smtplib
# from email.message import EmailMessage
# from email.utils import formataddr
# from subprocess import run, PIPE
# import logging

#@timeit
def run_hot_internal(config, op_vals_arr_ref, ep_vals_arr_ref, countTrees, tree_good_BranchLength, Branch):
    # cwd=WorkingDir replaces the shell "cd WorkingDir;" prefix;
    # stdout= replaces the ">> COS.std" shell redirect.
    HOT_COS_GUIDANCE2_cmd = f"{HOT_GUIDANCE3_PROGRAM} {config.dataset}_{countTrees} {config.HoT_MSA_Program}"

    if config.Seq_Type in ["AminoAcids", "Codons"]:
        HOT_COS_GUIDANCE2_cmd += " aa"
    elif config.Seq_Type == "Nucleotides":
        HOT_COS_GUIDANCE2_cmd += " nt"

    HOT_COS_GUIDANCE2_cmd += f" {config.codded_seq_fileName} . \"\" 0 {config.HoT_MSA_Program_path} {tree_good_BranchLength} {Branch}"

    if config.MSA_Program == "MAFFT":
        if config.PROGRAM == "GUIDANCE3":
            HOT_COS_GUIDANCE2_cmd += f" --- {config.align_param} --op {op_vals_arr_ref[countTrees]}"
        elif config.PROGRAM == "GUIDANCE3_HOT":
            HOT_COS_GUIDANCE2_cmd += f" --- {config.align_param} --op {op_vals_arr_ref[countTrees]} --ep {ep_vals_arr_ref[countTrees]}"
    elif config.MSA_Program == "PRANK":
        HOT_COS_GUIDANCE2_cmd += f" --- {config.align_param} -gaprate={op_vals_arr_ref[countTrees]}"
    elif config.MSA_Program == "CLUSTALO":
        # HOT_COS_GUIDANCE2_cmd += f" --- {config.align_param} -GAPOPEN={op_vals_arr_ref[countTrees]}"
        HOT_COS_GUIDANCE2_cmd += f" --- {config.align_param}"

    print(HOT_COS_GUIDANCE2_cmd)

    cos_std_path = os.path.join(config.WorkingDir, "COS.std")
    with open(cos_std_path, "a") as cos_std:
        result = subprocess.run(
            HOT_COS_GUIDANCE2_cmd,
            shell=True,
            cwd=config.WorkingDir,
            stdout=cos_std,
            stderr=cos_std,
        )

    if result.returncode != 0:
        raise RuntimeError(f"HOT command failed (rc={result.returncode}): {HOT_COS_GUIDANCE2_cmd}")

    return HOT_COS_GUIDANCE2_cmd

#@timeit
def run_hot_process_on_tree(config, epsilon, proc, RandomBranches,op_vals_arr_ref, ep_vals_arr_ref, Num_of_Aln_from_HoT_per_Run, lock):
    try:
        log_file = open(config.OutLogFile, "a")
    except OSError:
        print("run_guidance() could not open log file\n")
        sys.exit()

    bp_per_proc = (config.Bootstraps // config.proc_num) + 1

    for tree_num in range(bp_per_proc):

        convergence = 0

        countTrees = proc + config.proc_num * tree_num
        print(f"proc num {proc}\ttree num {tree_num} --> global tree index {countTrees}\n")

        if countTrees >= config.Bootstraps:
            break

        try:
            Branch = RandomBranches[countTrees]

            if config.MSA_Program == "MAFFT":
                if 'addfragments' in config.align_param:
                    tree = f"{config.prune_BootStrap_Dir}tree_{countTrees}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{countTrees}CORE.rooted"
                else:
                    tree = f"{config.BootStrap_Dir}tree_{countTrees}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{countTrees}.rooted"

            elif config.MSA_Program == "CLUSTALO":
                tree = f"{config.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{countTrees}.rooted"

            else: # PRANK
                tree = f"{config.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{countTrees}"

            if config.MSA_Program == "PRANK":
                if 'iterate' in config.align_param:
                    config.align_param = re.sub(r'-iterate=\d+', '', config.align_param).strip()
                    log_file.write("[WARNING] -iterate argument is ignored for the perturbed alignments stage\n")
                    print("[WARNING] -iterate argument is ignored when reconstructing the perturbed alignments\n")
                if 'once' not in config.align_param:
                    config.align_param += " -once"

            tree_good_BranchLength = f"{tree}.GoodBranchLength"
            reformat_trees_branch_length(tree, tree_good_BranchLength)

            HOT_COS_GUIDANCE2_cmd = run_hot_internal(config, op_vals_arr_ref, ep_vals_arr_ref, countTrees, tree_good_BranchLength, Branch)
            log_file.write(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")
            print(f"run_HOT_COS_GUIDANCE2: {HOT_COS_GUIDANCE2_cmd}\n")

            local_dataset = f"{config.dataset}_{countTrees}_cos_"
            pathname = os.path.join(
                f"{config.WorkingDir}", f"{local_dataset}{config.HoT_MSA_Program}/b[01]*.fasta")
            pertubed_aln = glob.glob(pathname)

            if not pertubed_aln:
                raise RuntimeError(f"HOT produced no output files for tree {countTrees} (pattern: {pathname})")

            n_to_use = min(len(pertubed_aln), Num_of_Aln_from_HoT_per_Run)
            if n_to_use < Num_of_Aln_from_HoT_per_Run:
                log_file.write(f"[WARNING] tree {countTrees}: HOT produced {len(pertubed_aln)} files, expected {Num_of_Aln_from_HoT_per_Run}, using {n_to_use}\n")
                print(f"[WARNING] tree {countTrees}: HOT produced {len(pertubed_aln)} files, expected {Num_of_Aln_from_HoT_per_Run}, using {n_to_use}")

            shuffled = random.sample(pertubed_aln, len(pertubed_aln))

            for j in range(n_to_use):
                aln = shuffled[j]
                base = os.path.basename(aln).split(".fasta")[0]

                if config.PROGRAM == "GUIDANCE3":
                    cp_from = os.path.join(config.WorkingDir, f"{config.dataset}_{countTrees}_cos_{config.HoT_MSA_Program}", f"{base}.fasta")
                    cp_to = os.path.join(config.Scoring_Alignments_Dir, f"{base}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_Split_{Branch}.fasta")
                elif config.PROGRAM == "GUIDANCE3_HOT":
                    cp_from = os.path.join(config.WorkingDir, f"{config.dataset}_{countTrees}_cos_{config.HoT_MSA_Program}", f"{base}.fasta")
                    cp_to = os.path.join(config.Scoring_Alignments_Dir, f"{base}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_EP_{ep_vals_arr_ref[countTrees]}_Split_{Branch}.fasta")
                else:
                    continue

                shutil.copy(cp_from, cp_to)

            if config.PROGRAM == "GUIDANCE3":
                mv_src = os.path.join(config.WorkingDir, f"{config.dataset}_{countTrees}_cos_{config.HoT_MSA_Program}")
                mv_dst = os.path.join(config.BootStrap_Dir, f"{config.dataset}_cos_{config.HoT_MSA_Program}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}_Split_{Branch}")
            elif config.PROGRAM == "GUIDANCE3_HOT":
                mv_src = os.path.join(config.WorkingDir, f"{config.dataset}_{countTrees}_cos_{config.HoT_MSA_Program}")
                mv_dst = os.path.join(config.BootStrap_Dir, f"{config.dataset}_cos_{config.HoT_MSA_Program}_tree_{countTrees}_OP_{op_vals_arr_ref[countTrees]}EP_{ep_vals_arr_ref[countTrees]}_Split_{Branch}")
            else:
                mv_src = mv_dst = None

            if mv_src and mv_dst:
                shutil.move(mv_src, mv_dst)

            # check convergence starting from the 20th tree (starting from 80 MSAs)
            if not getattr(config, 'disable_convergence', False) and countTrees >= 20:
                # check the convergence only for every nth tree
                if config.proc_num >= 2 or (config.proc_num == 1 and countTrees % 3 == 0):
                    try:
                        # Serialize only the FS read: calculate_sp_scores_convergence reads
                        # Scoring_Alignments_Dir while other processes are still copying into it.
                        # add_scores_to_dict and check_convergence manage their own thread safety.
                        with lock:
                            alt_msas = calculate_sp_scores_convergence(config, countTrees)
                        add_scores_to_dict(config, epsilon, countTrees, lock)
                        for scr_file in glob.glob(os.path.join(config.WorkingDir, f"{config.Output_Prefix}_tree_{countTrees}_*.scr")):
                            os.unlink(scr_file)
                        convergence = check_convergence(config, epsilon)
                        print(f"convergence of proc num {proc}\ttree num {tree_num} --> global tree index {countTrees} is {convergence} \n")
                    except Exception as e:
                        log_file.write(f"failed to calculate scores for convergence of proc num {proc}\ttree num {tree_num} \t# of alternative MSAs {alt_msas} error {e}\n")
                        print(f"failed to calculate scores for convergence of proc num {proc}\ttree num {tree_num} \t# of alternative MSAs {alt_msas} \n")
            if convergence == 1:
                with lock:
                    config.count_convergence.value += 1
                # With a single process there is no parallelism to coordinate,
                # so break on the first detection. With multiple processes, require
                # at least 2 to independently agree before stopping any worker.
                convergence_threshold = 1 if config.proc_num == 1 else 2
                if config.count_convergence.value >= convergence_threshold:
                    print(f'.done {proc}', flush=True)
                    break

        except Exception as e:
            log_file.write(f"[ERROR] proc {proc} tree {countTrees} failed, skipping: {e}\n")
            print(f"[ERROR] proc {proc} tree {countTrees} failed, skipping: {e}")
            continue

    # end of child // end of process
    return


#@timeit
def run_guidance(config):
    # Align
    ##############
    if config.form_user_MSA_File == None or config.form_user_MSA_File == "":  # align if user did not supplied alignment
        align(config)  # ADD ALIGN

    # handle the adjustdirection
    if "adjustdirection" in config.align_param or "adjustdirectionaccurately" in config.align_param:  # support the --adjustdirectionaccurately and --adjustdirection option in MAFFT
        if config.form_MSA_Program == "MAFFT" and config.form_Seq_Type == "Nucleotides" and not config.form_userMSA_File:
            # back up the original seq provided
            shutil.copy(f"{config.WorkingDir}{config.codded_seq_fileName}",
                        f"{config.WorkingDir}{config.codded_seq_fileName}.OrigDirection")
            shutil.copy(f"{config.WorkingDir}{config.Alignment_File}",
                        f"{config.WorkingDir}{config.Alignment_File}.OrigMAFFT")

            try:
                log_file = open(config.OutLogFile, "a")
            except OSError:
                print("run_guidance() could not open log file\n")
                sys.exit()
            # After the base alignment, take all the sequences out of the alignment and continue with them
            log_file.write("adjustdirection type option is in use; going to:\n")
            log_file.write(
                f"\t\t(1) copy {config.WorkingDir}{config.codded_seq_fileName} to {config.WorkingDir}{config.codded_seq_fileName}.OrigDirection\n")
            log_file.write(
                f"\t\t(2) copy {config.WorkingDir}{config.Alignment_File} to {config.WorkingDir}{config.Alignment_File}.OrigMAFFT\n")
            log_file.write(
                f"\t\t(3) Fix {config.WorkingDir}{config.codded_seq_fileName} to contain the sequences in the correct direction\n")
            log_file.write(
                f"\t\t(4) Fix {config.WorkingDir}{config.Alignment_File} to remove the '_R_' prefix added by MAFFT to the flipped seq names\n")

            try:
                MSA = open(f"{config.WorkingDir}{config.Alignment_File}.OrigMAFFT")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {config.WorkingDir}{config.Alignment_File}' for reading {e}\n", config)
            try:
                SEQ_FILE = open(f"{config.WorkingDir}{config.codded_seq_fileName}", "w")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {config.WorkingDir}{config.codded_seq_fileName}' for writing {e}\n", config)
            try:
                MSA_NEW = open(f"{config.WorkingDir}{config.Alignment_File}", "w")
            except Exception as e:
                exit_on_error("sys_error",
                              f"Can't open {config.WorkingDir}{config.Alignment_File} for writing {e}", config)

            FRAGMENT_FILE = None
            NEW_CORE_ALN = None
            if "addfragments" in config.align_param:  # alignment with addfragments
                log_file.write(
                    f"\t\t(5) Copy original core alignment {config.WorkingDir}{config.Core_Alignment_File} to {config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT")
                shutil.copy(f"{config.WorkingDir}{config.Core_Alignment_File}",
                            f"{config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT")
                log_file.write(
                    f"\t\t(6) Copy {config.WorkingDir}{config.fragments_file_name_seqName_coded} to {config.WorkingDir}{config.fragments_file_name_seqName_coded}.OrigDirection")
                shutil.copy(f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}",
                            f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}.OrigDirection")
                log_file.write(
                    f"\t\t(7) Fix {config.WorkingDir}{config.fragments_file_name_seqName_coded} to contain the fragments in the correct direction")
                log_file.write(
                    f"\t\t(8) Fix {config.WorkingDir}{config.Core_Alignment_File} to remove the '_R_' prefix added by MAFFT to the flipped seq names")

                try:
                    FRAGMENT_FILE = open(
                        f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}", "w")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {config.WorkingDir}{config.fragments_file_name_seqName_coded} for writing {e}", config)
                try:
                    CORE_ALN = open(f"{config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT",
                                    "r")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT for reading {e}", config)
                try:
                    NEW_CORE_ALN = open(f"{config.WorkingDir}{config.Core_Alignment_File}", "w")
                except Exception as e:
                    exit_on_error("sys_error",
                                  f"Can't open {config.WorkingDir}{config.Core_Alignment_File} for writing {e}", config)

                log_file.write(f"TOTAL SEQ IN CORE ALN: {config.NumOfSeq}")

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
                    if i <= config.NumOfSeq:
                        SEQ_FILE.write(f">{reverse_seq_name}\n")
                        log_file.write(f"[NOTICE] MAFFT REVERSE SEQUENCE {reverse_seq_name}\n")
                    elif "addfragments" in config.align_param:
                        FRAGMENT_FILE.write(f">{reverse_seq_name}\n")
                        log_file.write(f"[NOTICE] MAFFT REVERSE FRAGMENT SEQUENCE {reverse_seq_name}\n")
                    MSA_NEW.write(f">{reverse_seq_name}\n")
                elif match := re.match(r'^>', line):
                    i += 1
                    if i <= config.NumOfSeq:
                        SEQ_FILE.write(line + "\n")
                    elif "addfragments" in config.align_param:  # fragment
                        FRAGMENT_FILE.write(line + "\n")
                    MSA_NEW.write(line + "\n")
                else:
                    MSA_NEW.write(line + "\n")
                    line = line.replace("-", "")
                    if line != "" and i <= config.NumOfSeq:
                        SEQ_FILE.write(line + "\n")
                    elif line != "" and "addfragments" in config.align_param:  # fragment
                        FRAGMENT_FILE.write(line + "\n")

            MSA.close()
            SEQ_FILE.close()
            MSA_NEW.close()
            if "addfragments" in config.align_param:
                FRAGMENT_FILE.close()

            if IsReversedSeq == 0:
                shutil.copy(f"{config.WorkingDir}{config.Alignment_File}.OrigMAFFT",
                            f"{config.WorkingDir}{config.Alignment_File}")
                shutil.copy(f"{config.WorkingDir}{config.codded_seq_fileName}.OrigDirection",
                            f"{config.WorkingDir}{config.codded_seq_fileName}")
                os.unlink(f"{config.WorkingDir}{config.codded_seq_fileName}.OrigDirection")
                os.unlink(f"{config.WorkingDir}{config.Alignment_File}.OrigMAFFT")

                if '--addfragments' in config.align_param:
                    shutil.copy(f"{config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT",
                                f"{config.WorkingDir}{config.Core_Alignment_File}")
                    shutil.copy(
                        f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}.OrigDirection",
                        f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}")
                    os.unlink(
                        f"{config.WorkingDir}{config.fragments_file_name_seqName_coded}.OrigDirection")
                    os.unlink(f"{config.WorkingDir}{config.Core_Alignment_File}.OrigMAFFT")

                log_file.write("No flipped sequences were found... NOTHING OF THE ABOVE MENTIONED WAS NEEDED!!\n")

            log_file.close()

        ans = convert_names_of_align_with_seed(f"{config.WorkingDir}{config.Alignment_File}",
                                               f"{config.WorkingDir}{config.Alignment_File}_new")
        if ans[0] != "ok":
            exit_on_error("sys_error",
                          f"ConvertNamesOfAlignWithSeed({config.WorkingDir}{config.Alignment_File}, {config.WorkingDir}{config.Alignment_File}_new): {' '.join(ans)}\n", config)
        else:
            os.rename(f"{config.WorkingDir}{config.Alignment_File}_new",
                      f"{config.WorkingDir}{config.Alignment_File}")

#@timeit
def run_guidance3(config):
    config.Scoring_Alignments_Dir = config.GUIDANCE3_MSAs_Dir

    # INIT
    if config.MSA_Program in ["MAFFT", "MAFFT_LINSI"]:
        config.HoT_MSA_Program = "MFT"
        config.HoT_MSA_Program_path = config.mafft_prog
        check_mafft_profile = subprocess.getoutput(f"which {config.HoT_MSA_Program_path}-profile")
        if "Command not found" in check_mafft_profile:
            raise Exception(
                f"It seems that {config.HoT_MSA_Program_path}-profile is not properly installed or found in PATH config. "
                "Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n")

    if config.MSA_Program == "MUSCLE":
        exit_on_error('user_error', "GUIDANCE2 currently does not support MUSCLE, please run GUIDANCE<br>")
    if config.MSA_Program == "PAGAN":
        exit_on_error('user_error', "GUIDANCE2 currently does not support PAGAN, please run GUIDANCE<br>")
    elif config.MSA_Program == "CLUSTALO":
        config.HoT_MSA_Program = "CLO"
        config.HoT_MSA_Program_path = config.clustalw_prog
    elif config.MSA_Program == "PRANK":
        config.HoT_MSA_Program = "PRK"
        config.HoT_MSA_Program_path = config.prank_prog

    try:
        log_file = open(f'{config.OutLogFile}', "a")
    except Exception as e:
        print(f"Could not open Log file: {e}\n")
        sys.exit()
    # number of alignment to sample out of HoT (GUIDANCE2) perturbed alignment (per split,tree and op value, 4 or 8 HoT outputs are created)
    Num_of_Aln_from_HoT_per_Run = 4
    # HoT assumes that all the sequences are upper case
    log_file.write(f"convert_fs_to_upper_case({config.WorkingDir}{config.codded_seq_fileName})\n")
    convert_fs_to_upper_case(f"{config.WorkingDir}{config.codded_seq_fileName}")

    # Align
    ##############
    if config.userMSA_File == "" or config.userMSA_File is None:  # align if user did not supply alignment
        align(config)
        # HOT ASSUME THAT THE SEQUENCES ARE ALL UPPER CASE, SO WE CONVERT THE ALN TO UPPER CASE
        convert_fs_to_upper_case(
            f"{config.WorkingDir}{config.Alignment_File}")

    else:
        convert_fs_to_upper_case(
            f"{config.WorkingDir}{config.Alignment_File}")

    # TO DO: handle the adjustdirection

    # BootStrap Trees
    ##################
    Bootstrap_Trees(config)

    # pull out the trees
    ######################
    numUniqueTrees = ""
    numRepeats4UniqueTree = []

    if str(config.BBL).upper() == "YES":
        log_file.write(
            f"Guidance::pullOutBPtrees_BBL({config.WorkingDir}, {config.dataset}, {config.Bootstraps}, {config.MSA_Program});\n")
        ans = pull_out_bp_trees_bbl(config.WorkingDir, config.dataset, config.Bootstraps,
                                    config.MSA_Program)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::pullOutBPtrees_BBL: {' '.join(ans)}\n")
        if config.MSA_Program != "MAFFT":
            numUniqueTrees = ans[1]
            numRepeats4UniqueTree = ans[2]
    else:
        log_file.write(
            f"Guidance::pullOutBPtrees({config.WorkingDir}, {config.dataset}, {config.Bootstraps}, {config.MSA_Program});\n")
        ans = pull_out_bp_trees(config.WorkingDir, config.dataset, config.Bootstraps,
                                config.MSA_Program, config)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::pullOutBPtrees: {' '.join(ans)}\n", config)
        if config.MSA_Program != "MAFFT":
            numUniqueTrees = ans[1]
            numRepeats4UniqueTree = ans[2]

    # Convert trees to MAFFT format
    #################################
    if config.MSA_Program == "MAFFT" or config.MSA_Program == "MAFFT_LINSI":  # FOR MAFFT BUILDED ALIGNMENT FOR ALL TREES
        # if addfragments is used first prune the trees for the 'core' alignment
        # and afterward will use the created 'BP core' and the 'full tree' to create the full alignment
        if config.align_param.find('addfragments') != -1:
            config.prune_BootStrap_Dir = config.BootStrap_Dir + "PRUNE_BP_FOR_CORE_ALN/"
            # prune trees
            if not os.path.exists(config.prune_BootStrap_Dir):
                os.system(f"mkdir {config.prune_BootStrap_Dir}")
            for full_tree_dir in glob.glob(f"{config.BootStrap_Dir}tree*"):
                if os.path.isdir(full_tree_dir):
                    match = re.search(r'tree_([0-9]+)', full_tree_dir)
                    if match:
                        tree_num = match.group(1)
                        pruned_tree_dir = f"{config.prune_BootStrap_Dir}{tree_num}"
                        if not os.path.exists(pruned_tree_dir):
                            os.system(f"mkdir {pruned_tree_dir}")
                        full_tree = f"{full_tree_dir}/{config.dataset}.{config.MSA_Program}.iqtree.{tree_num}"
                        pruned_tree = f"{pruned_tree_dir}/{config.dataset}.{config.MSA_Program}.iqtree.{tree_num}CORE"
                        os.system(
                            f"{config.remove_taxa_prog} {full_tree} {config.WorkingDir}{config.fragments_codes} {pruned_tree}")
            # root pruned trees
            log_file.write(
                f"Guidance::root_BP_trees({config.prune_BootStrap_Dir},{config.dataset}, {config.MSA_Program}, {config.Bootstraps},'',{config.rooting_type});\n")
            ans = root_BP_trees(config.prune_BootStrap_Dir, config.dataset,
                                config.MSA_Program,
                                config.Bootstraps, "", config.rooting_type)
            if ans[0] != "ok":
                exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", config)
            if not os.path.exists(
                    f"{config.prune_BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}CORE.rooted") or os.path.getsize(
                f"{config.prune_BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}CORE.rooted") == 0:
                exit_on_error("sys_error",
                              f"{config.prune_BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}CORE.rooted does not exist/empty\n", config)  # TO DO: Consider to be numUniqueTrees instead Bootstraps
        # prepare the trees
        log_file.write(
            f"Guidance::root_BP_trees({config.BootStrap_Dir},{config.dataset}, {config.MSA_Program}, {config.Bootstraps},'',{config.rooting_type});\n")
        ans = root_BP_trees(config.BootStrap_Dir, config.dataset, config.MSA_Program,
                            config.Bootstraps,
                            "", config.rooting_type)
        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", config)
        if not os.path.exists(
                f"{config.BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted") or os.path.getsize(
            f"{config.BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted") == 0:
            exit_on_error("sys_error",
                          f"{config.BootStrap_Dir}tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted does not exist/empty\n", config)

    if config.MSA_Program == "CLUSTALO":
        # prepare the trees
        log_file.write(
            f"Guidance::root_BP_trees({config.BootStrap_Dir},{config.dataset}, {config.MSA_Program}, {config.Bootstraps},'',{config.rooting_type});\n")
        ans = root_BP_trees(f"{config.BootStrap_Dir}nonUniqueTrees/", config.dataset, config.MSA_Program,
                            config.Bootstraps,
                            "", config.rooting_type)

        # tree = f"{config.BootStrap_Dir}nonUniqueTrees/tree_{countTrees}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{countTrees}.rooted"

        if ans[0] != "ok":
            exit_on_error("sys_error", f"Guidance::root_BP_trees: {' '.join(ans)}\n", config)
        if not os.path.exists(
                f"{config.BootStrap_Dir}nonUniqueTrees/tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted") or os.path.getsize(
            f"{config.BootStrap_Dir}nonUniqueTrees/tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted") == 0:
            exit_on_error("sys_error",
                          f"{config.BootStrap_Dir}nonUniqueTrees/tree_{config.Bootstraps - 1}/{config.dataset}.{config.MSA_Program}.iqtree.tree_{config.Bootstraps - 1}.rooted does not exist/empty\n",
                          config)

    config.MSA_Depth = calculate_msa_depth(config.WorkingDir + config.Alignment_File, config)
    # Sample OP
    if (config.MSA_Program == "MAFFT" and config.PROGRAM == "GUIDANCE3"):
        OP_DistFile = MAFFT_OP_DIST
        # OP_DistFile = MAFFT_OP_DIST_0_25
    elif (config.MSA_Program == "MAFFT" and config.PROGRAM == "GUIDANCE3_HOT"):
        OP_DistFile = MAFFT_OP_DIST_0_25
        EP_DistFile = MAFFT_EP_DIST_0_25

    OutEP = os.path.join(config.WorkingDir, "SampledEPVals.log")
    OutOP = os.path.join(config.WorkingDir, "SampledOPVals.log")
    op_vals_arr_ref = ""
    ep_vals_arr_ref = ""

    if config.PROGRAM == "GUIDANCE3":
        if config.GapPenDist.upper() == "EMP":
            log_file.write(
                f"Sample op according to empiric distribution: Guidance::SampelFromEmpiricDistribution({OP_DistFile},{OutOP},{config.Bootstraps})\n")
            op_vals_arr_ref = sample_from_empirical_distribution(OP_DistFile, OutOP, config.Bootstraps)
        elif config.GapPenDist.upper() == "UNIF":
            if config.MSA_Program == "MAFFT":
                log_file.write(
                    f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(1,3,{OutOP},{config.Bootstraps})\n")
                # according to mafft web-site defaults: http://mafft.cbrc.jp/alignment/server/index.html
                op_vals_arr_ref = sample_from_uniform_dist(1, 3, OutOP, config.Bootstraps)
            elif config.MSA_Program == "PRANK":
                log_file.write(
                    f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(0,0.5,{OutOP},{config.Bootstraps})\n")
                op_vals_arr_ref = sample_from_uniform_dist(0, 0.5, OutOP, config.Bootstraps)  # for prank v.140110 the defaults are: dna 0.025 / prot 0.005
            elif config.MSA_Program == "CLUSTALO":
                log_file.write(
                    f"Sample gap opening panelty according to uniform distribution: Guidance::SampleFromUniformDist(4,16,{OutOP},{config.Bootstraps})\n")
                op_vals_arr_ref = sample_from_uniform_dist(4, 16, OutOP, config.Bootstraps)

    if config.PROGRAM == "GUIDANCE3_HOT":
        if config.GapPenDist.upper() == "EMP":
            log_file.write(
                f"Sample op according to empiric distribution: Guidance::SampleFromEmpiricDistribution({OP_DistFile},{OutOP},{config.Bootstraps})\n")
            op_vals_arr_ref = sample_from_empirical_distribution(OP_DistFile, OutOP, config.Bootstraps)
            log_file.write(
                f"Sample ep according to empiric distribution: Guidance::SampleFromEmpiricDistribution({EP_DistFile},{OutEP},{config.Bootstraps})\n")
            ep_vals_arr_ref = sample_from_empirical_distribution(EP_DistFile, OutEP, config.FORM['Bootstraps'])
        if config.GapPenDist.upper() == "UNIF":
            log_file.write(
                f"Sample op according to uniform distribution: Guidance::SampleFromUniformDist(0,6,{OutOP},{config.Bootstraps})\n")
            op_vals_arr_ref = sample_from_uniform_dist(0, 6, OutOP, config.Bootstraps)
            log_file.write(
                f"Sample ep according to uniform distribution: Guidance::SampleFromUnifomDist(0,4,{OutEP},{config.Bootstraps})\n")
            ep_vals_arr_ref = sample_from_uniform_dist(0, 4, OutEP, config.Bootstraps)

    if config.isServer == 1:
        config.status_file = config.WorkingDir + "MSA_STATUS.txt"
        update_progress(f"{config.WorkingDir}{config.progress_report}", "Started generating alternative alignments")

        with open(f"{config.server_output}", "a") as OUTPUT:
            OUTPUT.write(
                "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t$fil = fopen('MSA_STATUS.txt', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

    # Get random branches for HoT
    NumOfBranches = (2 * config.MSA_Depth) - 3
    RandomBranches = [random.randint(0, NumOfBranches - 1) for _ in range(config.Bootstraps)]

    # CREATE THE PERTURBED ALN DIR
    os.mkdir(config.Scoring_Alignments_Dir)

    #Run HoT to create alternative MSAs
    countTrees = 0
    epsilon = 0.0006
    manager = Manager()
    lock = manager.Lock()
    config.mean_res_pair_score = manager.list()
    config.mean_col_score = manager.list()
    config.count_convergence = manager.Value('i', 0)
    # this config.convergence value might be unnecessary at the end and can be removed
    config.convergence = config.Bootstraps * Num_of_Aln_from_HoT_per_Run

    # Running parallel processes using multiprocessing, each will run an equal share of the BP alignments (?)
    processes = [Process(target=run_hot_process_on_tree, args=(config, epsilon, proc, RandomBranches, op_vals_arr_ref, ep_vals_arr_ref, Num_of_Aln_from_HoT_per_Run, lock)) for proc in range(config.proc_num)]
    for process in processes:
        process.start()

    # Use a generous timeout to catch truly hung processes (e.g. deadlock) without
    # killing legitimate long runs on large datasets. Configurable via config.hot_timeout.
    timeout_per_worker = getattr(config, 'hot_timeout', 48 * 3600)  # default: 48 hours
    for process in processes:
        process.join(timeout=timeout_per_worker)
        if process.is_alive():
            print(f"[WARNING] Worker process {process.pid} timed out after {timeout_per_worker}s, terminating")
            log_file.write(f"[WARNING] Worker process {process.pid} timed out after {timeout_per_worker}s, terminating\n")
            process.terminate()
            process.join(timeout=10)
            if process.is_alive():
                print(f"[WARNING] Worker process {process.pid} did not terminate, killing")
                log_file.write(f"[WARNING] Worker process {process.pid} did not terminate, killing\n")
                process.kill()
                process.join()

    # this is the number of alternative MSAs produced at the end
    alt_msas = len(os.listdir(config.Scoring_Alignments_Dir))
    config.convergence = alt_msas
    log_file.write(f"run_HOT_COS_GUIDANCE2 converged at tree #{alt_msas}\n")
    log_file.close()


    # Ksenia removed this part
    if config.isServer == 1:
        update_progress(f"{config.WorkingDir}{config.progress_report}",
                        f"Finished generating {alt_msas} alternative alignments")

    # To validate all alns were created
    # aln_count = len(os.listdir(config.Scoring_Alignments_Dir))
    # expected_count = Num_of_Aln_from_HoT_per_Run * config.Bootstraps
    # expected_count = (config.convergence) * Num_of_Aln_from_HoT_per_Run
    print(f"the convergence final number is {config.convergence}")
    # if aln_count < expected_count:
    #     exit_on_error("sys_error",
    #                   f"run_Guidance2: Only {aln_count} alignments were created on {config.Scoring_Alignments_Dir} while expecting {expected_count}\n", config)
    # else:
    #     print("\nSUCCESS!\n")
    # print(config.mean_res_pair_score)
    # print(config.mean_col_score)
    print("\nSUCCESS!\n")

#@timeit
def run_hot(config):
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
        with open(config.OutLogFile, "a") as log_file:
            if config.MSA_Program in ["MAFFT", "MAFFT_LINSI"]:
                config.HoT_MSA_Program = "MFT"
                config.HoT_MSA_Program_path = config.mafft_prog
                check_mafft_profile = subprocess.getoutput(f"which {config.HoT_MSA_Program_path} -profile")
                if "Command not found" in check_mafft_profile:
                    raise ValueError(
                        f"It seems that {config.HoT_MSA_Program_path}-profile is not properly installed or found in PATH config. Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n")

            elif config.MSA_Program == "MUSCLE":
                exit_on_error('user_error', "HoT currently does not support MUSCLE, please run GUIDANCE<br>", config)

            elif config.MSA_Program == "PAGAN":
                exit_on_error('user_error', "HoT currently does not support PAGAN, please run GUIDANCE<br>",config)

            elif config.MSA_Program == "CLUSTALO":
                config.HoT_MSA_Program = "CLO"
                config.HoT_MSA_Program_path = config.clustalw_prog

            elif config.MSA_Program == "PRANK":
                config.HoT_MSA_Program = "PRK"
                config.HoT_MSA_Program_path = config.prank_prog

            cmd = f"cd {config.WorkingDir}; {HOT_PROGRAM} {config.dataset} {config.HoT_MSA_Program}"

            if config.Seq_Type in ["AminoAcids", "Codons"]:
                cmd += " aa"
            elif config.Seq_Type == "Nucleotides":
                cmd += " nt"

            log_file.write(f"convert_fs_to_upper_case({config.WorkingDir}{config.codded_seq_fileName})")
            convert_fs_to_upper_case(
                f"{config.WorkingDir}{config.codded_seq_fileName}")  # HoT assumes that all the sequences are upper case

            if config.align_param == "":
                cmd += f" {config.codded_seq_fileName} . {config.WorkingDir}MSA_STATUS.txt 0 {config.HoT_MSA_Program_path} \"\" all > COS.std"
            else:
                cmd += f" {config.codded_seq_fileName} . {config.WorkingDir}MSA_STATUS.txt 0 {config.HoT_MSA_Program_path} \"\" all --- {config.align_param} > COS.std"

            log_file.write(f"run_HoT: {cmd}\n")

            if config.isServer == 1:
                status_file = f"{config.WorkingDir}MSA_STATUS.txt"
                with open(status_file, "w") as STATUS:
                    STATUS.write("\n")

                with open(f"{status_file}.0", "w") as STATUS0:
                    STATUS0.write("\n")

                with open(f"{config.server_output}", "a") as OUTPUT:
                    OUTPUT.write(
                    "<?php\n\tif (file_exists('MSA_STATUS.txt.0'))\n\t{\n\t\t$fil =fopen('MSA_STATUS.txt.0', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt.0'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

                    OUTPUT.write(
                    "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t$fil = fopen('MSA_STATUS.txt', r);\n\t\t$dat = fread($fil, filesize('MSA_STATUS.txt'));\n\t\techo \"$dat\";\n\tfclose($fil);\n\t}\n?>\n")

            os.system(cmd)

            if not os.path.exists(f"{config.WorkingDir}{config.Alignment_File}"):
                cmd = f"cp {config.WorkingDir}{config.dataset}_cos_{config.HoT_MSA_Program}/hot_H.fasta {config.WorkingDir}{config.Alignment_File}"
                log_file.write(f"run_HoT: {cmd}\n")  # Copy Alignment
                os.system(cmd)

            if config.Align_Order == "as_input":
                print(f"MSA_parser::sort_alignment({config.WorkingDir}{config.Alignment_File},fasta);\n")
                ans = sort_alignment(f"{config.WorkingDir}{config.Alignment_File}", "fasta")
                print("".join(ans))
                config.Alignment_File_NOT_SORTED = config.Alignment_File
                config.Alignment_File = config.Alignment_File + ".Sorted"

            os.system(f"mkdir {config.WorkingDir}{config.HoT_MSAs_Dir}")

            if config.NumOfSeq > 2:
                cmd = f"cp {config.WorkingDir}{config.dataset}_cos_{config.HoT_MSA_Program}/b[01]*.fasta {config.WorkingDir}{config.HoT_MSAs_Dir}"
            else:
                cmd = f"cp {config.WorkingDir}{config.dataset}_cos_{config.HoT_MSA_Program}/hot*.fasta {config.WorkingDir}{config.HoT_MSAs_Dir}"

            log_file.write(f"run_HoT: {cmd}\n")
            os.system(cmd)
            # HoT assumes that all the sequences are upper case
            if os.path.exists(f"{config.WorkingDir}{config.Alignment_File}.WithCodesName"):
                convert_fs_to_upper_case(f"{config.WorkingDir}{config.Alignment_File}.WithCodesName")
            # if os.path.exists(f"{config.WorkingDir}{config.Alignment_File}"):
            #     convert_fs_to_upper_case(f"{config.WorkingDir}{config.Alignment_File}")
            config.Scoring_Alignments_Dir = f"{config.WorkingDir}{config.HoT_MSAs_Dir}"
    except Exception as e:
        sys.exit("run_hot() Error: " + str(e) + "\n")

def prepare_rerun_parameters(config):
    mafft_max_iterates = {
        0: '0',
        1: '1',
        2: '2',
        5: '3',
        10: '4',
        20: '5',
        50: '6',
        80: '7',
        100: '8',
        1000: '9',
    }

    with open(f"{config.WorkingDir}/rerun_param", "w") as param:
        param.write("<?php\n")

        if config.PROGRAM == "GUIDANCE":
            param.write("$PROGRAM=1;\n")
            if config.MSA_Program == "MAFFT":
                param.write("$MSA_Prog=0;\n")
            elif config.MSA_Program == "PRANK":
                param.write("$MSA_Prog=1;\n")
            elif config.MSA_Program == "CLUSTALW":
                param.write("$MSA_Prog=2;\n")
            elif config.MSA_Program == "MUSCLE":
                param.write("$MSA_Prog=3;\n")
            elif config.MSA_Program == "PAGAN":
                param.write("$MSA_Prog=4;\n")

        elif config.PROGRAM == "HoT":
            param.write("$PROGRAM=2;\n")
            if config.MSA_Program == "MAFFT":
                param.write("$MSA_Prog=0;\n")
            elif config.MSA_Program == "PRANK":
                param.write("$MSA_Prog=1;\n")
            elif config.MSA_Program == "CLUSTALW":
                param.write("$MSA_Prog=2;\n")
            elif config.MSA_Program == "CLUSTALO":
                param.write("$MSA_Prog=2;\n")

        elif config.PROGRAM == "GUIDANCE3":
            param.write("$PROGRAM=0;\n")
            if config.MSA_Program == "MAFFT":
                param.write("$MSA_Prog=0;\n")
            elif config.MSA_Program == "PRANK":
                param.write("$MSA_Prog=1;\n")
            elif config.MSA_Program == "CLUSTALW":
                param.write("$MSA_Prog=2;\n")
            elif config.MSA_Program == "CLUSTALO":
                param.write("$MSA_Prog=2;\n")

        if config.MSA_Program == "MAFFT":
            if config.MAFFT_maxiterate == "":
                param.write("$MAFFT_MAX_ITERATES=0;\n")
            else:
                param.write(f"$MAFFT_MAX_ITERATES={mafft_max_iterates[int(config.MAFFT_maxiterate)]};\n")

            if config.MAFFT_refinement == "":
                param.write("$MAFFT_REFINE=0;\n")
            elif config.MAFFT_refinement == "localpair":
                param.write("$MAFFT_REFINE=1;\n")
            elif config.MAFFT_refinement == "genafpair":
                param.write("$MAFFT_REFINE=2;\n")
            elif config.MAFFT_refinement == "globalpair":
                param.write("$MAFFT_REFINE=3;\n")
        else:
            param.write("$MAFFT_MAX_ITERATES=0;\n")
            param.write("$MAFFT_REFINE=0;\n")

        if config.MSA_Program == "PRANK":
            if config.PRANK_F == "+F":
                param.write("$PRANK_INSERTION=0;\n")
            elif config.PRANK_F == "-F":
                param.write("$PRANK_INSERTION=1;\n")
        else:
            param.write("$PRANK_INSERTION=0;\n")

        param.write(f"$Bootstraps={config.Bootstraps};\n")
        param.write(f"$SP_COL_CUTOFF={config.SP_COL_CUTOFF};\n")
        param.write(f"$SP_SEQ_CUTOFF={config.SP_SEQ_CUTOFF};\n")

        if config.Align_Order == "as_input":
            param.write("$Align_Order=0;\n")
        elif config.Align_Order == "aligned":
            param.write("$Align_Order=1;\n")
        else:
            param.write("$Align_Order=1;\n")  # default

        if config.Seq_Type == "AminoAcids":
            param.write("$Seq_Type=0;\n")
            param.write("$CodonTable=0;\n")  # DEFAULT VALUE
        elif config.Seq_Type == "Nucleotides":
            param.write("$Seq_Type=1;\n")
            param.write("$CodonTable=0;\n")  # DEFAULT VALUE
        elif config.Seq_Type == "Codons":
            param.write("$Seq_Type=2;\n")
            if config.CodonTable == 1:
                param.write("$CodonTable=0;\n")
            elif config.CodonTable == 15:
                param.write("$CodonTable=1;\n")
            elif config.CodonTable == 6:
                param.write("$CodonTable=2;\n")
            elif config.CodonTable == 10:
                param.write("$CodonTable=3;\n")
            elif config.CodonTable == 2:
                param.write("$CodonTable=4;\n")
            elif config.CodonTable == 5:
                param.write("$CodonTable=5;\n")
            elif config.CodonTable == 3:
                param.write("$CodonTable=6;\n")
            elif config.CodonTable == 13:
                param.write("$CodonTable=7;\n")
            elif config.CodonTable == 9:
                param.write("$CodonTable=8;\n")
            elif config.CodonTable == 14:
                param.write("$CodonTable=9;\n")
            elif config.CodonTable == 4:
                param.write("$CodonTable=10;\n")

        param.write("?>\n")


# def send_finish_email_to_user(config):
#     # Set up logging
#     logging.basicConfig(filename=f"{config.WorkingDir}/log.txt", level=logging.INFO)
#
#     email_subject = ""
#     # http_path = f"http://guidance-dev.tau.ac.il/results/{vars['run_number']}"
#     base_http_path = "http://guidance-dev.tau.ac.il/"
#     http_path = base_http_path + f"/guidance/results/{config.run_number}"
#
#     if config.JOB_TITLE:
#         email_subject = f"Your Guidance results for {config.JOB_TITLE} are ready"
#     elif config.usrSeq_File:
#         email_subject = f"Your Guidance results for {config.usrSeq_File} are ready"
#     else:
#         email_subject = f"Your Guidance results for run number {config.run_number} are ready"
#
#     email_message = f"""Hello,
#
# The results for your Guidance run are ready at:
# {http_path}
#
# Running Parameters:
# Job Title: {config.JOB_TITLE}
# Sequences File: {config.usrSeq_File}
# MSA Algorithm: {config.MSA_Program}
# Number of Bootstraps: {config.Bootstraps}
# Scoring Method: {config.PROGRAM}
#
# Please note: the results will be kept on the server for three months.
#
# Thanks,
# GUIDANCE Team"""
#
#     # Set up email
#     msg = EmailMessage()
#     msg.set_content(email_message)
#     msg['Subject'] = email_subject
#     msg['From'] = formataddr(('GUIDANCE Team', 'admin@example.com'))  # Replace with ADMIN_EMAIL
#     msg['To'] = config.user_mail
#
#     # Send email
#     try:
#         with smtplib.SMTP(config.smtp_server) as server:
#             server.login(config.userName, config.userPass)
#             server.send_message(msg)
#             logging.info("Email sent successfully.")
#     except Exception as e:
#         logging.error(f"Failed to send email: {e}")
#
#     # Log the email message and command
#     log_msg = f"MESSAGE: {email_message}\n"
#     logging.info(log_msg)
#
#     # Run external command (if needed)
#     email_command = [
#         'perl', f'{CONST.BIN_DIR}/perl/sendEmail.pl',
#         '-f', 'admin@example.com',  # Replace with GENERAL_CONSTANTS::ADMIN_EMAIL
#         '-t', config.user_mail,
#         '-u', email_subject,
#         '-xu', config.userName,
#         '-xp', config.userPass,
#         '-s', config.smtp_server,
#         '-m', email_message
#     ]
#
#     result = run(email_command, stdout=PIPE, stderr=PIPE, text=True)
#     if "successfully" not in result.stdout:
#         logging.error(f"send_mail: The message was not sent successfully. System returned: {result.stdout}")


