import glob
import os.path
import re
import shutil
import sys
import os
import subprocess
from Bio import AlignIO, Phylo

from guidance_common_functions import print_message_to_output, exit_on_error, update_progress
from time_decorator import timeit

Bin = os.path.dirname(sys.argv[0])

# NEWIC2MAFFT = os.path.join(Bin, 'exec', 'newick2mafft.rb')
# newick2mafft = os.path.join(Bin, 'exec', 'newick2mafft.rb')
MSA_SET_SCORE = os.path.join(Bin, 'programs', 'msa_set_score', 'msa_set_score')
HOT_PROGRAM = os.path.join(Bin, 'hot_cos_main.py')
MAFFT_OP_DIST = os.path.join(Bin, 'balibase.mafft_7123_mafft.op.Dist20bins.txt')
MAFFT_OP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
MAFFT_EP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')
HOT_GUIDANCE2_PROGRAM = os.path.join(Bin, 'hot_cos_main.py')
MIDPOINT_ROOTING_R = os.path.join(Bin, 'programs', 'MidPoint_Rooting.R')

# MAFFT_OP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
# MAFFT_EP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')


# MSA_Score_CSS = "http://guidance.tau.ac.il/MSA_Colored.NEW.css"
MSA_Score_CSS = "https://taux.evolseq.net/guidance/static/css/MSA_Colored.NEW.EM.css"
MidPoint_Rooting_R = os.path.join(Bin, 'programs', 'MidPoint_Rooting.R')
# phylonet_prog = os.path.join(Bin, 'exec', 'phylonet_v1_7', 'phylonet_v1_7.jar')
isEqualTopologyProg = os.path.join(Bin, 'programs', 'isEqualTree', 'isEqualTree')

@timeit
def calculate_msa_depth(inMSA, args_library):
    try:
        with open(inMSA, "r") as inMSA_file:
            aln = AlignIO.read(handle = inMSA_file, format = "fasta")
            msa_depth = len(aln)
    except Exception as e:
        exit_on_error("sys_error", f"MSA_Depth: Can't read MSA: '{inMSA}' - {e}\n", args_library)


    # JS ? - print for flask version
    msa_depth_file = os.path.join(args_library.WorkingDir, 'MSA_DEPTH')
    with open(msa_depth_file, "w") as length_file:
        length_file.write(f"{msa_depth}\n")

    return msa_depth


# NEED TO UPDATE FOR IQTREE
@timeit
def pull_out_bp_trees_bbl(no_bp_dir, dataset, bp_repeats, aln_prog):
    if not no_bp_dir.endswith("/"):
        no_bp_dir += "/"

    semphy_log_file = f"{no_bp_dir}BP/{dataset}.{aln_prog}.semphy.out"
    bp_dir = f"{no_bp_dir}BP/"

    if os.path.exists(f"{no_bp_dir}{dataset}.{aln_prog}.semphy.tree"):
        os.remove(f"{no_bp_dir}{dataset}.{aln_prog}.semphy.tree")

    if not os.path.exists(bp_dir):
        os.mkdir(bp_dir)

    non_unique_trees_dir = ""
    if aln_prog != "MAFFT":
        non_unique_trees_dir = f"{bp_dir}nonUniqueTrees/"
        if not os.path.exists(non_unique_trees_dir):
            os.mkdir(non_unique_trees_dir)

    make_unique_trees = "yes" if aln_prog != "MAFFT" else "no"

    tree_line = ""
    read_reconstructed_tree = False

    with open(semphy_log_file, 'r') as log_file:
        count_trees = 0
        count_unique_trees = 0
        num_repeats = []

        if aln_prog != "MAFFT":
            make_unique_trees = "yes"
        else:
            make_unique_trees = "no"

        for line in log_file:
            if line.startswith("# Finished tree reconstruction."):
                _ = next(log_file)
                _ = next(log_file)
                _ = next(log_file)
                tree_line = next(log_file)
                tree_file = f"{no_bp_dir}{dataset}.{aln_prog}.semphy.tree"

                with open(tree_file, 'w') as out_file:
                    out_file.write(tree_line)

                tree_line = ""
                read_reconstructed_tree = True
            elif ((" # Tree after BBL." in line or "The reconsructed tree:" in line) and read_reconstructed_tree):
                _ = next(log_file) if " # Tree after BBL." in line else None
                tree_line = next(log_file)

                if make_unique_trees == "no":
                    tree_dir = f"{bp_dir}/tree_{count_trees}/"
                else:
                    tree_dir = f"{non_unique_trees_dir}/tree_{count_trees}/"

                if not os.path.exists(tree_dir):
                    os.mkdir(tree_dir)

                tree_file = f"{tree_dir}{dataset}.{aln_prog}.semphy.tree_{count_trees}"
                count_trees += 1

                with open(tree_file, 'w') as out_file:
                    out_file.write(tree_line)

                if make_unique_trees == "yes":
                    for i in range(count_unique_trees):
                        unique_tree_file = f"{bp_dir}tree_{i}/{dataset}.{aln_prog}.semphy.tree_{i}"
                        is_equal_topology_res_file = f"{tree_dir}isEqualTopology.{i}.std"
                        is_equal_topology_command = f"{isEqualTopologyProg} {tree_file} {unique_tree_file}"
                        is_equal_topology = os.system(is_equal_topology_command)

                        with open(is_equal_topology_res_file, 'w') as out_equal_top_file:
                            out_equal_top_file.write(str(is_equal_topology))

                        if is_equal_topology == 1:
                            num_repeats[i] += 1
                            break
                        elif is_equal_topology == 2:
                            print(f"Skipping ERROR in isEqualTopology of {tree_file} and {unique_tree_file}")
                            continue

                    else:
                        num_repeats.append(1)
                        unique_tree_dir = f"{bp_dir}tree_{count_unique_trees}/"
                        if not os.path.exists(unique_tree_dir):
                            os.mkdir(unique_tree_dir)

                        unique_tree_file = f"{unique_tree_dir}{dataset}.{aln_prog}.semphy.tree_{count_unique_trees}"
                        shutil.copy(tree_file, unique_tree_file)
                        count_unique_trees += 1

        if count_trees != bp_repeats:
            return f"ERROR: dataset: {dataset} \t count_trees: {count_trees} while it should be {bp_repeats}\n"

        if make_unique_trees == "yes":
            num_repeats_file = f"{bp_dir}numRepeats"
            with open(num_repeats_file, 'w') as out_num_repeats:
                out_num_repeats.write(" ".join(map(str, num_repeats)))

            return "ok", count_unique_trees, num_repeats
        else:
            return ["ok"]

@timeit
def pull_out_bp_trees(no_bp_dir, dataset, bp_repeats, aln_prog, args_library):
    ####################################################################################################################
    # pull out all the BP trees into the BP directory
    # pull out the original tree (that was done on the complete MSA file)
    ####################################################################################################################
    make_unique = ""

    if not no_bp_dir.endswith("/"):
        no_bp_dir += "/"

    bp_dir = f"{no_bp_dir}BP/"
    if not os.path.exists(bp_dir):
        os.mkdir(bp_dir)

    non_unique_trees_dir = ""
    if aln_prog != "MAFFT":  # BUILT ALIGNMENT ONLY FOR UNIQUE TREES
        non_unique_trees_dir = f"{bp_dir}nonUniqueTrees/"
        if not os.path.exists(non_unique_trees_dir):
            os.mkdir(non_unique_trees_dir)
        make_unique = "yes"
    else:
        make_unique = "no"

    iqtree_boottrees_file = f"{bp_dir}{args_library.Alignment_File}.boottrees"
    # iqtree_boottrees_file = f"{bp_dir}{dataset}.{aln_prog}.aln.boottrees"
    print(f"iqtree boottrees file: {iqtree_boottrees_file}\n")

    with (open(iqtree_boottrees_file, 'r') as boottrees_file):

        count = 0
        count_unqique = 0
        num_repeats = []

        for my_tree in boottrees_file:
            if aln_prog == "MAFFT":
                tree_dir = f"{bp_dir}/tree_{count}/"
            else:  # CHECK UNIQUE TREE ONLY NOT FOR MAFFT
                tree_dir = f"{non_unique_trees_dir}/tree_{count}/"

            if not os.path.exists(tree_dir):
                os.mkdir(tree_dir)

            tree_file = f"{tree_dir}{dataset}.{aln_prog}.iqtree.tree_{count}"

            try:
                with open(tree_file, 'w') as out_file:
                    out_file.write(my_tree)
            except Exception as e:
                return f"can't open file {tree_file}"
            count += 1

            if make_unique == "yes":
                i=0
                while i < count_unqique:
                    unique_tree_file = f"{bp_dir}/tree_{i}/{dataset}.{aln_prog}.iqtree.tree_{i}"
                    is_equal_topology_res_file = tree_dir + "isEqualTopology." + str(i) + ".std"
                    is_equal_topology_command = isEqualTopologyProg + " " + tree_file + " " + unique_tree_file
                    is_equal_topology = os.system(is_equal_topology_command)
                    # with open(is_equal_topology_res_file, 'w') as out_equal_top:
                    #     out_equal_top.write(str(is_equal_topology) + "\n")
                    if is_equal_topology == 1:
                        num_repeats[i] += 1
                        break
                    if is_equal_topology == 2:
                        print("skipping ERROR in isEqualTopology of", tree_file, "and", unique_tree_file)
                        continue
                    i += 1
                if i == count_unqique:
                    num_repeats.append(1)
                    unique_trees_dir = f"{bp_dir}/tree_{count_unqique}/"
                    if not os.path.exists(unique_trees_dir):
                        os.system("mkdir " + unique_trees_dir)
                    unique_tree_file = unique_trees_dir + dataset + "." + aln_prog + ".iqtree.tree_" + str(count_unqique)
                    os.system("cp " + tree_file + " " + unique_tree_file)
                    count_unqique += 1


    if count != bp_repeats:
        return f"ERROR: dataset: {dataset} \t count_trees: {count} while it should be {bp_repeats}\n"

    if aln_prog != "MAFFT":
        num_repeats_file = bp_dir + "/"+ "numRepeats"
        with open(num_repeats_file, 'w') as out_num_repeats:
            out_num_repeats.write(" ".join(map(str, num_repeats)))
        return "ok", count_unqique, num_repeats

    else:
        return ["ok"]

@timeit
def root_BP_trees(bsDir, dataset, orig_prog, bp_repeats, suffix=None, rooting_type="BioPerl"):
    ####################################################################################################################
    # Root all trees on BP dir
    ####################################################################################################################

    # Default values if not defined
    if suffix is None:
        suffix = ""
    if not bsDir.endswith("/"):
        bsDir += "/"

    for countTrees in range(bp_repeats):
        tree_file = f"{bsDir}tree_{countTrees}/{dataset}.{orig_prog}.iqtree.tree_{countTrees}{suffix}"

        # Check if the tree file exists
        if os.path.exists(tree_file):
            rooted_tree_file = f"{tree_file}.rooted"

            # Rooting based on the specified method (BioPerl or MidPoint)
            if rooting_type.upper() == "BIOPERL":
                root_tree(tree_file, rooted_tree_file)
            elif rooting_type.upper() == "MIDPOINT":
                subprocess.run(["R", "--slave", "--no-save", "--no-restore", "--no-environ", "--silent",
                                "--args", tree_file, rooted_tree_file, "<", MidPoint_Rooting_R], shell=True)

            # Reading the rooted tree and processing it
            with open(rooted_tree_file, "r") as infile:
                newick = infile.read()

            # Additional processing for branch lengths
            if ":-" in newick:
                rooted_tree_file_with_minus_lengths = f"{rooted_tree_file}.withMinusLengths"
                shutil.move(rooted_tree_file, rooted_tree_file_with_minus_lengths)
                newick = newick.replace(":-", ":")

                with open(rooted_tree_file, "w") as outfile:
                    outfile.write(newick)

            if re.search(r":\d+[^\.\d+]", newick):
                rooted_tree_file_bad_branch_length = f"{rooted_tree_file}.badBranchLength"
                shutil.move(rooted_tree_file, rooted_tree_file_bad_branch_length)
                newick = re.sub(r":(\d+)([^\.\d+])", r":\1.0\2", newick)

                with open(rooted_tree_file, "w") as outfile:
                    outfile.write(newick)
        else:
            print(f"File does not exist: {tree_file}")

    return ["ok"]

@timeit
def root_tree(in_tree, out_tree):
    """
    The input tree (in_tree) must be an unrooted tree, i.e., the root node has at least 3 sons.
    The output tree (out_tree) will have the root with 2 sons, and all direct sons of the root will be made bifurcating.
    The rest of the tree is left untouched.

    Parameters:
    - in_tree (str): Input tree file in Newick format.
    - out_tree (str): Output tree file in Newick format.
    """

    with open(in_tree, 'r') as infile:
        trees = Phylo.parse(infile, 'newick')

        with open(out_tree, 'w') as outfile:
            for tree in trees:
                root = tree.clade
                sons = root.clades.copy()

                # Remove edges between root-sons
                for son in sons:
                    root.clades.remove(son)

                # Iteratively add
                curr_father = root
                while len(sons) > 2:
                    son = sons.pop(0)
                    curr_father.clades.append(son)
                    mid_node = Phylo.BaseTree.Clade()
                    curr_father.clades.append(mid_node)
                    mid_node.branch_length = 0
                    # if sons[0].branch_length == 0:
                    #     sons[0].branch_length = 1e-6
                    curr_father = mid_node

                curr_father.clades.append(sons[0])
                curr_father.clades.append(sons[1])
                sons[0].branch_length = 0
                sons[1].branch_length = 0
                # tree.rooted = True
                # root = tree.clade
                # root.branch_length = None
                Phylo.write(tree, outfile, 'newick', format_branch_length='%0.6f')

    with open(out_tree, 'r') as infile:
        newstr = infile.read()
        # Replace :0.000000; at the end with ;
        newstr = re.sub(r':0\.000000;', ';', newstr)

    with open(out_tree, 'w') as outfile:
        outfile.write(newstr)

@timeit
def reformat_trees_branch_length(in_tree, out_tree):
    """
    Reformat tree branch lengths and write the tree in newick format.

    Parameters:
    - in_tree (str): Input tree file in newick format.
    - out_tree (str): Output tree file in newick format.

    Returns:
    - None
    """

    # Read input tree
    with open(in_tree, 'r') as infile:
        # input_tree = Phylo.read(infile, 'newick')
        trees = list(Phylo.parse(infile, 'newick'))

    for input_tree in trees:
        # Iterate through nodes and update branch lengths
        input_tree.clade.branch_length = None
        for node in input_tree.find_clades():
            if node.branch_length is not None:
                node.branch_length = float("{:.6f}".format(node.branch_length))

        # Write output tree
        with open(out_tree, 'w') as outfile:
            Phylo.write([input_tree], outfile, 'newick', format_branch_length='%0.6f')

    with open(out_tree, 'r') as infile:
        newstr = infile.read()
        # Replace :0.000000; at the end with ;
        newstr = re.sub(r':0\.000000;', ';', newstr)

    with open(out_tree, 'w') as outfile:
        outfile.write(newstr)

@timeit
def fix_mafft_rough_tree(tree_file):
    """
    Fix MAFFT RoughTree by removing underscores and extra symbols in node labels.

    Parameters:
    - tree_file (str): Input tree file in newick format.

    Returns:
    - None
    """
    tree_orig = tree_file + ".orig"
    shutil.copy(tree_file, tree_orig)

    with open(tree_file, 'r') as tree_file_handle:
        tree = tree_file_handle.readline().rstrip()

    if not tree.endswith(';'):
        tree += ';'

    # Remove underscores and extra symbols in node labels
    tree = tree.replace('_', '').replace(';', '')

    with open(tree_file, 'w') as tree_file_handle:
        tree_file_handle.write(f'{tree}\n')

@timeit
def Bootstrap_Trees(args_library):
    # ---------------------------------------------
    if args_library.isServer == 1:
        # print_message_to_output("Constructing bootstrap guide-trees", args_library)
        update_progress(f"{args_library.WorkingDir}{args_library.progress_report}", "Constructing bootstrap guide-trees")

    os.makedirs(args_library.BootStrap_Dir)
    # args_library.Tree_File = f"{args_library.dataset}.{args_library.MSA_Program}.aln.treefile"
    # args_library.Iqtree_LogFile = f"{args_library.dataset}.{args_library.MSA_Program}.aln.log"
    # args_library.Iqtree_Boottrees = f"{args_library.dataset}.{args_library.MSA_Program}.aln.boottrees"
    args_library.Tree_File = f"{args_library.Alignment_File}.treefile"
    args_library.Iqtree_LogFile = f"{args_library.Alignment_File}.log"
    args_library.Iqtree_Boottrees = f"{args_library.Alignment_File}.boottrees"

    cmd = ""
    msa_depth = calculate_msa_depth(f"{args_library.WorkingDir}{args_library.Alignment_File}", args_library)
    verbose_level = 8

    if args_library.BBL.upper() == "YES":
        args_library.semphy_prog = f"{args_library.semphy_prog} -n "
        verbose_level = 1

    if args_library.Seq_Type in ["AminoAcids", "Codons"]:
        if msa_depth > 150:  # use JC for distance estimation
            cmd = (f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m JTT -bo {args_library.Bootstraps} -nt {args_library.proc_num} -v -st AA -n 0 -fast")
            print(cmd + "\n")

        else:  # use JTT for distance estimation
            # cmd = (
            #     f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m JTT -bo {args_library.Bootstraps} -nt {args_library.proc_num} -st AA -n 0 -fast")
            cmd = (
                f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m JTT -bo {args_library.Bootstraps} -st AA -n 0 -fast")
            print(cmd + "\n")

    elif args_library.Seq_Type == "Nucleotides":
        if msa_depth > 150:  # use JC for distance estimation
            cmd = (
                f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m JC -bo {args_library.Bootstraps} -nt {args_library.proc_num} -st DNA -n 0 -fast")
        else:  # use HKY for distance estimation
            # cmd = (
            #     f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m HKY -bo {args_library.Bootstraps} -nt {args_library.proc_num} -st DNA -n 0 -fast")
            cmd = (
                f"{args_library.iqtree_prog} -s {args_library.WorkingDir}{args_library.Alignment_File} -m HKY -bo {args_library.Bootstraps} -st DNA -n 0 -fast")

    with open(f'{args_library.OutLogFile}', "a") as log_file:
        log_file.write(f"Bootstrap_Trees: {cmd}\n")
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # tree_files = os.path.join(f"{args_library.WorkingDir}", f"{args_library.dataset}.{args_library.MSA_Program}.aln.*")
    tree_files = os.path.join(f"{args_library.WorkingDir}", f"{args_library.Alignment_File}.*")
    for file in glob.glob(tree_files):
        if not file.endswith("ORIG"):
            shutil.move(file, f'{args_library.BootStrap_Dir}')

    if os.path.getsize(f"{args_library.BootStrap_Dir}{args_library.Iqtree_Boottrees}") == 0 or not os.path.exists(
            f"{args_library.BootStrap_Dir}{args_library.Iqtree_Boottrees}"):
        exit_on_error("sys_error",
                      f"Bootstrap_Trees: '{args_library.BootStrap_Dir}{args_library.Iqtree_Boottrees}' is empty/does not exist\n", args_library)

