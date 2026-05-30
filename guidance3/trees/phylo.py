import glob
import os.path
import re
import shutil
import sys
import os
import subprocess
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Clade
from io import StringIO
import random
from guidance3.utils.common import print_message_to_output, exit_on_error, update_progress
from guidance3.utils.timing import timeit
from guidance3.constants import (
    MSA_SET_SCORE, MAFFT_OP_DIST, MAFFT_OP_DIST_0_25, MAFFT_EP_DIST_0_25,
    MIDPOINT_ROOTING_R, HOT_GUIDANCE3_PROGRAM,
)

HOT_PROGRAM = HOT_GUIDANCE3_PROGRAM

from guidance3.constants import MidPoint_Rooting_R, isEqualTopologyProg


#@timeit
def calculate_msa_depth(inMSA, config):
    """Calculates depth of MSA.

      Args:
        inMSA: path to the MSA file
        config: the config object which stores all relevant paths, script arguments and constants

      Returns:
        Depth value of the MSA
      """
    try:
        with open(inMSA, "r") as inMSA_file:
            aln = AlignIO.read(handle = inMSA_file, format = "fasta")
            msa_depth = len(aln)
    except Exception as e:
        exit_on_error("sys_error", f"MSA_Depth: Can't read MSA: '{inMSA}' - {e}\n", config)


    # JS ? - print for flask version
    msa_depth_file = os.path.join(config.WorkingDir, 'MSA_DEPTH')
    with open(msa_depth_file, "w") as length_file:
        length_file.write(f"{msa_depth}\n")

    return msa_depth


# NEED TO UPDATE FOR IQTREE
#@timeit
def pull_out_bp_trees_bbl(no_bp_dir, dataset, bp_repeats, aln_prog):
    """Pulls out all the Bootstrap trees (after Binary Branch Lengths (BBL) optimization) into the BP directory.
    Pulls out the original tree (that was done on the complete MSA file)


          Args:
            no_bp_dir: the folder in which BP dir is located, this is a working directory
            dataset: name of the dataset, default name is "MSA"
            bp_repeats: number of bootstrap trees produces
            aln_prog: multiple sequence alignment program

          Returns:
            ["ok"] in case of success when align program is MAFFT
            or a tuple "ok", count_unique_trees, num_repeats in case when align program is not MAFFT
            an error in case if number of trees produced is not equal to a number of bootstrap trees (parameter) requested
          """

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

#@timeit
def pull_out_bp_trees(no_bp_dir, dataset, bp_repeats, aln_prog, config):
    """Pulls out all the Bootstrap trees into the BP directory.
    Pulls out the original tree (that was done on the complete MSA file)

          Args:
            no_bp_dir: the folder in which BP dir is located, this is a working directory
            dataset: name of the dataset, default name is "MSA"
            bp_repeats: number of bootstrap trees produces
            aln_prog: multiple sequence alignement program
            config: the object with the program paths, arguments and constants

          Returns:
            ["ok"] in case of success when align program is MAFFT
            or a tuple "ok", count_unique_trees, num_repeats in case when align program is not MAFFT
            an error in case if number of trees produces is not equal to a number of bootstrap trees (parameter) requested
          """

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

    iqtree_boottrees_file = f"{bp_dir}{config.Alignment_File}.boottrees"
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
                i = 0
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

#@timeit
def root_BP_trees(bsDir, dataset, orig_prog, bp_repeats, suffix=None, rooting_type="BioPerl"):
    """Roots all the bootstrap trees in the BP directory."""

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
                # remove_node_support(tree_file)  # TODO testing this for FastTree
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

#@timeit
def root_tree(in_tree, out_tree):
    """
    The input tree (in_tree) must be an unrooted tree, i.e., the root node has at least 3 sons.
    The output tree (out_tree) will have the root with 2 sons, and all direct sons of the root will be made bifurcating.
    The rest of the tree is left untouched.

    Args:
    - in_tree (str): Input tree file in Newick format.
    - out_tree (str): Output tree file in Newick format.

    Returns:
    - None
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

def remove_node_support(tree_file, min_branch_length=1e-6): #TODO testing this for FastTree tree
    backup_file = f"{tree_file}.backup"
    shutil.copy(tree_file, backup_file)
    # Step 1: Read and clean support values
    with open(tree_file) as f:
        tree_str = f.read()

    # Remove internal node support values (e.g. )0.95:)
    tree_str_cleaned = re.sub(r'\)([0-9.eE\-]+):', r'):', tree_str)

    # Step 2: Load tree into Biopython for branch length fixing
    tree = Phylo.read(StringIO(tree_str_cleaned), "newick")

    # Set root branch length to None (avoids ":0.000000;" at the end)
    if tree.root.branch_length is not None:
        tree.root.branch_length = None

    # Step 3: Update small branch lengths
    for node in tree.find_clades():
        if node is tree.root:
            continue  # Skip root
        if node.branch_length is not None:
            node.branch_length = float("{:.6f}".format(max(node.branch_length, min_branch_length)))

    # Step 4: Output cleaned tree back to file
    with open(tree_file, "w") as f:
        Phylo.write(tree, f, "newick", format_branch_length='%0.6f')

    with open(tree_file) as f:
        tree_output = f.read()

    # Remove ":0.000;" at end if present
    tree_output = re.sub(r':0\.0+;$', ';', tree_output)

    with open(tree_file, "w") as f:
        f.write(tree_output)

#@timeit
def reformat_trees_branch_length(in_tree, out_tree):
    """
    Reformat tree branch lengths and write the tree in newick format.

    Args:
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

#@timeit
def fix_mafft_rough_tree(tree_file):
    """
    Fix MAFFT RoughTree by removing underscores and extra symbols in node labels.

    Args:
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

#@timeit
def Bootstrap_Trees(config):
    """
    Runs IQ-Tree bootstrapping on the existing MSA file with -fast parameter for either nucleotides (JC or HKY model used), or amino acid sequences (JTT model).

    Args:
    - config: an object with program paths, arguments and constants.

    Returns:
    - None

    Produces 3 files with the following endings: .treefile, .log, .boottrees.
    .boottrees file holds all the bootstrap trees produced by IQ-Tree
    """
    if config.isServer == 1:
        update_progress(f"{config.WorkingDir}{config.progress_report}", "Constructing bootstrap guide-trees")

    os.makedirs(config.BootStrap_Dir, exist_ok=True) #TODO: change on the server , add exists_ok=True
    config.Tree_File = f"{config.Alignment_File}.treefile"
    config.Iqtree_LogFile = f"{config.Alignment_File}.log"
    config.Iqtree_Boottrees = f"{config.Alignment_File}.boottrees"

    cmd = ""
    msa_depth = calculate_msa_depth(f"{config.WorkingDir}{config.Alignment_File}", config)
    verbose_level = 8

    # ### START: TESTING WITH FASTTREE ###
    # fasttree = 1
    #
    # if fasttree == 1:
    #     #bootstrap alignment
    #
    #     # Load alignment
    #     alignment = AlignIO.read(f"{config.WorkingDir}{config.Alignment_File}", 'fasta')
    #     alignment_length = alignment.get_alignment_length()
    #     output_prefix = f"{config.WorkingDir}/bootstrap_msa"
    #
    #     with open(f"{config.WorkingDir}{config.Iqtree_Boottrees}", "w") as tree_out:
    #         for i in range(config.Bootstraps):
    #             indices = [random.randint(0, alignment_length - 1) for _ in range(alignment_length)]
    #             bootstrapped_records = []
    #             for record in alignment:
    #                 bootstrapped_seq = "".join(record.seq[j] for j in indices)
    #                 record.seq = record.seq.__class__(bootstrapped_seq)
    #                 bootstrapped_records.append(record)
    #
    #             bootstrap_alignment = MultipleSeqAlignment(bootstrapped_records)
    #
    #             tmp_align_file = f"{output_prefix}_{i + 1:03d}.fasta"
    #             AlignIO.write(bootstrap_alignment, tmp_align_file, "fasta")
    #             with open(f'{config.OutLogFile}', "a") as log_file:
    #                 log_file.write(f"FastTree Boostrapped alignment #{i}: {tmp_align_file} saved\n")
    #
    #             # for each bootstrapped alignment run fasttree
    #             cmd = f"FastTree {tmp_align_file}"
    #             with open(f'{config.OutLogFile}', "a") as log_file:
    #                 log_file.write(f"FastTree Bootstrap_Tree {i}: {cmd}\n")
    #             # subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    #             result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    #             tree = result.stdout.strip()
    #
    #             if not tree.endswith(";"):
    #                 tree += ";"
    #             tree_out.write(tree + "\n")
    #
    #         tree_files = os.path.join(f"{config.WorkingDir}", f"{config.Alignment_File}.*")
    #         for file in glob.glob(tree_files):
    #             if not file.endswith("ORIG"):
    #                 shutil.move(file, f'{config.BootStrap_Dir}')
    #
    #         bootstrap_alignment_files = os.path.join(f"{config.WorkingDir}", f"bootstrap_msa_*")
    #         for boot_msa_file in glob.glob(bootstrap_alignment_files):
    #             if os.path.isfile(boot_msa_file):
    #                 os.remove(boot_msa_file)
    #         with open(f'{config.OutLogFile}', "a") as log_file:
    #             log_file.write("All bootstrap_msa_* files deleted\n")
    #
    #         if os.path.getsize(
    #                 f"{config.BootStrap_Dir}{config.Iqtree_Boottrees}") == 0 or not os.path.exists(
    #                 f"{config.BootStrap_Dir}{config.Iqtree_Boottrees}"):
    #             exit_on_error("sys_error",
    #                           f"Bootstrap_Trees: '{config.BootStrap_Dir}{config.Iqtree_Boottrees}' is empty/does not exist\n",
    #                           config)
    #         return
    # else:
    #     return
    #
    # ### END : TESTING WITH FASTTREE ###
    #
    if config.BBL.upper() == "YES":
        config.semphy_prog = f"{config.semphy_prog} -n "
        verbose_level = 1

    if config.Seq_Type in ["AminoAcids", "Codons"]:
        if msa_depth > 150:  # use JC for distance estimation
            cmd = (f"{config.iqtree_prog} -s {config.WorkingDir}{config.Alignment_File} -m JTT -bo {config.Bootstraps} -nt {config.proc_num} -v -st AA -n 0 -fast")
            print(cmd + "\n")

        else:  # use JTT for distance estimation
            cmd = (
                f"{config.iqtree_prog} -s {config.WorkingDir}{config.Alignment_File} -m JTT -bo {config.Bootstraps} -st AA -n 0 -fast")
            print(cmd + "\n")

    elif config.Seq_Type == "Nucleotides":
        if msa_depth > 150:  # use JC for distance estimation
            cmd = (
                f"{config.iqtree_prog} -s {config.WorkingDir}{config.Alignment_File} -m JC -bo {config.Bootstraps} -nt {config.proc_num} -st DNA -n 0 -fast")
        else:  # use HKY for distance estimation
            cmd = (
                f"{config.iqtree_prog} -s {config.WorkingDir}{config.Alignment_File} -m HKY -bo {config.Bootstraps} -st DNA -n 0 -fast")

    with open(f'{config.OutLogFile}', "a") as log_file:
        log_file.write(f"Bootstrap_Trees: {cmd}\n")
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    tree_files = os.path.join(f"{config.WorkingDir}", f"{config.Alignment_File}.*")
    for file in glob.glob(tree_files):
        if not file.endswith("ORIG"):
            shutil.move(file, f'{config.BootStrap_Dir}')

    if os.path.getsize(f"{config.BootStrap_Dir}{config.Iqtree_Boottrees}") == 0 or not os.path.exists(
            f"{config.BootStrap_Dir}{config.Iqtree_Boottrees}"):
        exit_on_error("sys_error",
                      f"Bootstrap_Trees: '{config.BootStrap_Dir}{config.Iqtree_Boottrees}' is empty/does not exist\n", config)

