from .logger import log_print, write_to_file, cleanup, debug, print_usage, run_command_line, handle_termination_signal
from .sequencing_method import SequencingMethod
from .sequence import Sequence
from .tree import Tree_
from .file_handler import FileHandler
from .mafft import prepare_user_tree_mafft
import sys
import os
import shutil
import re
from time import time, localtime, strftime
from datetime import datetime


def main(argv=None):
    if argv is None:
        argv = sys.argv

    args = argv
    if len(args) < 4:
        print_usage()
        sys.exit()

    (idstr, method, seqtype_str, input_file, output_base_dir, status_file, tgz, msa_program_path, user_treefile,
     user_split_number, *param) = argv[1:]

    split_method_input = re.match(r'^(.{3})(.?)', method)
    method = split_method_input.group(1)
    run_only_hot = split_method_input.group(2)
    output_dir = f"{idstr}_cos_{method}"

    file_handler = FileHandler(input_file, output_base_dir, status_file, output_dir, tgz)
    handle_termination_signal(file_handler)
    sequencing_method = SequencingMethod(method, msa_program_path, param)
    sequence = Sequence(seqtype_str, file_handler.input_file)
    tree = Tree_(user_treefile, user_split_number)

    log_entry = (
        f"--- init {file_handler.current_script_file} pid={os.getpid()} {datetime.now()}\n"
        f"met={sequencing_method.name}\ninfile={file_handler.input_file}\n"
        f"seqtype={sequence.sequence_type}\noutdir={file_handler.output_base_dir}\n"
        f"tgz={file_handler.should_create_archive}\ndebug={debug}\n"
        f"basedir={file_handler.basedir}\noutput_id={file_handler.output_dir}\n---\n"
    )

    os.makedirs(file_handler.output_dir, exist_ok=True)

    if os.path.isfile(f'{tree.tree_file}'):
        shutil.copy(tree.tree_file, file_handler.output_dir)
        if sequencing_method.method_to_cmd_mapping[sequencing_method.name] == "MAF":
            prepare_user_tree_mafft(os.path.join(file_handler.output_dir, os.path.basename(tree.tree_file)))

    os.chdir(file_handler.output_dir)

    with open(file_handler.log_file, 'a') as log_file_handler:
        old_stdout = sys.stdout
        sys.stdout = log_file_handler
        try:
            log_print(0, 1, log_entry, file_handler)
        finally:
            sys.stdout = old_stdout
        log_file_handler.flush()

    with open(file_handler.time_file, 'a') as time_file_handler:
        time_file_handler.write(log_entry)
        time_file_handler.flush()

    run_command_line("find . -size 0 -delete", file_handler)

    sequence.construct_bidirectional_sequence_manager(file_handler)

    if sequence.bidirectional_sequences_manager['number_of_sequences'] < 2:
        print(f"\nERROR: {sequence.heads_sequence_file_name} has less than 2 sequences....\n")
        sys.exit()

    sequencing_method.met_init(sequence.sequence_type, file_handler)
    string = ''.join([str(item) for item in sequence.bidirectional_sequences_manager['fasta'][0]])
    write_to_file('input.fasta', string)

    if os.path.exists(tree.tree_file):
        log_print(1, 1, f"-{tree.tree_file} exists", file_handler)
    else:
        log_print(0, 1, "-Making guide tree ...\n", file_handler)
        if file_handler.status_file != "":
            with open(f"{file_handler.status_file}.0", "w") as status_file_handler:
                status_file_handler.write("<ul><li>Making guide tree</li></ul>\n")
        sequencing_method.make_guide_tree('input.fasta', tree.tree_file, sequence, file_handler)

    tree.tree2split(tree.tree_file, sequencing_method.command, sequence.bidirectional_sequences_manager['number_of_sequences'], file_handler)

    treefile = "in.dnd"
    write_to_file(treefile, tree.subtrees['tree'])

    outfile, msar = None, None

    for i in range(2):
        outfile = f"hot_{sequence.hot_direction[i].upper()}"
        if (msar := sequence.reverse_sequence(outfile + sequence.file_extensions[0], 1, file_handler)):
            log_print(1, 1, f"-{outfile}{sequence.file_extensions[0]} exists", file_handler)
            continue

        if i == 1 and (msar := sequence.reverse_sequence(outfile + sequence.file_extensions[1], 1, file_handler)):
            log_print(1, 1, f"-{outfile}{sequence.file_extensions[1]} exists", file_handler)
            write_to_file(f"{outfile}{sequence.file_extensions[0]}", msar)
            continue

        infile = f"in{sequence.file_extensions[i]}"
        write_to_file(infile, str(''.join(sequence.bidirectional_sequences_manager['fasta'][i])))
        log_print(0, 1, f"-Making {outfile} ...\n", file_handler)

        if file_handler.status_file != "":
            with open(f"{file_handler.status_file}.0", "a") as status_handler:
                status_handler.write(f"<ul><li>Making {outfile}</li></ul>\n")

        msar = sequencing_method.align_sequences(infile, treefile, f"{outfile}{sequence.file_extensions[i]}", sequence, file_handler)

        if i == 1:
            write_to_file(f"{outfile}{sequence.file_extensions[0]}", msar)

    if run_only_hot.lower().startswith('h'):
        cleanup(0, file_handler)
        exit(0)

    pfiles = [["", ""], ["", ""]]
    tfiles = [""] * 3

    if tree.tree_branch_split_number == "ALL" or tree.tree_branch_split_number == "all":
        SplitsToCheck = list(range(tree.subtrees["nbr"]))
    else:
        SplitsToCheck = [int(tree.tree_branch_split_number)]

    for i in SplitsToCheck:
        if i < tree.subtrees["notu"]:
            isprof = 0
        else:
            isprof = 1

        skip = 1
        for j in range(isprof + 1):
            for j1 in range(2):
                for k in range(2):
                    outfile = f"{tree.subtrees['br'][i][0]['name']}_{sequence.hot_direction[j]}{sequence.hot_direction[j1]}{sequence.hot_direction[k].upper()}"

                    if os.path.exists(f"{outfile}{sequence.file_extensions[0]}"):
                        log_print(1, 1, f"-{outfile}{sequence.file_extensions[0]} exists", file_handler)
                        continue

                    if k == 1 and os.path.exists(f"{outfile}{sequence.file_extensions[1]}"):
                        log_print(1, 1, f"-{outfile}{sequence.file_extensions[1]} exists", file_handler)
                        write_to_file(f"{outfile}{sequence.file_extensions[0]}",
                                      sequence.reverse_sequence(outfile + sequence.file_extensions[1], 1, file_handler))
                        continue

                    skip = 0

        if skip:
            log_print(0, 1, f"-Skipping branch {i + 1} of {tree.subtrees['nbr']} ...\n", file_handler)
            continue

        print(f"-Making branch {i + 1} of {tree.subtrees['nbr']} ...\n")
        if file_handler.status_file != "":
            with open(file_handler.status_file, "w") as MSA_STATUS:
                MSA_STATUS.write(f"<ul><li>Making branch {i + 1} of {tree.subtrees['nbr']}</li></ul>\n")

        for j in range(2):
            prof = dict(tree.subtrees["br"][i][j])

            if j == 0 and not isprof:
                tfiles[j] = ''

                for k in range(2):
                    pfiles[j][k] = f"prof{i}_{j}"
                    write_to_file(pfiles[j][k] + sequence.file_extensions[k],
                                  '\n'.join([sequence.bidirectional_sequences_manager['fasta'][k][m] for m in prof['otu']]))

            else:
                tfiles[j] = f"tree_{i}_{j}.dnd"
                write_to_file(tfiles[j], prof['tree'])

                for k in range(2):
                    pfiles[j][k] = f"prof{i}_{j}{k}"
                    outfile = f"{pfiles[j][k]}{sequence.file_extensions[k]}"
                    outfiler = f"{pfiles[j][k]}{sequence.file_extensions[1 - k]}"

                    if os.path.exists(outfile):
                        log_print(1, 1, f"-{outfile} exists", file_handler)
                        if not os.path.exists(outfiler):
                            write_to_file(outfiler, sequence.reverse_sequence(outfile, 1, file_handler))
                        continue

                    infile = f"in{sequence.file_extensions[k]}"
                    write_to_file(infile, ''.join([sequence.bidirectional_sequences_manager['fasta'][k][m] for m in prof['otu']]))

                    msar = sequencing_method.align_sequences(infile, tfiles[j], outfile, sequence, file_handler)
                    write_to_file(outfiler, msar)

        tfiles[2] = f"tree_{i}.dnd"
        write_to_file(tfiles[2], tree.subtrees['br'][i][2]['tree'])

        for j in range(isprof + 1):
            for j1 in range(2):
                for k in range(2):
                    outfile = f"{tree.subtrees['br'][i][0]['name']}_{sequence.hot_direction[j]}{sequence.hot_direction[j1]}{sequence.hot_direction[k].upper()}"

                    if os.path.exists(f"{outfile}{sequence.file_extensions[0]}"):
                        log_print(1, 1, f"-{outfile}{sequence.file_extensions[0]} exists", file_handler)
                        continue

                    if k == 1 and os.path.exists(f"{outfile}{sequence.file_extensions[1]}"):
                        log_print(1, 1, f"-{outfile}{sequence.file_extensions[1]} exists", file_handler)
                        write_to_file(outfile + sequence.file_extensions[0],
                                      sequence.reverse_sequence(outfile + sequence.file_extensions[1], 1, file_handler))
                        continue

                    msar = sequencing_method.align_profiles(
                        f"{pfiles[0][j]}{sequence.file_extensions[k]}", f"{pfiles[1][j1]}{sequence.file_extensions[k]}",
                        tfiles[0], tfiles[1],
                        tfiles[2],
                        f"{outfile}{sequence.file_extensions[k]}", sequence, file_handler
                    )

                    if k == 1:
                        write_to_file(f"{outfile}{sequence.file_extensions[0]}", msar)

        if debug < 3:
            run_command_line(f"rm -f prof* tr* in* *{sequence.file_extensions[1]}", file_handler)

    cleanup(0, file_handler)
    sys.exit(0)


if __name__ == "__main__":
    main()
