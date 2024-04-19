import os
import shutil
import re
import sys
from hot_cos_logger import log_print, cleanup, write_to_file, run_command_line, debug


def met_init_MAF(seqtype, sequencing_method, file_handler):
    """Initialize MAFFT: prepare the command line with the relevant version and parameters"""
    if seqtype == 2:
        print("ERROR: Codon models not supported by MAFFT, use NT instead.")
        exit()

    sequencing_method.version = f"{sequencing_method.path}-profile"
    rc = run_command_line(f"which {sequencing_method.version} 2>&1", file_handler)
    if "which: no" in rc:
        print(f"ERROR: Could not find {sequencing_method.version}, please make sure that {sequencing_method.version} is on your path and try again.")
        exit()

    sequencing_method.version = sequencing_method.path
    rc = run_command_line(f"which {sequencing_method.version} 2>&1", file_handler)
    if "which: no" in rc:
        print(f"ERROR: Could not find {sequencing_method.version}, please make sure that {sequencing_method.version} is on your path and try again.")
        exit()

    if sequencing_method.name[2:3] == 'T':
        sequencing_method.version = sequencing_method.version + ' ' + sequencing_method.parameters

    if sequencing_method.name[2:3] == 'M':
        sequencing_method.version = sequencing_method.version + ' --localpair --maxiterate 1000'

    if seqtype == 0:
        sequencing_method.version = sequencing_method.version + ' --amino'
    else:
        sequencing_method.version = sequencing_method.version + ' --nuc'
    # sequencing_method.version = sequencing_method.version + (' --amino' if seqtype == 0 else ' --nuc')

    # debug = int(os.environ.get('DEBUG_LEVEL', 0))
    if debug < 2:
        sequencing_method.version = sequencing_method.version + ' --quiet'
    log_print(1, 0, "metvar: " + sequencing_method.version, file_handler)
    return sequencing_method.version

def align_sequences_MAF(infile, treefile, outfile, sequencing_method_version, sequence, file_handler):
    """Run MAFFT sequence alignment with the given configuration and input sequences and save the aligned sequences to the output file"""

    cmdstr = f"({sequencing_method_version} --treein {treefile}  {infile} | tr '[:lower:]' '[:upper:]' | sed 's/>SEQ/>seq/' > {outfile} 2>&1 ) 2>&1"
    rc = run_command_line(cmdstr, file_handler)

    if "err" in rc.lower():
        log_print(0, 2, f"ERROR: {file_handler.current_script_file} $$ : mafft error:\n{cmdstr}\n---\n{rc}\n---\n", file_handler)
        cleanup(1, file_handler)

    rtime = [float(x) for x in re.findall(r"\n(?:real|user|sys) ([0-9\.]+)", rc)]
    print(file_handler.time_file, f"{outfile} {','.join(map(str, rtime))}\n")

    # return MSA_check2tails(outfile, 0)
    return sequence.reverse_sequence(outfile, 0, file_handler)


def align_profiles_MAF(pfile1, pfile2, ofile, sequencing_method_path, sequence, file_handler):
    """Run MAFFT profile alignment with the given configuration and input profiles and save the aligned sequences to the output file"""

    cmdstr = f"({sequencing_method_path}-profile {pfile2} {pfile1} | tr '[:lower:]' '[:upper:]' | sed 's/>SEQ/>seq/' 2>&1 >{ofile}  )2>&1"
    rc = run_command_line(cmdstr, file_handler)
    pid = os.getpid()
    if "err" in rc.lower():
        log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : mafft error:\n{cmdstr}\n---\n{rc}\n---\n", file_handler)
        cleanup(1, file_handler)
    rtime = re.findall(r"\nreal ([0-9.]+).*\nuser ([0-9.]+).*\nsys ([0-9.]+)", rc)
    with open(file_handler.time_file,'a') as time_file_handler:
        time_file_handler.write(f"{ofile} {','.join(map(str, rtime))}\n")
    # return MSA_check2tails(ofile, 0)
    return sequence.reverse_sequence(ofile, 0, file_handler)


# create a guide tree (named <infile>.tree) by mafft and align the original sequences and produce hot_H.fasta
# Rename <infile>.tree into <treefile> and delete <infile>.tree
# <infile> is original sequences file, <treefile> is the name to give to the guide tree created
def make_guide_tree_MAF(infile, treefile, sequencing_method_version, sequence, file_handler):
    """Make a guide tree"""
    command_line = f"({sequencing_method_version} --treeout {infile} | tr '[:lower:]' '[:upper:]' | sed 's/>SEQ/>seq/' > hot_H{sequence.file_extensions[0]} 2>&1 ) 2>&1 "
    # example:
    # (mafft  --reorder --op 2.13613016216243   --amino --quiet --treeout /Users/kpolonsky/PycharmProjects/simulationMSA2/Guidance_Scores/test/MSA_25_cos_MFT/input.fasta | tr '[:lower:]' '[:upper:]' | sed 's/>SEQ/>seq/' > /Users/kpolonsky/PycharmProjects/simulationMSA2/Guidance_Scores/test/MSA_25_cos_MFT/hot_H.fasta 2>&1 ) 2>&1

    rc = run_command_line(command_line, file_handler)  # produces input.fasta.tree and hot_H.fasta
    if 'err' in rc.lower():
        pid = os.getpid()
        log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : mafft error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
        cleanup(1, file_handler)

    rtime = re.search(r"\nreal ([0-9.]+).*\nuser ([0-9.]+).*\nsys ([0-9.]+)", rc)
    if rtime:
        real_time, user_time, sys_time = map(float, rtime.groups())
        with open(file_handler.time_file, "a") as timefile_handler:
            timefile_handler.write(f"{treefile} {','.join(map(str, rtime))}\n")
    else:
        print("Time values not found in the output")

    sequence.reverse_sequence(f"hot_H{sequence.file_extensions[0]}", 0, file_handler)
    # to do text substitutions using sed, to redirect the modified content to {treefile}, and to delete the original {infile}.tree
    command_line = f"sed 's/_*//g;s/[0-9]*s/s/g' {infile}.tree > {treefile};rm {infile}.tree"
    rc = run_command_line(command_line, file_handler)  # makes guide_tree.nwk
    if 'err' in rc.lower():
        pid = os.getpid()
        log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : mafft guide tree sed error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
        cleanup(1, file_handler)

def prepare_user_tree_mafft(tree_file):
    """Take a user provided guide tree and prepare it for further use removing assigned sequence names"""
    try:
        # Create a backup of the original file
        backup_file = f"{tree_file}.ORIG"
        shutil.copy(tree_file, backup_file)

        # Read and process the tree file
        with open(tree_file, 'r') as f_in:
            lines = f_in.readlines()
            tree = ''.join(lines).rstrip(';')
            tree = tree.replace(r'(seq[0-9]{4})', r'\n\1\n')

        # Write the processed tree back to the file
        with open(tree_file, 'w') as f_out:
            f_out.write(f'{tree}\n')

    except IOError as e:
        print(f"Error processing tree file '{tree_file}': {e}")

