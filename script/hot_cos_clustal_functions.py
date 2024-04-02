import os
import shutil
import re
import sys
from hot_cos_logger import log_print, cleanup, write_to_file, run_command_line, debug


def met_init_CLO(seqtype, sequencing_method, file_handler):
    if seqtype == 2:
        print("ERROR: Codon models not supported by CLUSTAL OMEGA, use NT instead.")
        sys.exit()

        # Assuming $msa_program_path is defined elsewhere
    sequencing_method.version = f"{sequencing_method.path}"
    rc = run_command_line(f"which {sequencing_method.version} 2>&1", file_handler)
    if 'which: no' in rc:
        print(f"ERROR: Could not find {sequencing_method.version}, please make sure that {sequencing_method.version} is on your path and try again.")
        sys.exit()

    sequencing_method.version += ' ' + sequencing_method.parameters
    sequencing_method.version += " --seqtype=Protein" if seqtype == 0 else " --seqtype=DNA"

    if debug >= 2:
        sequencing_method.version += " -v --force"

    log_print(1, 0, "metvar: " + sequencing_method.version, file_handler)
    return


def align_sequences_CLO(infile, treefile, outfile, sequencing_method, sequence, file_handler):
    """Run Clustal Omega sequence alignment with the given configuration and input sequences and save the aligned sequences to the output file"""

    command_line = f"({sequencing_method.version} --infile={infile} --guidetree-in={treefile} --outfile={outfile} --output-order=input-order --outfmt=fasta 2>&1) 2>&1"
    rc = run_command_line(command_line, file_handler)
    if 'err' in rc.lower():
        log_print(0, 2, f"ERROR: clustal omega error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
        cleanup(1, file_handler)

    rtime = [float(x) for x in re.findall(r"\n(?:real|user|sys) ([0-9\.]+)", rc)]
    print(file_handler.time_file, f"{outfile} {','.join(map(str, rtime))}\n")

    rmsa = sequence.reverse_sequence(outfile, 0, file_handler)
    return rmsa

def align_profiles_CLO(pfile1, pfile2, tfile1, tfile2, tfile3, ofile, sequencing_method, sequence, file_handler):
    """Run Clustal Omega profile alignment with the given configuration and input profiles and save the aligned sequences to the output file"""

    command_line = f"({sequencing_method.version} --is-profile --p1={pfile1} --p2={pfile2} --outfmt=fasta --outfile={ofile} 2>&1) 2>&1"
    # command_line = f"({sequencing_method.version} --is-profile --p1={pfile1} --p2={pfile2} --outfmt=fasta --outfile={ofile} --guidetree-in={tfile3} 2>&1) 2>&1"
    rc = run_command_line(command_line, file_handler)

    if 'err' in rc.lower():
        log_print(0, 2, f"ERROR: clustal omega error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
        cleanup(1,file_handler)

    rtime = re.findall(r"\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)", rc)
    with open(file_handler.time_file, "a") as timefile:
        timefile.write(f"{ofile} {','.join(rtime)}\n")

    rmsa = sequence.reverse_sequence(ofile, 0, file_handler)

    return rmsa


def make_guide_tree_CLO(infile, treefile, sequencing_method, sequence, file_handler):
    """Make a guide tree"""
    command_line = f"({sequencing_method.version} --infile={infile} --guidetree-out={treefile} 2>&1) 2>&1"
    rc = run_command_line(command_line, file_handler)

    if 'err' in rc.lower():
        log_print(0, 2, f"ERROR: clustal omega error:\n{command_line}\n---\n{rc}\n---\n")
        cleanup(1, file_handler)

    rtime = re.findall(r"\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)", rc)
    with open(file_handler.time_file, "a") as timefile:
        timefile.write(f"{treefile} {','.join(rtime)}\n")

    rmsa = sequence.reverse_sequence(f"hot_H{sequence.file_extensions[0]}", 0, file_handler)

    return rmsa