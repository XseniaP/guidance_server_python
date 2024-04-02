import os
import shutil
import re
import sys
from hot_cos_logger import log_print, cleanup, write_to_file, run_command_line, debug


# adjusted for prank v.170427: prank interface change, output file naming assumed as per this version
def met_init_MCL(seqtype, sequencing_method, file_handler):

    sequencing_method.version = f"{sequencing_method.path}"
    rc = run_command_line(f"which {sequencing_method.version} 2>&1", file_handler)
    if "which: no" in rc:
        print(f"ERROR: Could not find {sequencing_method.version}, please make sure that {sequencing_method.version} is on your path and try again.")
        exit()

    # rc = run_command_line(f"{sequencing_method.version} 2>&1", file_handler)
    #
    # sequencing_method.prkver = 1 if 'showtree' in rc else 0
    sequencing_method.version += ' ' + sequencing_method.parameters
    # if sequencing_method.prkver == 0:
    #     sequencing_method.version += " -noxml -nopost"
    # if debug < 2:
    #     sequencing_method.version += " -quiet"
    # if seqtype == 2:
    #     sequencing_method.version += " -codon"
    log_print(1, 0, "metver: muscle version is " + sequencing_method.version, file_handler)
    # print("metver: " + sequencing_method.version + "\nprank version is " + str(sequencing_method.prkver))
    return sequencing_method.version


def align_sequences_MCL(infile, treefile, outfile, sequencing_method, sequence, file_handler):
    """Run MUSCLE sequence alignment with the given configuration and input sequences and save the aligned sequences to the output file"""
    pass

    command_line = f"({sequencing_method.version} -align {infile} -usetree {treefile} -output muscle_{outfile} 2>&1) 2>&1; mv muscle_{outfile} {outfile}"
    # if sequencing_method.prkver == 1:
    #     command_line = f"({sequencing_method.version} -d={infile} -t={treefile} -o=prank_{outfile} 2>&1) 2>&1; mv prank_{outfile}.best.fas {outfile}"
    #
    # rc = run_command_line(command_line, file_handler)
    # if 'err' in rc.lower():
    #     log_print(0, 2, f"ERROR: prank error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
    #     cleanup(1, file_handler)
    #
    # rtime = re.findall(r"\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)", rc)
    # with open(file_handler.time_file, "a") as timefile:
    #     timefile.write(f"{outfile} {','.join(rtime)}\n")
    #
    # rmsa = sequence.reverse_sequence(outfile, 0, file_handler)
    # if debug < 3:
    #     os.system("rm -f prank_*")
    #
    # return rmsa

def align_profiles_MCL(pfile1, pfile2, tfile3, ofile, sequencing_method, sequence, file_handler):
    """Run MUSCLE profile alignment with the given configuration and input profiles and save the aligned sequences to the output file"""
    pass

    # command_line = f"sed 's/^>.*\$/& group_a/' {pfile1} >prank_{ofile}_inp; sed 's/^>.*\$/& group_b/' {pfile2} >>prank_{ofile}_inp"
    # rc = run_command_line(command_line, file_handler)  # makes prank profile input
    # if 'err' in rc.lower():
    #     log_print(0, 2, f"ERROR: prank input sed error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
    #     cleanup(1, file_handler)
    #
    # command_line = f"({sequencing_method.version} -partaligned -d=prank_{ofile}_inp -t={tfile3} -o=prank_{ofile} -notree 2>&1) 2>&1; mv prank_{ofile}.0.fas {ofile}"
    # if sequencing_method.prkver == 1:
    #     command_line = f"({sequencing_method.version} -partaligned -d=prank_{ofile}_inp -t={tfile3} -o=prank_{ofile} 2>&1) 2>&1; mv prank_{ofile}.fas {ofile}"
    #
    # rc = run_command_line(command_line, file_handler)
    # if 'err' in rc.lower():
    #     log_print(0, 2, f"ERROR: prank error:\n{command_line}\n---\n{rc}\n---\n", file_handler)
    #     cleanup(1, file_handler)
    #
    # rtime = re.findall(r"\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)", rc)
    # with open(file_handler.time_file, "a") as timefile:
    #     timefile.write(f"{ofile} {','.join(rtime)}\n")
    #
    # rmsa = sequence.reverse_sequence(ofile, 0, file_handler)
    # if debug < 3:
    #     os.system("rm -f prank_*")
    #
    # return rmsa

# # create a guide tree (named <infile>.tree) by mafft and align the original sequences and produce hot_H.fasta
# # Rename <infile>.tree into <treefile> and delete <infile>.tree
# # <infile> is original sequences file, <treefile> is the name to give to the guide tree created
def make_guide_tree_MCL(infile, treefile, sequencing_method, sequence, file_handler):
    """Make a guide tree"""
    pass

    # command_line = f"({sequencing_method.version} -d={infile} -o=prank_gdt 2>&1) 2>&1; mv prank_gdt.best.fas hot_H{sequence.file_extensions[0]}; mv prank_gdt.best.dnd {treefile}"
    # if sequencing_method.prkver == 1:
    #     command_line = f"({sequencing_method.version} -d={infile} -o=prank_gdt -showtree 2>&1) 2>&1; mv prank_gdt.best.fas hot_H{sequence.file_extensions[0]}; mv prank_gdt.best.dnd {treefile}"
    # rc = run_command_line(command_line, file_handler)  # produces input.fasta.tree and hot_H.fasta
    # if 'err' in rc.lower():
    #     log_print(0, 2, f"ERROR: prank error:\n{command_line}\n---\n{rc}\n---\n")
    #     cleanup(1, file_handler)
    # rtime = re.findall(r"\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)", rc)
    # with open(file_handler.time_file, "a") as timefile:
    #     timefile.write(f"{treefile} {','.join(rtime)}\n")
    #
    # rmsa = sequence.reverse_sequence(f"hot_H{sequence.file_extensions[0]}", 0, file_handler)
    # # MSA_check2tails(f"hot_H{sequence.file_extensions[0]}", 0)
    # if debug < 3:
    #     os.system("rm -f prank_*")
    #
    # return rmsa