import Bio.Data.CodonTable
import os
import os.path
import re
import shutil
import sys
from Bio.Seq import Seq
from Bio.Data import CodonTable
import subprocess
from Bio import SeqIO, AlignIO
from guidance_common_functions import print_message_to_output, exit_on_error, update_progress
from guidance_tree_functions import fix_mafft_rough_tree
import tarfile
import glob
from time_decorator import timeit

# script_dir = os.path.dirname(os.path.realpath(__file__))
Bin = os.path.dirname(sys.argv[0])

# NEWIC2MAFFT = os.path.join(Bin, 'exec', 'newick2mafft.rb')
MSA_SET_SCORE = os.path.join(Bin, 'programs', 'msa_set_score', 'msa_set_score')
# HOT_PROGRAM = os.path.join(Bin, 'exec', 'HoT', 'COS.pl')
HOT_PROGRAM = os.path.join(Bin, 'hot_cos_main.py')
MAFFT_OP_DIST = os.path.join(Bin, 'balibase.mafft_7123_mafft.op.Dist20bins.txt')
MAFFT_OP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
MAFFT_EP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')
# HOT_GUIDANCE2_PROGRAM = os.path.join(Bin, 'exec', 'HoT_COS_GUIDANCE2.pl')
# HOT_GUIDANCE2_PROGRAM = os.path.join(Bin, 'HoT_COS_GUIDANCE2.pl')
HOT_GUIDANCE2_PROGRAM = os.path.join(Bin, 'hot_cos_main.py')
MIDPOINT_ROOTING_R = os.path.join(Bin, 'programs', 'MidPoint_Rooting.R')

# MAFFT_OP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
# MAFFT_EP_DIST_0_25 = os.path.join(Bin, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')

# newick2mafft = os.path.join(Bin, 'exec', 'newick2mafft.rb')

# MSA_Score_CSS = "http://guidance.tau.ac.il/MSA_Colored.NEW.css"
MSA_Score_CSS = "https://taux.evolseq.net/guidance/static/css/MSA_Colored.NEW.EM.css"
MidPoint_Rooting_R = os.path.join(Bin, 'programs', 'MidPoint_Rooting.R')
# phylonet_prog = os.path.join(Bin, 'exec', 'phylonet_v1_7', 'phylonet_v1_7.jar')
isEqualTopologyProg = os.path.join(Bin, 'programs', 'isEqualTree', 'isEqualTree')


def trim(line):
    line = line.lstrip()
    line = line.rstrip()
    return line


def convertNewline(input_file, working_directory):
    if os.name == 'posix':  # Unix-like system
        newline_char = '\n'
    elif os.name == 'nt':  # Windows system
        newline_char = '\r\n'
    else:
        newline_char = '\n'  # Default to newline for other systems

    file_path = os.path.join(working_directory, input_file)

    with open(file_path, 'r', newline='') as file:
        content = file.read()

    with open(file_path, 'w', newline=newline_char) as file:
        file.write(content)

def removeEndLineExtraChars(input_file, working_directory):
    working_directory = working_directory.rstrip('/') + '/'
    # lines = []

    file_path = os.path.join(working_directory, input_file)

    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as new_file:
        for line in lines:
            line = line.rstrip()
            new_file.write(f"{line}\n")

def is_amino_acid(char):
    amino_acids = "ARNDCQEGHILKMFPSTWYVXZ"
    return char.upper() in amino_acids


def sort_alignment(in_msa_file, format):
    """
    Sorts the alignment alphabetically and saves the sorted alignment in Fasta format.

    Parameters:
    - in_msa_file (str): Input alignment file.
    - format (str): Format of the input alignment file.

    Returns:
    - str: Success or error message.
    """
    out_file = f"{in_msa_file}.Sorted"

    try:
        with open(out_file, "w") as out:
            alignment = AlignIO.read(in_msa_file, format)
            alignment.sort()

            for record in alignment:
                sequence = str(record.seq)
                sequence = '\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)])
                out.write(f">{record.id}\n{sequence}\n")

        return "Success"
    except Exception as e:
        return f"Error: {str(e)}"


def convert_msa_format(infile, infile_format, outfile, outfile_format):
    """
    Reads an MSA, writes it in a different format, and removes trailing /1-129 in the sequence names.

    Parameters:
    - infile (str): Full path for the input MSA.
    - infile_format (str): Format of the input MSA.
    - outfile (str): Full path for the output MSA.
    - outfile_format (str): Format of the output MSA.

    Returns:
    - str: Success or error message.
    """
    try:
        with open(outfile, "w") as out_file:
            alignment = AlignIO.read(infile, infile_format)
            AlignIO.write(alignment, out_file, outfile_format)

        # Read the file, remove trailing /1-129 and unnecessary lines
        with open(outfile, "r") as msa_file:
            lines = msa_file.readlines()

        with open(outfile, "w") as msa_file:
            for line in lines:
                if not line.startswith(" "):
                    # Remove trailing /1-129
                    line = line.split()[0]
                    msa_file.write(line + "\n")

        return ["OK"]
    except Exception as e:
        return f"Error: {str(e)}"

# building a hash that will hold for each sequence name, a pair of DNA-AA seq put in an array
# if a seqName was not given, give a tentative name.
# create an AA file that will be submitted to MUSCLE
def translate_DNA_to_AA(input_file, output_file, codon_table_index, x_codonfile, output_html, ref_DNA_AA_seq, www_dir="", ref_seq_name=None, out_name_format="NUM"):
    counter = 0
    xFlag = "no"
    ter_mark_found = ""
    print_last_stop_codon = "no"

    codonTable = Bio.Data.CodonTable.unambiguous_dna_by_id[codon_table_index]

    with open(output_file, "w") as out_AA, open(x_codonfile, "w") as x_codon:
        x_codon.write("<html><table width=50%>\n<tr><td>Sequence Name<\/td><td>Codon Position<\/td><td>Codon<\/td><\/tr>\n")

        for seqObj in SeqIO.parse(input_file, "fasta"):
            counter += 1
            DNASequence = str(seqObj.seq)

            if DNASequence.endswith('*'):
                DNASequence = DNASequence[:-1]
                ter_mark_found = "yes"

            DNASequenceName = seqObj.id
            if not DNASequenceName:
                DNASequenceName = f"Seq_{counter}"

            legal_DNA = check_DNA_seq(DNASequence, DNASequenceName, ter_mark_found, codon_table_index, "no")
            if legal_DNA[0] != "yes":
                if legal_DNA[0] == "no":
                    return "user", legal_DNA[1]
                elif legal_DNA[0] == "fix":
                    if print_last_stop_codon == "no":
                        with open(output_html, "a") as html_out:
                            html_out.write(f"\n<p><ul><li><font color='red'><b>Please note:</b></font> {legal_DNA[1]}. <br>\nThe calculation continues nevertheless.</li></ul></p>\n")
                        print_last_stop_codon = "yes"
                    DNASequence = DNASequence[:-3]

            if '-' in DNASequence:
                return "user", "You have chosen to upload a DNA file which is not codon-aligned. Despite that, the sign '-' (which stands for a gap) was found in your input file in sequence \"{DNASequenceName}\".<br>\nIf your file is codon-aligned, please use the second box for DNA codon-aligned sequences.<br>\nOtherwise, please remove '-' signs from your input file and resubmit your query."

            xFlag, AASeq = translate_sequence(DNASequence, DNASequenceName, codon_table_index, xFlag, x_codonfile)

            if AASeq.endswith('*'):
                AASeq = AASeq[:-1]

            if out_name_format == "SEQNUM":
                seq_num = f"seq{counter - 1:04d}"
                ref_DNA_AA_seq[seq_num] = [DNASequence, AASeq]
            else:
                ref_DNA_AA_seq[counter] = [DNASequence, AASeq]

            if ref_seq_name is not None:
                for i in range(1, len(ref_seq_name) + 1):
                    if ref_seq_name[i] is not None and DNASequenceName == ref_seq_name[i]:
                        return "user", f"The Sequence name \"{DNASequenceName}\" appears in your DNA input file more than once. Please make sure that each sequence name in your input files is unique and re-submit your query.\n"

                ref_seq_name.append(DNASequenceName)

            if out_name_format == "SEQNUM":
                seq_num = f"seq{counter - 1:04d}"
                out_AA.write(f">{seq_num}\n{AASeq}\n")
            else:
                out_AA.write(f">{counter}\n{AASeq}\n")

    with open(x_codonfile, "a") as x_codon:
        x_codon.write("<\/table><\/html>\n")

    if xFlag == "yes":
        os.chmod(x_codonfile, 0o644)
        with open(output_html, "a") as html_out:
            html_out.write(f"\n<p><ul><li><font color='red'><b>Warning:</b></font> Unknown codons were found in your query file and were translated to 'X'. Please look <a href=\"{www_dir}xCodons.html\" target=\"xcodon\">here</a> for details. Please note that many 'X' signs in your translated DNA sequence have an impact\n")
            html_out.write("on the alignment's quality. <br>\nThe calculation continues nevertheless.</li></ul></p>\n")
    else:
        os.unlink(x_codonfile)

    return ["ok"]



# Check if DNA sequence is legal:
# 1. it is divisible by 3
# 2. it has no stop codon or * sign in its middle   STOP CODONS are detected according to the chosen codon convert table
# input: DNA sequence
# output: "yes" if all tests are OK, otherwise - a string that describes the input problem

def check_DNA_seq(input_DNA, DNA_input_name, ter_mark, table_codon_index, is_dna_aligned):
    ans = ("yes", "yes")
    codon_table_obj = CodonTable.unambiguous_dna_by_id[table_codon_index]
    seq_length = len(input_DNA)

    if seq_length % 3 != 0:
        if ter_mark == "yes":     # in case an earlier * sign was cut from the sequence, we inform the user
            ter_mark = "(without the last * sign)"
        ans = ("no", f"The sequence {DNA_input_name} {ter_mark} is of length {seq_length}, which is not divisible by 3.")
        return ans

    i = 0
    while i < seq_length - 2:
        codon = input_DNA[i:i + 3]

        if (not codon_table_obj.is_unknown_codon(codon) and codon_table_obj.is_ter_codon(codon)) or '*' in codon:
            if i <= seq_length - 6:
                ans = ("no", f"A Stop codon, \"{codon}\", was found in sequence {DNA_input_name} in position {i + 1}. Please verify that there are no internal stop-codons in your sequences.")
                return ans

        elif i == seq_length - 3 and '-' in codon:
            pass
        # in case the DNA input file was aligned, we ask the user to remove stop codons from the end of the sequences, as some of his seuqneces might have stop codons and some are not - and we don't want to delete it for him (to decide for him whether to remove, or to put a gap etc.)
        elif is_dna_aligned == "yes" and i == seq_length - 3 and codon_table_obj.is_ter_codon(codon):
            ans = ("no", f"Please remove the Stop Codon \"{codon}\" from your sequence {DNA_input_name}.")
            return ans

        i += 3

    return ans

def translate_sequence(DNA_sequence, DNA_sequence_name, codon_table_index, x_flag, x_codonfile):
    codon_table_obj = CodonTable.unambiguous_dna_by_id[codon_table_index]
    seq_length = len(DNA_sequence)
    i = 0
    AA_seq = ""

    while i < seq_length - 2:
        codon = DNA_sequence[i:i + 3]
        if codon == '---':
            AA = '-'
        else:
            AA = str(Seq(codon).translate(table=codon_table_obj))
            # if the AA is X we print the codon to a file and later inform the user
            if AA == "X":
                x_flag = "yes"
                with open(x_codonfile, "a") as x_codon:
                    x_codon.write(f"<tr><td>{DNA_sequence_name}</td><td>{i + 1}</td><td>{codon}</td></tr>\n")

        AA_seq += AA
        i += 3

    return x_flag, AA_seq

def validate_seqs(working_dir, input_file, seq_type, msa, codon_table):
    try:
        with open(os.path.join(working_dir, input_file), 'r') as infile, open(
                os.path.join(working_dir, input_file) + ".FIXED", 'w') as outfile:
            seq = ""
            seq_name = ""
            seq_length = 0
            warning = ""
            counter = 0
            errors = ""

            for line in infile:
                line = line.strip()

                if line != "" and not line.startswith('>'):
                    seq += line

                elif line.startswith('>'):
                    # validate prev seq
                    if seq == "" and seq_name != "":
                        # return f"The sequence named '{seq_name}' is missing<br>"
                        errors += f"The sequence named '{seq_name}' is missing\n"

                    if seq != "" and seq_name != "":
                        # validate seq according to type
                        if msa == "Yes":  # Make sure alignment length equal
                            seq_length = len(seq) if seq_length == 0 else seq_length
                            if len(seq) != seq_length:
                                # return f"The sequences of the provided MSA are not properly aligned. For example, the sequence: '{seq_name}' does not align with others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                                errors += f"The sequences of the provided MSA are not properly aligned. For example, the sequence: '{seq_name}' does not align with others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only\n"

                            if seq_type == "Codons":
                                ans = validate_seq_in_codon_align(seq, seq_name, codon_table)
                                if ans != "OK":
                                    # return ans
                                    errors += ans
                        if msa == "No":
                            if seq.endswith('-'):
                                seq = seq.rstrip('-')
                                warning = "Gap characters (-) were removed from the end of the sequences"
                            if '-' in seq:
                                # return f"Sequence named '{seq_name}' contains a gap character '-' which is illegal when sequences are submitted to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>"
                                errors += f"A sequence named '{seq_name}' contains a gap character '-', which is illegal when sequences are submitted to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option\n"

                        if seq.endswith('*'):
                            warning = "Star character (*) were removed from the end of the sequences\n"

                        ans = validate_single_seq(seq_name, seq, seq_type)
                        if ans[0] == "OK":
                            outfile.write(f">{seq_name}\n")
                            outfile.write(f"{seq}\n")
                            counter += 1
                        else:
                            # return ans
                            errors += ans

                    # Start new seq
                    if re.match(r'^>(.*)', line):
                        seq_name = re.match(r'^>(.*)', line).group(1)
                        seq_name = seq_name.strip()

                        if seq_name == "":
                            seq_num = counter + 1
                            # return f"Seq number {seq_num} has no sequence name; Please fix and resubmit<br>"
                            errors += f"Seq number {seq_num} has no sequence name; Please fix and resubmit\n"
                        else:
                            # seq_name = re.match(r'^>(.*)', line).group(1)  # Commented out since it's already assigned above
                            seq = ""

            # validate last sequence
            if seq == "" and seq_name != "":
                # return f"The sequence named '{seq_name}' is missing<br>"
                errors += f"The sequence named '{seq_name}' is missing<br>"
            else:
                if msa == "Yes":
                    seq_length = len(seq) if seq_length == 0 else seq_length
                    if len(seq) != seq_length:
                        # return f"The sequences of the provided MSA are not properly aligned. For example, the sequence: '{seq_name}' does not align with others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                        errors += f"The sequences of the provided MSA are not properly aligned. For example, the sequence: '{seq_name}' does not align with others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only\n"
                    if seq_type == "Codons":
                        ans = validate_seq_in_codon_align(seq, seq_name, codon_table)
                        if ans != "OK":
                            # return ans
                            errors += ans
                if msa == "No":
                    if seq.endswith('-'):
                        seq = seq.rstrip('-')
                        warning = "Gap characters (-) were removed from the end of the sequences"
                    if '-' in seq:
                        # return f"Sequence named '{seq_name}' contains a gap character '-' which is illegal when sequences are submitted to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>"
                        errors += f"A sequence named '{seq_name}' contains a gap character '-', which is illegal when sequences are submitted to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option\n"

                ans = validate_single_seq(seq_name, seq, seq_type)
                if ans[0] == "OK":
                    outfile.write(f">{seq_name}\n")
                    outfile.write(f"{seq}\n")
                    counter += 1
                else:
                    # return ans
                    errors += ans
            # outfile.close()
            # infile.close()
            if errors != "":
                try:
                    f = open(f'{working_dir}/errors.txt', "w")
                    f.write(errors.replace("<br>", "\n"))
                    f.close()
                except:
                    error = f"validate_seqs:Can't open {f} for writing\n"
                    return 'sys_error', error
                return errors
            return "OK", warning, input_file + ".FIXED", str(counter)
    except Exception as e:
        # raise RuntimeError(f"validate_seqs:Can't open {working_dir}{infile} : {e}\n")
        error = f"validate_seqs:Can't open {working_dir}{infile} : {e}\n"
        return 'sys_error', error


def validate_single_seq(seq_name, seq, seq_type):
    if not seq or not seq_name:
        return f"Seq: '{seq_name}' is not valid<br>"

    if not re.search(r'[ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz]+', seq) and seq_type == "AminoAcids":
        return f"Seq: '{seq_name}' is empty<br>"
    elif not re.search(r'[ACTGUNactgun]+', seq) and seq_type != "AminoAcids":
        return f"Seq: '{seq_name}' is empty<br>"

    if re.search(r'[^ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz-]', seq) and seq_type == "AminoAcids":
        return f"Seq: '{seq_name}' contained the character '{re.search(r'[^ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz-]', seq).group(0)}', which is not a standard Amino Acid\n"

    # ----------- Amit -------------

    if re.search(r'[^ACGTRYWSMKHBVDNUXacgtrywsmkhbvdnux-]', seq) and seq_type != "AminoAcids":
        wrong_char = re.search(r'[^ACGTRYWSMKHBVDNUXacgtrywsmkhbvdnux-]', seq).group(0)
        if re.search(r'[Uu]', seq) and seq_type == "Nucleotides":
            return f"Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit. <br> In addition, seq: '{seq_name}' contained the character '{wrong_char}', which is not a standard Nucleotide \n"
        return f"Seq: '{seq_name}' contained the character '{wrong_char}', which is not a standard Nucleotide\n"
    if re.search(r'[Uu]', seq) and seq_type == "Nucleotides":
        return "Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit.\n"

    # ----------- Amit -------------
    return ["OK"]


def validate_seq_in_codon_align(dna_sequence, dna_sequence_name, codon_table_index):
    # stopCodon_Found = "NO"
    aa_sequence = ""
    # codonTable_obj = CodonTable.unambiguous_dna_by_id[codon_table_index]
    # codonTable_obj = Seq.IUPAC.unambiguous_dna
    dna_sequence = dna_sequence.rstrip('\n')
    seq_length = len(dna_sequence)
    i = 0

    if seq_length % 3 > 0:
        return f"Sequence '{dna_sequence_name}' is not a valid coding sequence: the sequence is of length {seq_length} which is not divided by 3\n"

    while i < seq_length - 2:
        codon = dna_sequence[i:i + 3]

        if codon == '---':
            amino_acid = '-'
        else:
            amino_acid = str(Seq(codon).translate(codon_table_index))
            # AA = Seq(codon, IUPAC.unambiguous_dna).translate(codonTable_obj)
            # AA = codon_table.get(codon, '')

        aa_sequence += amino_acid
        i += 3

    if '*' in aa_sequence:
        return f"Sequence '{dna_sequence_name}' contains a stop codon, please remove all stop codons (from all sequences) and submit to GUIDANCE again\n"

    return "OK"


def name2code_fasta_from1(in_file_name, code_file_name, out_file_name, counter_offset=None, out_name_format='num'):
    # Open files
    in_file = SeqIO.parse(in_file_name, 'fasta')
    try:
        code_file = open(code_file_name, 'w')
    except IOError as e:
        print(f"Could not open {code_file_name}: {str(e)}")
        sys.exit()
    try:
        out_file = open(out_file_name, 'w')
    except IOError as e:
        print(f"Could not open {out_file_name}: {str(e)}")
        sys.exit()

    # Set default values
    out_name_format = "num" if out_name_format is None else out_name_format
    out_name_format = "num" if out_name_format == "" else out_name_format
    counter_offset = 1 if counter_offset is None else counter_offset
    counter_offset = 1 if counter_offset == 0 and out_name_format != 'seqNum' else counter_offset
    counter = counter_offset

    for seq_record in in_file:
        name = seq_record.id
        if seq_record.description and seq_record.description != seq_record.id:
            name += " " + seq_record.description

        if out_name_format == 'seqNum':
            sn = f'seq{counter:04d}'
            code_file.write(f"{name}\t{sn}\n")
            out_file.write(f">{sn}\n")
        else:
            code_file.write(f"{name}\t{counter}\n")
            out_file.write(f">{counter}\n")

        seq = seq_record.seq
        for i in range(0, len(seq), 60):
            out_file.write(str(seq[i:i + 60]) + "\n")

        out_file.write("\n")
        counter += 1

    # Close files
    out_file.close()
    code_file.close()

    return ["ok"]


def name2code_fasta_without_codded_out(in_file_name, code_file_name, counter_offset=None):
    """
    Convert the names in a fasta file to numbers and create a code file with the names and codes (running number).

    :param in_file_name: Input fasta file name.
    :param code_file_name: Output code file name.
    :param counter_offset: Optional offset for the counter. Default is None
    :return: Tuple with status and counter value.
    """
    counter_offset = 1 if counter_offset is None else counter_offset
    counter_offset = 1 if counter_offset == 0 else counter_offset
    counter = counter_offset
    try:
        with open(in_file_name, "r") as in_file, open(code_file_name, "a") as code_file:
            for seq_record in SeqIO.parse(in_file, "fasta"):
                name = seq_record.id
                if seq_record.description:
                    name += " " + seq_record.description
                code_file.write(f"{name}\t{counter}\n")
                counter += 1
        return "ok", counter
    except EOFError as e:
        print(f"Caught EOFError: {e}")
        sys.exit()


def convert_names_of_align_with_seed(aln, out):
    try:
        with open(aln, 'r') as infile, open(out, 'w') as outfile:
            seed_counter = 0
            for line in infile:
                if line.startswith('>_seed_'):
                    seed_counter += 1
                    line = re.sub(r'>_seed_(.*)', rf'>{seed_counter}', line)
                    outfile.write(line)
                else:
                    outfile.write(line)
        return "ok"
    except Exception as e:
        return f"Guidance::ConvertNamesOfAlignWithSeed: can't open {str(e)}"


def convert_to_codons_numbering(score_file, score_codons_file):
    try:
        with open(score_file, "r") as infile, open(score_codons_file, "w") as outfile:
            for line in infile:
                line = line.strip()
                line = trim(line)
                line = line.split()
                if line[0].isdigit():
                    line_joined = '\t'.join(line[1:])
                    outfile.write(f"{int(line[0]) * 3 - 2}\t{line_joined}\n")
                    outfile.write(f"{int(line[0]) * 3 - 1}\t{line_joined}\n")
                    outfile.write(f"{int(line[0]) * 3}\t{line_joined}\n")
                else:
                    outfile.write('\t'.join(line) + '\n')

        return ["OK"]
    except Exception as e:
        return f"Can't open the file: {str(e)}\n"


def check_sequence_licit(sequence):
    no_regular_format_char = ""
    no_regular_note = ""

    sequence = sequence.rstrip()
    for char in sequence:
        if char not in "ACDEFGHIKLMNPQRSTVWXY-":
            if char in "BJOUZ":
                if f" \"{char}\"" not in no_regular_note:
                    no_regular_note += f" \"{char}\", "
            elif char == '*':
                if '\"*\"' not in no_regular_format_char:
                    no_regular_format_char += ' \"*\", '
            elif char == '$':
                if '\"$\"' not in no_regular_format_char:
                    no_regular_format_char += ' \"$\", '
            elif char == '.':
                if '\".\"' not in no_regular_format_char:
                    no_regular_format_char += ' \".\", '
            elif char == '?':
                if '\"?\"' not in no_regular_format_char:
                    no_regular_format_char += ' \"?\", '
            elif char == '|':
                if '\"|\"' not in no_regular_format_char:
                    no_regular_format_char += ' \"|\", '
            elif char == '\\':
                no_regular_format_char = "\"\\\\\", "
            elif char:
                if f"\"{char}\"" not in no_regular_format_char:
                    no_regular_format_char += f"\"{char}\", "

    if no_regular_format_char:
        no_regular_format_char = no_regular_format_char.rstrip(', ')

    return no_regular_format_char


def check_msa_licit_and_size(msa_file, msa_format, check_sequence_licit_bool="yes"):
    # read MSA according to mode, look for ilegal characers
    # supported format by bio:seqIO are: fasta, gcg, pir
    # supported format by Bio::AlignIO are: clustalw, fasta, msf (gcg), nexus,
    #
    # If illegal chars were found, returns an array with 3 cells:
    # ('user_error',"SEQ_NAME: <the sequence name where the error was found>", "IRR_CHAR: <a string with the iregular chars>")
    # If an exception found - the SeqIO and AlighIO may through an exception, returns an array with 2 cells:
    # ('user_error',"exception")
    # Any other falut:
    # ('user_error',"could not read msa")
    # ----------------------------------
    read_seq = 0
    ans = ""
    try:
        if msa_format in ["fasta", "gcg", "pir"]:
            # with SeqIO.parse(msa_file, msa_format) as msa_fh:
            # for seq_record in msa_fh:
            for seq_record in SeqIO.parse(msa_file, msa_format):
                read_seq += 1
                if check_sequence_licit_bool == "yes":
                    ans = check_sequence_licit(seq_record.seq)
                if ans != "":
                    # break
                    return 'user_error', f"SEQ_NAME: {seq_record.id}", f"IRR_CHAR: {ans}"

        elif msa_format in ["clustalo", "nexus"]:
            # with AlignIO.read(msa_file, msa_format) as msa_fh:
            #     for seq in msa_fh:
            #       for seq_record in seq:
            for aln in SeqIO.parse(msa_file, msa_format):
                for seq_record in aln.seq:
                    read_seq += 1
                    if check_sequence_licit_bool == "yes":
                        ans = check_sequence_licit(seq_record.seq)
                    if ans != "":
                        # break
                        return 'user_error', f"SEQ_NAME: {seq_record.id}", f"IRR_CHAR: {ans}"
        else:
            return 'user_error', "wrong msa format\n"
    except Exception as e:
        return 'user_error', f"could not read msa, exception: {e}"
    if read_seq == 0:
        return 'user_error', "could not read msa\n"
    else:
        return "OK", read_seq


def set_displayname_flat(alignment):
    for record in alignment:
        record.description = record.id


def codes2name_scoresFile_NEW(Score_File, Codes_File, MSA_File, Out):
    MSA_row_Num_to_Seq_Name = {}
    Code_Names = {}

    # Read codes
    with open(Codes_File, 'r') as codes:
        for line in codes:
            Seq_name, Code = line.strip().split("\t")
            Code_Names[Code] = Seq_name

    # Read MSA to see which seq in which row and assign the correct code name
    # MSA_Depth = None
    aln = AlignIO.read(MSA_File, 'fasta')
    aln.verbose = True
    set_displayname_flat(aln)
    ans = check_msa_licit_and_size(MSA_File, "fasta", "no")
    if ans[0] == "OK":
        MSA_Depth = ans[1]
    else:
        return ": " + " ".join(ans)

    for i in range(1, int(MSA_Depth) + 1):
        seq = aln.alignment.sequences[i - 1]
        # seq = aln[:, i - 1]
        Seq_Name = Code_Names[seq.id]
        MSA_row_Num_to_Seq_Name[i] = Seq_Name

    # Add names to score file
    try:
        with open(Out, 'w') as output, open(Score_File, 'r') as scores:
            header = scores.readline()  # Header
            output.write("SEQUENCE_NAME\tSEQUENCE_SCORE\n")
            for line in scores:
                line = line.strip()
                match = re.match(r'([0-9]+)\s+([0-9.]+)', line)
                if match is not None and "#END" not in line:
                    MSA_row_Num, score = match.group(1), match.group(2)
                    MSA_row_Num = int(MSA_row_Num)
                    if MSA_row_Num in MSA_row_Num_to_Seq_Name:
                        output.write(f"{MSA_row_Num_to_Seq_Name[MSA_row_Num]}\t{score}\n")
                else:
                    output.write(line)
        return ["OK"]
    except Exception as e:
        return f"codes2name_scoresFile_NEW: Can't open file: {e}\n"


def codes2name_fasta_from1(aln_with_codes, codes_file, aln_with_names):
    codes = {}
    with open(codes_file, 'r') as codes_file_handle:
        for line in codes_file_handle:
            seq_name, code = line.strip().split("\t")
            codes[code] = seq_name

    with open(aln_with_codes, 'r') as in_file, open(aln_with_names, 'w') as out_file:
        for line in in_file:
            line = line.strip()
            if line.startswith('>'):
                match = re.match(r'>(seq[0-9]+)$', line)
                if match:
                    out_file.write(f'>{codes[match.group(1)]}\n')
                else:
                    match = re.match(r'^>([0-9]+)', line)
                    if match:
                        out_file.write(f'>{codes[match.group(1)]}\n')
                    else:
                        out_file.write(f'{line}\n')
            else:
                out_file.write(f'{line}\n')

    return ["OK"]

@timeit
def add_original_seq_names_to_the_MSA(args_library):
    args_library.Alignment_File_With_Names = args_library.Alignment_File + ".With_Names"
    ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                 f"{args_library.WorkingDir}{args_library.code_fileName}",
                                 f"{args_library.WorkingDir}{args_library.Alignment_File_With_Names}")
    if ans[0] != "OK":
        exit_on_error("sys_error",
                      f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1({args_library.WorkingDir}{args_library.Alignment_File},{args_library.WorkingDir}{args_library.code_fileName},{args_library.WorkingDir}{args_library.Alignment_File_With_Names}) failed:" +
                      ''.join(ans) + "\n", args_library)


def extract_seq_from_MSA(in_msa, seq_file, args_library):
    line_num = 0
    try:
        with open(in_msa, 'r') as infile, open(seq_file, 'w') as outfile:
            seq_num = 1
            for line in infile:
                line_num += 1
                if '>' not in line:
                    line = line.rstrip('\n')
                    line = line.replace('-', '')
                    outfile.write(line)
                else:
                    if seq_num > 1:
                        outfile.write('\n' + line)
                    else:
                        outfile.write(line)
                    seq_num += 1
    except OSError as e:
        exit_on_error('sys_error', f"extract_seq_from_MSA:Can't open {e.filename}: {str(e)}", args_library)


def convert_fs_to_lower_case(file_path):
    try:
        with open(file_path, 'r') as file:
            file_content = file.readlines()

        with open(file_path, 'w') as out_file:
            for line in file_content:
                if line.startswith('>'):
                    out_file.write(line)
                else:
                    out_file.write(line.lower())
        return None  # Success
    except OSError as e:
        return f"convert_fs_to_lower_case: Fail to open {e.filename} : {str(e)}"

@timeit
def convert_fs_to_upper_case(file_path):
    try:
        with open(file_path, 'r') as file:
            file_content = file.readlines()

        with open(file_path, 'w') as out_file:
            for line in file_content:
                if line.startswith('>'):
                    out_file.write(line)
                else:
                    out_file.write(line.upper())

        return None  # Success
    except OSError as e:
        return f"convert_fs_to_upper_case: Failed to open {e.filename} : {str(e)}"


def names_according_cos(file_path):
    try:
        with open(file_path, 'r') as file:
            msa_content = file.readlines()

        with open(file_path, 'w') as msa_file:
            for line in msa_content:
                match = re.match(r'^>([0-9]+)', line)
                # if line.startswith('>'):
                if match:
                    seq_num = int(match.group(1))
                    # seq_num = int(line[1:])
                    seq_num -= 1
                    seq_prefix = "seq"
                    if seq_num < 10:
                        seq_prefix += "000"
                    elif seq_num < 100:
                        seq_prefix += "00"
                    elif seq_num < 1000:
                        seq_prefix += "0"
                    new_line = f">{seq_prefix}{seq_num}\n"
                    msa_file.write(new_line)
                else:
                    msa_file.write(line)
        return None  # Success
    except OSError as e:
        return f"names_according_cos: Failed to open {e.filename} - {str(e)}"

@timeit
def align(args_library):
    # ---------------------------------------------
    if args_library.isServer == 1:
        # print_message_to_output("Generating the base alignment\n", args_library)
        update_progress(f"{args_library.WorkingDir}{args_library.progress_report}", "Generating the base alignment")

    args_library.Alignment_File = f"{args_library.dataset}.{args_library.MSA_Program}.aln"
    args_library.Core_Alignment_File = f"{args_library.dataset}.{args_library.MSA_Program}.CORE.aln"

    if args_library.MSA_Program == "MAFFT":
        # align with mafft
        # create the core alignment
        cmd = ""
        if 'addfragments' in args_library.align_param:
            aln_param_for_core_msa = args_library.align_param
            tmp = aln_param_for_core_msa.split('--')
            tmp = [param for param in tmp if
                   'addfragments' not in param and 'multipair' not in param and '6merpair' not in param]
            aln_param_for_core_msa = '---'.join(tmp)

            if args_library.Seq_Type in ["AminoAcids", "Codons"]:
                cmd = f"{args_library.mafft_prog} {aln_param_for_core_msa} --amino --quiet --thread {args_library.proc_num} {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Core_Alignment_File}"
            elif args_library.Seq_Type == "Nucleotides":
                cmd = f"{args_library.mafft_prog} {aln_param_for_core_msa} --nuc --quiet --thread {args_library.proc_num} {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Core_Alignment_File}"

            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write(f"Core Align: {cmd}/n")
            subprocess.run(cmd, shell=True, check=True)

            if not os.path.exists(
                    f"{args_library.WorkingDir}{args_library.Core_Alignment_File}") or os.path.getsize(
                    f"{args_library.WorkingDir}{args_library.Core_Alignment_File}") == 0:
                exit_on_error("sys_error",
                              f"Align: '{args_library.WorkingDir}{args_library.Core_Alignment_File}' is empty/does not exist", args_library)

            # create the full alignment with fragments
            if args_library.Seq_Type in ["AminoAcids", "Codons"]:
                cmd = f"{args_library.mafft_prog} {args_library.align_param} --amino --quiet {args_library.WorkingDir}{args_library.Core_Alignment_File} > {args_library.WorkingDir}{args_library.Alignment_File}"
            elif args_library.Seq_Type == "Nucleotides":
                cmd = f"{args_library.mafft_prog} {args_library.align_param} --nuc --quiet {args_library.WorkingDir}{args_library.Core_Alignment_File} > {args_library.WorkingDir}{args_library.Alignment_File}"

        else:
            if args_library.Align_Order == "aligned" and 'reorder' not in args_library.align_param:
                args_library.align_param += " --reorder"

            if args_library.Seq_Type in ["AminoAcids", "Codons"]:
                cmd = f"{args_library.mafft_prog} {args_library.align_param} --amino --quiet {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Alignment_File}"
            elif args_library.Seq_Type == "Nucleotides":
                cmd = f"{args_library.mafft_prog} {args_library.align_param} --nuc --quiet {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Alignment_File}"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Align: {cmd}\n")
        subprocess.run(cmd, shell=True, check=True)

    elif args_library.MSA_Program == "MAFFT_LINSI":
        # align with mafft
        cmd = ""
        if args_library.Seq_Type in ["AminoAcids", "Codons"]:
            cmd = f"{args_library.mafft_prog} --localpair --maxiterate 1000 --amino --quiet {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Alignment_File}"
        elif args_library.Seq_Type == "Nucleotides":
            cmd = f"{args_library.mafft_prog} --localpair --maxiterate 1000 --nuc --quiet {args_library.WorkingDir}{args_library.codded_seq_fileName} > {args_library.WorkingDir}{args_library.Alignment_File}"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Align: {cmd}\n")
        subprocess.run(cmd, shell=True, check=True)

    elif args_library.MSA_Program == "PAGAN":
        cmd = ""
        RoughTree_MAFFT = "RoughTree_MAFFT_globalpair.tree"

        # build estimated tree with mafft
        subprocess.run(
            f"{args_library.mafft_prog} --retree 0 --treeout --globalpair --reorder {args_library.WorkingDir}{args_library.codded_seq_fileName} >/dev/null",
            shell=True, check=True)

        # fix the tree (add semicolon and remove the ___)
        subprocess.run(
            f"mv {args_library.WorkingDir}{args_library.codded_seq_fileName}.tree {args_library.WorkingDir}{RoughTree_MAFFT}",
            shell=True, check=True)
        fix_mafft_rough_tree(f"{args_library.WorkingDir}{RoughTree_MAFFT}")

        # pagan cmd
        cmd = f"{args_library.pagan_prog} --seqfile {args_library.WorkingDir}{args_library.codded_seq_fileName} --treefile {args_library.WorkingDir}{RoughTree_MAFFT} --outfile {args_library.WorkingDir}{args_library.Alignment_File}"
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Align: {cmd}\n")
        subprocess.run(cmd, shell=True, check=True)

        os.rename(f"{args_library.WorkingDir}{args_library.Alignment_File}.fas",
                  f"{args_library.WorkingDir}{args_library.Alignment_File}")

    if args_library.MSA_Program == "PRANK":

        PRANK_VERSION = 0

        # if PRANK find out its version

        prank_help = subprocess.getoutput(args_library.prank_prog)

        if 'prunedata' in prank_help:

            PRANK_VERSION = "121218"

        elif 'showanc' in prank_help:

            PRANK_VERSION = "120626"

        # print(f"PRANK VERSION:>={PRANK_VERSION}\n")

        cmd = ""

        if args_library.Seq_Type in ["AminoAcids", "Nucleotides", "Codons"]:

            if PRANK_VERSION == 0:

                cmd = f"{args_library.prank_prog} {args_library.align_param} -quiet -d={args_library.WorkingDir}{args_library.codded_seq_fileName} -o={args_library.WorkingDir}{args_library.Alignment_File} -noxml > {args_library.WorkingDir}{args_library.Alignment_File}.std"

            else:  # version >=120626

                cmd = f"{args_library.prank_prog} {args_library.align_param} -quiet -d={args_library.WorkingDir}{args_library.codded_seq_fileName} -o={args_library.WorkingDir}{args_library.Alignment_File} > {args_library.WorkingDir}{args_library.Alignment_File}.std"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Align: {cmd}\n")

        if not os.path.exists(
                f"{args_library.WorkingDir}{args_library.Alignment_File}.2.fas"):  # Just to save time

            subprocess.run(cmd, shell=True, check=True)

        if PRANK_VERSION == "121218":  # The output file now named best.fas

            shutil.copy(f"{args_library.WorkingDir}{args_library.Alignment_File}.best.fas",
                        f"{args_library.WorkingDir}{args_library.Alignment_File}")

        else:

            shutil.copy(f"{args_library.WorkingDir}{args_library.Alignment_File}.2.fas",
                        f"{args_library.WorkingDir}{args_library.Alignment_File}")

        if args_library.Align_Order == "as_input":
            sort_alignment(f"{args_library.WorkingDir}{args_library.Alignment_File}", "fasta")

            args_library.Alignment_File_NOT_SORTED = args_library.Alignment_File

            args_library.Alignment_File = args_library.Alignment_File + ".Sorted"

    if args_library.MSA_Program == "CLUSTALO":
        cmd = ""

        if args_library.Seq_Type in ["AminoAcids", "Codons"]:
            cmd = f"{args_library.clustalw_prog} --infile={args_library.WorkingDir}{args_library.codded_seq_fileName} --outfile={args_library.WorkingDir}{args_library.Alignment_File} --seqtype=Protein > {args_library.WorkingDir}{args_library.Alignment_File}.std"
        elif args_library.Seq_Type == "Nucleotides":
            cmd = f"{args_library.clustalw_prog} --infile={args_library.WorkingDir}{args_library.codded_seq_fileName} --outfile={args_library.WorkingDir}{args_library.Alignment_File} --seqtype=DNA > {args_library.WorkingDir}{args_library.Alignment_File}.std"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Align: {cmd}\n")

        subprocess.run(cmd, shell=True, check=True)

        # convert_msa_format(f"{args_library.WorkingDir}{args_library.Alignment_File}", "clustal",
        #                               f"{args_library.WorkingDir}{args_library.Alignment_File}.fs", "fasta")
        # shutil.move(f"{args_library.WorkingDir}{args_library.Alignment_File}",
        #             f"{args_library.WorkingDir}{args_library.Alignment_File}.orig")
        # shutil.move(f"{args_library.WorkingDir}{args_library.Alignment_File}.fs",
        #             f"{args_library.WorkingDir}{args_library.Alignment_File}")

        if args_library.Align_Order == "as_input":
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write(f"MSA_parser.sort_alignment({args_library.WorkingDir}{args_library.Alignment_File}, fasta);\n")
            ans = sort_alignment(f"{args_library.WorkingDir}{args_library.Alignment_File}", "fasta")
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write("".join(ans) + "\n")
            args_library.Alignment_File_NOT_SORTED = args_library.Alignment_File
            args_library.Alignment_File = args_library.Alignment_File + ".Sorted"

    if args_library.MSA_Program == "MUSCLE":
        cmd = ""

        if args_library.Seq_Type in ["AminoAcids", "Codons"]:
            cmd = f"{args_library.muscle_prog} -quiet -in {args_library.WorkingDir}{args_library.codded_seq_fileName} -out {args_library.WorkingDir}{args_library.Alignment_File} -seqtype protein > {args_library.WorkingDir}{args_library.Alignment_File}.std"
        elif args_library.Seq_Type == "Nucleotides":
            cmd = f"{args_library.muscle_prog} -quiet -in {args_library.WorkingDir}{args_library.codded_seq_fileName} -out {args_library.WorkingDir}{args_library.Alignment_File} -seqtype dna > {args_library.WorkingDir}{args_library.Alignment_File}.std"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"MUSCLE Align: {cmd}\n")
        subprocess.run(cmd, shell=True, check=True)

        if args_library.Align_Order == "as_input":
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write(f"MSA_parser::sort_alignment({args_library.WorkingDir}{args_library.Alignment_File}, fasta);\n")
            ans = sort_alignment(f"{args_library.WorkingDir}{args_library.Alignment_File}", "fasta")
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write("".join(ans))
            args_library.Alignment_File_NOT_SORTED = args_library.Alignment_File
            args_library.Alignment_File = args_library.Alignment_File + ".Sorted"

    if not os.path.exists(f"{args_library.WorkingDir}{args_library.Alignment_File}") or os.path.getsize(
            f"{args_library.WorkingDir}{args_library.Alignment_File}") == 0:
        exit_on_error("sys_error",
                      f"Align: '{args_library.WorkingDir}{args_library.Alignment_File}' is empty/does not exist\n", args_library)


def codon_alignment_to_aminoacids_alignment(codon_aln, AA_aln, codon_table, xcodon_MSAfile, args_library):
    try:
        with open(codon_aln, "r") as codon_file, open(AA_aln, "w") as aa_file:
            seq = ""
            seq_name = ""
            warning = ""
            stop_codon_warning = ""
            # x_flag = ""
            for line in codon_file:
                line = line.strip()
                if not line:
                    continue
                if not line.startswith(">") and line != "":
                    seq += line
                elif line.startswith(">"):
                    # Translate previous sequence
                    if seq and seq_name:
                        # x_flag, AA_Seq, stop_codon_in_seq = translate_sequence(seq, seq_name, codon_table, "no",
                        #                                                        f"{args_library.WorkingDir}{xcodon_MSAfile}")

                        x_flag, AA_Seq = translate_sequence(seq, seq_name, codon_table, "no",
                                                                               f"{args_library.WorkingDir}{xcodon_MSAfile}")
                        if AA_Seq:
                            aa_file.write(f">{seq_name}\n{AA_Seq}\n")
                        else:
                            return "Seq:{} is empty or without any legal codons".format(seq_name)
                        if x_flag == "yes":
                            warning = "Unknown codons on the alignment file were treated as 'X' see more details <A href=\"{}\">here</A>\n".format(
                            f"{args_library.WorkingDir}{xcodon_MSAfile}")
                    # Start new sequence
                    seq_name = line[1:]
                    seq_name = seq_name.strip()
                    seq = ""

            # Validate the last sequence
            if seq and seq_name:
                # x_flag, AA_Seq, stop_codon_in_seq = translate_sequence(seq, seq_name, codon_table, "no",
                #                                                        f"{args_library.WorkingDir}{xcodon_MSAfile}")
                x_flag, AA_Seq = translate_sequence(seq, seq_name, codon_table, "no",
                                                                       f"{args_library.WorkingDir}{xcodon_MSAfile}")
                if AA_Seq:
                    aa_file.write(f">{seq_name}\n{AA_Seq}\n")
                else:
                    return "Seq:{} is empty or without any legal codons".format(seq_name)
                if x_flag == "yes":
                    warning = "Unknown codons on the alignment file were treated as 'X' see more details <A href=\"{}\">here</A>\n".format(
                    f"{args_library.WorkingDir}{xcodon_MSAfile}")
                # stop_codon_warning += ",{}".format(stop_codon_in_seq) if stop_codon_in_seq
                # stop_codon_warning = "Stop codons were removed from all the sequences" if stop_codon_in_seq

        # Close files
        # codon_file.close()
        # aa_file.close()

        warning = warning + stop_codon_warning
        return "OK", warning

    except Exception as e:
        return "sys_error", str(e)

def AA_to_DNA_aligned(input_AA_file, output_DNA_file, ref_DNA_AA_seq):
    try:
        with open(output_DNA_file, "w") as output_file_aligned:
            with open(input_AA_file, "r") as muscle_file:
                DNASequenceName = None
                AA_seq_pointer = 0
                line = 0

                for row in muscle_file:
                    row = row.strip()
                    if row.startswith('>'):
                        DNASequenceName = row[1:]
                        output_file_aligned.write(row + "\n")
                        AA_seq_pointer = 0
                        line = 0
                    else:
                        line += 1
                        amino_acids_muscle = list(row)
                        for amino_acid in amino_acids_muscle:
                            if amino_acid == '-':
                                output_file_aligned.write("---")
                            elif amino_acid.isalnum():
                                ref_aa = ref_DNA_AA_seq[DNASequenceName][1][AA_seq_pointer]
                                if amino_acid == ref_aa:
                                    output_file_aligned.write(ref_DNA_AA_seq[DNASequenceName][0][:3])
                                    ref_DNA_AA_seq[DNASequenceName][0] = ref_DNA_AA_seq[DNASequenceName][0][3:]
                                    AA_seq_pointer += 1
                                else:
                                    return "sys", f"\nAA_to_DNA_aligned: In seq name: {DNASequenceName} " \
                                                  f"Read from muscle file char: {amino_acid} at index: {AA_seq_pointer} line: {line}. " \
                                                  f"In hash found {ref_aa}\n"
                        output_file_aligned.write("\n")

        return ["ok"]
    except Exception as e:
        return "sys", f"\nAA_to_DNA_aligned: can't open file, {str(e)}\n"
    
    
def create_tar_archives(args_library):
    if args_library.PROGRAM == "GUIDANCE2":
        alt_msas_dir = os.path.join(args_library.GUIDANCE2_MSAs_Dir, 'wNames')
        os.makedirs(alt_msas_dir, exist_ok=True)
        guidance2_msas_name_pattern = os.path.join(args_library.GUIDANCE2_MSAs_Dir, '*.fasta')

        for alt_msa in glob.glob(guidance2_msas_name_pattern):
            alt_msa_no_path = os.path.basename(alt_msa)

            if args_library.Seq_Type == "Codons":
                DNA_AA_cp = {key: args_library.DNA_AA[key][:] for key in args_library.DNA_AA}
                AA_to_DNA_aligned(alt_msa, f"{alt_msa}.AA", DNA_AA_cp)
                codes2name_fasta_from1(f"{alt_msa}.AA", os.path.join(args_library.WorkingDir, args_library.code_fileName),
                                     os.path.join(alt_msas_dir, alt_msa_no_path))
                os.remove(f"{alt_msa}.AA")
            else:
                codes2name_fasta_from1(alt_msa, os.path.join(args_library.WorkingDir, args_library.code_fileName),
                                     os.path.join(alt_msas_dir, alt_msa_no_path))

        # Create a tar.gz file for alternative MSAs
        with tarfile.open(os.path.join(args_library.WorkingDir, f"{args_library.Output_Prefix}_AlternativeMSA.tar.gz"),
                          "w:gz") as tar:
            for alt_msa in glob.glob(os.path.join(alt_msas_dir, '*.fasta')):
                tar.add(alt_msa, arcname=os.path.basename(alt_msa))

        # for the webserver leave this directory open so we can create SuperMSA
        os.system(f"cp -r {alt_msas_dir} {args_library.WorkingDir}{args_library.Output_Prefix}_AlternativeMSA/")
        # list the default and alternative MSAs
        args_library.List_Of_Alternative_MSAs = f"{args_library.WorkingDir}List_Of_Default_and_AltMSAs.txt"
        os.system(f"ls -1 {args_library.WorkingDir}{args_library.Alignment_File_With_Names} {args_library.WorkingDir}{args_library.Output_Prefix}_AlternativeMSA/*.fasta > {args_library.List_Of_Alternative_MSAs}")

    if args_library.PROGRAM in ["GUIDANCE", "GUIDANCE2", "GUIDANCE3"]:
        # Tar and remove the BP dir
        cmd = f"cd {args_library.WorkingDir}; tar -czf {args_library.Output_Prefix}_BP_Dir.tar.gz ./BP"
        os.system(cmd)
        shutil.rmtree(args_library.BootStrap_Dir, ignore_errors=True)

    if args_library.PROGRAM == "HoT":
        # Tar and remove the HoT MSAs dir
        cmd = f"cd {args_library.WorkingDir}; tar -czf {args_library.dataset}.{args_library.MSA_Program}_HoT_Dir.tar.gz ./{args_library.HoT_MSAs_Dir} ./{args_library.dataset}_cos_{args_library.HoT_MSA_Program}"
        os.system(cmd)
        shutil.rmtree(os.path.join(args_library.WorkingDir, f"{args_library.dataset}_cos_{args_library.HoT_MSA_Program}"), ignore_errors = True)
        shutil.rmtree(os.path.join(args_library.WorkingDir, args_library.HoT_MSAs_Dir), ignore_errors=True)
