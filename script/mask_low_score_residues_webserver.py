#!/usr/bin/env python
import json
import math
import sys
import os
import re
import Bio.SeqIO as SeqIO

from Bio import AlignIO
# from Bio.Align.AlignInfo import GapInfo
from Bio.Seq import Seq
from warnings import warn

from Bio.SeqRecord import SeqRecord

import SharedConsts

# Load Storable module
# try:
#     import pickle as pickle
# except ImportError:
#     import pickle

# Function to print a message to the output
def print_message_to_output(msg):
    with open(output_page, "a") as output_file:
        output_file.write(f"\n<ul><li>{msg}</li></ul>\n")

# Main program
if len(sys.argv) < 4:
    sys.exit(f"USAGE: {sys.argv[0]} STORED_DATA_FILE ALPHABET CUTOFF\n"
             "ALPHABET can be either aa or nuc\n")

stored_data_file, alphabet, cutoff = sys.argv[1:]

# Load variables from the stored data file
# var_path = os.path.join(stored_data_file)
with open(stored_data_file, 'r') as vars_file:
    json_string = vars_file.read()
    VARS = json.loads(json_string)


# Set file paths
msa_file = os.path.join(VARS["WorkingDir"], VARS["Alignment_File"])
score_file = os.path.join(VARS["WorkingDir"], f"{VARS['Output_Prefix']}_res_pair_res.scr")
out_website = f"Mask_Residues_Res_{cutoff}.aln"
out_file = os.path.join(VARS["WorkingDir"], out_website)
seq_names_index = os.path.join(VARS["WorkingDir"], VARS["code_fileName"])
output_page = os.path.join(VARS["WorkingDir"], VARS["output_page"])

# Set missing data character based on alphabet
if alphabet == "aa":
    missing_data_char = "X"
elif alphabet == "nuc":
    missing_data_char = "N"
else:
    sys.exit("ALPHABET must be either 'aa' or 'nuc'\n")

# Read sequences from MSA file
seqs = []
ids = []
with open(msa_file, "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seqs.append(list(str(seq_record.seq)))
        ids.append(seq_record.id)

# Read column scores from score file and mask unreliable residues
with open(score_file, "r") as in_file:
    for line in in_file:
        line = line.strip()
        if line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) == 3:
            # col, row, score = map(int, cols)
            col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
            # if score != float('nan') and score < float(cutoff):
            if score != float('nan') and math.isnan(score) != True and score < float(cutoff):
                seqs[row-1][col-1] = missing_data_char
        else:
            warn(f"WARNING: failed to parse line: '{line}'\n")

# Read sequence names index
id_names = {}
with open(seq_names_index, "r") as seq_index:
    for line in seq_index:
        name, seq_id = line.strip().split("\t")
        id_names[seq_id] = name

# Write masked sequences to output file
with open(out_file, "w") as out:
    for i, seq_chars in enumerate(seqs):
        # seq = "".join(seq_chars)
        # seq_id = ids[i]
        # seq_name = id_names.get(seq_id, seq_id)
        # out.write(f">{seq_name}\n")
        # out.write(f"{seq}\n")
        seq = "".join(seq_chars)
        seq_record = SeqRecord(Seq(seq), id=ids[i], description='')
        SeqIO.write(seq_record, out, "fasta")

os.chmod(out_file, 0o664)

# Update the output page
with open(output_page, "r") as output:
    out_lines = output.readlines()

with open(output_page, "w") as output:
    for line in out_lines:
        if "Mask specific residues below a certain cutoff:" in line:
            print('I am here')
            output.write(line)
            os.system(f"chmod +r {out_file}")
            print_message_to_output(f"<A HREF={out_website} TARGET=_blank>The MSA after masking unreliable residues (below {cutoff})</A>")
        else:
            output.write(line)
