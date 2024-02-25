import os
import sys
import random
from os.path import basename

def print_message_to_output(msg, output):
    output.write(f"\n<ul><ul><li>{msg}</li></ul></ul>\n")

def update_msa_hash(msa_file, msa_ref, first_msa):
    if first_msa is None or first_msa is "":
        first_msa = "No"
    msa_order = []
    with open(msa_file, 'r') as file:
        fasta_line = file.readline()
        while fasta_line:
            header = fasta_line[1:].strip()
            # fasta_line = file.readline().strip()
            fasta_line = file.readline()
            seq = ""
            while fasta_line and fasta_line[0] != ">":
                seq += fasta_line.strip()
                # fasta_line = file.readline().strip()
                fasta_line = file.readline()
            if header in msa_ref:
                msa_ref[header] += seq
            else:
                if first_msa != "Yes":
                    sys.exit("[ERROR] sequence named {} is missing on some of the MSAs....! Notice that this script assume no missing data...\n".format(header))
                else:
                    msa_ref[header] = seq
                    msa_order.append(header)
    if first_msa == "Yes":
        return msa_ref, msa_order
    else:
        return msa_ref

if len(sys.argv) < 3:
    sys.exit("USAGE: <MSA_LIST> <OUT_ALN> <?Num_Of_Aln_to_concat [Number of MSAs in LIST]> <?Shuffle {YES[NO]}> <?isWebServerMode> <?outHTML>\n")

msa_list = sys.argv[1]
out_msa = sys.argv[2]
num_of_aln_to_concat = int(sys.argv[3]) if len(sys.argv) > 3 else None
shuffle = sys.argv[4].upper() if len(sys.argv) > 4 else "NO"
is_web_server = sys.argv[5].upper() if len(sys.argv) > 5 else "NO"
out_html = sys.argv[6] if len(sys.argv) > 6 and is_web_server == "YES" else ""

new_msa = {}
msa_order = []

with open(msa_list, 'r') as file:
    msa_files = file.readlines()
    msa_files = [line.strip() for line in msa_files]

if shuffle == "YES":
    random.shuffle(msa_files)

if num_of_aln_to_concat is None:
    num_of_aln_to_concat = len(msa_files)

for i in range(num_of_aln_to_concat):
    if i == 0:
        new_msa, msa_order = update_msa_hash(msa_files[i], new_msa, "Yes")
    else:
        new_msa = update_msa_hash(msa_files[i], new_msa, None)

with open(out_msa, 'w') as output_msa:
    for seq_id in msa_order:
        output_msa.write(f">{seq_id}\n")
        output_msa.write(f"{new_msa[seq_id]}\n")

if is_web_server == "YES":
    out_msa_no_path = basename(out_msa)
    os.chmod(out_msa, 0o664)
    with open(out_html, 'r') as output:
        out_lines = output.readlines()
    with open(out_html, 'w') as output:
        super_msa_section = False
        for line in out_lines:
            if "SuperMSA_results_section" in line:
                super_msa_section = True
                output.write(line)
            elif super_msa_section:
                output.write(line)
                num_of_alt = num_of_aln_to_concat - 1
                print_message_to_output(f"<A HREF='{out_msa_no_path}' TARGET=_blank>The SuperMSA composed of the base MSA and {num_of_alt} alternative MSAs</A><br>", output)
                super_msa_section = False
            else:
                output.write(line)
