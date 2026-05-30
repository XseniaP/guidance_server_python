"""SuperMSA concatenation utility — concatenate a list of MSAs into one."""
import os
import sys
import random
from os.path import basename
from Bio import SeqIO


def _print_to_output(msg, output_file):
    output_file.write(f"\n<ul><ul><li>{msg}</li></ul></ul>\n")


def update_msa_hash(msa_file, msa_ref, first_msa):
    if first_msa is None or first_msa == "":
        first_msa = "No"
    msa_order = []
    for rec in SeqIO.parse(msa_file, "fasta"):
        header = rec.id
        seq = str(rec.seq)
        if header in msa_ref:
            msa_ref[header] += seq
        else:
            if first_msa != "Yes":
                sys.exit(f"[ERROR] sequence '{header}' is missing in some MSAs — no missing data allowed.\n")
            msa_ref[header] = seq
            msa_order.append(header)
    if first_msa == "Yes":
        return msa_ref, msa_order
    return msa_ref


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 3:
        sys.exit("USAGE: <MSA_LIST> <OUT_ALN> [Num_Of_Aln] [Shuffle YES|NO] [isWebServer YES|NO] [outHTML]\n")

    msa_list = argv[1]
    out_msa = argv[2]
    num_of_aln_to_concat = int(argv[3]) if len(argv) > 3 else None
    shuffle = argv[4].upper() if len(argv) > 4 else "NO"
    is_web_server = argv[5].upper() if len(argv) > 5 else "NO"
    out_html = argv[6] if len(argv) > 6 and is_web_server == "YES" else ""

    with open(msa_list) as f:
        msa_files = [line.strip() for line in f.readlines()]

    if shuffle == "YES":
        random.shuffle(msa_files)

    if num_of_aln_to_concat is None:
        num_of_aln_to_concat = len(msa_files)

    new_msa, msa_order = {}, []
    for i in range(num_of_aln_to_concat):
        if i == 0:
            new_msa, msa_order = update_msa_hash(msa_files[i], new_msa, "Yes")
        else:
            new_msa = update_msa_hash(msa_files[i], new_msa, None)

    with open(out_msa, "w") as f:
        for seq_id in msa_order:
            f.write(f">{seq_id}\n{new_msa[seq_id]}\n")

    if is_web_server == "YES":
        out_msa_no_path = basename(out_msa)
        os.chmod(out_msa, 0o664)
        with open(out_html) as f:
            out_lines = f.readlines()
        with open(out_html, "w") as f:
            super_msa_section = False
            for line in out_lines:
                if "SuperMSA_results_section" in line:
                    super_msa_section = True
                    f.write(line)
                elif super_msa_section:
                    f.write(line)
                    num_of_alt = num_of_aln_to_concat - 1
                    _print_to_output(
                        f"<A HREF='{out_msa_no_path}' TARGET=_blank>The SuperMSA composed of the base MSA "
                        f"and {num_of_alt} alternative MSAs</A><br>", f)
                    super_msa_section = False
                else:
                    f.write(line)


if __name__ == "__main__":
    main()
