"""SuperMSA concatenation utility — concatenate a list of MSAs into one."""
import argparse
import glob
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
        header = rec.description  # full header after ">", matching original manual reader
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

    # New argparse interface: guidance3-concat-msa --msas-dir <dir> --n <N> --out-file <file>
    if len(argv) > 1 and argv[1].startswith("--"):
        p = argparse.ArgumentParser(
            description="Concatenate N randomly sampled MSAs from a folder into a super-MSA."
        )
        p.add_argument("--msas-dir", required=True, help="Folder containing alternative MSA .fasta files.")
        p.add_argument("--n",        required=True, type=int, help="Number of MSAs to concatenate (sampled at random).")
        p.add_argument("--out-dir",  required=True, help="Output directory. Writes concatenated_msa.fasta inside it.")
        args = p.parse_args(argv[1:])

        fasta_files = sorted(glob.glob(os.path.join(args.msas_dir, "*.fasta")))
        if not fasta_files:
            sys.exit(f"No .fasta files found in: {args.msas_dir}")
        if args.n > len(fasta_files):
            sys.exit(f"--n ({args.n}) exceeds number of available MSAs ({len(fasta_files)})")

        os.makedirs(args.out_dir, exist_ok=True)
        msa_files = random.sample(fasta_files, args.n)
        out_msa = os.path.join(args.out_dir, "concatenated_msa.fasta")
        num_of_aln_to_concat = args.n
        is_web_server = "NO"
        out_html = ""

    # Legacy positional interface used by the web server:
    # <MSA_LIST> <OUT_ALN> [Num_Of_Aln] [Shuffle YES|NO] [isWebServer YES|NO] [outHTML]
    else:
        if len(argv) < 3:
            sys.exit("USAGE: guidance3-concat-msa --msas-dir <dir> --n <N> --out-file <file>\n"
                     "       (legacy) <MSA_LIST> <OUT_ALN> [Num_Of_Aln] [Shuffle YES|NO]\n")

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
