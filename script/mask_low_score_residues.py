import math
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv) < 6:
    sys.exit(f"USAGE: {sys.argv[0]} CODON_MSA_FILE GUIDANCE_RESIDUE_SCORES_FILE OUT_FILE CUTOFF ALPHABET\n ALPHABET can be either aa or nuc\n")

msa_file, score_file, out_file, cutoff, alphabet = sys.argv[1:]

if alphabet == "aa":
    missing_data_char = "X"
elif alphabet == "nuc":
    missing_data_char = "N"
else:
    sys.exit("ALPHABET must be either 'aa' or 'nuc'\n")

seqs = []
ids = []

with open(msa_file, "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seqs.append(list(str(seq_record.seq)))
        ids.append(seq_record.id)

with open(score_file, "r") as in_file:
    for line in in_file:
        line = line.strip()
        if line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) == 3:
            col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
            if score != float('nan') and math.isnan(score) != True and score < float(cutoff):
                seqs[row-1][col-1] = missing_data_char

with open(out_file, "w") as out:
    for i, seq_chars in enumerate(seqs):
        seq = "".join(seq_chars)
        seq_record = SeqRecord(Seq(seq), id=ids[i], description='')
        SeqIO.write(seq_record, out, "fasta")
        # out.write(">{}\n".format(ids[i]))
        # out.write("{}\n".format(seq))
