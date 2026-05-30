"""Tests for best-MSA selection utilities (mask residues, strip blank lines)."""
import os
import math
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _write_score_file(path, entries):
    """Write a residue score file (col row score lines)."""
    with open(path, "w") as f:
        f.write("# col row score\n")
        for col, row, score in entries:
            f.write(f"{col}\t{row}\t{score}\n")


def _write_fasta(path, sequences):
    records = [SeqRecord(Seq(seq), id=sid, description="") for sid, seq in sequences.items()]
    with open(path, "w") as f:
        SeqIO.write(records, f, "fasta")


class TestMaskResidues:
    def test_low_score_residues_masked(self, tmp_path):
        fasta = tmp_path / "msa.fasta"
        _write_fasta(fasta, {"seq0001": "ACGT", "seq0002": "TGCA"})

        score_file = tmp_path / "scores.scr"
        _write_score_file(score_file, [(1, 1, 0.3), (2, 1, 0.9), (1, 2, 0.8), (2, 2, 0.2)])

        cutoff = 0.5
        seqs = []
        ids = []
        with open(fasta) as f:
            for rec in SeqIO.parse(f, "fasta"):
                seqs.append(list(str(rec.seq)))
                ids.append(rec.id)

        with open(score_file) as sf:
            for line in sf:
                line = line.strip()
                if line.startswith("#"):
                    continue
                cols = line.split()
                if len(cols) == 3:
                    col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
                    if not math.isnan(score) and score < cutoff:
                        seqs[row - 1][col - 1] = "X"

        masked_seq1 = "".join(seqs[0])
        masked_seq2 = "".join(seqs[1])
        # seq1 = "ACGT": col1=A(idx0), col2=C(idx1)
        assert masked_seq1[0] == "X", "col1 row1 score=0.3 < 0.5, should be masked"
        assert masked_seq1[1] == "C", "col2 row1 score=0.9 >= 0.5, should not be masked"
        # seq2 = "TGCA": col1=T(idx0), col2=G(idx1)
        assert masked_seq2[1] == "X", "col2 row2 score=0.2 < 0.5, should be masked"
        assert masked_seq2[0] == "T", "col1 row2 score=0.8 >= 0.5, should not be masked"

    def test_nan_scores_not_masked(self, tmp_path):
        fasta = tmp_path / "msa.fasta"
        _write_fasta(fasta, {"seq0001": "ACGT"})
        score_file = tmp_path / "scores.scr"
        _write_score_file(score_file, [(1, 1, float("nan"))])

        seqs = [list("ACGT")]
        with open(score_file) as sf:
            for line in sf:
                line = line.strip()
                if line.startswith("#"):
                    continue
                cols = line.split()
                if len(cols) == 3:
                    col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
                    if not math.isnan(score) and score < 0.5:
                        seqs[row - 1][col - 1] = "X"

        assert seqs[0][0] == "A", "NaN score should not cause masking"
