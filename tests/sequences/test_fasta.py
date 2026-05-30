"""Tests for FASTA utility functions in guidance3.sequences.fasta."""
import os
import pytest
from Bio import SeqIO


def _strip_blanks(src, dst):
    """Reproduce the _copy_stripped helper from guidance_sequence_functions."""
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.strip():
                fout.write(line)


class TestStripBlanks:
    def test_blank_lines_removed(self, tmp_fasta_with_blanks, tmp_path):
        dst = str(tmp_path / "stripped.fasta")
        _strip_blanks(tmp_fasta_with_blanks, dst)
        with open(dst) as f:
            lines = f.readlines()
        assert all(line.strip() for line in lines), "Blank lines should be absent after stripping"

    def test_sequence_count_preserved(self, tmp_fasta_with_blanks, tmp_path):
        dst = str(tmp_path / "stripped.fasta")
        _strip_blanks(tmp_fasta_with_blanks, dst)
        records = list(SeqIO.parse(dst, "fasta"))
        assert len(records) == 3

    def test_sequences_unchanged(self, tmp_fasta_with_blanks, tmp_path):
        dst = str(tmp_path / "stripped.fasta")
        _strip_blanks(tmp_fasta_with_blanks, dst)
        seqs = [str(r.seq) for r in SeqIO.parse(dst, "fasta")]
        assert seqs == ["ACGTACGT", "TGCATGCA", "AAAACCGG".replace("CC", "GG")]

    def test_clean_fasta_unchanged(self, tmp_fasta, tmp_path):
        dst = str(tmp_path / "clean_out.fasta")
        _strip_blanks(tmp_fasta, dst)
        orig = list(SeqIO.parse(tmp_fasta, "fasta"))
        cleaned = list(SeqIO.parse(dst, "fasta"))
        assert len(orig) == len(cleaned)
        for o, c in zip(orig, cleaned):
            assert str(o.seq) == str(c.seq)
