import os
import pytest


@pytest.fixture
def tmp_fasta(tmp_path):
    """Write a minimal FASTA file with no blank lines and return its path."""
    content = ">seq1\nACGTACGT\n>seq2\nTGCATGCA\n>seq3\nAAAAGGGG\n"
    fasta = tmp_path / "test.fasta"
    fasta.write_text(content)
    return str(fasta)


@pytest.fixture
def tmp_fasta_with_blanks(tmp_path):
    """FASTA with empty lines between sequences (common user upload issue)."""
    content = ">seq1\nACGTACGT\n\n>seq2\nTGCATGCA\n\n>seq3\nAAAAGGGG\n\n"
    fasta = tmp_path / "blanks.fasta"
    fasta.write_text(content)
    return str(fasta)
