"""Smoke tests for guidance3.hot_cos subpackage."""
import pytest
from guidance3.hot_cos.file_handler import FileHandler
from guidance3.hot_cos.sequence import Sequence
from guidance3.hot_cos.tree import Tree_
from guidance3.hot_cos.sequencing_method import SequencingMethod


class TestSequencingMethod:
    def test_known_method_accepted(self):
        sm = SequencingMethod("MFT", "mafft", [])
        assert sm.command == "MAF"

    def test_unknown_method_raises(self):
        with pytest.raises(SystemExit):
            SequencingMethod("ZZZ", "mafft", [])

    def test_parameters_joined(self):
        sm = SequencingMethod("MFT", "mafft", ["--auto"])
        assert "--auto" in sm.parameters


class TestSequence:
    def test_nucleotide_type(self):
        seq = Sequence("nt", "dummy.fasta")
        assert seq.sequence_type == 1

    def test_amino_acid_type(self):
        seq = Sequence("aa", "dummy.fasta")
        assert seq.sequence_type == 0

    def test_bidirectional_manager_initialized(self):
        seq = Sequence("nt", "dummy.fasta")
        assert "assigned_ids" in seq.bidirectional_sequences_manager
