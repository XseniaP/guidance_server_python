"""Tests for scoring functions in guidance3.pipeline.scoring."""
import re
import pytest
from guidance3.pipeline.scoring import (
    _strip_all_gap_columns,
    remove_low_sp_sites_no_bioperl,
    remove_low_sp_seq,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path, seqs: dict):
    """Write {name: seq} to a FASTA file (all seqs must be the same length)."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _read_fasta(path) -> dict:
    """Return OrderedDict {name: seq} from a FASTA file (skips blank lines)."""
    seqs = {}
    order = []
    with open(path) as f:
        header = None
        parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header = line[1:]
                parts = []
                order.append(header)
            else:
                parts.append(line)
        if header is not None:
            seqs[header] = "".join(parts)
    return {k: seqs[k] for k in order}


def _write_col_score_file(path, scores: list[tuple[int, float]]):
    """Write a column score file (1-indexed col, score)."""
    with open(path, "w") as f:
        f.write("#COL_NUMBER\t#RES_PAIR_COLUMN_SCORE\n")
        for col, score in scores:
            f.write(f"{col}\t{score:.3f}\n")


def _write_seq_score_file(path, scores: list[tuple[int, float]]):
    """Write a sequence score file (1-indexed row, score)."""
    with open(path, "w") as f:
        f.write("#ROW_NUMBER\t#RES_PAIR_SEQUENCE_SCORE\n")
        for row, score in scores:
            f.write(f"{row}\t{score:.3f}\n")


def _parse_removed_col_positions(path) -> list[int]:
    """Return sorted list of 1-indexed column positions from a Removed_Col file."""
    positions = []
    with open(path) as f:
        for line in f:
            m = re.match(r"Remove Pos:\s*(\d+)", line)
            if m:
                positions.append(int(m.group(1)))
    return sorted(positions)


# ---------------------------------------------------------------------------
# _strip_all_gap_columns
# ---------------------------------------------------------------------------

class TestStripAllGapColumns:

    def test_all_gap_column_is_removed(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        _write_fasta(aln, {"s1": "A-G", "s2": "T-C"})
        n = _strip_all_gap_columns(str(aln))
        seqs = _read_fasta(str(aln))
        assert n == 1
        assert seqs["s1"] == "AG"
        assert seqs["s2"] == "TC"

    def test_returns_count_of_removed_columns(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        # cols 0,2,4 are all-gap; cols 1,3 have residues
        _write_fasta(aln, {"s1": "-A-G-", "s2": "-T-C-"})
        n = _strip_all_gap_columns(str(aln))
        assert n == 3

    def test_no_all_gap_columns_returns_zero(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        original = {"s1": "ACG-T", "s2": "TG-CA"}
        _write_fasta(aln, original)
        n = _strip_all_gap_columns(str(aln))
        seqs = _read_fasta(str(aln))
        assert n == 0
        assert seqs == original

    def test_residue_content_unchanged_after_strip(self, tmp_path):
        """Ungapped sequences must be identical before and after stripping."""
        aln = tmp_path / "aln.fasta"
        _write_fasta(aln, {"s1": "AC--GT", "s2": "TG--CA", "s3": "AA--GG"})
        _strip_all_gap_columns(str(aln))
        seqs = _read_fasta(str(aln))
        assert seqs["s1"].replace("-", "") == "ACGT"
        assert seqs["s2"].replace("-", "") == "TGCA"
        assert seqs["s3"].replace("-", "") == "AAGG"

    def test_partial_gap_column_is_kept(self, tmp_path):
        """A column with at least one residue must not be stripped."""
        aln = tmp_path / "aln.fasta"
        _write_fasta(aln, {"s1": "A-G", "s2": "TAC"})  # col 1: '-' and 'A' → keep
        n = _strip_all_gap_columns(str(aln))
        assert n == 0

    def test_blank_lines_in_input_tolerated(self, tmp_path):
        """Input with blank lines between records should not crash."""
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nA-G\n\n>s2\nT-C\n\n")
        n = _strip_all_gap_columns(str(aln))
        seqs = _read_fasta(str(aln))
        assert n == 1
        assert seqs["s1"] == "AG"
        assert seqs["s2"] == "TC"

    def test_output_column_count(self, tmp_path):
        """Original columns minus removed count equals output column count."""
        aln = tmp_path / "aln.fasta"
        orig = {"s1": "ACGT--", "s2": "TGCA--", "s3": "AAGG--"}
        _write_fasta(aln, orig)
        n = _strip_all_gap_columns(str(aln))
        seqs = _read_fasta(str(aln))
        orig_len = len(next(iter(orig.values())))
        assert len(next(iter(seqs.values()))) == orig_len - n


# ---------------------------------------------------------------------------
# remove_low_sp_sites_no_bioperl — column masking integrity
# ---------------------------------------------------------------------------

class TestRemoveLowSpSites:
    """Verify that the masked output + removed-column list are consistent with the original."""

    def _run(self, tmp_path, seqs, col_scores, cutoff):
        msa_file = str(tmp_path / "msa.fasta")
        sp_file = str(tmp_path / "col.scr")
        out_file = str(tmp_path / "filtered.fasta")
        removed_file = str(tmp_path / "removed_cols.txt")
        _write_fasta(msa_file, seqs)
        _write_col_score_file(sp_file, col_scores)
        ans = remove_low_sp_sites_no_bioperl(msa_file, sp_file, out_file, cutoff, removed_file)
        return ans, out_file, removed_file

    def test_return_ok(self, tmp_path):
        seqs = {"s1": "ACGT", "s2": "TGCA"}
        scores = [(1, 1.0), (2, 0.9), (3, 0.2), (4, 1.0)]
        ans, _, _ = self._run(tmp_path, seqs, scores, 0.5)
        assert ans[0] == "OK"

    def test_removed_positions_file_matches_low_score_columns(self, tmp_path):
        """The removed-col file must list exactly the columns whose score < cutoff."""
        seqs = {"s1": "ACGTA", "s2": "TGCAT"}
        scores = [(1, 1.0), (2, 0.3), (3, 0.9), (4, 0.1), (5, 0.8)]
        _, _, removed_file = self._run(tmp_path, seqs, scores, 0.5)
        removed_positions = _parse_removed_col_positions(removed_file)
        assert removed_positions == [2, 4]

    def test_filtered_output_equals_original_kept_columns(self, tmp_path):
        """Output MSA[:, kept] must equal original MSA[:, kept_cols]."""
        seqs = {"s1": "ACGTA", "s2": "TGCAT"}
        scores = [(1, 1.0), (2, 0.3), (3, 0.9), (4, 0.1), (5, 0.8)]
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)

        removed = set(_parse_removed_col_positions(removed_file))  # 1-indexed
        kept_indices = [i for i in range(1, 6) if i not in removed]

        filtered_seqs = _read_fasta(out_file)
        for name, orig_seq in seqs.items():
            expected = "".join(orig_seq[i - 1] for i in kept_indices)
            assert filtered_seqs[name] == expected, (
                f"{name}: expected '{expected}', got '{filtered_seqs[name]}'"
            )

    def test_output_column_count(self, tmp_path):
        seqs = {"s1": "ACGTA", "s2": "TGCAT"}
        scores = [(1, 1.0), (2, 0.3), (3, 0.9), (4, 0.1), (5, 0.8)]
        ans, out_file, _ = self._run(tmp_path, seqs, scores, 0.5)
        n_removed = ans[1]
        out_seqs = _read_fasta(out_file)
        orig_len = len(next(iter(seqs.values())))
        assert len(next(iter(out_seqs.values()))) == orig_len - n_removed

    def test_high_cutoff_removes_all_columns(self, tmp_path):
        seqs = {"s1": "ACG", "s2": "TGC"}
        scores = [(1, 0.1), (2, 0.2), (3, 0.3)]
        ans, out_file, removed_file = self._run(tmp_path, seqs, scores, 1.0)
        assert ans[1] == 3
        out_seqs = _read_fasta(out_file)
        for seq in out_seqs.values():
            assert seq == ""

    def test_no_columns_below_cutoff_leaves_msa_unchanged(self, tmp_path):
        seqs = {"s1": "ACGT", "s2": "TGCA"}
        scores = [(1, 0.9), (2, 0.8), (3, 0.95), (4, 1.0)]
        ans, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)
        assert ans[1] == 0
        assert _read_fasta(out_file) == seqs
        assert _parse_removed_col_positions(removed_file) == []


# ---------------------------------------------------------------------------
# remove_low_sp_seq — sequence removal partition integrity
# ---------------------------------------------------------------------------

class TestRemoveLowSpSeq:
    """Verify that kept + removed sequences form the complete original MSA."""

    def _run(self, tmp_path, seqs, seq_scores, cutoff):
        msa_file = str(tmp_path / "msa.fasta")
        score_file = str(tmp_path / "seq.scr")
        out_file = str(tmp_path / "kept.fasta")
        removed_file = str(tmp_path / "removed_seqs.fasta")
        _write_fasta(msa_file, seqs)
        _write_seq_score_file(score_file, seq_scores)
        ans = remove_low_sp_seq(msa_file, score_file, out_file, cutoff, removed_file)
        return ans, out_file, removed_file

    def _degap(self, seq):
        return seq.replace("-", "")

    def test_return_ok(self, tmp_path):
        seqs = {"s1": "ACGT", "s2": "TGCA", "s3": "AAGG"}
        scores = [(1, 0.9), (2, 0.3), (3, 0.8)]
        ans, _, _ = self._run(tmp_path, seqs, scores, 0.5)
        assert ans == ["OK"]

    def test_kept_plus_removed_equals_original(self, tmp_path):
        """Every original sequence must appear in exactly one of the output files."""
        seqs = {"s1": "ACGT", "s2": "TGCA", "s3": "AAGG", "s4": "CCTT"}
        scores = [(1, 0.9), (2, 0.3), (3, 0.8), (4, 0.2)]
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)

        kept_seqs = _read_fasta(out_file)
        removed_seqs = _read_fasta(removed_file)

        # Degapped sequences from both outputs
        all_output = {self._degap(s) for s in kept_seqs.values()} | \
                     {self._degap(s) for s in removed_seqs.values()}
        all_original = {self._degap(s) for s in seqs.values()}
        assert all_output == all_original

    def test_no_sequence_is_duplicated(self, tmp_path):
        seqs = {"s1": "ACGT", "s2": "TGCA", "s3": "AAGG"}
        scores = [(1, 0.9), (2, 0.3), (3, 0.8)]
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)

        kept = [self._degap(s) for s in _read_fasta(out_file).values()]
        removed = [self._degap(s) for s in _read_fasta(removed_file).values()]
        combined = kept + removed
        assert len(combined) == len(set(combined)), "A sequence appears more than once across kept and removed"

    def test_sequences_partitioned_correctly_by_cutoff(self, tmp_path):
        """Sequences with score >= cutoff go to kept; those below go to removed."""
        seqs = {"s1": "ACGT", "s2": "TGCA", "s3": "AAGG"}
        cutoff = 0.5
        scores = [(1, 0.9), (2, 0.3), (3, 0.8)]  # s2 (row 2) is below cutoff
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, cutoff)

        kept = {self._degap(s) for s in _read_fasta(out_file).values()}
        removed = {self._degap(s) for s in _read_fasta(removed_file).values()}

        assert self._degap(seqs["s1"]) in kept
        assert self._degap(seqs["s3"]) in kept
        assert self._degap(seqs["s2"]) in removed

    def test_all_sequences_kept_when_all_above_cutoff(self, tmp_path):
        seqs = {"s1": "ACGT", "s2": "TGCA", "s3": "AAGG"}
        scores = [(1, 0.9), (2, 0.8), (3, 0.7)]
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)

        kept = _read_fasta(out_file)
        removed = _read_fasta(removed_file)
        assert len(removed) == 0
        assert {self._degap(s) for s in kept.values()} == {self._degap(s) for s in seqs.values()}

    def test_residues_are_preserved_in_output(self, tmp_path):
        """Degapped output sequences must match original degapped sequences exactly."""
        seqs = {"s1": "AC-GT", "s2": "TG-CA", "s3": "AA-GG"}
        scores = [(1, 0.9), (2, 0.3), (3, 0.8)]
        _, out_file, removed_file = self._run(tmp_path, seqs, scores, 0.5)

        for path in (out_file, removed_file):
            for name, seq in _read_fasta(path).items():
                # recover original name (remove_low_sp_seq strips gaps and may wrap)
                # match by degapped content
                orig_degapped = {self._degap(v): self._degap(v) for v in seqs.values()}
                assert self._degap(seq) in orig_degapped, (
                    f"Sequence '{name}' in {path} has residues not in the original MSA"
                )
