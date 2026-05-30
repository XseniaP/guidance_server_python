"""Tests for the progress-report update logic in guidance3.utils.common."""
import pytest
from guidance3.utils.common import update_progress


class TestUpdateProgress:
    def _write(self, path, content):
        with open(path, "w") as f:
            f.write(content)

    def test_in_progress_becomes_finished(self, tmp_path):
        report = tmp_path / "progress.html"
        self._write(report, '<ul class="in_progress"><li>Calculating GUIDANCE2 scores</li></ul>\n')
        update_progress(str(report), "Calculating GUIDANCE2 scores")
        content = report.read_text()
        assert "finished" in content
        assert "in_progress" not in content

    def test_unrelated_line_untouched(self, tmp_path):
        report = tmp_path / "progress.html"
        original = '<ul class="in_progress"><li>Generating the base alignment</li></ul>\n'
        self._write(report, original)
        update_progress(str(report), "Calculating GUIDANCE2 scores")
        assert report.read_text() == original

    def test_running_model_line_updated(self, tmp_path):
        report = tmp_path / "progress.html"
        self._write(report, '<ul class="in_progress"><li>Running the model and selecting the best MSA</li></ul>\n')
        update_progress(str(report), "Running the model and selecting the best MSA")
        content = report.read_text()
        assert "finished" in content

    def test_multiple_lines(self, tmp_path):
        report = tmp_path / "progress.html"
        lines = (
            '<ul class="in_progress"><li>Generating the base alignment</li></ul>\n'
            '<ul class="in_progress"><li>Calculating GUIDANCE2 scores</li></ul>\n'
        )
        self._write(report, lines)
        update_progress(str(report), "Calculating GUIDANCE2 scores")
        content = report.read_text()
        assert content.count("in_progress") == 1  # only first line still in_progress
        assert content.count("finished") == 1
