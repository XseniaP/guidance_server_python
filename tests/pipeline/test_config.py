"""Tests for guidance3.config.RunConfig."""
import pytest
from guidance3.config import RunConfig


class TestRunConfig:
    def test_instantiation(self):
        config = RunConfig()
        assert config is not None

    def test_defaults(self):
        config = RunConfig()
        assert config.PROGRAM == "GUIDANCE3"
        assert config.Bootstraps == 100
        assert config.disable_convergence is False
        assert config.proc_num == 2

    def test_disable_convergence_default_false(self):
        config = RunConfig()
        assert not config.disable_convergence

    def test_codon_table_default(self):
        config = RunConfig()
        assert config.CodonTable == 1  # Nuclear standard code
