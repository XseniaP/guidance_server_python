"""Thin wrapper — delegates to guidance3.sequences.filters.mask_residues_main."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from guidance3.sequences.filters import mask_residues_main
if __name__ == "__main__":
    mask_residues_main()
