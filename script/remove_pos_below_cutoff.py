"""Thin wrapper — delegates to guidance3.sequences.filters.remove_pos_main."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from guidance3.sequences.filters import remove_pos_main
if __name__ == "__main__":
    remove_pos_main()
