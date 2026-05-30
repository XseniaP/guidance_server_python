"""Thin wrapper — delegates to guidance3.sequences.concat.main."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from guidance3.sequences.concat import main
if __name__ == "__main__":
    main()
