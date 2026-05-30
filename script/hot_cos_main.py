"""Thin wrapper — delegates to guidance3.hot_cos.main so the Flask app's HOT_PROGRAM
path remains stable while the canonical implementation lives in the guidance3 package."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from guidance3.hot_cos.main import main

if __name__ == "__main__":
    main()
