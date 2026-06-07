"""Thin wrapper — delegates to guidance3.pipeline.main.
Called by the web server as: python3 script/guidance_main.py VARS.json FORM.json
"""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from guidance3.pipeline.main import main

if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    print("Time: ", elapsed)