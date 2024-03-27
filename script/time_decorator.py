import re
from functools import wraps
import time
import os
import sys

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)
RESULTS = os.path.join(BIN_DIR, "results/Guidance")

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        # job_id = 0
        # if match := re.search("results/Guidance/([\d]*)", f"{args}"):
        #     job_id = match.group(1)
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        # if job_id != 0:
        #     time_file = os.path.join(RESULTS, f"{job_id}", "timestamps.txt")
        # else:
        #     time_file = os.path.join(RESULTS, "timestamps.txt")
        time_file = os.path.join(RESULTS, "timestamps.txt")
        with open(time_file, "a") as file:
            file.write(f'{func.__name__} \t {args} \t {kwargs} \t {total_time:.4f} \n')
        return result
    return timeit_wrapper