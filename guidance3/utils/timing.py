import time
import os
import sys
from functools import wraps

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)
RESULTS = os.path.join(BIN_DIR, "results/Guidance")


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        time_file = os.path.join(RESULTS, "timestamps.txt")
        mode = "a" if os.path.exists(time_file) else "w"
        with open(time_file, mode) as f:
            f.write(f'{func.__name__} \t {args} \t {kwargs} \t {total_time:.4f} \n')
        return result
    return timeit_wrapper
