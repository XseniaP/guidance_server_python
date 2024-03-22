from functools import wraps
import time


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        with open("/Users/kpolonsky/Downloads/timestamps.txt","a") as file:
            file.write(f'{func.__name__} \t {args} \t {kwargs} \t {total_time:.4f} \n')
            # file.write(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds\n')
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds\n')
        return result
    return timeit_wrapper