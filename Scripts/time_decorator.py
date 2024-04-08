import time
import functools


def timefunc(f):

    @functools.wraps(f)
    def tc(*args, **kwargs):
        time_start = time.perf_counter()
        res = f(*args, **kwargs)
        print(f"Function: {f.__name__}, Time: {time.perf_counter()-time_start}")
        return res

    return tc