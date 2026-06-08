import time
import functools
import contextlib


def timefunc(f):

    @functools.wraps(f)
    def tc(*args, **kwargs):
        time_start = time.perf_counter()
        res = f(*args, **kwargs)
        print(f"Function: {f.__name__}, Time: {time.perf_counter()-time_start}")
        return res

    return tc


@contextlib.contextmanager
def timed(label):
    """Wall-clock timer for a stage, logged with a greppable [STAGE] prefix. Works across
    multiprocessing (unlike cProfile, which only sees the main process)."""
    t0 = time.perf_counter()
    print(f"[STAGE] {label}: start", flush=True)
    try:
        yield
    finally:
        print(f"[STAGE] {label}: done in {time.perf_counter()-t0:.3f}s", flush=True)