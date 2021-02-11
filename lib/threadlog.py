
import threading
import sys
import os

thr_local = threading.local()

def open_logger(base_dir):

    me = threading.current_thread().name
    out_fh = open(base_dir / f"{me}.out", "w", 1)
    thr_local.out_fh = out_fh

    return out_fh

def get_logger():
    if hasattr(thr_local, "out_fh"):
        return thr_local.out_fh
    else:
        return sys.stdout
