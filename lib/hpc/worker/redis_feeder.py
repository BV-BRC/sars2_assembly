
import os
import sys
import threading
import pickle
import threadlog

def worker(aff, redis_conn, output_queue, output_path):
    """Redis feeder worker

    Block on a pop from the redis service. If it ever returns empty, exit the thread.

    For each data item received, unpickle to generate a SraSample.
    Push the sample onto our output queue.
    """
    me = threading.current_thread().name
    out_fh = sys.stdout

    if output_path:
        out_fh = threadlog.open_logger(output_path)

    if aff:
        print(f"{me} starting with affinity {aff}", file=out_fh)
        os.sched_setaffinity(0, aff)
    while True:

        pitem = redis_conn.rpop("sra")
        if pitem is None:
            break
        item = pickle.loads(pitem)
            
        if item is None:
            break

        print(f"{me}: got item {item.id}", file=out_fh)

        output_queue.put(item)

