
import os
import threading
import pickle

def worker(aff, redis_conn, output_queue):
    """Redis feeder worker

    Block on a pop from the redis service. If it ever returns empty, exit the thread.

    For each data item received, unpickle to generate a SraSample.
    Push the sample onto our output queue.
    """
    me = threading.current_thread().name

    if aff:
        print(f"{me} starting with affinity {aff}")
        os.sched_setaffinity(0, aff)
    while True:

        pitem = redis_conn.rpop("sra")
        if pitem is None:
            break
        item = pickle.loads(pitem)
            
        if item is None:
            break

        print(f"{me}: got item ", item)

        output_queue.put(item)

