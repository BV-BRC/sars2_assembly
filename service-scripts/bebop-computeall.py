#
# Pipeline runner variant that runs all the pieces in a single worker.
#
# We feed using a redis list.
#

import sys
import os
import queue
import threading
import subprocess
import glob
import time
import socket
import threading
import json
import argparse
import redis
import re
import pickle
from pathlib import Path

import sra_sample
from hpc.worker import redis_feeder, compute_all

#
# Set up for redis.
#
# If we are node 0, we create the redis process, create a
# client, and push the sra defs to it.
# If we are not node 0, we create a client connecting to the node 0 hostname.
#

def redis_setup(args, sra_defs):

    nodeid = int(os.getenv("SLURM_NODEID"))

    redis_host = args.redis_host

    print(f"Redis host is {redis_host}")

    conn = redis.Redis(host=redis_host)

    if nodeid == 0:
        conn.delete("sra")
        for d in sra_defs:
            conn.lpush("sra", pickle.dumps(d))

    return conn
        
def compute_affinity_knl(cpu):
    aff = []
    for x in range(0,4):
        aff.append(cpu + x * 64 + 1)
    return aff

def compute_affinity_bdw(cpu):
    aff = []
    aff.append(cpu)
    return aff

def main():

    #
    # commandline processing
    #

    parser = argparse.ArgumentParser("Run pipelined SRA annotation and assembly")
    parser.add_argument('redis_host', type=str, help='Use redis at the given host for job distribution. Implies non-task-array mode')
    parser.add_argument('sra_def_file', type=str, help='File of SRA identifiers to process')
    parser.add_argument('output_dir', type=str, help='Output directory base')
    parser.add_argument('--compute-queue-size', type=int, help='Size of compute queue backlog', default=4)
    parser.add_argument('--n-computes', type=int, help='Number of compute threads', default=4)
    parser.add_argument('--n-app-threads', type=int, help='Number of threads for apps', default=4)
    parser.add_argument('--knl', action='store_true', help='Running on KNL node')
    parser.add_argument('--scratch', type=str, help='Scratch directory', default='/scratch')
    parser.add_argument('--metadata-cache', type=str, help='Cache of SRA metadata files', default='/home/olson/sra-output/data-files')
    parser.add_argument('--log-output', type=str, help='Directory to write per-thread outputs')
    parser.add_argument('--fastq-temp', type=str, help='fastq temp dir')

    args = parser.parse_args()

    sra_sample.SraSample.md_cache = Path(args.metadata_cache)
    if args.fastq_temp:
        sra_sample.SraSample.fastq_tmp = args.fastq_temp

    output_path = None
    slurm_job = os.getenv("SLURM_JOB_ID")
    if args.log_output and slurm_job:
        output_path = Path(args.log_output) / slurm_job
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"logging thread output to {output_path}")

    output = args.output_dir
    sra_defs = sra_sample.read_defs_from_file(args.sra_def_file, output)

    #
    # Find our NCBI config file and from there the SRA scratch folder. We
    # want to keep that cleared out.
    #

    ncbi_dir = None
    ncbi_config = Path.home().joinpath(".ncbi/user-settings.mkfg")
    if ncbi_config.exists():
        with ncbi_config.open() as fh:
            txt = fh.read();
            m = re.search(r'/repository/user/default-path\s*=\s*"(.*)"', txt)
            if m:
                ncbi_dir = m.group(1)

    scratch = os.getenv("SCRATCH_DIR")
    if scratch is None:
        scratch = args.scratch

    #
    # If we are in redis mode and are on node 0 of the nodelist,
    # and push the contents of the job file into the job list.
    # (We are going to have the sbatch start up script create the
    # redis server and pass the hostname here)
    #
    # In redis mode, the downloaders will block on a redis list read.
    #

    redis_conn = redis_setup(args, sra_defs)

    #
    # Our queues. This app just needs one, from the downloader to compute.
    # We limit its size so we don't pull down the entire redis queue to this host.
    #

    compute_queue = queue.Queue(args.compute_queue_size)

    N_compute = args.n_computes
    app_threads = args.n_app_threads

    cpu = 0
    if args.knl:
        compute_affinity = compute_affinity_knl
    else:
        compute_affinity = compute_affinity_bdw
        
    compute_threads = []

    for i in range(N_compute):
        
        aff = compute_affinity(cpu)
        cpu += 1
            
        t = threading.Thread(target = compute_all.worker, name=f"compute-{i}", args=[aff, app_threads, compute_queue, output_path])
        t.start()
        compute_threads.append(t)

    #
    # We run the downloader in this thread.
    #
    redis_feeder.worker([0], redis_conn, compute_queue, output_path)

    #
    # Clean up and wait.
    #
    compute_queue.join()
    print("computes done")
    for i in range(N_compute):
        compute_queue.put(None)
    for t in compute_threads:
        t.join()
    print("computes joined")

if __name__ == "__main__":
    main()
