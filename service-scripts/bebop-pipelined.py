#
# Pipelined download/execute script
#
# Input queue is initialized with [id, sra] pairs.
# Download threads read that and write to the compute queue.
# Compute threads read the compute queue.
# 

import sys
import os
import queue
import threading
import subprocess
import glob
import time
import socket
import json
import argparse
import redis
import re
import pickle
from pathlib import Path

def compute_out_dir(out_dir_base, sra):
    path = f"{out_dir_base}/{sra[0:7]}/{sra}"
    return path
    
def create_out_dir(out_dir_base, sra):
    path = compute_out_dir(out_dir_base, sra)
    os.makedirs(path, exist_ok=True)
    return path

def create_fq_dir(scratch, task):
    path = f"{scratch}/task-{task}"
    os.makedirs(path, exist_ok=True)
    return path

def dl_worker(aff, out_dir_base, scratch_dir, ncbi_dir, redis_conn, input_queue, output_queue):
    me = threading.current_thread().name

    if aff:
        print(f"{me} starting with affinity {aff}")
        os.sched_setaffinity(0, aff)
    while True:

        if redis_conn is not None:
            pitem = redis_conn.rpop("sra")
            if pitem is None:
                break
            item = pickle.loads(pitem)
        else:
            item = input_queue.get()
            
        if item is None:
            break
        print("got item ", item)

        id, sra = item
        out_dir = create_out_dir(out_dir_base, sra)
        fq_dir = create_fq_dir(scratch_dir, id)
        fail_file = f"{out_dir}/download.failure"

        if os.path.exists(fail_file):
            os.unlink(fail_file);

        ret = subprocess.run(["p3-sra", "--id", sra, "--out", fq_dir,
                              "--metadata-file", f"{out_dir}/{sra}.json",
                              "--sra-metadata-file", f"{out_dir}/{sra}.xml",
                              ],
                             stdout=open(f"{out_dir}/download.stdout", "w"),
                             stderr=open(f"{out_dir}/download.stderr", "w"))

        if ret.returncode != 0:
            print(f"Nonzero returncode {ret.returncode} from p3-sra download of {sra}", file=sys.stderr)
            fh = open(fail_file, "w")
            print(f"Nonzero returncode {ret.returncode} from p3-sra download of {sra}", file=fh)
            fh.close()

        else:
        
            fq_files = glob.glob(f"{fq_dir}/*.fastq")

            if len(fq_files) == 3:
                fq_files = glob.glob(f"{fq_dir}/*_[12].fastq")

            fq_files.sort();
            print(f"found fq files {fq_files}", file=sys.stderr)

            if False:
                handles = []
                for i in range(1, len(fq_files)):
                    h = subprocess.Popen(["gzip", "-v", "-f", fq_files[i]])
                    handles.append(h)
                subprocess.run(["gzip", "-v", "-f", fq_files[0]])
                for h in handles:
                    h.communicate()

                fq_files = glob.glob(f"{fq_dir}/*.fastq.gz")
                fq_files.sort();

            output_queue.put([id, sra, fq_files, out_dir])

        if input_queue:
            input_queue.task_done()

        if ncbi_dir:
            path = f"{ncbi_dir}/sra/{sra}.sra"
            if os.path.exists(path):
                os.unlink(path)

def compute_worker(aff,threads, input_queue, output_queue):
    me = threading.current_thread().name
    if aff:
        print(f"{me} starting with affinity {aff}")
        os.sched_setaffinity(0, aff)
    while True:
        print(f"{me} waiting")
        item = input_queue.get()
        if item is None:
            print(me, " got none")
            break
        print(f"{me} got {item}")

        id, sra, fq_files, out_dir = item

        start = time.time()
        cmd = ["sars2-onecodex"]
        cmd.extend(fq_files)
        cmd.extend([sra, out_dir, "--threads", str(threads), "--delete-reads"])

        print(cmd, file=sys.stderr)
        ret = subprocess.run(cmd,
                             stdout=open(f"{out_dir}/assemble.stdout", "w"),
                             stderr=open(f"{out_dir}/assemble.stderr", "w"))
        end = time.time()

        elapsed = end - start

        md = {
            "sra": sra,
            "run_index": id,
            "start": start,
            "end": end,
            "elapsed": elapsed,
            "host": socket.gethostname(),
            "slurm_task": os.getenv("SLURM_ARRAY_TASK_ID"),
            "slurm_job": os.getenv("SLURM_JOB_ID"),
            "slurm_cluster": os.getenv("SLURM_CLUSTER_NAME")
            }
        labels = "/.singularity.d/labels.json"
        if os.path.exists(labels):
            with open(labels) as f:
                label = json.load(f)
                md["container_metadata"] = label
        with open(f"{out_dir}/meta.json", "w") as f:
            json.dump(md, f, indent=2)
        with open(f"{out_dir}/RUNTIME", "w") as f:
            print(f"{start}\t{end}\t{elapsed}", file=f)
        for fq in fq_files:
            if os.path.exists(fq):
                os.unlink(fq)
        input_queue.task_done()

        if ret.returncode == 0:
            output_queue.put([id, sra, md, out_dir])
        else:
            print(f"Nonzero returncode {ret.returncode} from assembly of {sra}", file=sys.stderr)
            fh = open(f"{out_dir}/assembly.failure", "w")
            print(f"Nonzero returncode {ret.returncode} from assembly of {sra}", file=fh)
            fh.close()
            

def annotate_worker(aff, input_queue):
    me = threading.current_thread().name
    if aff:
        print(f"{me} starting with affinity {aff}")
        os.sched_setaffinity(0, aff)
    while True:
        print(f"{me} waiting")
        item = input_queue.get()
        if item is None:
            print(me, " got none")
            break
        print(f"{me} got {item}")

        id, sra, md, out_dir = item

        start = time.time()

        cmd = ["p3x-create-sars-gto",
               f"{out_dir}/{sra}.fasta",
               f"{out_dir}/{sra}.json",
               f"{out_dir}/{sra}.raw.gto"];

        print(cmd, file=sys.stderr)
        subprocess.run(cmd,
                       stdout=open(f"{out_dir}/annotate.stdout", "w"),
                       stderr=open(f"{out_dir}/annotate.stderr", "w"))

        cmd = ["p3x-annotate-vigor4",
               "-i", f"{out_dir}/{sra}.raw.gto",
               "-o", f"{out_dir}/{sra}.gto"];

        print(cmd, file=sys.stderr)
        subprocess.run(cmd,
                       cwd=out_dir,
                       stdout=open(f"{out_dir}/annotate.stdout", "a"),
                       stderr=open(f"{out_dir}/annotate.stderr", "a"))
        end = time.time()

        elapsed = end - start

        md["annotation_elapsed"] = elapsed

        with open(f"{out_dir}/meta.json", "w") as f:
            json.dump(md, f, indent=2)
        with open(f"{out_dir}/RUNTIME_ANNO", "w") as f:
            print(f"{start}\t{end}\t{elapsed}", file=f)

        input_queue.task_done()

def read_sra_defs(sra_defs, output_dir):
    #
    # Read our SRA defs to find the inputs needed. Stuff them into input queue.
    #
    defs = []
    with open(sra_defs) as fh:
        idx = 1
        for line in fh:
            cols = line.rstrip().split("\t")
            sra = cols[0]

            dir = compute_out_dir(output_dir, sra)
            if not os.path.exists(f"{dir}/{sra}.gto"):
                defs.append([idx, sra])
            idx += 1

    return defs


#
# Set up for redis.
#
# If we are node 0, we create the redis process, create a
# client, and push the sra defs to it.
# If we are not node 0, we create a client connecting to the node 0 hostname.
#
def redis_setup(args, sra_defs):

    nodeid = int(os.getenv("SLURM_NODEID"))

    host = socket.gethostname()
    if args.hostlist:
        with open(args.hostlist) as f:
            line = f.readline()
            nodes = line.rstrip().split('\n')
    else:
        proc = subprocess.run(["scontrol", "show", "hostname", os.getenv("SLURM_NODELIST")], capture_output=True)
        if proc.returncode != 0:
            print("Cannot determine nodelist", file=sys.stderr)
            sys.exit(1);
        nodes = proc.stdout.decode().rstrip().split('\n')

    redis_host = nodes[0]

    print(f"Redis host is {redis_host}")

    redis_proc = None

    redis_server = "/opt/patric-common/runtime/bin/redis-server"
    
    redis_proc = None

    if nodeid == 0:
        print(f"Starting redis on {host}")
        redis_proc = subprocess.Popen([redis_server, "--protected-mode", "no"])
    else:
        print(f"Sleep on {host} to wait for redis to start on {redis_host}")
        time.sleep(10)

    conn = redis.Redis(host=redis_host)

    if nodeid == 0:
        time.sleep(3)
        for d in sra_defs:
            conn.lpush("sra", pickle.dumps(d))

    return (redis_host, conn, redis_proc)
        
def compute_affinity_knl(cpu):
    aff = []
    for x in range(0,4):
        aff.append(cpu + x * 64)
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
    parser.add_argument('job_offset', type=int, help='Value to add to the Slurm task id to determine the work number for this execution')
    parser.add_argument('entries_per_job', type=int, help='Number of entries in the data file to be processed per job')
    parser.add_argument('sra_def_file', type=str, help='File of SRA identifiers to process')
    parser.add_argument('output_dir', type=str, help='Output directory base')
    parser.add_argument('--sra-output-queue-size', type=int, help='Size of SRA downloader output queue', default=3)
    parser.add_argument('--n-downloaders', type=int, help='Number of downloader threads', default=4)
    parser.add_argument('--n-assemblers', type=int, help='Number of assembler threads', default=18)
    parser.add_argument('--n-annotators', type=int, help='Number of annotator threads', default=18)
    parser.add_argument('--n-app-threads', type=int, help='Number of threads for apps', default=4)
    parser.add_argument('--redis', action='store_true', help='Use redis for job distribution. Implies non-task-array mode')
    parser.add_argument('--knl', action='store_true', help='Running on KNL node')
    parser.add_argument('--hostlist', type=str, help='Slurm hostlist')
    parser.add_argument('--scratch', type=str, help='Scratch directory', default='/scratch')
    
    args = parser.parse_args()

    job_offset = args.job_offset
    entries_per_job = args.entries_per_job
    output = args.output_dir
    sra_defs = read_sra_defs(args.sra_def_file, output)

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
    # start the redis server and push the contents of the job file
    # into the job list.
    #
    # In redis mode, the downloaders will block on a redis list read.
    #
    # In non-redis mode, they will block on the input queue.
    #

    redis_info = redis_conn = None
    input_queue = None
    if args.redis:
        redis_info = redis_setup(args, sra_defs)
        print(redis_info)
        redis_host, redis_conn, redis_proc = redis_info

    else:
        
        start = job_offset + (int(os.getenv("SLURM_ARRAY_TASK_ID")) - 1) * entries_per_job + 1
        end = start + entries_per_job - 1
        input_queue = queue.Queue()

        for i in range(start, end+1):
            input_queue.put(sra_defs[i])

    #
    # Our queues
    #

    compute_queue = queue.Queue(args.sra_output_queue_size)
    annotate_queue = queue.Queue()

    N_download = args.n_downloaders
    N_compute = args.n_assemblers
    N_annotate = args.n_annotators
    app_threads = args.n_app_threads

    cpu = 0
    if args.knl:
        compute_affinity = compute_affinity_knl
    else:
        compute_affinity = compute_affinity_bdw
        

    download_threads = []
    for i in range(N_download):

        if args.knl:
            aff = compute_affinity(cpu)
        else:
            aff = [cpu, cpu + 1]
            cpu += 1
        cpu += 1

        t = threading.Thread(target = dl_worker, name = f"dl-{i}", args=[aff, output, scratch, ncbi_dir, redis_conn, input_queue, compute_queue])
        t.start()
        download_threads.append(t)

    compute_threads = []
    for i in range(N_compute):
        
        aff = compute_affinity(cpu)
        cpu += 1
            
        t = threading.Thread(target = compute_worker, name=f"compute-{i}", args=[aff, app_threads, compute_queue, annotate_queue])
        t.start()
        compute_threads.append(t)

    annotate_threads = []
    for i in range(N_annotate):

        aff = compute_affinity(cpu)
        cpu += 1

        t = threading.Thread(target = annotate_worker, name=f"annotate-{i}", args=[aff,annotate_queue])
        t.start()
        annotate_threads.append(t)

    if input_queue:
        input_queue.join()
        print("inputs done")
        for i in range(N_download):
            input_queue.put(None)
    for t in download_threads:
        t.join()
    compute_queue.join()
    print("computes done")
    for i in range(N_compute):
        compute_queue.put(None)
    for t in compute_threads:
        t.join()
    print("computes joined")
    annotate_queue.join()
    print("Annotates done")
    for i in range(N_annotate):
        annotate_queue.put(None)
    for t in annotate_threads:
        t.join()

if __name__ == "__main__":
    main()
