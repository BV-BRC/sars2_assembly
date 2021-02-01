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

def create_out_dir(out_dir_base, sra):
    path = f"{out_dir_base}/{sra[0:7]}/{sra}"
    os.makedirs(path, exist_ok=True)
    return path

def create_fq_dir(scratch, task):
    path = f"{scratch}/task-{task}"
    os.makedirs(path, exist_ok=True)
    return path

def dl_worker(out_dir_base, scratch_dir, ncbi_dir, input_queue, output_queue):
    print(threading.current_thread().name)
    while True:
        item = input_queue.get()
        if item is None:
            break
        print("got item ", item)

        id, sra = item
        out_dir = create_out_dir(out_dir_base, sra)
        fq_dir = create_fq_dir(scratch_dir, id)

        subprocess.run(["p3-sra", "--id", sra, "--out", fq_dir,
                        "--metadata-file", f"{out_dir}/{sra}.json",
                        "--sra-metadata-file", f"{out_dir}/{sra}.xml",
                        ],
                       stdout=open(f"{out_dir}/download.stdout", "w"),
                       stderr=open(f"{out_dir}/download.stderr", "w"))

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
        
        input_queue.task_done()

        output_queue.put([id, sra, fq_files, out_dir])

        t = f"{ncbi_dir}/sra/{sra}.sra"
        if os.path.exists(t):
            os.unlink(t)

def compute_worker(threads, input_queue, output_queue):
    me = threading.current_thread().name
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
        subprocess.run(cmd,
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
        output_queue.put([id, sra, md, out_dir])

def annotate_worker(input_queue):
    me = threading.current_thread().name
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

def main():

    #
    # commandline processing
    #

    if len(sys.argv) != 5:
        print(f"usage: {sys.argv[0]} job-offset entries-per-job sra-def-file output-dir", file=sys.stderr)
        sys.exit(1)

    job_offset = int(sys.argv[1])
    entries_per_job = int(sys.argv[2])
    sra_defs = sys.argv[3]
    output = sys.argv[4]

    ncbi_dir = "/scratch/olson-ncbi"

    scratch = os.getenv("SCRATCH_DIR")
    if scratch is None:
        scratch = "/scratch"


    start = job_offset + (int(os.getenv("SLURM_ARRAY_TASK_ID")) - 1) * entries_per_job + 1
    end = start + entries_per_job - 1

    #
    # Our queues
    #

    input_queue = queue.Queue()
    compute_queue = queue.Queue(3)
    annotate_queue = queue.Queue()

    #
    # Read our SRA defs to find the inputs needed. Stuff them into input queue.
    #
    fh = open(sra_defs)
    idx = 1
    print("Running the following libraries:")
    for line in fh:
        cols = line.split("\t")
        sra = cols[0]
        if idx >= start and idx <= end:
            print(f"{idx}\t{sra}")
            input_queue.put([idx, sra])
        idx += 1

    N_download = 4
    N_compute = 18
    N_annotate = N_compute
    app_threads = 4

    download_threads = []
    for i in range(N_download):
        t = threading.Thread(target = dl_worker, name = f"dl-{i}", args=[output, scratch, ncbi_dir, input_queue, compute_queue])
        t.start()
        download_threads.append(t)

    compute_threads = []
    for i in range(N_compute):
        t = threading.Thread(target = compute_worker, name=f"compute-{i}", args=[app_threads, compute_queue, annotate_queue])
        t.start()
        compute_threads.append(t)

    annotate_threads = []
    for i in range(N_annotate):
        t = threading.Thread(target = annotate_worker, name=f"annotate-{i}", args=[annotate_queue])
        t.start()
        annotate_threads.append(t)

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
