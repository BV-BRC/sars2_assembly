import os
import socket
import subprocess
import threading
import time
import sys
import json
import shutil

def worker(aff, threads, input_queue, output_path):
    """ Worker that runs both assembly and annotation

    The items we receive from the input_queue are SraSample instances.
    We invoke the download method to either find existing fastq files,
    use fasterq-dump or fastq-dump to produce fastq from a .sra file, or
    download from SRA. We prefer not to download from SRA so that we can
    have more consistent runtimes and not risk throttling from SRA.
    """

    me = threading.current_thread().name

    if output_path:
        print (f"{me} log to {output_path}")
        sys.stdout = open(output_path / f"{me}.stdout", "w", 1)
        sys.stderr = open(output_path / f"{me}.stderr", "w", 1)

    if aff:
        print(f"{me} starting with affinity {aff}")
        os.sched_setaffinity(0, aff)

    while True:
        print(f"{me} waiting")
        item = input_queue.get()
        if item is None:
            print(me, " got none")
            break
        print(f"{me} got {item.id}")

        sra = item.id
        out_dir = item.path

        #
        # Download
        #

        dl_output = item.download()
        if dl_output is None:
            return
        fq_files, delete_reads = dl_output

        print(f"{me} has {fq_files} delete={delete_reads}")

        #
        # Assemble
        #
        
        start = time.time()
        cmd = ["sars2-onecodex", "--max-depth", "8000"]
        cmd.extend(fq_files)
        cmd.extend([sra, out_dir, "--threads", str(threads)])
        if delete_reads:
            cmd.append("--delete-reads")

        print(cmd, file=sys.stderr)
        ret = subprocess.run(cmd,
                             stdout=open(f"{out_dir}/assemble.stdout", "w"),
                             stderr=open(f"{out_dir}/assemble.stderr", "w"))
        end = time.time()

        asm_elapsed = end - start
        with open(f"{out_dir}/RUNTIME", "w") as f:
            print(f"{start}\t{end}\t{asm_elapsed}", file=f)

        anno_elapsed = 0

        if ret.returncode != 0:
            print(f"Nonzero returncode {ret.returncode} from assembly of {sra}", file=sys.stderr)
            with open(f"{out_dir}/assembly.failure", "w") as fh:
                print(f"Nonzero returncode {ret.returncode} from assembly of {sra}", file=fh)
        else:

            #
            # Annotate
            # 

            start = time.time()

            md_file = item.metadata_file()
            if not md_file.exists():
                md_file = "/dev/null"

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
            ret = subprocess.run(cmd,
                                 cwd=out_dir,
                                 stdout=open(f"{out_dir}/annotate.stdout", "a"),
                                 stderr=open(f"{out_dir}/annotate.stderr", "a"))
            end = time.time()

            anno_elapsed = end - start

            with open(f"{out_dir}/ANNO_RUNTIME", "w") as f:
                print(f"{start}\t{end}\t{anno_elapsed}", file=f)

            if ret.returncode != 0:
                print(f"Nonzero returncode {ret.returncode} from annotation of {sra}", file=sys.stderr)
                with open(f"{out_dir}/annotation.failure", "w") as fh:
                    print(f"Nonzero returncode {ret.returncode} from annotation of {sra}", file=fh)
                #
                # Copy the raw GTO to the output gto. Best we can dow.
                #
                shutil.copyfile(f"{out_dir}/{sra}.raw.gto", f"{out_dir}/{sra}.gto")

        #
        # Create metadata to save based on this run and on the
        # container information if we are running in a container.
        #

        md = {
            "sra": sra,
            "run_index": item.idx,
            "start": start,
            "end": end,
            "elapsed": asm_elapsed,
            "annotation_elapsed": anno_elapsed,
            "host": socket.gethostname(),
            "slurm_task": os.getenv("SLURM_ARRAY_TASK_ID"),
            "slurm_job": os.getenv("SLURM_JOB_ID"),
            "slurm_cluster": os.getenv("SLURM_CLUSTER_NAME")
            }

        print(md)
        labels = "/.singularity.d/labels.json"
        if os.path.exists(labels):
            with open(labels) as f:
                label = json.load(f)
                md["container_metadata"] = label
        with open(f"{out_dir}/meta.json", "w") as f:
            json.dump(md, f, indent=2)
        for fq in fq_files:
            if os.path.exists(fq):
                os.unlink(fq)
        input_queue.task_done()

