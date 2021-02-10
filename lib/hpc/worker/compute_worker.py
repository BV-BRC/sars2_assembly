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
            

