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

