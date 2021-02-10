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

