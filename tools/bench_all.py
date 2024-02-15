#!/usr/bin/env python3

import sys
import os
import shutil
import re
import argparse
import random
import datetime
import argparse
import subprocess
import time
import functools
from multiprocessing import Pool

aparser = argparse.ArgumentParser()
aparser.add_argument("--timeout", type=int, default=60*60*3, help="max seconds to allow each job to run")
aparser.add_argument("--r", type=int, default=512, help="parameter r for SRG generation")
aparser.add_argument("--rerun", action="store_true", help="force rerun of existing benchmarks")
aparser.add_argument("--hpc", action="store_true", help="run benchmarks on HPC cluster")
aparser.add_argument("--number_of_jobs", type=int, default=4, help="number of jobs to keep running on HPC")
aparser.add_argument("--sort", action="store_true", help="run the benches in fast to slow order")
aparser.add_argument("--tcp_base_dir", default=".", help="tcp base directory to pass to comfa")

args = aparser.parse_args()

def create_dir_if_not_exists(path):
    if not os.path.isdir(path):
        print(f"Creating directory {path}")
        os.mkdir(path)

def bench_single_hpc(input_key, run_key):
    timeout_secs  = args.timeout
    timeout_hours = timeout_secs // (60*60)
    timeout_mins  = max(0, (timeout_secs - timeout_hours*60*60)) // 60

    create_dir_if_not_exists("tools/bench/tmp")
    out_path = "tools/bench/tmp/bench_{}-{}.sh".format(input_key, run_key)
    with open(out_path, "w") as f:
        f.write("#!/bin/sh\n")
        f.write("### General options\n")
        f.write("#BSUB -q hpc\n")
        f.write("#BSUB -J ComFA_Bench_{}-{}\n".format(input_key, run_key))
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -R \"span[hosts=1]\"\n")
        f.write("#BSUB -R \"rusage[mem=64GB]\"\n")
        f.write("#BSUB -M 64GB\n")
        f.write("#BSUB -R \"select[model == XeonGold6226R]\"\n")
        f.write("#BSUB -W {:02}:{:02}\n".format(timeout_hours, timeout_mins))
        f.write("#BSUB -u mhrpe@dtu.dk\n")
        f.write("#BSUB -B\n")
        f.write("#BSUB -N\n")
        f.write("#BSUB -o hpc_logs/Output_%J_{}_{}.out\n".format(input_key, run_key))
        f.write("#BSUB -e hpc_logs/Error_%J_{}_{}.err\n".format(input_key, run_key))
        f.write("\n")
        f.write("SCRATCH=\"{}\"\n".format(args.tcp_base_dir))
        f.write("RESULT_DIR=\"data/runs/bench\"\n")
        f.write("\n")
        f.write("./comfa {} --tcp_r {} --tcp_base_dir \"$SCRATCH\" --bench_fill \"$RESULT_DIR\" {} {} > hpc_logs/bench_{}-{}.txt\n".format(
            "--bench_rerun" if args.rerun else "",
            args.r, input_key, run_key, input_key, run_key))

    job_file = open(out_path)
    res = subprocess.run(["bsub"], stdin=job_file)
    job_file.close()

def bench_single_local(input_key, run_key):
    run_args = [
        "./comfa",
        "--tcp_r", f"{args.r}",
        "--tcp_base_dir", args.tcp_base_dir,
        "--bench_fill", "data/runs/bench",
        f"{input_key}", f"{run_key}",
    ]
    if args.rerun:
        run_args.append("--bench_rerun")

    print("Running {} {}".format(input_key, run_key))
    res = subprocess.run(run_args, capture_output=True, text=True, timeout=args.timeout)
    if res.returncode != 0:
        print("[ERROR] Failed command: {}\n".format(" ".join(run_args)))
        print("Out: \n{}\n\n".format(res.stdout))
        print("Err: \n{}\n\n".format(res.stderr))
        exit(res.returncode)

input_big_num = 100000000

def get_input_key(target_name, size):
    if target_name == "snort":
        base = 1 * input_big_num
    elif target_name == "suricata":
        base = 2 * input_big_num
    elif target_name == "zeek":
        base = 3 * input_big_num
    else:
        print("[ERROR] Unrecognized target '{}'".format(target_name))
        exit(1)

    assert(size < input_big_num)
    return base + size

def compare_bench(a, b):
    (a_input, a_run) = a
    (b_input, b_run) = b
    (a_source, a_size) = (a_input / input_big_num, a_input % input_big_num)
    (b_source, b_size) = (b_input / input_big_num, b_input % input_big_num)

    if a_size != b_size:
        return a_size - b_size
    elif a_source != b_source:
        return a_source - b_source
    elif a_run != b_run:
        return a_run - b_run
    else:
        return 0

# Execute all benchmarks.
benches = [] # (input_key, run_key)
with open("data/all_dfas.txt", "r") as f:
    for line in f:
        parts = line.split(" ")
        target_name = parts[1]
        actual_size = int(parts[3])
        input_key = get_input_key(target_name, actual_size)

        for run_key in [0, 1, 2, 3, 4, 5, 6, 7]:
            benches.append((input_key, run_key))

if args.sort:
    benches.sort(key=functools.cmp_to_key(compare_bench))

if not args.hpc:
    for (input_key, run_key) in benches:
        bench_single_local(input_key, run_key)
else:
    left = [b for b in benches]
    t0 = time.time()
    print("Starting benchmark of {} measurements.".format(len(left)))
    while left:
        time.sleep(1)
        res = subprocess.run(["bstat"], capture_output=True, text=True)
        if res.returncode != 0:
            print("[Error] 'bstat' returned {}\tErr:".format(res.returncode, res.stderr))
        else:
            out = res.stdout
            running_jobs = len(out.split("\n")) - 1

            print("[{:.2f}s] ActiveJobs={}/{} Benches={}/{}".format(
                time.time() - t0,
                running_jobs, args.number_of_jobs,
                len(left), len(benches)))

            if running_jobs >= args.number_of_jobs:
                continue;

            (input_key, run_key) = left[0]
            left.pop(0)
            print("  Starting ({}, {})".format(input_key, run_key))
            bench_single_hpc(input_key, run_key)
