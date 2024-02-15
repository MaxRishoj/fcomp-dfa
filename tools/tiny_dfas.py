#!/usr/bin/env python3

import os
import sys
import re
import random
import datetime
import argparse
import subprocess
import shutil
import time
from multiprocessing import Pool

aparser = argparse.ArgumentParser()
aparser.add_argument("--max_size", default=10000, type=int, help="maximum size of a DFA")
aparser.add_argument("--n_cores", default=1, type=int, help="number of cores to utilize")
aparser.add_argument("in_dir",  type=str, help="directory to read rules from")
aparser.add_argument("out_dir", type=str, help="directory to output sizes to")
args = aparser.parse_args()

size_re  = re.compile("Size (\d+)")

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
tmp_dir = "tmp"
base_run_args = [
    "becchi/regex",
    "--no-compress",
    "--max-dfa-size", "{}".format(args.max_size),
]

timeout = 60 + 2 * int(args.max_size / 1000)

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)

def create_dir_if_not_exists(path):
    if not os.path.isdir(path):
        print(f"Creating directory {path}")
        os.mkdir(path)

create_dir_if_not_exists(args.out_dir)
create_dir_if_not_exists(tmp_dir)

def build(params):
    (rule, i, base_name) = params

    tmp_path = f"{tmp_dir}/tiny_{base_name}_{i}.re"
    with open(tmp_path, "w") as f:
        f.write(f"{rule}\n")

    out_dir  = f"{args.out_dir}/{base_name}"
    out_path = f"{out_dir}/r{i}.dfa"
    run_args = base_run_args + ["-p", tmp_path, "-e", out_path]

    success = False
    try:
        res = subprocess.run(run_args, capture_output=True, text=True, timeout=timeout)
        if res.returncode == 0:
            match = size_re.search(res.stdout)
            if match:
                success = True
    except subprocess.TimeoutExpired as e:
        pass

    remove_if_exists(tmp_path)
    if not success:
        remove_if_exists(out_path)
    else:
        # Tag the DFA file with a comment so we know the rule it came from
        with open(out_path, "a") as f:
            f.write("\n# {} : {}".format(base_name, rule))

    return success

print(f"TmpDir='{tmp_dir}' MaxSize={args.max_size} Timeout={timeout}")
total_param_list = []
for target in os.listdir(args.in_dir):
    base_name = os.path.splitext(target)[0]

    in_path = f"{args.in_dir}/{target}"
    out_dir  = f"{args.out_dir}/{base_name}"
    create_dir_if_not_exists(out_dir)

    print("Processing rules from '{}' into folder '{}'".format(in_path, out_dir))
    with open(in_path) as f:
        rules = f.readlines()

    param_list = zip(rules, range(len(rules)))
    param_list = [(rules, i, base_name) for (rules, i) in param_list]
    total_param_list += param_list

with Pool(args.n_cores) as p:
    results = p.map(build, total_param_list)

total_count = len(total_param_list)
success_count = sum(results)
print("Built {} / {} DFAs ({}%)".format(success_count, total_count, 100.0 * success_count / total_count))
