#!/usr/bin/env python3
#
import os
import sys
import re
import random
import argparse
import subprocess
import shutil
import time

aparser = argparse.ArgumentParser()
aparser.add_argument("in_dir",  help="input directory")
aparser.add_argument("out_dir", help="output directory")
aparser.add_argument("result_file_path", help="path to file that will contain an index of all the files")
aparser.add_argument("--single_max_size", default=1000, type=int, help="maximum size for a single DFA of a rule")
aparser.add_argument("--total_max_size",  default=1000*1000, type=int, help="target size for a the largest DFA")

args = aparser.parse_args()

run_args_without_args = [
    "./comfa",
    "--gen_single_max_size",
    "{}".format(args.single_max_size),
    "--gen",
]

targets_to_short_names = {
    "snort_pcre":    "snort",
    "suricata_pcre": "suricata",
    "zeek_payload":  "zeek",
}

def compute_out_path(target, size):
    assert(target in targets_to_short_names)
    target_name = targets_to_short_names[target]
    return f"{args.out_dir}/{target_name}_{size//1000}k.dfa"

size_re = re.compile("Saving DFA of size (\d+)")
max_size = args.total_max_size

result_file = open(args.result_file_path, "w")

for target in ["snort_pcre", "suricata_pcre", "zeek_payload"]:
    size = 1000
    while size < max_size:
        in_path = "{}/{}".format(args.in_dir, target)
        out_path = compute_out_path(target, size)
        run_args = run_args_without_args + [
            in_path, out_path, f"{size}",
        ]

        try:
            t0 = time.time()
            res = subprocess.run(run_args, capture_output=True, text=True)
            dt = time.time() - t0
            if res.returncode == 42:
                sys.exit("Exceeded maximum DFA size estimate")
            elif res.returncode == 0:
                match = size_re.search(res.stdout)
                if match:
                    found_size = int(match.group(1))
                    print("Constructed DFA '{}' ({} states) in {:.2f}s".format(out_path, found_size, dt))

                    out_base_name = out_path[(len(args.out_dir) + 1):]
                    result_file.write("{} {} {} {}\n".format(out_path, targets_to_short_names[target], size, found_size))

                    size *= 2
                else:
                    print("Failed to fetch size. Aborting!")
                    print("  >> StdOut: {}".format(res.stdout))
                    print("  >> StdErr: {}".format(res.stderr))
                    exit(1)
            else:
                print("Unexpected return code: {}".format(res.returncode))
                print("  >> StdOut: {}".format(res.stdout))
                print("  >> StdErr: {}".format(res.stderr))
                sys.exit(1)
        except subprocess.TimeoutExpired as e:
            sys.exit("Timed out. Probably an issue in regex parsing.")

result_file.close()
