#!/usr/bin/env python3

import sys
import re
import argparse
from collections import defaultdict

aparser = argparse.ArgumentParser()
aparser.add_argument("in_file", type=str, help="input rules file")
aparser.add_argument("out_dir", type=str, help="directory to output extracted rules into")
aparser.add_argument("base_name", type=str, help="base name of each data file")
aparser.add_argument("--zeek", action='store_true', help="parse Zeek instead of Snort/Suricata")

args = aparser.parse_args()

pathIn = args.in_file
try:
    fileIn = open(pathIn, "r")
except:
    print(f"Error: Failed to read input file '{pathIn}'")
    exit(1)

payloadLines = [
    "payload",
    "http-request",
    "http-request-header",
    "http-request-body",
    "http-reply-header",
    "http-reply-body",
    "ftp",
    "finger",
    "file-magic",
]

payloadRe = re.compile("({})".format("|".join(map(lambda l: l + " ", payloadLines))))
pcreRe = re.compile("(pcre)")

outPerCategory = defaultdict(list)

def extractSnortOrSuricata(line):
    for match in pcreRe.finditer(line):
        i = match.end()
        assert(i >= 0)

        if line[i:i+3] != ":\"/":
            continue # We do not handle these for now. Could be ones starting with '!' (do not match expression)
        assert(line[i:i+3] == ":\"/")
        i += 3 # Skip past initial '/'

        # Scan to matching unescaped quote.
        j = i
        while not (line[j] == '"' and line[j-1] != '\\'):
            if j >= len(line) - 1:
                print(line)
            j += 1

        # Retract end to exlude tailing '/xyz'
        while line[j] != '/':
            j -= 1

        out = line[i:j]
        cat = match.group(1).strip()
        outPerCategory[cat].append(out)

def extractZeek(line):
    line = line.lstrip()
    match = payloadRe.match(line)
    if not match:
        return

    i = 0
    while line[i] != ' ':
        i += 1

    i += 1
    assert(line[i] == '/')
    i += 1 # Skip past initial '/'

    # Scan to matching unescaped quote.
    j = i
    while not (line[j] == '/' and line[j-1] != '\\'):
        if j >= len(line) - 1:
            assert(False)
        j += 1

    out = line[i:j]
    cat = match.group(1).strip()
    outPerCategory[cat].append(out)

for line in fileIn:
    if args.zeek:
        extractZeek(line)
    else:
        extractSnortOrSuricata(line)

for category in outPerCategory:
    outPath = f"{args.out_dir}/{args.base_name}_{category}.re"
    print(f"Writing to {outPath}")

    duplicateCount = 0
    seen = set()
    with open(outPath, "w") as f:
        for line in outPerCategory[category]:
            if line in seen:
                # print("  > Already saw: {}".format(line))
                duplicateCount += 1
                continue

            seen.add(line)
            f.write(line + "\n")

    nFull   = len(outPerCategory[category])
    nDup    = duplicateCount
    nOut    = len(seen)
    print("  Parsed {} rules, removed {} duplicates ({:.2f}%), total {} rules".format(nFull, nDup, 100.0 * nDup / nFull, nOut))
