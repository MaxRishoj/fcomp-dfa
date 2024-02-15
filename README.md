# Instructions

If you have the password for `data/dfa.zip`, unzip the file (`unzip data/dfas.zip`) and go to [Step 2](#step-2-benchmark):
Otherwise, proceed to [Step 1](#step-1-generate-dfas).

## Step 1) Generate DFAs

First use the tiny (single rule) DFAs to create larger DFAs:
``` shell
mkdir data/dfas
python3 tools/gen_increasing_sizes.py --total_max_size 2000000 data/tiny_dfas data/dfas data/all_dfas.txt
```
Now `data/all_dfas.txt` contains a file that lists the path and size all input DFAs, used in the next step.

## Step 2) Benchmark

Run benchmarks for all the DFAs:
``` shell
mkdir data/runs/bench -p
python3 tools/bench_all.py --sort
```
The result of each benchmark is its own file in `data/runs/bench`. The next step collects it.

## Step 3) Output

Collect the output of the previous step and output into a digestible format:
``` shell
./comfa --bench_collect data/runs/bench --tcp=x  > data/results.csv
```
The output `data/runs/results.csv` can be plotted with the Jupyter notebooks, or however.
