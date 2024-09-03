# Instructions

To run the experiments of this project follow the five steps below.

## Step 1) Download

The data files are compressed in the file `data/dfas.zip` which is hosted by Github LFS.
This is not automatically pulled by `git clone` so it is recommended to just download the entire repository as a zip-file using the green button on [Github](https://github.com/MaxRishoj/fcomp-dfa).

## Step 2) Extract

Unzip the file `data/dfas.zip` (e.g., `unzip data/dfas.zip`).
The input DFAs are then located under `data/dfas/`. 

## Step 3) Build

Run `make`.
This results in the binary `comfa`.

## Step 4) Run

Run benchmarks for all the DFAs:
``` shell
mkdir data/runs/bench -p
python3 tools/bench_all.py --sort
```
The result of each benchmark is its own file in `data/runs/bench`.
The next step gathers the results.

## Step 5) Output

Collect the output of the previous step and output into a digestible format:
``` shell
./comfa --bench_collect data/runs/bench --tcp=x  > data/results.csv
```
The output `data/runs/results.csv` can be plotted with the Jupyter notebooks, or however.
