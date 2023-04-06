#!/usr/bin/env python3
import subprocess
import csv
import sys 

###
# Obtain convergence results
###
ranks = [2, 2, 1]
base_res = 32
num_res = 5

for boundary in ["4,4,4,4,4,4", "4,4,3,3,3,3"]:
    for kernel in [1, 8, 9, 10, 13, 14]:

        args  = ["./flups_validation_a2a"]
        args += ["--center=1"]
        args += ["--nsolve=1"]
        args += ["--lda=1"]
        args += ["--dom=1,1,1"]
        args += ["--bc=" + boundary]
        args += [f"--kernel={kernel}"]
        args += [f"--nres={num_res}"] 
        args += [f"--np={ranks[0]},{ranks[1]},{ranks[2]}"]
        args += [f"--res={base_res},{base_res},{base_res}"]

        nranks = ranks[0]*ranks[1]*ranks[2]
        completed_process = subprocess.run(["mpirun", "-np", f"{nranks}"] + args, capture_output=True)