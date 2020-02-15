#!/usr/bin/python
from taskgen import StaffordRandFixedSum
import sys
import os

# @n: number of tasks in the task set
# @U: total utilization of n tasks.
# @lower: lower-bound of task utilization
# @upper: upper-bound of task utilization
# Return: a list of task utilizations
def generate_utils(n, U, lower, upper):
    # Generate list of utilizations in [0.0, 1.0].
    x = StaffordRandFixedSum(n, float(U-n*lower)/(upper-lower), 1)
    # Convert them into the range [lower, upper].
    X = x*(upper-lower) + lower
    
    return X[0]

# Inputs are the same as generate_utils() except the last one is the taskset number.
# Generate utils for the tasks in a task set and write to utils.txt file.
# This is aimed to generate heavy tasks only.
def gen_single_tset_utils(n, m, frac, lower, upper, i, path):
    path += "/taskset" + str(i)
    try:
        if not os.path.exists(path):
            os.mkdir(path)
        f = open(path + "/utils.txt", 'w+')
        utils = generate_utils(n, m * frac, lower, upper)
        for util in utils:
            f.write(str(util) + "\n")
        f.close()
    except OSError:
        print "ERROR: Could not create folder", path, "!!"

# Main func for generating set of heavy tasks.
def main_heavy_only():
    if len(sys.argv) != 8:
        print "Usage: Program <NumTasks> <NumCores> <TotalNormUtil> <LowerUtil> <UpperUtil> <NumTasksets> <OutFolder> !!"
        exit(1)
    n = int(sys.argv[1])
    m = int(sys.argv[2])
    frac = float(sys.argv[3])
    lower = float(sys.argv[4])
    upper = float(sys.argv[5])
    no_tasksets = int(sys.argv[6])
    path = sys.argv[7]
    for i in range(1, no_tasksets+1):
        gen_single_tset_utils(n, m, frac, lower, upper, i, path)
        
        
# Generate utils for both heavy and light tasks in a task set. Write the utils to
# files, one for heavy tasks and one for light tasks.
# @m: Total number of processors.
# @util_frac: Normalized total utilization.
# @heavy_ratio: The ratio of total utilization occupied by heavy tasks.
# @n_heavy: Number of heavy tasks.
# @min_heavy, @upper_heavy: Lower-bound and upper-bound for heavy task utilization.
# @n_light: Number of light tasks.
# @min_light, @max_light: Lower-bound and upper-bound for light task utilization.
# @i: ID of the task set being generated.
# @path: Path to the output folder.
# @gen_density: True if it is generating densities.
def gen_full_taskset(m, util_frac, heavy_ratio, n_heavy, min_heavy, max_heavy, n_light, min_light, max_light, i, path, gen_density):
    total_util = m * util_frac
    heavy_total_util = heavy_ratio * total_util
    light_total_util = total_util - heavy_total_util
    path += "/taskset" + str(i)
    try:
        if not os.path.exists(path):
            os.mkdir(path)
        heavy_filename = "heavy_utils.txt"
        if gen_density == True:
            heavy_filename = "heavy_densities.txt"
        heavy_f = open(path + "/" + heavy_filename, 'w+')
        if n_heavy > 0:
            heavy_utils = generate_utils(n_heavy, heavy_total_util, min_heavy, max_heavy)
            for util in heavy_utils:
                heavy_f.write(str(util) + "\n")
        heavy_f.close()

        light_filename = "light_utils.txt"
        if gen_density == True:
            light_filename = "light_densities.txt"
        light_f = open(path + "/" + light_filename, 'w+')
        if n_light > 0:
            light_utils = generate_utils(n_light, light_total_util, min_light, max_light)
            for util in light_utils:
                light_f.write(str(util) + "\n")
        light_f.close()
    except OSError:
        print "ERROR: Could not create folder", path, "!!"
        

# Main func for generating task sets with both heavy and light tasks.
# @gen_density: True if it is generating densities.
def main_heavy_and_light(gen_density):
    if len(sys.argv) != 12:
        print "Usage: Program <NumCores> <NormalizedUtil> <HeavyRatio> <NumHeavy> <MinHeavyUtil> \
<MaxHeavyUtil> <NumLight> <MinLightUtil> <MaxLightUtil> <NumTasksets> <OutFolder>!!"
        exit(1)
    m = int(sys.argv[1])
    frac = float(sys.argv[2])
    heavy_ratio = float(sys.argv[3])
    n_heavy = int(sys.argv[4])
    lower_heavy = float(sys.argv[5])
    upper_heavy = float(sys.argv[6])
    n_light = int(sys.argv[7])
    lower_light = float(sys.argv[8])
    upper_light = float(sys.argv[9])
    no_tasksets = int(sys.argv[10])
    path = sys.argv[11]
    for i in range(1, no_tasksets+1):
        gen_full_taskset(m, frac, heavy_ratio, n_heavy, lower_heavy, upper_heavy, n_light, lower_light, upper_light, i, path, gen_density)

if __name__ == "__main__":
    #main_heavy_only()
    main_heavy_and_light(True)
