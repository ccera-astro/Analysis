# prepare a set of pulsar transit scan analysis jobs  

import glob 
import json

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--start",default="2024-12-10",help="start date")
    parser.add_argument("--search",default="",help="search term")
    return parser.parse_args()

args = getArgs() 
start = args.start 
data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/" 

outLines = [] 
files = glob.glob("/mnt/c/Users/marlow/Documents/CCERA/data/{0:s}*.json".format(args.search))
print("Number of files={0:d}".format(len(files)))
for file in  files :
    base_name = file.split("/")[-1].strip(".json")
    print("Checking: {0:s}".format(base_name))
    tt = file.split("/")[-1].split(".")[0]
    if tt < start : continue 
    with open(file) as json_file : metadata = json.load(json_file)
    if not metadata['run_mode'] == 'pulsar' : continue 
    if not metadata['target']  == 'J0332+5434' : continue 
    if metadata['az'] > 1. and metadata['az'] < 359. : continue 
    if abs(metadata['alt'] - 80.8) > 1. : continue  
    print("   Good file: base_name={0:s}".format(base_name))
    cmd = "python anaPulsar.py -b {0:s} -m phase --no_plot --high_pass_off --no_roll --nPhaseBins 200".format(base_name)
    outLines.append(cmd +"\n")

open("runAnaPulsar.csh",'w').writelines(outLines)



