# 
# 1) loop through JSON files in data directory, selecting those that 
#    meet the desired criterion.
# 2) create a soft link in temp_soft_link directory
# 3) copy all files from the temp_soft_link directory to laptop 

import glob
import json 
import os 
import socket 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--start",default="2025-01-01",help="start date")
    parser.add_argument("--search",default="",help="search term")
    return parser.parse_args()

args = getArgs() 
start = args.start 
data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/"
if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"

search_term = "{0:s}{1:s}*.json".format(data_dir,args.search)
print("Searching: {0:s}".format(search_term))
files = glob.glob(search_term)
print("Total of {0:d} files found.".format(len(files)))

for file in files :
    print("file={0:s}".format(file))
    tt = file.split("/")[-1].split(".")[0]
    if tt < start : continue 
    print("  check metadata")
    with open(file) as json_file : metadata = json.load(json_file)
    if not metadata['run_mode'] == 'pulsar' : continue 
    if not metadata['target']  == 'J0332+5434' : continue 
    print("   good file")
    os.system("ln -s {0:s} {1:s}temp_soft_link/.".format(file,data_dir))
    os.system("ln -s {0:s} {1:s}temp_soft_link/.".format(file.replace(".json","_1.sum"),data_dir))
    os.system("ln -s {0:s} {1:s}temp_soft_link/.".format(file.replace(".json","_2.sum"),data_dir))

