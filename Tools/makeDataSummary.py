#!/usr/bin/env python3

# read all JSON files and make a spreadsheet summarizing the results

import numpy as np
import json 
import socket
import glob 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    return parser.parse_args()

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

# begin execution

args = getArgs()

if args.data_dir :
    data_dir = args.data_dir 
else :
    data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Scan/" 
    if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"

files = glob.glob(data_dir + "*.json").sort() 

outLines = [] 
for file in files :
    line = file.split('/')[-1]
    metadata = getMetaData(file)
    line += metadata['run_mode'] + ", "
    line += metadata['target'] + ", "
    line += metadata['run_type'] + ", "
    line += metadata['fft_size'] + ", "
    line += "{0:.6f}, ".format(float(metadata['t_sample']))
    line += "{0:.1f}, ".format(1.e-6*float(metadata['freq']))
    line += "{0:.3f}, ".format(1.e-6*float(metadata['srate']))
    line += "\n"
    outLines.append(line)

with open("data_summary.csv",'w') as file : file.writelines(outLines)




