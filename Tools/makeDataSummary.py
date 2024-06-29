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
        #print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

# begin execution

args = getArgs()

if args.data_dir :
    data_dir = args.data_dir 
else :
    data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Scan/" 
    if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"

print("data_dir={0:s}".format(data_dir))
files = glob.glob(data_dir + "*.json")
print("files={0:s}".format(str(files)))

outLines = [" file, mode, target, type, fft_size, t_sample, freq, srate \n"]
for file in files :
    print("Reading file={0:s}".format(file))
    line = file.split('/')[-1] + ", "
    metadata = getMetaData(file)
    try : line += metadata['run_mode'] + ", "
    except KeyError : line += "Unknown, "
    try : line += metadata['target'] + ", "
    except KeyError : line += "Unknown, "
    line += metadata['run_type'] + ", "
    line += "{0:d}, ".format(metadata['fft_size']) 
    line += "{0:.6f}, ".format(float(metadata['t_sample']))
    line += "{0:.1f}, ".format(1.e-6*float(metadata['freq']))
    line += "{0:.3f} ".format(1.e-6*float(metadata['srate']))
    line += "\n"
    outLines.append(line)

with open("data_summary.csv",'w') as file : file.writelines(outLines)




