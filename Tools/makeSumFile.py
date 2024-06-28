#!/usr/bin/env python3

# convert .raw data format to .sum data format by summing over columns in the 
# data array to make a time series of single values 

import numpy as np
import json 
import socket

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-25-2100",help="File(s) to be analyzed.")
    return parser.parse_args()

def getFileName(args) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Scan/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
    return data_dir + args.base_name 

def getData(file,fft_size) :
    print("Reading from file: {0:s}".format(file))
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size)
    vals = np.reshape(vals, (rows,cols))   
    return vals, rows, cols

# begin execution

args = getArgs()
base_name = getFileName(args)
with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
n_decimate = metadata['decimation_factor']
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 

for chan in [1,2] :
    file = base_name + "_{0:d}.raw".format(chan)
    data, nRows, nCols = getData(file,fft_size)
    print("Read {0:d} {1:d}-channel spectra from {2:s}".format(nRows,nCols,file))
    power = np.sum(data,1)
    file_out = base_name + "_{0:d}.sum".format(chan)
    with open(file_out,'w') as file : power.tofile(file)



