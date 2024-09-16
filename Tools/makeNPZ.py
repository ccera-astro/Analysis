#!/usr/bin/env python3

# convert .sum data format to .npz data format that is compatible 
# with PHY312 lab writeup

import numpy as np
import json 
import socket
import astropy 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-25-2100",help="File(s) to be analyzed.")
    parser.add_argument("--down_sample",type=int,default=1,help="Down sample (averaging) factor")
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

def downSample(x,n):
    nIn = len(x)
    nOut = int(nIn/n)
    print("In downsample: nIn={0:d} nOut={1:d} nOut*n={2:d}".format(nIn,nOut,n*nOut))
    y = np.zeros(nOut)  
    for j in range(nOut) : y[j] = np.sum(x[j*n:(j+1)*n])
    return y

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

# sum the two channels 
p1 = np.fromfile(base_name+"_1.sum", dtype=np.float32)
p2 = np.fromfile(base_name+"_2.sum", dtype=np.float32)
print("len(p1)={0:d} len(p2)={1:d}".format(len(p1),len(p2)))

# if the lengths are different choose the shorter one
if len(p1) > len(p2)   : p1 = p1[:len(p2)]
elif len(p1) < len(p2) : p2 = p2[:len(p1)]
p = p1 + p2 

# down sample, if requested 
if args.down_sample > 1 : p = downSample(p,args.down_sample)
    
nSamples = len(p)
t0 = metadata['t_start']
run_time_1 = metadata['t_sample']*nSamples*args.down_sample  
run_time_2 = metadata['run_time']
print("t0={0:f} nSamples={1:d} run_time_1={2:.3f} run_time_2={3:.3f}".format(t0,nSamples,run_time_1,run_time_2))

firstMJD = (t0 / 86400) + 40587 
lastMJD = ((t0 + run_time_1)/86400) + 40587 
mjd = np.linspace(firstMJD,lastMJD,nSamples)
#print("mjd={0:s}".format(str(mjd)))

dict = {"a":"A", "b":"B"}
outFile = base_name + ".npz" 
print("Saving file to: {0:s}".format(outFile))
np.savez(outFile,mjd=mjd,data=p,fileInfo=dict)

exit()


