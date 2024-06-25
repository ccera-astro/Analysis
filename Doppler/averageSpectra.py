# average the spectra in a raw data file 

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt 
import socket 
import json 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    return parser.parse_args()

def getFileName(args) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Doppler/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
    return data_dir + args.base_name 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

def getMetaData(file) :
    import json
    print("file={0:s}".format(file))
    with open(file) as json_file:
        dict = json.load(json_file)
        print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

# begin execution here

args = getArgs() 
base_name = getFileName(args)
#with open(base_name + ".json") as json_file : metadata = json.load(json_file)

metadata = getMetaData(base_name + ".json")
f_sample = metadata['srate']
fft_size = metadata['fft_size']

for chan in [1,2] :
    data, rows, cols = getData(base_name + "_{0:d}.raw".format(chan),fft_size)
    print("After getData Chan {0:d}: rows={1:d} cols={2:d}".format(chan,rows,cols))

    # reshape array into a series of row 
    data = np.reshape(data, (rows,cols))   
    print("len(data[0])={0:d}".format(len(data[0])))

    avg_data = np.mean(data,0)
    print("len(avg_data)={0:d}".format(len(avg_data)))

    avg_data.tofile(base_name + "_{0:d}.avg".format(chan))








