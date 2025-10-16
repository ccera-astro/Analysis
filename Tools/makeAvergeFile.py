# average the spectra in a raw data file 

import numpy as np
import socket 
import glob

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    return parser.parse_args()

def getFileName(args) :
    base_names = [] 
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Doppler/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"

    glob_string = data_dir + args.base_name  
    print("glob_string={0:s}".format(glob_string))
    files = glob.glob(glob_string)
    
    print("files = {0:s}".format(str(files)))
    for file in files : 
        if "json" in file : base_names.append(file.strip(".json"))

    print("\n\n base_names = {0:s}".format(str(base_names)))
    return base_names  

def getData(file,fft_size) :
    print("Reading from file = {0:s}".format(file))
    vals = np.fromfile(file, dtype=np.float32)
    print("len(vals)={0:d}".format(len(vals)))
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
base_names = getFileName(args)

for base_name in base_names :
    print("\n** base_name={0:s}".format(base_name))
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








