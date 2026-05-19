# test the numpy roll function

import numpy as np
import glob 
import json

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--base_name",default=None,help="File basename to be analyzed, e.g. 2026-01-03-1124")
    parser.add_argument("--data_dir",default="./",help="data directory")
    parser.add_argument("--start_time",default=None,help="Start time, e.g. 2025-01-01")
    parser.add_argument("--stop_time", default=None,help="Stop time, e.g. 2025-01-01")
    parser.add_argument("--target",default="J0332+5434",help="Target pulsar")
    parser.add_argument("--az",type=float,default=0.0,help="Azimuth")
    parser.add_argument("--alt",type=float,default=80.75,help="Altitude")
    parser.add_argument("-r","--roll_mult",type=float,default=1.0,help='Roll multiplier.')
    return parser.parse_args()

def getData(file,fft_size) :
    print("Reading from file: {0:s}".format(file))
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size)
    vals = np.reshape(vals, (rows,cols))   
    return vals, rows, cols

def getDelays(meta_data) :
    f0 = 1.e-6*meta_data['freq']
    dF = 1.e-6*meta_data['srate']
    fMin = f0 - 0.5*dF
    fMax = f0 + 0.5*dF
    print("getDelays: fMin={0:e} fMax={1:e}".format(fMin,fMax))
    freqs = np.linspace(fMin,fMax,meta_data['fft_size'])
    DM = 26.7641          # pc cm^-3, approximate ATNF value for J0332+5434
    # Cold-plasma dispersion formula: Delta t(ms) = K_DM * DM * (1/f1^2 - 1/f2^2), with f in MHz
    K_DM = 4.148808e6     # ms MHz^2 pc^-1 cm^3
    delays  = 0.001* K_DM * DM * (1.0/f0**2 - 1.0/freqs**2)  # delay time in seconds 
    print("dts={0:s}".format(str(delays)))
    if False :
        import matplotlib.pyplot as plt 
        plt.plot(freqs,1.e6*delays,'b-')
        plt.xlabel("f (MHz)")
        plt.ylabel("Delay (us)")
        plt.show() 
    return delays 

def roll_columns(arr, k):
    """
    Roll each column of a 2D NumPy array independently.

    Parameters
    ----------
    arr : numpy.ndarray
        Input array with shape (M, N)

    k : array-like
        Length-N array of integer shifts.
        k[j] specifies how much to roll column j.
        Positive values roll downward;
        negative values roll upward.

    Returns
    -------
    numpy.ndarray
        Copy of the array with all columns rolled.
    """
    arr = np.asarray(arr)
    k = np.asarray(k)

    if arr.ndim != 2:
        raise ValueError("arr must be a 2D array")

    M, N = arr.shape

    if len(k) != N:
        raise ValueError("k must have length equal to number of columns")

    result = np.empty_like(arr)

    for j in range(N):
        result[:, j] = np.roll(arr[:, j], k[j])

    return result

def print_file(data,nRows,nCols) :
    for i in range(nRows) :
        line = "row={0:3d}".format(i)
        for j in range(nCols) :
            line += "{0:8.3f}".format(1000.*data[i][j])
        print(line)
    return

# Begin execution here 

args = getArgs()

# determine which file(s) to process 
if args.base_name :    
    files = [args.data_dir + args.base_name] 
else :
    files = glob.glob(args.data_dir + "*.json") 
    print("{0:d} files found in {1:s}".format(len(files),args.data_dir))
    # apply filters 
    filtered_files = [] 
    for file in files :
        with open(file) as json_file : metadata = json.load(json_file)
        file_time = file.split("/")[-1].strip(".json")
        if args.start_time and file_time < args.start_time : continue 
        if args.stop_time  and file_time > args.stop_time : continue
        if not metadata["target"] == args.target : continue 
        diff = abs(metadata["az"] - args.az)
        if not (diff < 1. or diff > 359.) : continue 
        filtered_files.append(file[:-5]) 
    print("{0:d} good files found.".format(len(filtered_files)))
    files = filtered_files[:]

print("files={0:s}".format(str(files)))

for file in files :
    with open(file+".json") as json_file : metadata = json.load(json_file)
    delays = getDelays(metadata)
    delays *= args.roll_mult 
    rolls = np.round(delays/metadata["t_sample"]).astype(int)
    print("rolls={0:s}".format(str(rolls)))
    print("t_sample={0:.3f} us".format(1.e6*metadata["t_sample"]))
    for chan in [1,2] :
        in_file = file + "_{0:d}.raw".format(chan)
        data, nRows, nCols = getData(in_file,metadata['fft_size'])
        print("Read {0:d} {1:d}-channel spectra from {2:s}".format(nRows,nCols,file))
        print("Before roll:")
        print_file(data,10,10)
        
        data = roll_columns(data, rolls)
        print("After roll:")
        print_file(data,10,10)

        out_file = in_file.replace(".raw",".dsp")
        power = np.sum(data,1)
        print("outfile={0:s}".format(out_file))
        #with open(out_file,'w') as ff : power.tofile(ff)
