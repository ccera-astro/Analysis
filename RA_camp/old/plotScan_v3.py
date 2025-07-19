#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import glob
from matplotlib.backends.backend_pdf import PdfPages

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--plot",action="store_true",help="Plot to screen.")
    parser.add_argument("--start_time",default="2025-06-28-00",help="Start time yyyy-mm-dd-hh")
    parser.add_argument("--stop_time", default="2025-06-29-00",help="Stop time yyyy-mm-dd-hh")
    return parser.parse_args()

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        #print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

def downSample(x,n):
    nIn = len(x)
    nOut = int(nIn/n)
    #print("In downsample(): nIn={0:d} nOut={1:d} nOut*n={2:d}".format(nIn,nOut,n*nOut))
    y = np.zeros(nOut)  
    for j in range(nOut) : y[j] = np.sum(x[j*n:(j+1)*n])
    y /= n 
    return y

def filterSpectrum(power,filter_on = True) :
    if filter_on :
        fft_size, slice_size, q = len(power), 128, 0.90
        for i in range(0,fft_size,slice_size) :
            x = np.linspace(0.,float(slice_size),slice_size)
            slice = power[i:i+slice_size]
            q_value = np.quantile(slice,q)
            indices_1 = np.where(slice > q_value)[0]
            filtered_slice = np.delete(slice, indices_1)
            filtered_x = np.delete(x,indices_1)
            
            # fit the filtered array 
            coeffs = np.polyfit(filtered_x, filtered_slice, 1)
            fit_fcn = coeffs[1] + coeffs[0]*filtered_x

            # calculate the residuals with respect to the fit
            resid = filtered_slice - fit_fcn
            rms = np.std(resid)

            # redo the filter 
            fit_fcn = coeffs[1] + coeffs[0]*x
            resid = np.abs((slice-fit_fcn)/rms)
            indices = np.where(resid > 3.0)[0]
            filtered_slice = np.delete(slice, indices)
            filtered_x = np.delete(x,indices)

            # refit 
            coeffs = np.polyfit(filtered_x, filtered_slice, 1)
            fit_fcn = coeffs[1] + coeffs[0]*filtered_x

            for index in indices :
                ff = coeffs[1] + coeffs[0]*index 
                slice[index] = np.random.normal(ff,rms)

            power[i:i+slice_size] = slice 

    power = downSample(power,16)

    return power

def getFiles(args,chan) :
    # Get a list of all files.  Make sure that they are ordered by time
    all_files = glob.glob("data/*.json")
    all_files.sort() 

    # select files that match time, channel, and input 
    files = []  
    for file in all_files :
        file_chan = int(file.split("/Ch")[1][:2])
        #print("file={0:s} file_chan={1:d}".format(file,file_chan))
        if file_chan == chan :
            file_time = file.split("_")[1][0:13] 
            if file_time >= args.start_time and file_time < args.stop_time :
                #print("Keeping file={0:s}".format(file))
                files.append(file)
    return files 

def makeKey(chan,inp) :
    return "Ch{0:02d}_{1:d}".format(chan,inp)

args = getArgs() 

# build dictionary of names
channel_lookup = {} 
for line in open("Channel_lookup.csv").readlines()[1:] :
    vals = line.strip().split(',')
    ch, inp, name = int(vals[0]), int(vals[1]), vals[2]
    channel_lookup[makeKey(ch,inp)] = vals[2]
print("channel_lookup={0:s}".format(str(channel_lookup)))

pdf = PdfPages("Scan_{0:s}.pdf".format(args.start_time))

for chan in range(2) :
    files = getFiles(args,chan)
    #print("files={0:s}".format(str(files)))
    base_name_0 =  files[0].strip(".json")
    meta_data = getMetaData(base_name_0 + ".json")

    for inp in [1,2] : 
        times, total_power = [], []
        for row, file in enumerate(files) :
            base_name = "./" + file.strip(".json")
            data_type = 'avg'
            if data_type == 'raw' :
                data, rows, cols = getData(base_name + "_{0:d}.raw".format(inp),meta_data['fft_size'])
                # reshape array into a series of row 
                data = np.reshape(data, (rows,cols))   
                power = (100./75000.)*np.mean(data,0)
                power = np.sum(data[0])
            else :
                power = np.fromfile(base_name +"_{0:d}.avg".format(inp), dtype=np.float32)
                power = 0.00094*np.sum(power)

            hour = int(base_name.split("-")[-1][0:2]) 
            min = int(base_name.split("-")[-1][2:4])
            #print("base_name={0:s} hour={1:d} min={2:d}".format(base_name,hour,min))
            times.append(hour + min/60.)
            
            total_power.append(np.sum(power))

        fig1 = plt.figure() 
        plt.plot(times,total_power,'r.')
        plt.xlabel("t (hrs UTC)")
        plt.ylabel("Total power (K)")
        ky = makeKey(chan,inp)
        plt.title("Power vs. time: {0:s} {1:s}\n  {2:s} to {3:s}".format(ky,channel_lookup[ky],args.start_time,args.stop_time))
        if args.plot : plt.show() 
        pdf.savefig(fig1)

pdf.close() 
