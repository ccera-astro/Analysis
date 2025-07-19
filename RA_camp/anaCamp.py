#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import time
import socket
import glob 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("-c","--channel",type=int,default=1,help="Radio channel (1 or 2)")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    parser.add_argument("-n","--name",default=None,help="Group name.")
    parser.add_argument("--start_time",default="2025-06-28-00",help="Start time yyyy-mm-dd-hh")
    #parser.add_argument("--stop_time", default="2025-06-29-00",help="Stop time yyyy-mm-dd-hh")
    parser.add_argument("-g","--gain_factor",type=float,default=66000.,help="Gain factor")
    #parser.add_argument("-r","--removeSpurs",action='store_true',help="Remove spurious peaks from spectrum.")
    return parser.parse_args()

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

def getFreqs(meta_data) :
    dF = 1.e-6*meta_data['srate']
    fMin = 1.e-6*meta_data['freq'] - 0.5*dF
    fMax = 1.e-6*meta_data['freq'] + 0.5*dF
    print("getFreqs: fMin={0:e} fMax={1:e}".format(fMin,fMax))
    freqs = np.linspace(fMin,fMax,meta_data['fft_size'])
    return freqs

def getVelocities(f) :
    f0, c = 1420.41, 3.0e5
    v = c*(f/f0 - 1.)
    return v

def getGain(fName) :
    gain = np.fromfile(fName,dtype=float) 
    return gain

def removeSpurs(series,x) :
    spurChannels = range(613,618)    
    x = np.delete(x,spurChannels)
    return x

def subtractBackground(vDoppler,power) :
    N = len(vDoppler) 
    v1, v2 = -200., 200.
    i1 = np.searchsorted(vDoppler,v1)
    i2 = np.searchsorted(vDoppler,v2) - 1 
    #print("In subtractBackground(): i1={0:d} i2={1:d} len(power)={2:d}".format(i1,i2,len(power)))
    p1, p2 = power[i1], power[i2]
    b = p1
    m = (p2-p1)/(i2-i1)
    #print("i1={0:d} i2={1:d} p1={2:f} p2={3:f} m={4:f}".format(i1,i2,p1,p2,m))
    x = np.array(range(N))
    background = b - m*(i1 - x) 
    power -= background
    return power

def getFileName(args) :
    inp = 1 
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/RA_camp/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
        elif "student" in socket.gethostname().lower() : data_dir = "/home/student/data/"
    print("In getFileName(): data_dir={0:s} base_name={1:s}".format(data_dir,args.base_name))
    if args.name :
        for line in open("Channel_lookup.csv",'r').readlines()[1:] :
            vals = line.strip().split(',')
            inp = int(vals[1])
            if vals[2].lower() == args.name.lower() :
                query = "{0:s}Ch{1:02d}*.json".format(data_dir,int(vals[0]))
                files = glob.glob(query)
                files.sort() 
                print("In getFileName(): len(files)={0:d}".format(len(files)))
                for file in files :
                    file_time = file.split("/")[-1].split("_")[1][0:13]
                    if file_time >= args.start_time :
                        return_name = data_dir + file.split("/")[-1].strip(".json")
                        return return_name, inp 
                print("Valid file not found for args.name={0:s}".format(args.name))
                return None 
    return data_dir + args.base_name , inp

def getFiles(args) :
    import glob
    return glob.glob('./data/{0:s}/{0:s}_???.dat'.format(args.series))

# fit spectrum to Chebyshev polynomial 
# restrict range of fit to |vDoppler| > vSignal 
def fitBackground(vDoppler,power,n,vSignal) :
    weights = np.ones_like(vDoppler)
    for i in range(len(vDoppler)) :
        if abs(vDoppler[i]) < vSignal : weights[i] = 1.e-6 
    series = np.polynomial.chebyshev.Chebyshev.fit(vDoppler, power, n, w=weights)
    background = series(vDoppler) 
    return background

def downSample(x,n):
    nIn = len(x)
    nOut = int(nIn/n)
    print("In downsample(): nIn={0:d} nOut={1:d} nOut*n={2:d}".format(nIn,nOut,n*nOut))
    y = np.zeros(nOut)  
    for j in range(nOut) : y[j] = np.sum(x[j*n:(j+1)*n])
    y /= n 
    return y

def filterSpectrum(power,freqs,vDoppler, filter_on = True) :
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
    vDoppler = downSample(vDoppler,16)
    freqs = downSample(freqs,16)

    return power, freqs, vDoppler 


# Begin execution here

args = getArgs()
chan = args.channel 
gain = args.gain_factor 
base_name, inp = getFileName(args)

print("base_name={0:s}".format(base_name))
meta_data = getMetaData(base_name + ".json")

data_type = 'avg'
if data_type == 'raw' :
    data, rows, cols = getData(base_name + "_{0:d}.raw".format(chan),meta_data['fft_size'])
    print("After getData Chan {0:d}: rows={1:d} cols={2:d}".format(chan,rows,cols))
    fft_size = cols 

    # reshape array into a series of row 
    data = np.reshape(data, (rows,cols))   
    print("len(data[0])={0:d}".format(len(data[0])))
    power = np.mean(data,0)
    print("len(power)={0:d}".format(len(power)))
else :
    power = np.fromfile(base_name +"_{0:d}.avg".format(inp), dtype=np.float32)

freqs = getFreqs(meta_data) 
print("freqs[0]={0:.2f} freqs[-1]={1:.2f}".format(freqs[0],freqs[-1]))
vDoppler = getVelocities(freqs)

if True :    
    #plt.plot(freqs,power,'b-',label="Chan 1")
    plt.plot(freqs,power,'b.')
    plt.title("Raw spectrum: {0:s}".format(base_name.split("/")[-1]))
    plt.xlabel("f (MHz)")
    plt.ylabel("PSD (AU)")
    #plt.ylim(0.,1.1*np.max(power)) 
    plt.show() 

power, freqs, vDoppler = filterSpectrum(power,freqs,vDoppler, filter_on = True)

if True : 
    plt.plot(freqs,power)
    plt.title("Filtered spectrum: {0:s}".format(base_name.split("/")[-1]))
    plt.xlabel("f (MHz)")
    plt.ylabel("PSD (AU)")
    plt.show() 

vMin, vMax = -300., 300.
i1 = np.searchsorted(vDoppler,vMin)
i2 = np.searchsorted(vDoppler,vMax)
print("i1={0:d} i2={1:d}".format(i1,i2))
freqs = freqs[i1:i2]
vDoppler = vDoppler[i1:i2]
power = power[i1:i2]
background = fitBackground(vDoppler,power,5,150.)

subtract_background = True 
if subtract_background : 
    plt.plot(vDoppler,power-background,'r.')
    plt.plot([-300.,300.],[0., 0.],'g-')
    yMax = np.max(power-background)
    plt.text(-200.,0.8*yMax,'Chan {0:d}'.format(chan))
else :
    plt.plot(vDoppler,power,'r.')
    plt.plot(vDoppler,background,'g-')
    yMax = np.max(power)
    plt.text(-200.,0.8*yMax,'Chan {0:d}'.format(chan))

plt.title("Spectrum with gain factor {0:s}".format(base_name.split("/")[-1]))
plt.xlabel('Doppler Velocity (km/s)')
plt.ylabel('Power (K)')
plt.legend() 
plt.show()