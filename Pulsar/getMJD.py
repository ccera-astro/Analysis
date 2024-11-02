# analyze raw gnuradio file 

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt 
import socket 
import json 
import runtempo as rt
import scipy.stats as stats

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    parser.add_argument("-r","--raw",action="store_true",help="Read raw data file.")
    return parser.parse_args()

def getFileName(args) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/J0332+5434/" 
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/" 
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
    print("In downsample(): nIn={0:d} nOut={1:d} nOut*n={2:d}".format(nIn,nOut,n*nOut))
    y = np.zeros(nOut)  
    for j in range(nOut) : y[j] = np.sum(x[j*n:(j+1)*n])
    return y

def highpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
    print("In   highpass(): cutoff={0:f} sample_rate={1:f}".format(cutoff,sample_rate))
    sos = scipy.signal.butter(poles, cutoff, 'highpass', fs=sample_rate, output='sos')
    filtered_data = scipy.signal.sosfiltfilt(sos, data)
    return filtered_data

# use lower tail to estimate mean and sigma and the remove samples
# that are more than 5 sigma above the mean 
def denoise(array) :
    p16 = np.percentile(array, 16)
    p50 = np.percentile(array, 50)
    sigma = p50 - p16
    threshold = p50 + 10.*sigma 
    indices = np.where(array > threshold)[0]
    filtered_array = np.delete(array, indices)
    if True : 
        la, lf = len(array), len(filtered_array)
        print("In    denoise(): p16={0:f} p50={1:f} sigma={2:f} fraction removed={3:.4f}%".format(
            p16,p50,sigma,100.*(la-lf)/float(la)))
    return filtered_array, indices 

def plotFilterHist(array, file, args=None) :
    p16 = np.percentile(array, 16)
    p50 = np.percentile(array, 50)
    sigma = p50 - p16
    threshold = p50 + 5.*sigma 
    hist, bins = np.histogram(array, bins=1000, range=(-1.,5.))
    plt.semilogy(bins[:-1], hist, label="All samples")
    plt.title("Histogram for {0:s}".format(file.split('/')[-1]))
    plt.xlabel("Value")
    plt.ylabel("Counts")
    filtered_array, dummy = denoise(array)
    hist, bins = np.histogram(filtered_array, bins=1000, range=(-1.,5.))
    plt.semilogy(bins[:-1], hist, 'r.', label="Filtered samples")
    plt.text(1.5,1000.,"Median={0:.3f} Sigma={1:.3f} Threshold={2:.3f}".format(p50,sigma,threshold))
    plt.legend()
    if True or args == None or not args.printPlot :
        plt.show()
    else  :  
        plotFileName = "./plots/FilterHist_{0:s}.png".format(getDay(args))
        plt.savefig(plotFileName,format='png')
        plt.clf()

def plotTimeSeries(times,array, file, args=None) :
    plt.subplot() 
    plt.plot(times,array,'r-',label="All samples")
    plt.title("Time Series for {0:s}".format(file.split('/')[-1]))
    plt.xlabel("Time (s)")
    plt.ylabel("Sample Value")
    filtered_array, bad_elements = denoise(array)
    filtered_times = np.delete(times,bad_elements)
    plt.plot(filtered_times,filtered_array,'b.', label="Filtered samples")
    plt.legend()
    if True or args == None or not args.printPlot :
        plt.show()
    else  :  
        plotFileName = "./plots/TimeSeries_{0:s}.png".format(getDay(args))
        plt.savefig(plotFileName,format='png')
        plt.clf()
    
def detectSignal(times, array, period, phase0 = 0., nBins = 128) :
    hist_array, count_array = np.zeros(nBins), np.zeros(nBins)
    times = times - phase0*period
    print("times[0:10]={0:s} \nperiod={1:f}".format(str(times[0:10]),period))
    for i, t in enumerate(times) :
        tm = t % period 
        j = int(nBins*tm/period) 
        hist_array[j] += array[i]
        count_array[j] += 1 
    hist_array = np.divide(hist_array,count_array)
    return hist_array 

def getSigmaArray(x) :
    iMax = np.argmax(x)
    # remove phase bins within +/- 3 of peak 
    iLow, iHigh = max(0,iMax-3), min(len(x),iMax+3)
    y = np.delete(x,range(iLow,iHigh))
    mean, sigma = np.mean(y), np.std(y)
    xx = (x-mean)/sigma 
    return xx 

def getDay(args) :
    f = args.fileName 
    #print("f={0:s}".format(f))
    date = f.strip(".raw")
    day = date[-4:]
    #print("day={0}".format(day))
    return day

def time2phase(time, best_coeff):	#Converts a set of times (mjd) into phases for the specified pulsar
    # these define the indices	
    TMID = 1
    RPHASE = 3
    F0 = 4
    coeff1 = 6
    coeff2 = 7
    coeff3 = 8
    coeff4 = 9

    dt = (time - best_coeff[TMID]) * 1440.0

    # Construct the polynomial coefficients	
    pcoeff = (best_coeff[RPHASE] + best_coeff[coeff1], 60.0 * best_coeff[F0] + best_coeff[coeff2], best_coeff[coeff3], best_coeff[coeff4]) 
    phase = np.polynomial.polynomial.polyval(dt, pcoeff)
    return phase

def bindata( phase, value, NumBins):	# Averages 'value' into 'NumBins' bins based on the 'phase' values
    mphase = phase - np.floor(phase)
    iphase = np.floor(mphase * NumBins)
    bdata = np.zeros((NumBins))
    bnum = np.zeros((NumBins))

    (bdata,d1,d2) = stats.binned_statistic(iphase, value, statistic='sum', bins=NumBins, range = [0,NumBins-1])
    (bnum,d1,d2) = stats.binned_statistic(iphase, value, statistic='count', bins=NumBins, range = [0,NumBins-1])

    return (bdata, bnum)

# begin execution here

args = getArgs() 

day = args.base_name 
base_name = getFileName(args)

with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
freq = 1.0e-6*metadata['freq']
try : n_decimate = metadata['decimation_factor']
except KeyError : n_decimate = 250 
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 

# set up time series 
pts_1 = 1000*np.fromfile(base_name + "_1.sum", dtype=np.float32)
pts_2 = 1000*np.fromfile(base_name + "_2.sum", dtype=np.float32)
if len(pts_1) > len(pts_2) : pts_1 = pts_1[:len(pts_2)]
if len(pts_2) > len(pts_1) : pts_2 = pts_2[:len(pts_1)]
power_time_series = pts_1 + pts_2 

nRows = len(power_time_series) 
times = np.linspace(0.,nRows*t_fft,nRows) 
MJD = (metadata['t_start'] / 86400.) + 40587. 
MJDs = MJD + times/86400. 
startMJD, stopMJD = MJD, MJDs[-1]
midMJD = 0.5*(startMJD + stopMJD)
print("Pulsar:{0:s}".format(metadata['target']))
print("MJD: start={0:f} stop={1:f} mid={2:f}".format(startMJD,stopMJD,midMJD))



