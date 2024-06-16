# analyze raw gnuradio file 

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
    #parser.add_argument("-f","--fileName",default="../data/J0332+5434/20240605.raw",help="File(s) to be analyzed.")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    parser.add_argument("-n","--name",default="J0332+5434",help="Pulsar to be analyzed.")
    parser.add_argument("--p0",type=float,default=0.7145196,help="Period in s.")
    parser.add_argument("--DM",type=float,default=26.74,help="Dispersion measure.")
    parser.add_argument("--eps",type=float,default=3.0e-5,help="Period Scan Range")
    parser.add_argument("--down_sample",type=int,default=1,help="Down sample (averaging) factor")
    parser.add_argument("--filter_plots",action="store_true",help="Enable plots related to noise filter.")
    parser.add_argument("--printPlot",action="store_true",help="Write plots to file.")
    parser.add_argument("--high_pass_off",action="store_true",help="Turn off high pass filter.")
    return parser.parse_args()

def getFileName(args) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/J0332+5434/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
    return data_dir + args.base_name 

def getData(file,fft_size) :
    print("Reading from file: {0:s}".format(file))
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size)
    vals = np.reshape(vals, (rows,cols))   
    return vals, rows, cols

# average every n elements 
def downSample(x,n): 
    m = len(x) % n  
    if m > 0 : x = x[:-m]
    return np.mean(x.reshape(-1,n), axis=1)

#high pass filter
def highpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
    print("In highpass(): cutoff={0:f} sample_rate={1:f}".format(cutoff,sample_rate))
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
        print("In denoise(): p16={0:f} p50={1:f} sigma={2:f} fraction removed={3:.4f}%".format(
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
    
def detectSignal(times, array, period, nBins = 128) :
    hist_array, count_array = np.zeros(nBins), np.zeros(nBins)
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

def getOffsetTimes(times,day) : 
    ts_times = {'0525':1716653760.24831, '0526':1716739924.41415, '0528':1716912253.08713,
                '0529': 1716998417.45714}
    # take the May 25 time as the zero reference 
    try :
        tOffset = ts_times[day] - ts_times['0525']
        offset_times = times + tOffset
        print("In getOffsetTimes: Day={0:s} tOffset={1:.6f}".format(day,tOffset))
    except KeyError :
        print("In genOffsetTimes(): **KeyError** Day={0:s}\nNo offset added.".format(day))        
        offset_times = times 
    return offset_times 

# begin execution here

args = getArgs() 

day = args.base_name 
base_name = getFileName(args)
with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
try : n_decimate = metadata['decimation_factor']
except KeyError : n_decimate = 250 
n_downsample = args.down_sample 
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 
c_rate /= n_downsample
plot_filter_hist, plot_time_series = args.filter_plots, args.filter_plots 

for i in [1,2] :
    file = base_name + "_{0:d}.raw".format(i)
    data, nRows, nCols = getData(file,fft_size)
    print("Read {0:d} {1:d} spectra from {2:s}".format(nRows,nCols,file))

    power_time_series = 1000.*np.sum(data,1) 
    power_time_series = downSample(power_time_series,n_downsample)

    if not args.high_pass_off :
        # apply a high pass filter to mitigate baseline wandering
        f_cutoff = 0.1
        power_time_series = highpass(power_time_series,f_cutoff,c_rate)

    times = np.linspace(0.,nRows*t_fft,nRows) 
    times = downSample(times,n_downsample)
    if plot_filter_hist : plotFilterHist(power_time_series, file, args=args)
    if plot_time_series : plotTimeSeries(times, power_time_series, file, args=args)

    power_time_series, bad_elements = denoise(power_time_series)
    times = np.delete(times,bad_elements)

    fig, axs = plt.subplots(4, 2, figsize=(8,7))
    fig.suptitle("{0:s} {1:s} Period Scan: eps={2:.3e}".format(args.name,file.split('/')[-1],float(args.eps)), fontsize=16)

    nBins = 100
    periods = np.linspace(args.p0*(1.-args.eps),args.p0*(1.+args.eps),8)

    sigma_mode, best_sigma = True, 0. 
    for i, period in enumerate(periods) :
        phase_hist = detectSignal(times, power_time_series, period, nBins = nBins)        
        sigma_array = getSigmaArray(phase_hist)
        nx, ny = i % 4 , int(i/4)    
        if not sigma_mode :
            mn, mx = np.amin(phase_hist), np.amax(phase_hist)   
            axs[nx,ny].plot(phase_hist)
        else :
            mn, mx = np.amin(sigma_array), np.amax(sigma_array)   
            axs[nx,ny].plot(sigma_array)
            this_sigma = np.max(sigma_array)
            best_string = ""
            if  this_sigma > best_sigma :
                best_period = period 
                best_sigma = np.max(sigma_array)
                best_sigma_array = 1.001*sigma_array 
                best_string = "**"
            print("i={0:d} period={1:.3f} best_period={2:.3f} sigma={3:.3f} best_sigma={4:.3f}{5:s}".format(
                    i,1000.*period,1000.*best_period,this_sigma,best_sigma,best_string))
        axs[nx,ny].text(0.45*nBins,0.2*mn+0.8*mx,"P={0:.3f}ms".format(1000.*period),fontsize=11)

    for ax in axs.flat: ax.set(xlabel='Phase Bin', ylabel='Power (A.U.)')
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat: ax.label_outer()

    if False and args.printPlot :  
        plotFileName = "./plots/PeriodScan_{0:s}.png".format(getDay(args))
        plt.savefig(plotFileName,format='png')
    else :
        plt.show()

    if sigma_mode :
        nPoints = len(best_sigma_array)
        x = np.linspace(0.,1.,nPoints)
        x_err = np.zeros_like(x)
        y_err = np.ones_like(best_sigma_array)
        fig2, axs2 = plt.subplots(figsize=(8,7))
        fig2.suptitle("{0:s} {1:s} Best Period={2:.3f} (ms)".format(args.name,file.split('/')[-1],1000.*best_period), fontsize=16)
        axs2.errorbar(x,best_sigma_array,xerr=x_err,yerr=y_err,color='red',ecolor='black',fmt='o')
        plt.xlabel("Phase",fontsize=16)
        plt.ylabel("Signal/Noise",fontsize=16)
        bs = max(best_sigma_array)
        plt.text(0.1,0.8*bs,'Best sigma={0:.2f}'.format(bs))
        axs2.plot([0.,1.],[0.,0.],'k--')
        axs2.set_xlim(0.,1.)
        if False and args.printPlot :  
            plotFileName = "./plots/BestSigma_{0:s}.png".format(getDay(args))
            plt.savefig(plotFileName,format="png")
        else :
            plt.show()
 








