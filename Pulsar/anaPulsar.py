# analyze raw gnuradio file 

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt 
import socket 
import json 
import runPolyco as rp 
import scipy.stats as stats
import time 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    parser.add_argument("-n","--name",default="J0332+5434",help="Pulsar to be analyzed.")
    parser.add_argument("-m","--mode",default="polyco",help="Analysis mode (polyco, scan, or f0)")
    parser.add_argument("--p0",type=float,default=0.,help="Period in s.")
    parser.add_argument("--f0",type=float,default=0.,help="Frequency in Hz.")
    parser.add_argument("--phase0",type=float,default=0.,help="Reference phase.")
    parser.add_argument("--nPhaseBins",type=int,default=50,help="Number of bins in phase plot")
    parser.add_argument("--nPeriods",type=int,default=100,help="Number of trial periods")
    parser.add_argument("--no_roll",action="store_true",help="Disable roll")
    parser.add_argument("--power_mode",action="store_true",help="Plot power values")
    parser.add_argument("--eps",type=float,default=2.0e-4,help="Period Scan Range")
    parser.add_argument("--down_sample",type=int,default=1,help="Down sample (averaging) factor")
    parser.add_argument("--filter_plots",action="store_true",help="Enable plots related to noise filter.")
    parser.add_argument("--printPlot",action="store_true",help="Write plots to file.")
    parser.add_argument("--high_pass_off",action="store_true",help="Turn off high pass filter.")
    parser.add_argument("--middle_scan",action="store_true",help="Limit data to middle of transit scan.")
    parser.add_argument("--gain",type=float,default=10.,help="Gain (K/counts)")
    parser.add_argument("--no_plot",action="store_true",help="Don't show plots.")
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
    y /= n 
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

def plotFilterHist(array, args=None) :
    p16 = np.percentile(array, 16)
    p50 = np.percentile(array, 50)
    sigma = p50 - p16
    threshold = p50 + 5.*sigma 
    hist, bins = np.histogram(array, bins=1000, range=(-1.,5.))
    plt.semilogy(bins[:-1], hist, label="All samples")
    #plt.title("Histogram for {0:s}".format(file.split('/')[-1]))
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

def plotTimeSeries(times,array, args=None) :
    plt.subplot() 
    plt.plot(times,array,'r-',label="All samples")
    #plt.title("Time Series for {0:s}".format(file.split('/')[-1]))
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
    
def getSigmaArray(x) :
    iMax = np.argmax(x)
    # remove phase bins within +/- 3 of peak 
    iLow, iHigh = max(0,iMax-3), min(len(x),iMax+3)
    y = np.delete(x,range(iLow,iHigh))
    mean, sigma = np.mean(y), np.std(y)
    xx = (x-mean)/sigma 
    return xx 

def getMeanSigma(x) :
    iMax = np.argmax(x)
    # remove phase bins within +/- 3 of peak 
    iLow, iHigh = max(0,iMax-3), min(len(x),iMax+3)
    y = np.delete(x,range(iLow,iHigh))
    mean, sigma = np.mean(y), np.std(y)
    return mean, sigma 

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

#   write results out to JSON file 
def writeResults(results,file_base_name) :
    file_name = './outputs/' + file_base_name.split('/')[-1] + '_out.json'
    with open(file_name, 'w') as fp :
        json.dump(results, fp)
    return

# begin execution here

args = getArgs() 

day = args.base_name 
base_name = getFileName(args)
mode = args.mode 
high_pass_off = args.high_pass_off
if args.power_mode : high_pass_off = True 
gain = args.gain

results = {}       # dictionary of results to be written out 

with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
freq = 1.0e-6*metadata['freq']
try : n_decimate = metadata['decimation_factor']
except KeyError : n_decimate = 250 
n_downsample = args.down_sample 
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 
c_rate /= n_downsample
plot_filter_hist, plot_time_series = args.filter_plots, args.filter_plots 

pname = metadata['target']
nBins = args.nPhaseBins 
sigma_mode = not args.power_mode  

if mode == "scan" : 
    n_downsample = 4 
    p0 = args.p0  
    if not args.p0 > 0. :
        if pname.lower() == 'j0332+5434' : p0 = 0.714519699726
        elif pname.lower() == 'j0953+0755' : p0 = 0.253065164948
    print("Pulsar={0:s} p0={1:.6f}".format(metadata['target'],p0))
    
# set up time series 
pts_1 = 1000*np.fromfile(base_name + "_1.sum", dtype=np.float32)
pts_2 = 1000*np.fromfile(base_name + "_2.sum", dtype=np.float32)
if len(pts_1) > len(pts_2) : pts_1 = pts_1[:len(pts_2)]
if len(pts_2) > len(pts_1) : pts_2 = pts_2[:len(pts_1)]
power_time_series = pts_1 + pts_2 
#power_time_series = pts_1

results['raw_mean'] = float(np.average(power_time_series))
results['raw_sigma'] = float(np.std(power_time_series))
results['noise_temperature'] = gain*results['raw_mean'] 

#n_downsample = int(0.10/t_fft) 
#print("t_fft={0:e} n_downsample={1:d}".format(t_fft,n_downsample))
if n_downsample > 1 : power_time_series = downSample(power_time_series,n_downsample)
nRows = len(power_time_series) 

if args.middle_scan :
    # limit the range of the time series 
    i1, i2 = int(nRows/4), int(3*nRows/4) 
    power_time_series = power_time_series[i1:i2]
    nRows = len(power_time_series)

if not high_pass_off :
    # apply a high pass filter to mitigate baseline wandering
    f_cutoff = 0.1
    power_time_series = highpass(power_time_series,f_cutoff,c_rate)

times = np.linspace(0.,nRows*t_fft*n_downsample,nRows) 

if plot_filter_hist : plotFilterHist(power_time_series, args=args)
if plot_time_series : plotTimeSeries(times, power_time_series, args=args)

if True :
    power_time_series, bad_elements = denoise(power_time_series)
    times = np.delete(times,bad_elements)


results['mode'] = mode 
MJD = (metadata['t_start']/ 86400.) + 40587.
results['MJD'] = MJD 

if mode == "f0" :
    f0 = args.f0
    if f0 < 0.01 :
        print("ERROR: f0 must be >0 in f0 mode.  Exiting.")
        exit()  
    p0 = 1./f0      
    phase = f0*times + args.phase0
    #phase = phase - np.floor(phase) 
    bdata, bnum = bindata(phase, power_time_series, nBins)
    phase_hist = np.divide(bdata,bnum)
    sigma_array = getSigmaArray(phase_hist)
    best_sigma = max(sigma_array)
    results["best_sigma"] = float(best_sigma)

elif mode == "polyco" :
    MJDs = MJD + times/86400. 
    pulsarName = pname 
    if pulsarName[0] == 'J' : pulsarName = pulsarName[1:]
    #coeff = rp.getpolycoeff(MJD, metadata, base_name, file_name="polyco.dat")
    coeff = rp.getpolycoeff(MJD, metadata, base_name)
    #print("MJD={0:f} coeff={1:s}".format(MJD,str(coeff)))
    p0 = 1./coeff[4]
    
    phase = time2phase(MJDs, coeff) + args.phase0 
    average_period = (times[-1]-times[0])/(phase[-1]-phase[0])
    results['period'] = average_period 
    #print("Average period={0:.4f} ms".format(1000.*average_period))
    bdata, bnum = bindata(phase, power_time_series, nBins)
    phase_hist = np.divide(bdata,bnum)
    mean, sigma = getMeanSigma(phase_hist)
    results["hist_baseline"], results["hist_sigma"] = mean, sigma 
    results["best_signal"] = 1000.*gain*(np.max(phase_hist)-mean)
    sigma_array = getSigmaArray(phase_hist)
    best_sigma = max(sigma_array)
    print("polyco mode: best_sigma={0:f}".format(best_sigma))
    results["best_sigma"] = float(best_sigma)
else : 
    print("In scan: p0={0:f}".format(p0))
    nPeriods, nBins = args.nPeriods, args.nPhaseBins 
    periods = np.linspace(p0*(1.-args.eps),p0*(1.+args.eps),nPeriods)
    mapData = np.zeros((nPeriods,nBins))   
    best_sigma = 0. 
    for i, period in enumerate(periods) :
        bdata, bnum = bindata(times/period, power_time_series, nBins)
        phase_hist = np.divide(bdata,bnum)
        sigma_array = getSigmaArray(phase_hist)
        if not sigma_mode :
            mapData[i] = phase_hist 
        else :
            mapData[i] = sigma_array 
        this_sigma = np.max(sigma_array)
        best_string = ""
        if  this_sigma > best_sigma :
            best_period = period 
            best_sigma = np.max(sigma_array)
            best_sigma_array = 1.001*sigma_array 
            mean, sigma = getMeanSigma(phase_hist)
            results["hist_baseline"], results["hist_sigma"] = mean, sigma 
            results["best_signal"] = 1000.*gain*(np.max(phase_hist)-mean)
            best_string = "**"
            print("i={0:d} period={1:.4f} best_period={2:.4f} sigma={3:.3f} best_sigma={4:.3f}{5:s}".format(
                i,1000.*period,1000.*best_period,this_sigma,best_sigma,best_string))

    results["best_sigma"] = float(best_sigma)
    results['period'] = best_period 
    fig = plt.figure(figsize=(10.,7.5))
    ax = fig.add_subplot(111)
    ax.set_title("S/N: {0:s}  Best Period={1:.4f} ms".format(metadata['target'],1000.*best_period))
    ax.patch.set_facecolor('white')
    ax.set_xlabel("Phase",size=14)
    ax.set_ylabel("Period (ms)",size=14)
    eps = float(args.eps)
    mapData = np.maximum(mapData,0.)
    im = ax.imshow(mapData,extent=[0,1.,1000.*periods[-1],1000.*periods[0]],aspect='auto')
    im.set_cmap('jet')
    plt.colorbar(im, use_gridspec=True)
    if not args.no_plot : plt.show() 

nPoints = len(phase_hist)
if mode in ['polyco','vlsr','f0'] :
    best_sigma_array = getSigmaArray(phase_hist)
    nPoints = len(best_sigma_array)
    best_period = p0

roll = 0 

if not args.no_roll :  
    roll = int(nPoints/2) - np.argmax(best_sigma_array) 
    best_sigma_array = np.roll(best_sigma_array,roll)
    phase_hist = np.roll(phase_hist,roll)
results["roll"] = int(roll)

x = np.linspace(0.,1.,nPoints)
x_err = np.zeros_like(x)
y_err = np.ones_like(best_sigma_array)
fig2, axs2 = plt.subplots(figsize=(8,7))
if sigma_mode :
    axs2.errorbar(x,best_sigma_array,xerr=x_err,yerr=y_err,color='red',ecolor='black',fmt='o')
    plt.ylabel("Signal/Noise",fontsize=16)
    #plt.ylim(bottom=-5., top=40.)
else :
    mean, sigma = getMeanSigma(phase_hist)
    phase_hist = 1000.*gain*(phase_hist-mean)
    axs2.errorbar(x,phase_hist,xerr=x_err,yerr=1000.*gain*sigma*y_err,color='red',ecolor='black',fmt='o')
    plt.ylabel("Signal (mK)",fontsize=16)

plt.xlabel("Phase",fontsize=16)
plt.title("{0:s}".format(args.base_name, fontsize=14))
bs = max(best_sigma_array)


if True : writeResults(results,base_name)

ymin, ymax = plt.ylim() 
dy = ymax - ymin 

if mode == "polyco" :
    plt.text(0.05,0.89*dy+ymin,'Polyco mode ' + metadata['target'] ,fontsize=14) 
    plt.text(0.05,0.83*dy+ymin,"Period={0:.4f} ms".format(1000.*average_period),fontsize=14)
    plt.text(0.05,0.77*dy+ymin,'Sigma={0:.2f} Roll={1:d}'.format(bs,roll),fontsize=14)

elif mode == "scan" :
    plt.text(0.05,0.89*dy+ymin,'Scan mode ' + metadata['target'] ,fontsize=14) 
    plt.text(0.05,0.83*dy+ymin,"Best period={0:.4f} ms".format(1000.*best_period),fontsize=14)
    plt.text(0.05,0.77*dy+ymin,'Best sigma={0:.2f} Roll={1:d}'.format(bs,roll),fontsize=14)

elif mode == "f0" :
    plt.text(0.05,0.89*dy+ymin,'f0 mode ' + metadata['target'] ,fontsize=14) 
    plt.text(0.05,0.83*dy+ymin,"Period={0:.4f} ms".format(1000.*best_period),fontsize=14)
    plt.text(0.05,0.77*dy+ymin,'Best sigma={0:.2f} Roll={1:d}'.format(bs,roll),fontsize=14)
    
plt.text(0.05,0.71*dy+ymin,'Az={0:.1f} Alt={1:.1f}'.format(metadata['az'],metadata['alt']),fontsize=14)

axs2.plot([0.,1.],[0.,0.],'k--')
axs2.set_xlim(0.,1.)
if not args.no_plot : plt.show() 



     






