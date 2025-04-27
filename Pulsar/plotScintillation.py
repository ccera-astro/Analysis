# analyze pulsar data using polyco method
# make a heatmap showing phase vs. time slice 
# this is derived from anaPulsar.py, but with many features removed 

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt 
import socket 
import json 
import runPolyco as rp 
import scipy.stats as stats

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-13-1858",help="File(s) to be analyzed.")
    parser.add_argument("-n","--name",default="J0332+5434",help="Pulsar to be analyzed.")
    parser.add_argument("--nPhaseBins",type=int,default=50,help="Number of bins in phase plot")
    parser.add_argument("--power_mode",action="store_true",help="Plot power values")
    parser.add_argument("--down_sample",type=int,default=1,help="Down sample (averaging) factor")
    parser.add_argument("--filter_plots",action="store_true",help="Enable plots related to noise filter.")
    parser.add_argument("--high_pass_off",action="store_true",help="Turn off high pass filter.")
    return parser.parse_args()

def getFileName(args, base_name) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/J0332+5434/" 
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
    return data_dir + base_name 

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

b = args.base_name 
if b == '2024-06-13-1858': 
    base_names = [getFileName(args,b)]
else : 
    base_names = [] 
    for b in ['2024-12-28-1831','2024-12-28-1858','2024-12-28-1925','2024-12-28-1952'] :
        base_names.append(getFileName(args,b)) 

slices_per_run = 8
n_time_slices = slices_per_run*len(base_names)
mapData = np.zeros((n_time_slices,args.nPhaseBins)) 
slice_time, signal = np.zeros(n_time_slices), np.zeros(n_time_slices)

for iRun, base_name in enumerate(base_names) :    
    with open(base_name + ".json") as json_file : metadata = json.load(json_file)
    if iRun == 0 : t0 = metadata['t_start'] 
    t_run = metadata['t_start'] - t0 
    print(" iRun={0:d} t_run={1:f} base_name={2:s}".format(iRun,t_run,base_name)) 
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

    pname = metadata['target']
    nBins = args.nPhaseBins 
    sigma_mode = not args.power_mode  

    # set up time series 
    pts_1 = 1000*np.fromfile(base_name + "_1.sum", dtype=np.float32)
    pts_2 = 1000*np.fromfile(base_name + "_2.sum", dtype=np.float32)
    if len(pts_1) > len(pts_2) : pts_1 = pts_1[:len(pts_2)]
    if len(pts_2) > len(pts_1) : pts_2 = pts_2[:len(pts_1)]
    power_time_series = pts_1 + pts_2 

    if n_downsample > 1 : power_time_series = downSample(power_time_series,n_downsample)
    nRows = len(power_time_series) 

    if not args.high_pass_off :
        # apply a high pass filter to mitigate baseline wandering
        f_cutoff = 0.1
        power_time_series = highpass(power_time_series,f_cutoff,c_rate)

    times = np.linspace(0.,nRows*t_fft*n_downsample,nRows) 
    run_time = times[-1] - times[0]
    middle_time = 0.5*(times[-1] + times[0])
    print("run_time={0:.2f}".format(run_time))

    power_time_series, bad_elements = denoise(power_time_series)
    times = np.delete(times,bad_elements)

    samples_per_slice = int(len(times)/slices_per_run)     

    pulsarName = pname 
    if pulsarName[0] == 'J' : pulsarName = pulsarName[1:]
    MJD = (metadata['t_start']/ 86400.) + 40587. 
    coeff = rp.getpolycoeff(MJD, metadata, base_name)
    print("MJD={0:f} coeff={1:s}".format(MJD,str(coeff)))
    p0 = 1./coeff[4]
    print("iRun={0:d} len(times)={1:d} samples_per_slice={2:d}".format(iRun,len(times),samples_per_slice))
    for i in range(slices_per_run) :
        j1 = i*samples_per_slice 
        j2 = j1 + samples_per_slice - 1
        print("   i={0:d}  j1={1:d}  j2={2:d}".format(i,j1,j2))
        MJDs = MJD + times[j1:j2]/86400. 
        phase = time2phase(MJDs, coeff)
        bdata, bnum = bindata(phase, power_time_series[j1:j2], nBins)
        phase_hist = np.divide(bdata,bnum)
        sigma_array = getSigmaArray(phase_hist)
        mapData[iRun*slices_per_run + i] = sigma_array
        slice_time[iRun*slices_per_run + i] = t_run + 0.5*(times[j1]+times[j2]) 
        signal[iRun*slices_per_run + i] = np.sum(sigma_array[35:38])

fig = plt.figure(figsize=(10.,7.5))
ax = fig.add_subplot(111)
ax.set_title("S/N: {0:s}  {1:s}  Period={2:.3f} ms".format(metadata['target'],base_name.split('/')[-1],1000.*p0))
ax.patch.set_facecolor('white')
ax.set_xlabel("Phase",size=14)
ax.set_ylabel("Time (s)",size=14)
mapData = np.maximum(mapData,0.)
im = ax.imshow(mapData,extent=[0,1.,slice_time[-1],0.],aspect='auto')
im.set_cmap('jet')
plt.colorbar(im, use_gridspec=True)
plt.show()

plt.plot(slice_time,signal,'ro')
plt.title("Signal/noise vs. time {0:s} {1:s}".format(pname,base_name.split('/')[-1]))
plt.xlabel("t (s)")
plt.ylabel("Signal/Noise")
plt.ylim(bottom=0.)
plt.show() 








