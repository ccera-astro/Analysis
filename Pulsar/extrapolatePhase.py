# extrapolate phase from one observing run to another 
import numpy as np 
import json 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import runPolyco as rp 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--ref",default="2024-12-24-0217",help="Reference run.")
    parser.add_argument("-t","--target",default="2024-12-28-1858",help="Target run.")
    parser.add_argument("--nPhaseBins",type=int,default=50,help="Number of bins in phase plot")
    parser.add_argument("-p","--phase_0",type=float,default=-1000.,help="Phase shift for reference data.")
    return parser.parse_args()  

def getData(base_name) :
    pts_1 = 1000*np.fromfile(base_name + "_1.sum", dtype=np.float32)
    pts_2 = 1000*np.fromfile(base_name + "_2.sum", dtype=np.float32)
    if len(pts_1) > len(pts_2) : pts_1 = pts_1[:len(pts_2)]
    if len(pts_2) > len(pts_1) : pts_2 = pts_2[:len(pts_1)]
    return pts_1 + pts_2 

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

def bindata( phase, value, NumBins):	# Averages 'value' into 'NumBins' bins based on the 'phase' values
    mphase = phase - np.floor(phase)
    iphase = np.floor(mphase * NumBins)
    bdata = np.zeros((NumBins))
    bnum = np.zeros((NumBins))

    (bdata,d1,d2) = stats.binned_statistic(iphase, value, statistic='sum', bins=NumBins, range = [0,NumBins-1])
    (bnum,d1,d2) = stats.binned_statistic(iphase, value, statistic='count', bins=NumBins, range = [0,NumBins-1])

    return (bdata, bnum)

def getSigmaArray(x) :
    iMax = np.argmax(x)
    # remove phase bins within +/- 3 of peak 
    iLow, iHigh = max(0,iMax-3), min(len(x),iMax+3)
    y = np.delete(x,range(iLow,iHigh))
    mean, sigma = np.mean(y), np.std(y)
    xx = (x-mean)/sigma 
    return xx 

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

# do a standard polyco analysis on this dataset 
def anaPolyco(args, base_name, data) :
    print("data={0:s}".format(str(data)))
    # read or calculate various run parameters 
    c_rate = data['srate']/data['fft_size']/data['decimation_factor'] 
    t_fft = 1./c_rate 
    nBins = args.nPhaseBins 

    power_series = getData(base_name)
    nRows = len(power_series)

    times = np.linspace(0.,nRows*t_fft,nRows)
    power_series, bad_elements = denoise(power_series)
    times = np.delete(times,bad_elements)

    MJD = data['t_start']/ 86400. + 40587.
    MJDs = MJD + times/86400. 
    pulsarName = data['target']
    if pulsarName[0] == 'J' : pulsarName = pulsarName[1:]
    coeff = rp.getpolycoeff(MJD, data, base_name)

    phase = time2phase(MJDs, coeff) 
    p = (times[-1]-times[0])/(phase[-1]-phase[0])
    bdata, bnum = bindata(phase, power_series, nBins)
    ph = np.divide(bdata,bnum)
    sa = getSigmaArray(ph)
    s = max(sa)

    return p, 1./p, s, sa

def anaf0(args, base_name, data, f, phase_offset) :
    # read or calculate various run parameters 
    c_rate = data['srate']/data['fft_size']/data['decimation_factor'] 
    t_fft = 1./c_rate 
    nBins = args.nPhaseBins 

    power_series = getData(base_name)
    nRows = len(power_series)

    times = np.linspace(0.,nRows*t_fft,nRows)
    power_series, bad_elements = denoise(power_series)
    times = np.delete(times,bad_elements)

    phase = f*times + phase_offset 
    bdata, bnum = bindata(phase, power_series, nBins)
    ph = np.divide(bdata,bnum)
    sa = getSigmaArray(ph)

    nc = np.argmax(sa)
    center_of_mass = (-sa[nc-1]+sa[nc+1])/np.sum(sa[nc-1:nc+2]) + nc + 0.5
    print("nc={0:d} center_of_mass={1:.3f}".format(nc,center_of_mass))
    phase_peak = center_of_mass/nBins - 0.5
    print("nc={0:d} phase_peak={1:.3f} f_ref={2:.7f}".format(nc,phase_peak,f_ref))

    return sa, phase_peak 

# begin execution here 
args = getArgs() 
nBins = args.nPhaseBins 
data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/"
ref_base_name, target_base_name = data_dir + args.ref, data_dir + args.target  
with open(ref_base_name + ".json") as json_file: ref_data = json.load(json_file)
if ref_data['t_start'] % 1. > 0.45 : ref_data['t_start'] += 1. 
with open(target_base_name + ".json") as json_file: target_data = json.load(json_file)
if target_data['t_start'] % 1. > 0.45 : target_data['t_start'] += 1.

# do a standard PINT/polyco analysis on both data sets 
p_ref, f_ref, sigma_ref, sigma_array_ref = anaPolyco(args,ref_base_name, ref_data)
p_target, f_target, sigma_target, sigma_array_target = anaPolyco(args,target_base_name, target_data)

# do a fixed frequency analysis on the reference data set 
sigma_array, ref_phase  = anaf0(args, ref_base_name, ref_data, f_ref, 0.)

#plt.plot(sigma_array,'bo')
#plt.show() 

# extrapolate to target 
f_avg = 0.5*(f_ref+f_target)
print("f_ref={0:.7f} f_target={1:.7f} avg={2:.7f}".format(f_ref,f_target,f_avg))
dt = target_data['t_start'] - ref_data['t_start']
target_phase = f_avg*dt - ref_phase
print("t_ref={0:.3f} t_target={1:.3f} dt={2:.3f}".format(ref_data['t_start'],target_data['t_start'],dt))
print("ref_phase={0:.3f} d_phi={1:.3f} phi_1={2:.3f}".format(ref_phase,f_avg*dt,target_phase))

# do a fixed frequency analysis on the target data set 
sigma_array, phase_peak  = anaf0(args, target_base_name, target_data, f_target, target_phase)

# plot the result 
x = np.linspace(-0.5,0.5,nBins)
x_err = np.zeros_like(x)
y_err = np.ones_like(sigma_array)
fig2, axs2 = plt.subplots(figsize=(8,7))
axs2.errorbar(x,sigma_array,xerr=x_err,yerr=y_err,color='red',ecolor='black',fmt='o')
plt.ylabel("Signal/Noise",fontsize=16)
plt.xlabel("Phase",fontsize=16)
plt.title('Extrapolation from {0:s} to {1:s}'.format(args.ref,args.target) ,fontsize=14) 
bs = max(sigma_array)

ymin, ymax = plt.ylim() 
dy = ymax - ymin 

plt.text(0.05,0.77*dy+ymin,'Best sigma={0:.2f}'.format(bs),fontsize=14)
plt.text(0.05,0.71*dy+ymin,'Az={0:.1f} Alt={1:.1f}'.format(target_data['az'],target_data['alt']),fontsize=14)

plt.arrow(phase_peak,ymax/6.,0.,-ymax/6.,head_width=0.02,head_length=ymax/30.,fc='k',length_includes_head=True)
plt.text(phase_peak+0.02,ymax/6.,r'$\phi={0:.3f}$'.format(phase_peak),fontsize=16)
axs2.plot([-0.5,0.5],[0.,0.],'k--')
axs2.set_xlim(-0.5,0.5)
plt.show() 


