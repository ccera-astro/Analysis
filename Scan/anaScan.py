import numpy as np
import matplotlib.pyplot as plt
import time 
from math import pi, factorial, degrees, radians, cos, sin, acos 
import ephem
from datetime import datetime , timezone
import json 
import socket 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("--data_dir",default=None,help="data directory")
    parser.add_argument("-b","--base_name",default="2024-06-27-0029",help="File(s) to be analyzed.")
    return parser.parse_args()

def getFileName(args) :
    if args.data_dir :
        data_dir = args.data_dir 
    else :
        data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/Scan/" 
        if "receiver" in socket.gethostname().lower() : data_dir = "/home/dmarlow/data/"
    return data_dir + args.base_name 

def getData(file,fft_size) :
    print("Reading from file: {0:s}".format(file))
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size)
    vals = np.reshape(vals, (rows,cols))   
    return vals, rows, cols

def getObserver() :
    princeton=ephem.Observer()
    princeton.lat='40.344892'
    princeton.lon='-74.651692'
    carp = ephem.Observer()
    carp.lat = '45.34998'
    carp.lon = '-76.05693'
    return carp 

def airy(mean_time, base_temp, peak_temp, width) : 
    times = np.linspace(mean_time-700.,mean_time+700.,200)
    airy_function = [] 
    for tt in times :
        r = (tt-mean_time)/width 
        t, I = 0.5*pi*r, 1.
        for k in range(1,20) :
            dI = (-1)**k * t**(2*k) / (factorial(k)*factorial(k+1))
            I = I + dI 
        I *= I
        airy_function.append(I) 
    return times, base_temp + peak_temp*np.array(airy_function)

def getDay(args) :
    #2024-06-27-0029
    args.base_name  
    date = f.strip(".raw")
    day = date[-4:]
    return day

# begin execution

args = getArgs()

# set up pyephem stuff
obs = getObserver()
CasA = ephem.readdb("CasA,f|J|F7,23:23:25.8,58:48:00,2.02,2000")
CasA.compute(obs)
CygA = ephem.readdb("CygA,f|J|F7,19:59:28.3,40:44:02.1,2000")
CygA.compute(obs) 

base_name = getFileName(args)
with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
n_decimate = metadata['decimation_factor']
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 

# analyze channel 2 only 

file = base_name + "_2.sum"
power = np.fromfile(file, dtype=np.float32)
nVals = len(power)
times = np.linspace(0.,nVals*t_fft,nVals) 

day = args.base_name[5:10]

# for x310 data, divide by gain function
x310 = day in ['06-25','06-27','06-29'] 
if x310 :
    gain_file = 'gain.dat'
    gain = np.fromfile(gain_file,dtype=np.float32)
    gain_len = len(gain)
    # trim longer array to length of shorter array 
    if nVals > gain_len : 
        power = power[:gain_len]
        times = times[:gain_len]
    else : 
        gain = gain[:nVals]
    power *= gain
else :
    calibration_factor = 17.5
    power *= calibration_factor 

if metadata['target'] == 'CygnusA' :
    mean_time = {'06-25':3620, '06-27':3622, '06-29':3620, '06-30':3610., '07-01':3610. }
    base_temp = {'06-25':89., '06-27':93.5, '06-29':94., '06-30':92.9, '07-01': 94.}
    peak_temp = {'06-25':20.6, '06-27':21.3, '06-29':22., '06-30':22., '07-01': 22.} 
    width = {'06-25':370., '06-27':370., '06-29':370., '06-30':370., '07-01': 370.} 
else :
    mean_time = {'07-02':3227. }
    base_temp = {'07-02': 94.}
    peak_temp = {'07-02': 15050.} 
    width = {'07-02': 310.} 

airy_times, airy_function = airy(mean_time[day], base_temp[day], peak_temp[day], width[day])
max_gmt = metadata['t_start'] + mean_time[day] 

# convert the fitted width to an angle 
dec = metadata['dec']
dPhi = (360./24./3600.)*cos(radians(dec))*width[day]

plt.figure(2)
if x310 : gain_inverse = np.divide(90.*np.ones(nVals,dtype=np.float32),gain)
plt.plot(times,power,'b-',label="Chan 2")
plt.plot(airy_times,airy_function,'r-') 
ts = datetime.fromtimestamp(max_gmt).strftime('%Y-%m-%d %H:%M:%S')
yMin, yMax = 80., np.max(power)
plt.text(100.,yMin + 0.9*(yMax-yMin),'Peak={0:s}'.format(ts))
plt.text(100,yMin+0.8*(yMax-yMin),'Width={0:.2f} deg.'.format(dPhi))
plt.title("Power vs. Time  {0:s}".format(base_name.split("/")[-1]))
plt.xlabel("t (s)")
plt.ylabel("Power (K)")
#plt.ylim(80.,120.)
plt.legend() 
plt.show()



# Airy function parameters
#phiMean, width = 0.0 , 1.17
#baseTemp, peakTemp, meanTime = 100., 20500., 17.010

if False :
    # computer Airy curve based on phiDot method
    plt.figure(2)
    plt.plot(dPhi,airyFun,'r-',label="Airy Function Width = {0:.2f} deg".format(width))
    phiDot = (360./24.)*cos(sun.dec)
    print("phiDot={0:e}".format(phiDot))
    times = (times - meanTime)
    dPhi = phiDot*times   
    plt.plot(dPhi,power,'bo',label='Observed Power (K)')
    txt = r'$\Delta\phi$'
    maxY = 1.25*peakTemp
    plt.title("Power vs. " + txt)
    plt.text(-4.8,0.94*maxY,"Sun scan: {0:d}-{1:02d}-{2:02d}".format(yy,mon,dd))
    plt.text(-4.8,0.88*maxY,"Peak time: {0:.3f} hrs".format(meanTime))
    plt.xlim(-5.,5.)
    plt.ylim(0.,1.25*peakTemp)
    plt.grid(True)
    plt.xlabel(txt + " (deg)")
    plt.ylabel("Power (K)")
    plt.legend(loc="upper right") 
    plt.show()


exit()
dec = 63.15
phiDot = cos(radians(dec))*15./3600.     # 
angles = phiDot*(times - tMean)
angleWidth = phiDot*width
print("Angle Width={0:.2f}".format(angleWidth))

plt.figure(1)
plt.plot(angles,power,'b.')
#plt.plot(angles,airyFun,'r-')
plt.title("Power vs. Time")
plt.xlabel("t (s)")
plt.ylabel("Power (arb. units)")
plt.show()
