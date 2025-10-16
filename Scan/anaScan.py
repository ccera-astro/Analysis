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
    parser.add_argument("-g","--gain",type=float,default=1.,help="Calibration factor")
    parser.add_argument("-c","--chan",type=int,default=1,help="Receiver channel")
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
    times = np.linspace(mean_time-1000.,mean_time+1000.,200)
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

chan = args.chan 

file = base_name + "_{0:d}.sum".format(chan)
if False :
    file = base_name + "_{0:d}.raw".format(chan)
    vals, rows, cols = getData(file,metadata['fft_size'])
    print("rows={0:d} cols={1:d} len(vals)={2:d}".format(rows,cols,len(vals)))
    power = np.zeros(rows)
    for i in range(rows) : 
        power[i] = 8.*np.sum(vals[i][1024:1280])    
    print("power[0:10]={0:s}".format(str(power[0:10])))  
else : 
    power = args.gain*np.fromfile(file, dtype=np.float32)[2:]

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
    calibration_factor = 17.5*1.5
    if metadata['target'] == 'M87' : calibration_factor = 18.65 
    if metadata['target'] == 'CasA' : calibration_factor = 45.1 
    if metadata['target'] == '3C353' : calibration_factor = 45.1 
    calibration_factor = 45.1 
    power *= calibration_factor 

if metadata['target'] == 'CygnusA' or metadata['target'] == 'test' :
    mean_time = {'06-25':3620, '06-27':3622, '06-29':3620, '06-30':3610., '07-01':3610., '07-23': 1822., '07-24': 1825., '07-26': 1800. }
    base_temp = {'06-25':89., '06-27':93.5, '06-29':94., '06-30':92.9, '07-01': 94., '07-23': 113., '07-24': 110., '07-26': 100.}
    peak_temp = {'06-25':20.6, '06-27':21.3, '06-29':22., '06-30':22., '07-01': 22., '07-23': 28., '07-24': 25.5, '07-26': 10.} 
    width = {'06-25':370., '06-27':370., '06-29':370., '06-30':370., '07-01': 370., '07-23': 366., '07-24': 370., '07-26': 370.} 
elif metadata['target'] == 'Sun' :  
    mean_time = {'07-02':3227., '07-03':1854., '07-04':1839., '07-05':1812. ,'07-06':1869. , '07-22' :1820., '07-23': 3945.}
    base_temp = {'07-02': 94., '07-03': 94., '07-04': 94., '07-05': 94., '07-06': 94., '07-22': 94., '07-23':94.}
    peak_temp = {'07-02': 15050., '07-03':13450., '07-04': 11580., '07-05': 9438., '07-06': 7410., '07-22': 20895., '07-23': 19010.} 
    width = {'07-02': 310., '07-03':310., '07-04':312., '07-05':312., '07-06': 312. , '07-22': 312., '07-23': 312.} 
elif metadata['target'] == 'M87' or metadata['target'] == 'dummy' or metadata['target'] == 'CasA':
    mean_time = {'07-12':1824., '07-13':1824. ,'07-16': 1807., '07-19' : 1820., '07-20': 1820.}
    base_temp = {'07-12': 90.1 , '07-13': 90.1, '07-16':92.8, '07-19': 91.5, '07-20':91.5}
    peak_temp = {'07-12': 2.35, '07-13':2.35, '07-16':3.52, '07-19': 20.4, '07-20':20.4}
    width = {'07-12': 270. ,'07-13': 270., '07-16': 328., '07-19': 525., '07-20':525.}
elif metadata['target'] == '3C353' : 
    mean_time = {'07-20': 1800. , '07-21': 1800.}
    base_temp = {'07-20': 88.2 , '07-21': 90.}
    peak_temp = {'07-20': 0.5, '07-21': 0.5}
    width = {'07-20': 310. , '07-21': 310. }
else :
    mean_time = {'07-07': 1800., '07-09':3600. }
    base_temp = {'07-07': 86., '07-09': 95. }
    peak_temp = {'07-07': 0.6, '07-09': 2.4}
    width = {'07-07': 310., '07-09': 310. }

airy_times, airy_function = airy(mean_time[day], base_temp[day], peak_temp[day], width[day])
max_gmt = metadata['t_start'] + mean_time[day] 

# convert the fitted width to an angle 
dec = metadata['dec']
dPhi = (360./24./3600.)*cos(radians(dec))*width[day]

plt.figure(2)
if x310 : gain_inverse = np.divide(90.*np.ones(nVals,dtype=np.float32),gain)
plt.plot(times,power,'b.',label="Chan {0:d}".format(chan))
plt.plot(airy_times,airy_function,'r-') 
ts = datetime.fromtimestamp(max_gmt).strftime('%Y-%m-%d %H:%M:%S')

plt.title("Power vs. Time  {0:s}".format(base_name.split("/")[-1]))
plt.xlabel("t (s)")
plt.ylabel("Power (K)")
if metadata['target'] == '3C273' : 
    plt.ylim(85.,105.)
    plt.text(1500.,100.0,'Peak={0:s}'.format(ts))
    plt.text(1500,98.0,'Width={0:.2f} deg.'.format(dPhi))
elif metadata['target'] == 'M87' :
    #plt.ylim(89.,93.)
    plt.text(100.,95.5,'Peak={0:s}'.format(ts))
    plt.text(100,95.0,'Width={0:.2f} deg.'.format(dPhi))
else :
    yMin, yMax = 80., np.max(power)
    plt.text(100.,yMin + 0.9*(yMax-yMin),'Peak={0:s}'.format(ts))
    plt.text(100,yMin+0.8*(yMax-yMin),'Width={0:.2f} deg.'.format(dPhi))
plt.legend() 
plt.show()

# plot spectrum
file = base_name + "_{0:d}.avg".format(chan)
vals, rows, cols = getData(file,metadata['fft_size']) 
print("rows={0:d} cols={1:d} len(vals[0])={2:d}".format(rows,cols,len(vals[0])))
plt.figure(3)
plt.plot(vals[0])
plt.show() 
