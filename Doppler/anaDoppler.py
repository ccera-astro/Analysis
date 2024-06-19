#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import time
import earthVelocity as eV
import earthVelocityLookup as eVLU
import galCoords 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("-b","--base_name",default="./data/2024-06-15-1335",help="Base name of file to be analyzed.")
    #parser.add_argument("--freqs",action='store_true',help="Use frequency for x axis.")
    parser.add_argument("-g","--gain_factor",type=float,default=66000.,help="Gain factor")
    #parser.add_argument("--keepSpurs",action='store_true',help="Don't remove spurious peaks from spectrum.")
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
    #print("background={0:s}".format(str(background)))
    return background

# Begin execution here

args = getArgs()
#gain = getGain(series) 

# get the first file to establish the parameters of the plot
base_name = args.base_name 
meta_data = getMetaData(base_name + ".json")
freqs = getFreqs(meta_data) 
print("freqs[0]={0:.2f} freqs[-1]={1:.2f}".format(freqs[0],freqs[-1]))
vDoppler = getVelocities(freqs)

scale2 = 0.6
if True :
    power1 = np.fromfile(base_name +"_1.avg", dtype=np.float32)
    power2 = scale2*np.fromfile(base_name +"_2.avg", dtype=np.float32)
    plt.plot(freqs,power1,'b-',label="Chan 1")
    plt.plot(freqs,power2,'r-',label="Chan 2 x {0:.3f}".format(scale2))
    plt.title("Raw spectrum no corrections: {0:s}".format(base_name))
    plt.xlabel("f (MHz)")
    plt.ylabel("PSD (AU)")
    plt.ylim(0.,1.1*np.max(power1)) 
    plt.legend() 
    plt.show() 

power = scale2*args.gain_factor*np.fromfile(base_name +"_2.avg", dtype=np.float32)

vMin, vMax = -300., 300.
i1 = np.searchsorted(vDoppler,vMin)
i2 = np.searchsorted(vDoppler,vMax)
print("i1={0:d} i2={1:d}".format(i1,i2))
freqs = freqs[i1:i2]
vDoppler = vDoppler[i1:i2]
power = power[i1:i2]
background = fitBackground(vDoppler,power,3,200.)
#gain = gain[i1:i2]



#plt.plot(vDoppler,power,'r.',label="Power")
#plt.plot(vDoppler,background,'g-',label="3rd Order Chebyshev")
plt.plot(vDoppler,power-background,'r.')
plt.plot([-300.,300.],[0., 0.],'g-')
plt.title("Spectrum with gain factor {0:s}".format(base_name))
plt.xlabel('Doppler Velocity (km/s)')
plt.ylabel('Power (K)')
plt.legend() 
plt.show()