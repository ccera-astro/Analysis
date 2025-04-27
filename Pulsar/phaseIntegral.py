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
    parser.add_argument("-r","--ref",default="2024-12-24-0217",help="Reference run.")
    parser.add_argument("-t","--target",default="2024-12-28-1858",help="Target run.")
    parser.add_argument("--nPhaseBins",type=int,default=50,help="Number of bins in phase plot")
    parser.add_argument("-p","--phase_0",type=float,default=-1000.,help="Phase shift for reference data.")
    return parser.parse_args()  
 
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

def time2freq(time, best_coeff):	#Converts a set of times (mjd) into phases for the specified pulsar
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
    pcoeff = (60.0 * best_coeff[F0] + best_coeff[coeff2], best_coeff[coeff3], best_coeff[coeff4]) 
    #FREQ(Hz) = F0 + (1/60)*(COEFF(2) + 2*DT*COEFF(3) + 3*DT^2*COEFF(4) + ....)
    print("time={0:f} best_coeff[TMID]={1:f}".format(time,best_coeff[TMID]))
    freq = best_coeff[F0] + (best_coeff[coeff2] + 2.*dt*best_coeff[coeff3] + 3.*dt*dt*best_coeff[coeff4])/60. 
    #freq = np.polynomial.polynomial.polyval(dt, pcoeff)
    return freq

# begin execution here

args = getArgs() 

data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/"
ref_base_name, target_base_name = data_dir + args.ref, data_dir + args.target  
with open(ref_base_name + ".json") as json_file : ref_data = json.load(json_file)
with open(target_base_name + ".json") as json_file : target_data = json.load(json_file)

MJD_ref = ref_data['t_start']/ 86400. + 40587.
MJD_target = target_data['t_start']/ 86400. + 40587.

MJDs =np.linspace(MJD_ref,MJD_target,50)

pulsarName = ref_data['target']
if pulsarName[0] == 'J' : pulsarName = pulsarName[1:]

freqs = []
for MJD in MJDs :
    coeff = rp.getpolycoeff(MJD, ref_data, ref_base_name)
    freq = time2freq(MJD, coeff)
    freqs.append(freq)
    print("MJD={0:.13f} freq={1:.7f}".format(MJD,freq))

f_avg = np.average(np.array(freqs))
print("f_avg={0:.7f}".format(f_avg))


plt.plot(MJDs,freqs,'bo')
plt.plot([MJDs[0],MJDs[-1]],[f_avg,f_avg],'r-')
plt.text(MJDs[5],f_avg + 0.5e-7,'f_bar={0:.7f} Hz'.format(f_avg),color='r')
plt.xlabel('MJD')
plt.ylabel('f (Hz)')
plt.show() 







