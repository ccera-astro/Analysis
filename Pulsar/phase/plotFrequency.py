# plot properties of polycos 
import runPolyco as rp
import numpy as np
import matplotlib.pyplot as plt 

def getpolycoeff(mjd, pulsar_name, rfreq) :
    file_name = 'polyco.dat'
    print("Entering getplotcoeff mjd={0:f} file_name={1:s} rfreq={2:f}".format(mjd,file_name,rfreq))
    fp = open(file_name,'r')
    best_coeff, ind, num, tfreq = rp.readpolycoeff(mjd,file_name)
    file_valid = best_coeff[0] == pulsar_name[1:] and ind != num and (abs(tfreq-rfreq)/tfreq < 0.01)

    if not file_valid  :
        print("********* Unable to make valid polyco file ********")
        print("best_coeff[0]={0:s} name={1:s} ind={2:d} num={3:d} tfreq={4:f} rfreq={5:f}".format(
             best_coeff[0],pulsar_name[1:],ind,num,tfreq,rfreq))
        
        exit() 
    
    print("Leaving getplotcoeff()")
    return best_coeff 

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

# begin execution here 
times = np.linspace(0.,3600.,1000)
pulsar_name, rfreq = 'J0332+5434', 1420.4 
coeff, phases = {}, {} 
hours = np.linspace(0.,64.,25)
freqs = np.zeros_like(hours)
for i, hour in enumerate(hours) :
    mjd = 60679.5 + hour/24. 
    coeff = getpolycoeff(mjd, pulsar_name, rfreq) 
    MJD0 = coeff[1]
    MJDs = MJD0  + times/86400. 
    phases = time2phase(MJDs, coeff) - coeff[3]
    average_frequency = (phases[-1]-phases[0])/(times[-1]-times[0])
    freqs[i] = average_frequency
    print("MJD={0:f} hour={1:.3f} average_frequency={2:.12f}".format(mjd,hour,average_frequency))
    
dfdt = 24.*(freqs[-5]-freqs[2])/(hours[-5]-hours[2])
print("dfdt={0:e} Hz/day  df/dt/f={1:.2e}".format(dfdt,dfdt/freqs[2]))
print("phase difference per day={0:.2f}".format(dfdt*86400.))

plt.plot(hours,freqs,'ro')
plt.plot([hours[2],hours[-5]],[freqs[2],freqs[-5]],'b--')
plt.text(hours[7],freqs[4],'df/dt={0:.2e} Hz/day'.format(dfdt))
plt.xlabel("Hour")
plt.ylabel("Pulsar frequency (Hz)")
plt.show() 


