# plot properties of polycos 
import runPolyco as rp
import numpy as np
import matplotlib.pyplot as plt 

def getpolycoeff(mjd, pulsar_name, rfreq) :
    file_name = 'polyco_{0:d}.dat'.format(int(mjd))
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
times = np.linspace(-24*3600.,24*3600.,5000)
pulsar_name, rfreq = 'J0332+5434', 1420.4 
coeff, phases = {}, {} 
for MJD in [60679, 60680] :
    mjd = MJD + 0.6875
    coeff[MJD] = getpolycoeff(mjd, pulsar_name, rfreq) 
    MJD0 = coeff[60679][1]
    MJDs = MJD0  + times/86400. 
    p0 = 1./coeff[MJD][4]
    print("p0={0:f}".format(p0))
    phases[MJD] = time2phase(MJDs, coeff[MJD]) - coeff[60679][3]
    
c = coeff[60679]
print("f={0:e} + {1:e} + {2:e}".format(c[4],c[7]/60.,2.0*c[8]))
plt.plot(times/3600.,phases[60680]-phases[60679],'b-')
plt.show() 

