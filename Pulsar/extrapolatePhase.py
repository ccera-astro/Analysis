# extrapolate phase from one observing run to another 
import numpy as np 
import json 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import runPolyco as rp 
import fitPhase  
from math import sqrt 
import random 
from EarthOrbit import earth_state_xy_km 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--ref",default="2024-12-24-0217",help="Reference run.")
    parser.add_argument("-t","--tgt",default="2024-12-25-0214",help="Target run.")
    parser.add_argument("--nPhaseBins",type=int,default=200,help="Number of bins in phase plot")
    parser.add_argument("-p","--phase_0",type=float,default=-1000.,help="Phase shift for reference data.")
    parser.add_argument("--plot",action="store_true",help="Plot results")
    parser.add_argument("--run_all",action="store_true",help="Plot results")
    parser.add_argument("--dither",action="store_true",help="dither fitted phase by 0.1")
    parser.add_argument("-m","--mode",default="polyco",help="Method used for initial freq:  polyco, H1, pass1, pass2")
    parser.add_argument("--nSigma",type=float,default=5.0,help="Minimum SNR")
    parser.add_argument("--n_files",type=int,default=999999,help="Maximum number of files to analyze")
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
    if False : 
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
    #print("data={0:s}".format(str(data)))
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
    #print("f={0:f} phase_offset={1:f}".format(f,phase_offset))
    phase = f*times + phase_offset 
    bdata, bnum = bindata(phase, power_series, nBins)
    ph = np.divide(bdata,bnum)
    sa = getSigmaArray(ph)
     
    amp, mean, sigma, mean_err  = fitPhase.fitPhasePlot(sa) 

    return sa, mean, mean_err 

def orbit_function_3(MJD, MJD0, e, coslat, perihelion_angle):
    t_year = 365.24
    MJD1 = MJD0 + (perihelion_angle/360.)*t_year
    x, y, vx, vy = earth_state_xy_km(MJD-MJD1,e=e,omega=perihelion_angle)
    return vx*coslat  

def orbit_function_H1_DRM(MJD, offset, e, coslat, perihelion_angle):
    MJD0 = 60822.75
    t_year = 365.24
    MJD1 = MJD0 + (perihelion_angle/360.)*t_year
    x, y, vx, vy = earth_state_xy_km(MJD-MJD1,e=e,omega=perihelion_angle)
    return vx*coslat + offset

def getScanFreq(args, base_name, data) :
    if args.mode.lower() == "polyco" :
        p, f, sigma, sigma_array = anaPolyco(args,base_name,data)
        return f
    
    if args.mode.lower() in ['h1','pass1','pass2'] :
        MJD00 = 60729.95 
        PEPOCH = 46473.0
        F0 = 1.399541538720
        F1 = -4.011970E-15     
        f0 = F0 - 86400.*(PEPOCH-MJD00)*F1
        c = 299792.458
        MJD = data['t_start']/ 86400. + 40587.
    
    if args.mode.lower() == 'pass1' :
        MJD0, e, coslat, perihelion = 6.08239457e+04, 2.19226138e-02, 8.19508363e-01, 2.41564142e+02
        vx = orbit_function_3(MJD, MJD0, e, coslat, perihelion)
    elif args.mode.lower() == 'h1' :
        offset, MJD0, e, coslat, perihelion = 1.31850274e-01, 6.08203485e+04, 1.64021937e-02, 8.27603619e-01, 2.07466866e+02
        vx = orbit_function_H1_DRM(MJD, offset, e, coslat, perihelion)
    elif args.mode.lower() == 'pass2' :
        MJD0, e, coslat, perihelion = 6.08227742e+04, 1.67813192e-02, 8.26257985e-01, 2.18123231e+02
        vx = orbit_function_3(MJD, MJD0, e, coslat, perihelion)
    else :
        print("args.mode ={0:s} not implemented . . . exiting.".format(args.mode))
        exit() 

    f = f0*(1. - vx/c)
    return f 


def extrapolatePhase(args) :

    nBins = args.nPhaseBins 
    data_dir = "/mnt/c/Users/marlow/Documents/CCERA/data/"
    ref_base_name, tgt_base_name = data_dir + args.ref, data_dir + args.tgt
    with open(ref_base_name + ".json") as json_file: ref_data = json.load(json_file)
    if ref_data['t_start'] % 1. > 0.45 and ref_data['t_start'] < 1745604855. : 
        ref_data['t_start'] += 1.
    with open(tgt_base_name + ".json") as json_file: tgt_data = json.load(json_file)
    if tgt_data['t_start'] % 1. > 0.45 and tgt_data['t_start'] < 1745604855. : 
        tgt_data['t_start'] += 1.

    print("ref_base_name={0:s} \ntgt_base_name={1:s}".format(ref_base_name,tgt_base_name))
    
    f_ref = getScanFreq(args, ref_base_name, ref_data)
    f_tgt = getScanFreq(args, tgt_base_name, tgt_data)
        
    sigma_array, ref_phase, phase_err_ref  = anaf0(args, ref_base_name, ref_data, f_ref, 0.)

    # extrapolate to target 
    print("In extrapolatePhase(): f_ref={0:.9f} f_tgt={1:.9f}".format(f_ref,f_tgt))
    f_avg = 0.5*(f_ref+f_tgt) 
    dt = tgt_data['t_start'] - ref_data['t_start']
    tgt_phase = f_avg*dt - ref_phase

    sigma_array, phase_peak, phase_err_tgt  = anaf0(args, tgt_base_name, tgt_data, f_tgt, tgt_phase)
    phase_err = sqrt(phase_err_ref*phase_err_ref + phase_err_tgt*phase_err_tgt)
    
    if args.plot : 
        # plot the result 
        x = np.linspace(-0.5,0.5,nBins)
        x_err = np.zeros_like(x)
        y_err = np.ones_like(sigma_array)
        fig2, axs2 = plt.subplots(figsize=(8,7))
        axs2.errorbar(x,sigma_array,xerr=x_err,yerr=y_err,color='red',ecolor='black',fmt='o')
        plt.ylabel("Signal/Noise",fontsize=16)
        plt.xlabel("Phase",fontsize=16)
        plt.title('Extrapolation from {0:s} to {1:s}'.format(args.ref,args.tgt) ,fontsize=14) 
        bs = max(sigma_array)

        ymin, ymax = plt.ylim() 
        dy = ymax - ymin 

        plt.text(0.05,0.77*dy+ymin,'Best sigma={0:.2f}'.format(bs),fontsize=14)
        plt.text(0.05,0.71*dy+ymin,'Az={0:.1f} Alt={1:.1f}'.format(tgt_data['az'],tgt_data['alt']),fontsize=14)

        plt.arrow(phase_peak,ymax/6.,0.,-ymax/6.,head_width=0.02,head_length=ymax/30.,fc='k',length_includes_head=True)
        plt.text(phase_peak+0.02,ymax/6.,r'$\phi={0:.3f}$'.format(phase_peak),fontsize=16)
        axs2.plot([-0.5,0.5],[0.,0.],'k--')
        axs2.set_xlim(-0.5,0.5)
        plt.show()
    print("\n") 
    return phase_peak, phase_err, f_avg   

def MJD_to_unix(MJD) : return 86400.*(MJD-40587.)

def unix_to_MJD(t) : return (t / 86400.) + 40587.

if __name__ == "__main__":
    args = getArgs() 
    if args.run_all :
        import glob 
        files = glob.glob("./outputs/*.json")
        print("len(files)={0:d}".format(len(files)))
           
        lhead = "   reference        target       dPhi      f_avg     f_corr (nHz)   f_err (nHz) days\n"
        outlines, MJD_lines = [lhead], []
        ref_file = files[0]
        args.ref = ref_file.split("/")[-1].replace("_out.json","")
        nDays = 1  # Maximum number of days between reference and target
        if args.mode.lower() == "h1" : nDays = 2  
        if args.mode.lower() in ["polyco","pass2"] : nDays = 5
        for i in range(1,min(args.n_files,len(files))) :
            tgt_file = files[i]
            #print("i={0:3d} ref_file={1:s} tgt_file={2:s}".format(i,ref_file,tgt_file))
            if tgt_file == "./outputs/2025-01-28-2356_out.json" : continue 
            if tgt_file == "./outputs/2025-01-28-0000_out.json" : continue 
            if tgt_file == "./outputs/2025-02-27-2158_out.json" : continue  
            with open(ref_file) as json_file : ref_data = json.load(json_file)
            t_ref = ref_data["t_start"] 
            args.ref = ref_file.split("/")[-1].replace("_out.json","")
            with open(tgt_file) as json_file : tgt_data = json.load(json_file)
            args.tgt = tgt_file.split("/")[-1].replace("_out.json","")
            t_tgt = tgt_data["t_start"] 
            #print("ref: {0:s} tgt: {1:s}".format(args.ref,args.tgt))
            print("best_sigma: ref={0:7.2f} tgt={1:7.2f}  dT={2:7.2f}".format(
                ref_data["best_sigma"],tgt_data["best_sigma"],(t_tgt - t_ref)/86400.))
            # if the target signal is too weak, move on to the next one, keeping the same reference 
            if tgt_data["best_sigma"] < args.nSigma : continue 
            if t_tgt - t_ref < 86400.*nDays :
                dPhi, phi_err, f_avg = extrapolatePhase(args)
                if args.dither :
                    dPhi = dPhi + random.gauss(0.,0.1)
                    if dPhi < -0.5 : dPhi += 1.
                    if dPhi >  0.5 : dPhi -= 1. 
                f_corr = dPhi/(t_tgt-t_ref)
                ref_timing_error = 1.e-6*(1260. - 58.3*ref_data["best_sigma"])
                tgt_timing_error = 1.e-6*(1260. - 58.3*tgt_data["best_sigma"])
                timing_error = sqrt(ref_timing_error**2 + tgt_timing_error**2)   
                f_err = (timing_error/(t_tgt-t_ref))*f_avg
                t_avg = 0.5*(t_ref + t_tgt)
                MJD_avg = unix_to_MJD(t_avg)
                MJD_lines += "{0:15.6f} {1:13.10f} {2:13.10f}\n".format(MJD_avg,f_avg-f_corr,f_err)
                ll = "{0:s} {1:s} {2:6.3f} {3:12.9f} {4:12.1f} {5:12.1f} {6:6.2f}\n".format(
                    args.ref,args.tgt,dPhi, f_avg, 1.e9*f_corr, 1.e9*f_err, (t_tgt-t_ref)/86400.)
                outlines += ll 
                print(lhead + ll)
            # use existing target as new reference 
            ref_file = files[i]

        open("out.txt","w").writelines(outlines) 
        open("MJDs_freqs_temp.txt","w").writelines(MJD_lines)
    else :
        ref_json = "./outputs/{0:s}_out.json".format(args.ref)
        with open(ref_json) as json_file : ref_data = json.load(json_file)
        tgt_json = "./outputs/{0:s}_out.json".format(args.tgt)
        with open(tgt_json) as json_file : tgt_data = json.load(json_file)
        print("  best_sigma: ref={0:6.2f} tgt={1:6.2f}".format(ref_data["best_sigma"],tgt_data["best_sigma"]))
        dPhi, phi_err, f_avg = extrapolatePhase(args)
        print("ref={0:s} tgt={1:s} dPhi={2:.5f} +/- {3:.5f} f_avg={4:.8f}".format(
            args.ref,args.tgt,dPhi,phi_err,f_avg))
