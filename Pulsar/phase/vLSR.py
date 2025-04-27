# loop over all JSON files in the specified directory
# make a new set of JSON files that include vLSR information

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import glob 
import json 
import matplotlib.pyplot as plt 
import numpy as np 
import runPolyco as rp

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("-i","--in_dir",default="./data/",help="input directory")
    parser.add_argument("-o","--out_dir",default="./new/",help="output directory")
    parser.add_argument("-b","--base_name",default=None,help="base_name to convert")
    return parser.parse_args()

def vlsr(metadata,loc,verbose=False):
    """Compute the line of sight radial velocity

    psrc: SkyCoord object or source
    loc: EarthLocation object of observer
    t: Time object
    """
    tra, tdec = metadata['RA'], metadata['dec']
    psrc = SkyCoord(ra = tra, dec = tdec ,frame = "icrs", unit = (u.deg,u.deg)) 
    t = Time(metadata['t_start'],scale="utc",format="unix")

    # Direction of motion of the Sun. Info from
    psun = SkyCoord(ra = "18:03:50.29", dec = "+30:00:16.8",frame = "icrs",unit = (u.hourangle,u.deg))
    vsun = -20.0*u.km/u.s
    vsun = 0.*u.km/u.s

    # Radial velocity correction to solar system barycenter
    vsrc = psrc.radial_velocity_correction(obstime = t, location = loc)

    # Projection of solar velocity towards the source
    vsun_proj = psrc.cartesian.dot(psun.cartesian)*vsun

    if verbose:
        print("Barycentric radial velocity: {0:+8.3f}".format(vsrc.to(u.km/u.s)))
        print("Projected solar velocity:    {0:+8.3f}".format(vsun_proj.to(u.km/u.s)))
    
    return vsun_proj-vsrc

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
    # DT = (T-TMID)*1440
    # PHASE = RPHASE + DT*60*F0 + COEFF(1) + DT*COEFF(2) + DT^2*COEFF(3) + ....
    #FREQ(Hz) = F0 + (1/60)*(COEFF(2) + 2*DT*COEFF(3) + 3*DT^2*COEFF(4) + 
    pcoeff = (best_coeff[RPHASE] + best_coeff[coeff1], 60.0 * best_coeff[F0] + best_coeff[coeff2], best_coeff[coeff3], best_coeff[coeff4]) 
    phase = np.polynomial.polynomial.polyval(dt, pcoeff)
    fcoeff = ( best_coeff[F0] + (1./60.)*best_coeff[coeff2], (2./60.)*best_coeff[coeff3], (3./60.)*best_coeff[coeff4])
    freq = np.polynomial.polynomial.polyval(dt, fcoeff)
    return phase, freq 

# Observatory location 
lat, lon = 45.3491, -76.0413 
geo_loc = EarthLocation.from_geodetic(lat = lat*u.deg, lon = lon*u.deg, height = 100*u.m)

#in_dir, out_dir = "./data36/", "./new/"
args = getArgs() 
in_dir, out_dir, base = args.in_dir, args.out_dir, args.base_name 

if not base :
    files = glob.glob(in_dir + "*.json")
else :
    files = ['{0:s}/{1:s}'.format(in_dir,base) + '.json']

print("In add_vLSR.py: in_dir={0:s} files={1:s}".format(in_dir,str(files)))
f0 = 1.399541538720
c = 3.0e5
days, f_polycos, f_corrs = [], [], [] 
for file in files :
    month = int(file.split("/")[-1].split("-")[1]) 
    day = int(file.split("/")[-1].split("-")[2])
    if month == 1 : day += 31
    days.append(day)
    with open(file) as json_file: metadata  = json.load(json_file)
    results_name = "./outputs/" + file.split('/')[-1].strip(".json") + "_out.json"
    with open(results_name) as results_json : results = json.load(results_json)
 
    # set up polyco stuff
    
    MJD = (metadata['t_start']/ 86400.) + 40587.
    base_name = file.split("/")[-1].strip(".json")
    coeff = rp.getpolycoeff(MJD, metadata, base_name) 
    times = np.linspace(0.,960.,10000)   # a typical run is 16 minutes 
    MJDs = MJD + times/86400. 
    phase, freq = time2phase(MJDs, coeff)
    f_polyco_2 = freq[5000]
 
    f_polyco = 1./results['period']
    vLSR =  -vlsr(metadata,geo_loc) 
    f_corr = f0*(1. + (vLSR.value)/c)
    f_corrs.append(f_corr)
    f_polycos.append(f_polyco)
    f_diff = f_polyco - f_polyco_2
    print("file={0:s} vLSR={1:f} f_corr={2:f} f_polyco={3:f} f_diff={4:e} f_polyco_2={5:f}".format(
        file,float(vLSR.value),f_corr,f_polyco,f_diff,f_polyco_2)) 
    #base_name = file.split("/")[-1]
    #with open(out_dir + base_name, 'w') as fp : json.dump(metadata, fp)

plt.plot(days,f_corrs,'bo',label='f_corr')
plt.plot(days,f_polycos,'ro',label='f_polyco')
plt.xlabel('Day')
plt.ylabel('f (Hz)')
plt.legend()
plt.show()