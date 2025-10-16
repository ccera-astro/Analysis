# loop over all JSON files in the specified directory
# make a new set of JSON files that include vLSR information

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import glob 
import json 

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

    # Radial velocity correction to solar system barycenter
    vsrc = psrc.radial_velocity_correction(obstime = t, location = loc)

    # Projection of solar velocity towards the source
    vsun_proj = psrc.cartesian.dot(psun.cartesian)*vsun

    if verbose:
        print("Barycentric radial velocity: {0:+8.3f}".format(vsrc.to(u.km/u.s)))
        print("Projected solar velocity:    {0:+8.3f}".format(vsun_proj.to(u.km/u.s)))
    
    return vsun_proj-vsrc

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
for file in files :
    with open(file) as json_file: metadata  = json.load(json_file)
    vLSR =  -vlsr(metadata,geo_loc) 
    metadata['vlsr'] = vLSR.value 
    base_name = file.split("/")[-1]
    with open(out_dir + base_name, 'w') as fp : json.dump(metadata, fp)
