#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from matplotlib.backends.backend_pdf import PdfPages 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--start_time",default="2025-06-28-00",help="Start time yyyy-mm-dd-hh")
    parser.add_argument("--stop_time", default="2025-06-29-00",help="Stop time yyyy-mm-dd-hh")
    parser.add_argument("-p","--plot",action="store_true",help="Plot to screen.")
    parser.add_argument("-n","--num_chan",type=int,default=4,help="Number of channels to read")
    return parser.parse_args()

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        #print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

def getFreqs(metadata) :
    dF = 1.e-6*metadata['srate']
    fMin = 1.e-6*metadata['freq'] - 0.5*dF
    fMax = 1.e-6*metadata['freq'] + 0.5*dF
    #print("getFreqs: fMin={0:e} fMax={1:e}".format(fMin,fMax))
    freqs = np.linspace(fMin,fMax,metadata['fft_size'])
    return freqs

def getVelocities(f) :
    f0, c = 1420.41, 3.0e5
    v = c*(f/f0 - 1.)
    return v

# fit spectrum to Chebyshev polynomial 
# restrict range of fit to |vDoppler| > vSignal 
def fitBackground(vDoppler,power,n,vSignal) :
    weights = np.ones_like(vDoppler)
    for i in range(len(vDoppler)) :
        if abs(vDoppler[i]) < vSignal : weights[i] = 1.e-6 
    series = np.polynomial.chebyshev.Chebyshev.fit(vDoppler, power, n, w=weights)
    background = series(vDoppler) 
    return background

def downSample(x,n):
    nIn = len(x)
    nOut = int(nIn/n)
    #print("In downsample(): nIn={0:d} nOut={1:d} nOut*n={2:d}".format(nIn,nOut,n*nOut))
    y = np.zeros(nOut)  
    for j in range(nOut) : y[j] = np.sum(x[j*n:(j+1)*n])
    y /= n 
    return y

def filterSpectrum(power,freqs,vDoppler) :
    fft_size, slice_size, q = len(power), 128, 0.90
    for i in range(0,fft_size,slice_size) :
        x = np.linspace(0.,float(slice_size),slice_size)
        slice = power[i:i+slice_size]
        q_value = np.quantile(slice,q)
        indices_1 = np.where(slice > q_value)[0]
        filtered_slice = np.delete(slice, indices_1)
        filtered_x = np.delete(x,indices_1)
        
        # fit the filtered array 
        coeffs = np.polyfit(filtered_x, filtered_slice, 1)
        fit_fcn = coeffs[1] + coeffs[0]*filtered_x

        # calculate the residuals with respect to the fit
        resid = filtered_slice - fit_fcn
        rms = np.std(resid)

        # redo the filter 
        fit_fcn = coeffs[1] + coeffs[0]*x
        resid = np.abs((slice-fit_fcn)/rms)
        indices = np.where(resid > 3.0)[0]
        filtered_slice = np.delete(slice, indices)
        filtered_x = np.delete(x,indices)

        # refit 
        coeffs = np.polyfit(filtered_x, filtered_slice, 1)
        fit_fcn = coeffs[1] + coeffs[0]*filtered_x

        for index in indices :
            ff = coeffs[1] + coeffs[0]*index 
            slice[index] = np.random.normal(ff,rms)

        power[i:i+slice_size] = slice 

    power = downSample(power,16)
    vDoppler = downSample(vDoppler,16)
    freqs = downSample(freqs,16)

    return power, freqs, vDoppler 

def anaSpectrum(base_name, inp) :
    #print("In anaSpectrum() base_name={0:s}".format(base_name))
    metadata = getMetaData(base_name + ".json")
    fft_size = metadata['fft_size']

    data_type = 'avg'
    if data_type == 'raw' :
        data, rows, cols = getData(base_name + "_1.raw" , metadata['fft_size'])
        #print("After getData Chan 1: rows={0:d} cols={1:d}".format(rows,cols))

        # reshape array into a series of row 
        data = np.reshape(data, (rows,cols))   
        power = np.mean(data,0)
    else :
        power = np.fromfile(base_name +"_{0:d}.avg".format(inp), dtype=np.float32)

    freqs = getFreqs(metadata) 
    vDoppler = getVelocities(freqs)

    power, freqs, vDoppler = filterSpectrum(power,freqs,vDoppler)

    vMin, vMax = -300., 300.
    i1 = np.searchsorted(vDoppler,vMin)
    i2 = np.searchsorted(vDoppler,vMax)
    #print("i1={0:d} i2={1:d}".format(i1,i2))
    freqs = freqs[i1:i2]
    vDoppler = vDoppler[i1:i2]
    power = power[i1:i2]

    background = fitBackground(vDoppler,power,5,200.)
  
    return vDoppler, power-background 

def getFiles(args,chan) :
    # Get a list of all files.  Make sure that they are ordered by time
    all_files = glob.glob("data/*.json")
    all_files.sort() 
    #print("In getFiles(): len(all_files)={0:d}".format(len(all_files)))

    # select files that match time, channel, and input 
    files = []  
    for file in all_files :
        file_chan = int(file.split("/Ch")[1][:2])
        #print("file={0:s} file_chan={1:d}".format(file,file_chan))
        if file_chan == chan :
            file_time = file.split("_")[1][0:13] 
            if file_time >= args.start_time and file_time < args.stop_time :
                #print("Keeping file={0:s}".format(file))
                files.append(file)
    return files 

#def vlsr(t,loc,psrc,verbose=False):
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

def makeKey(chan,inp) :
    return "Ch{0:02d}_{1:d}".format(chan,inp)

# begin execution here 

args = getArgs() 

# build dictionary of names
channel_lookup = {} 
for line in open("Channel_lookup.csv").readlines()[1:] :
    vals = line.strip().split(',')
    ch, inp, name = int(vals[0]), int(vals[1]), vals[2]
    channel_lookup[makeKey(ch,inp)] = vals[2]
print("channel_lookup={0:s}".format(str(channel_lookup)))

pdf = PdfPages("DriftScan_{0:s}.pdf".format(args.start_time))

for chan in range(args.num_chan+1) :
    files = getFiles(args,chan)
    #print("chan={0:d} len(files)={1:d}".format(chan,len(files)))
    if len(files) == 0 : break  
    base_name_0 =  files[0].strip(".json")
    meta_data = getMetaData(base_name_0 + ".json")

    for inp in [1,2] : 

        # analyse the first file to establish the parameters of the plot
        vDoppler, power = anaSpectrum(base_name_0,inp)
        nRows, nCols = len(files), len(vDoppler)
        print("Create mapData nRows={0:d} nCols={1:d}".format(nRows,nCols))
        mapData = 0.3*np.zeros((nRows,nCols))

        # Observatory location 
        lat, lon = 45.3491, -76.0413 
        geo_loc = EarthLocation.from_geodetic(lat = lat*u.deg, lon = lon*u.deg, height = 100*u.m)

        # Scan time information 
        firstMeta, lastMeta = getMetaData("./" + files[0]), getMetaData("./" + files[-1])
        scanDuration = (lastMeta['t_start'] - firstMeta['t_start'])/3600.

        times, vals, gain = [], [], 50./1.6 
        title = files[0].strip(".json") + "_{0:d}".format(inp)
        for row, file in enumerate(files) :
            base_name = "./" + file.strip(".json")
            vDoppler, power = anaSpectrum(base_name,inp)
            power *= gain 

            mapData[row] = np.maximum(power,0.) 
            if False and row % 4 == 0 :
                metadata = getMetaData(file)
                #vLSR = -vlsr(metadata,geo_loc,verbose=False)
                #idx = np.searchsorted(vDoppler,vLSR)
                times.append((metadata['t_start']-firstMeta['t_start'])/3600.)
                #vals.append(vDoppler[idx])
                vals.append(metadata['vlsr'])

        ky = makeKey(chan,inp) 
        fig = plt.figure(figsize=(10.5,8.))
        ax = fig.add_subplot(111)
        #ax.set_title("Drift Scan {0:s}".format(files[0].split('/')[-1].strip(".json")))
        ax.set_title("Drift Scan: {0:s} {1:s}\n {2:s} ".format(ky,channel_lookup[ky],args.start_time))
        ax.set_xlabel("Approach velocity (km/s)")
        ax.set_ylabel("Time (hours)")

        ax.patch.set_facecolor('white')

        im = ax.imshow(mapData,extent=[-300.,300.,scanDuration,0.],aspect='auto')
        im.set_cmap('jet')
        plt.colorbar(im, use_gridspec=True)
        plt.plot(vals,times,'w-')
        if args.plot : plt.show() 
        pdf.savefig(fig)

pdf.close() 