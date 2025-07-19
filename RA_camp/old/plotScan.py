#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

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

# Get a list of files.  Make sure that they are ordered by time
files = glob.glob("data/Ch00_2025-06-26*.json")
files.sort() 
print("files={0:s}".format(str(files)))
base_name_0 =  files[0].strip(".json")
meta_data = getMetaData(base_name_0 + ".json")

times, total_power = [], []
for row, file in enumerate(files) :
    base_name = "./" + file.strip(".json")
    data_type = 'avg'
    if data_type == 'raw' :
        data, rows, cols = getData(base_name + "_1.raw",meta_data['fft_size'])
        # reshape array into a series of row 
        data = np.reshape(data, (rows,cols))   
        power = (100./75000.)*np.mean(data,0)
        power = np.sum(data[0])
    else :
        power = np.fromfile(base_name +"_1.avg", dtype=np.float32)
        power = 0.00094*np.sum(power)

    hour = int(base_name.split("-")[-1][0:2]) 
    min = int(base_name.split("-")[-1][2:4])
    print("base_name={0:s} hour={1:d} min={2:d}".format(base_name,hour,min))
    times.append(hour + min/60.)
    
    total_power.append(np.sum(power))

plt.plot(times,total_power,'r.')
#plt.ylim(0.,80.)
plt.xlabel("t (hrs UTC)")
plt.ylabel("Total power (K)")
plt.title("Power vs. time: {0:s}".format(base_name_0))
plt.show() 
