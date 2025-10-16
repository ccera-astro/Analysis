import numpy as np
from pint.models import get_model 
import polycos 
from pint.pulsar_mjd import data2longdouble 
import pint.toa as toa
import pint.logging 

def tofloat(str):	#Converts to float allowing 'double' notation 
	str = str.replace('D', 'E')
	return float(str)

def readpoly(filename):	# Reads the polynomial file output by tempo and returns it as a list of tuples,
	#print("Entering readpoly")
	f = open(filename, 'r')	# each tuple containing the coefficients for a give time period

	poly = list()
	while True:
		l0_items = f.readline().split()
		l1_items = f.readline().split()
		l2_items = f.readline().split()
		l3_items = f.readline().split()
		#print("In readpoly l0_items={0:s}".format(str(l0_items)))
		if len(l0_items) == 0:
			break
		#              name(0)              TMID(1)             DM(2)                 RPHASE(3)
		poly.append((l0_items[0], tofloat(l0_items[3]), tofloat(l0_items[4]), tofloat(l1_items[0]),
		#               F0(4)             Ncoeff(5)             coeff1(6)             coeff2(7)
			tofloat(l1_items[1]), int(l1_items[4]), tofloat(l2_items[0]), tofloat(l2_items[1]),
		#               coeff3(8)             coeff4(9)
			tofloat(l2_items[2]), tofloat(l3_items[0])))
		#print(l1_items[5])
		tfreq = tofloat(l1_items[5])
	f.close()
	return  poly, tfreq

def readpolycoeff(mjd,file_name): # Read the 'polyco.dat' file from tempo and finds the entry with the closest Tmid
    poly, tfreq   = readpoly(file_name)
    min_diff = 1e10
    ind = 0
    num = 0
    for coeff in poly:
        num += 1
        diff  = coeff[1] -  mjd
        if abs(diff) < min_diff:
            min_diff = abs(diff)
            best_coeff = coeff
            ind += 1

    #print('Leaving readpolycoeff: Best fit at ' + str(ind) + ' / ' + str(num))
    #print("best_coeff={0:s} ind={1:d} num={2:d} tfreq={3}".format(str(best_coeff), ind, num, tfreq))
    return ( best_coeff, ind, num, tfreq )


def getpolycoeff(mjd, metadata, base_name, file_name = None) :
    if not file_name : file_name = "./polycos/polyco_{0:s}_{1:s}.dat".format(base_name.split("/")[-1][:-5],metadata['target'])
    #print("Entering getpolycoeff mjd={0:f} file_name={1:s}".format(mjd,file_name))
    rfreq = 1.e-6*(metadata['freq'])
    pint.logging.setup(level="ERROR") 
    try :
        fp = open(file_name,'r')
        best_coeff, ind, num, tfreq = readpolycoeff(mjd,file_name)
        con1 = best_coeff[0] == metadata['target'][1:]
        con2 = ind != 0
        con3 = ind != num
        con4 = (abs(tfreq-rfreq)/tfreq < 0.01)
        #print("con1={0} con2={1} con3={2} con4={3}".format(con1,con2,con3,con4))
        file_valid = best_coeff[0] == metadata['target'][1:] and ind != 0 and ind != num and (abs(tfreq-rfreq)/tfreq < 0.01)
    except :
        par_file = "{0:s}.par".format(metadata['target']) 
        lines = open(par_file,'r').readlines() 
        #print("lines={0:s}".format(str(lines)))
        lines = lines[:-1]
        lines.append("TZRMJD          {0:f}\n".format(mjd)) 
        #print("lines={0:s}".format(str(lines))) 
        open(par_file,'w').writelines(lines)
        model = get_model(par_file)
        #print("Making new polyco.dat file: model={0:s}".format(str(model)))
        start, stop = (mjd - 0.5) // 1., (mjd + 1.5) // 1. 
        #print("   Call to Polycos: start={0} stop={1}\n XXXmodel={2:s}\nXXX".format(start,stop,str(model)))
        p = polycos.Polycos.generate_polycos(model, start, stop, "CCERA", 1440/4, 4, 1420.4)
        polycos.tempo_polyco_table_writer(p,file_name) 
        #print("In runPolyco.py: calling readpolycoeff({0},{1})".format(mjd,file_name))
        best_coeff, ind, num, tfreq = readpolycoeff(mjd,file_name)
        file_valid = best_coeff[0] == metadata['target'][1:] and ind != 0 and ind != num and (abs(tfreq-rfreq)/tfreq < 0.01)

    if not file_valid  :
        print("********* Unable to make valid polyco file ********")
        print("best_coeff[0]={0:s} name={1:s} ind={2:d} num={3:d} tfreq={4:f} rfreq={3:f}".format(
             best_coeff[0],metadata['target'][1:],ind,num,tfreq,rfreq))
        exit() 

    return best_coeff 

    

    


