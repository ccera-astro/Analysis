# prepare a set of pulsar transit scan analysis jobs  

import glob 

outLines = [] 
files = glob.glob("./inputs/*.json")
print("files={0:s}".format(str(files)))
for file in  files :
    base_name = file.split("/")[-1].strip(".json")
    print("base_name={0:s}".format(base_name))
    cmd = "python anaPulsar.py -b {0:s} -m phase --no_plot --high_pass_off --no_roll ".format(base_name)
    print("cmd={0:s}".format(cmd))
    outLines.append(cmd +"\n")

open("runAnaPulsar.csh",'w').writelines(outLines)



