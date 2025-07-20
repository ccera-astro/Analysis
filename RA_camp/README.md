# CCERA RA Camp Analysis Software

## `anaCamp.py`

This script analyses a single 21 cm spectrum.    It shows: i) a raw spectrum; ii) a filtered spectrum; and iii) a fully-processed background-subtracted spectrum.

### Usage 

There are two ways to access files: by base name and by group name.   For base name

`python anaCampy.py -b Chnn_yyyy-mm-dd-hhhhhh`

`python anaCamp.py -n group_name --start_time yyyy-mm-dd-hhhh`

where `group_name` is the name of the telescope group--e.g., Mercury, Venus, etc.
The first spectrum after the value of `start_time` will be analyzed. 


## `sunScan.py`

This script plots total power versus time over a specified time interval.   It is aimed 
at displaying Sun scan data. 

### Usage 

`python plotScan.py --start_time yyyy-mm-dd-hhhh --stop_time yyyy-mm-dd-hhhh -d n -p`

where `-d 3` indicates that the power time series is downsampled (averaged) by `n` and 
the `-p` option indicates that a plot should be made.   The script looks through all 
files in the time range and sorts them by group.   A PDF file for each group is printed.

## `driftScan.py`

This script is similar to `sunScan.py` but makes a heatmap showing power spectral density 
as a function of time over a specified time interval.    

### Usage 

`python driftScan.py --start_time yyyy-mm-dd-hhhh --stop_time yyyy-mm-dd-hhhh -p`

where  the `-p` option indicates that a plot should be made.   The script looks through all 
files in the time range and sorts them by group.   A PDF file containing the results for 
all groups is printed. 

