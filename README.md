# Analysis Tools

**WARNING:**  the analysis scripts in this repo are provided on an "as is" (or, "as was," since they are frequently being updated) basis.
They are constantly in flux and there may be bugs.   Users should view them  
of examples of the sorts of things that can be done, but should make no assumptions 
about their validity.

## Data format 
Data are stored in binary format.   This has the advantage of compact files that can be read quickly, but it also means that one needs to know the format in advance.    Information on the various file formats can be found in the `./DataFormats/' folder of this repository.

## Scripts of potential interest 
### `./Doppler/anaDoppler.py` 
plots the Doppler spectrum for a 21 cm observation.   A sample `.avg` data file has been included in the directory.   To view it enter 

`python anaDoppler.py -b 2026-02-19-2218 --data_dir "./"` 

### `./Pulsar/anaPulsar.py`
does various different analyses on pulsar data.  A sample analysis can be done with the command 

`python anaPulsar.py -b 2026-01-19-0036 -m scan --no_roll --data_dir ../../data/`

Note that the data files corresponding to 2026-01-19-0036 must be placed in a data directory somewhere on 
the local machine.    These files are available on the CCERA Google Drive.

