#!/usr/bin/env python

import os
import glob
import aspicker
from multiprocessing import Pool


#===================================================#
# Parameters
#
eventpath       = "/home/weijias/disk1/archive_data/Tibet/Tibet/YA2013MIT_teleseismic/processedSeismograms/Event_2004_02_20_05_58_45/*BHZ*.sac"
wavepath        = "aspicking/waveforms"
unit_evdp       = "m"
phase           = "P"
filt            = False
freqmin         = 0.05
freqmax         = 1
sample_rate     = 40
snr_threshold   = 2.0

#
aqpath          = "aspicking/aqfiles"
refpath         = "aspicking/ref"
min_station     = 10
#
asname          = "./tcas"

# PGPLOT parameters
pgplotname      = "./aqplot"
figurepath      = "aspicking/figures"
device          = "/vps"
width           = 5.0
height_scale    = 1.0
cheight         = 1.0
line_thickness  = 1
label           = 0

#===================================================#

# initial processing
files = glob.glob(eventpath)
files.sort()
aspicker.init_process(sacfiles=files)

# converting SAC to AQ
files = glob.glob("aspicking/waveforms/20040220055845/*.sac")
files.sort()
aspicker.sac2aq(sacfiles=files, outpath=aqpath, outpath2=refpath)

# adaptive stacking
files = glob.glob(aqpath+"/rts*.aq")
for f in files:
    aspicker.run_as_fortran(aqfile=f, binname=asname)


# plotting .aq files
files = glob.glob(aqpath+"/*.aq")
for f in files:
    aspicker.pgplot(aqfile=f, outpath=figurepath, binname=pgplotname,
                     device=device, width=width, scale=height_scale,
                     cheight=cheight, line_thickness=line_thickness, label=label)

# assign picks
# for f in glob.glob(aqpath+"/*.ttr"):
#     aspicker.picks(ttrfile=f, wavepath=wavepath, refpath=refpath)



# path = "picks/20040220055845/"
# ttrfile = "ts1101600.ttr"
# sac2aq.picks(ttrfile, path)

