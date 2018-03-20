#!/usr/bin/env python

import os, glob
from obspy import read
from obspy import UTCDateTime
from obspy import Stream
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
from scipy.ndimage.interpolation import shift
import numpy as np
import subprocess


def reftime_from_trace(tr):

    year = tr.stats.sac.nzyear
    jday = tr.stats.sac.nzjday
    hour = tr.stats.sac.nzhour
    min  = tr.stats.sac.nzmin
    sec  = tr.stats.sac.nzsec
    msec = tr.stats.sac.nzmsec * 1000
    reftime = UTCDateTime(year=year, julday=jday, hour=hour, minute=min, second=sec, microsecond=msec)
    # print reftime
    return reftime

def eventid_from_sac(tr):

    """
    Get event ID from a given trace of obspy.

    Parameters
    ----------
    tr : type
         obpsy trace

    Returns
    -------
    eventid : str
              e.g., "20161208173846"

    """

    year = tr.stats.sac.nzyear
    jday = tr.stats.sac.nzjday
    hour = tr.stats.sac.nzhour
    min  = tr.stats.sac.nzmin
    sec  = tr.stats.sac.nzsec

    o = UTCDateTime(year=year, julday=jday, hour=hour, minute=min, second=sec) + tr.stats.sac.o
    o = str(o)
    o = o.replace("-", "")
    o = o.replace("T", "")
    o = o.replace(":", "")
    eventid = o[:14]

    return eventid

def get_travel_time(evdp, gcarc, phase):

    """
    Calculation traveltime for a given phase of an event with depth and epicentral distance.
    """

    model = TauPyModel(model='ak135')

    phase_list = [phase]
    arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=phase_list)
    arr = arrivals[0]
    return arr.name, arr.time

def init_process(sacfiles, outpath="aspicking/waveforms", unit_evdp="m", phase="P",
                 filt=False, freqmin=0.05, freqmax=1, sample_rate=40,
                 snr_threshold=2.0):

    for file in sacfiles:

        # read files in SAC format and calculation the ak135 traveltime
        tr = read(file)[0]
        evdp = tr.stats.sac.evdp
        if unit_evdp=="m":
            evdp /= 1000.0

        gcarc = tr.stats.sac.gcarc
        name, time = get_travel_time(evdp, gcarc, phase)

        # pre-processing data including detreand, filter, resampling
        tr1 = tr.copy()

        try:
            tr1.detrend(type="linear")
            tr1.detrend(type="demean")
        except:
            continue

        if filt:
            tr1.taper(type="cosine", max_percentage=0.05, side="both")
            tr1.filter(type="bandpass", freqmin=freqmin, freqmax=freqmax, zerophase=True)

        if tr.stats.sampling_rate!=sample_rate:
            tr1.resample(sampling_rate=sample_rate)


        reftime = reftime_from_trace(tr1)
        eventid = eventid_from_sac(tr1)

        # discard some bad traces in accordance with SNR
        ttak135 = reftime + tr.stats.sac.o + time
        tnoise1 = ttak135 - 50
        tnoise2 = ttak135 - 20
        tsigna1 = ttak135 - 10
        tsigna2 = ttak135 + 30

        tr2 = tr1.copy()
        tr3 = tr1.copy()
        tr2.normalize()
        tr3.normalize()
        trnoi = tr2.trim(starttime=tnoise1, endtime=tnoise2, fill_value=0)
        trsig = tr3.trim(starttime=tsigna1, endtime=tsigna2, fill_value=0)

        noi = trnoi.data
        sig = trsig.data

        sig = np.sqrt(np.sum(sig*sig))/len(sig)
        noi = np.sqrt(np.sum(noi*noi))/len(noi)

        if noi>0:
            snr = sig/noi
        else:
            continue

        if snr<snr_threshold:
            print file, snr
            continue

        t1 = ttak135 - 100
        t2 = ttak135 + 100

        tr1.trim(starttime=t1, endtime=t2, fill_value=0)
        tr1.stats.sac.t1 = tr.stats.sac.o + time
        tr1.stats.sac.kt1 = name
        tr1.stats.sac.user1 = snr
        if unit_evdp=="m":
            tr1.stats.sac.evdp /= 1000.0


        mypath = outpath + "/" + eventid
        if not os.path.exists(mypath):
            os.makedirs(mypath)

        bn = os.path.basename(file)
        fn = mypath + "/" + bn

        tr1.write(filename=fn, format="SAC")



def sac2aq(sacfiles, min_station=10, outpath="aspicking/aqfiles", outpath2="aspicking/ref"):

    try:
        os.makedirs(outpath)
    except:
        pass

    try:
        os.makedirs(outpath2)
    except:
        pass

    if len(sacfiles)<min_station:
        return

    # choose the reference trace
    sacfiles.sort()
    saclst = []
    tak135lst = []
    for sacfile in sacfiles:
        tr = read(sacfile)[0]
        saclst.append(tr)
        reftime = reftime_from_trace(tr)
        a = reftime.timestamp
        tak135lst.append(a + tr.stats.sac.t1)

    st = Stream(traces=saclst)
    # print st
    # print min(tak135lst), max(tak135lst)

    mid = 0.5 * (max(tak135lst) + min(tak135lst))
    tt = abs(tak135lst-mid)
    aa = tt.tolist()
    index = aa.index(min(tt))

    # min and max time differences
    tle = min(tak135lst) - tak135lst[index]
    tri = max(tak135lst) - tak135lst[index]
    # print tle, tri

    # write reference trace
    trref = st[index]
    eventid = eventid_from_sac(trref)
    aqfile = outpath + "/rts" + eventid + ".aq"
    fn = outpath2 + "/" + trref.id + "." + eventid + ".sac"
    trref.write(fn, format="SAC")

    # prepare for AQ waveform files
    ttak135ref = reftime_from_trace(trref) + trref.stats.sac.o + trref.stats.sac.t1

    tlen = 30
    t1 = ttak135ref - tlen

    starttime = t1
    year = starttime.year
    month = starttime.month
    day = starttime.day
    hour = starttime.hour
    minute = starttime.minute
    sec = starttime.second
    msec = starttime.microsecond

    second = sec + msec * 1.0e-6
    print starttime


    delta = trref.stats.sac.delta

    evla = trref.stats.sac.evla
    evlo = trref.stats.sac.evlo
    evdp = trref.stats.sac.evdp
    phase = trref.stats.sac.kt1

    tdiff = t1 - reftime_from_trace(trref) - trref.stats.sac.o

    fp = open(aqfile, "w")
    # number of stations
    fp.write("%d\n" % len(sacfiles))
    # lat, lon, evdp
    fp.write("%f\t%f\t%f\n" % (evla, evlo, evdp))
    # date
    fp.write("%d\t%d\t%d\n" % (year, month, day))
    # time and diff
    fp.write("%d\t%d\t%f\t%f\n" % (hour, minute, second, tdiff))
    # delta and phase
    fp.write("%f\t%s\n" % (delta, phase))

    for file in sacfiles:

        tr0 = read(file)[0]
        ttak135 = reftime_from_trace(tr0) + tr0.stats.sac.o + tr0.stats.sac.t1
        tshift = ttak135 - ttak135ref
        kstnm = tr0.stats.sac.kstnm
        print file, tshift

        # cut data ranging 10s before and 30s after the ak135 of each station
        t10 = ttak135 - tlen
        t20 = ttak135 + tlen
        tr0.trim(starttime=t10, endtime=t20, fill_value=0)

        npts = len(tr0.data)

        try:
            if tr0.stats.sac.t9 != -12345.0:
                flag = 0
        except:
            flag = 1

        tr0.normalize()
        data = tr0.data

        # make time shift to be compatable with AQ
        nshift = -tshift/delta
        data = shift(data, nshift)

        fp.write("%d\t%d\t%f\t%s\n" % (flag, npts, tshift, kstnm))
        for d in data:
            fp.write("%g " % (d))
        fp.write("\n")
    fp.close()

def run_as_fortran(aqfile, binname="./tcas"):

    niter = 10
    index = 3.0
    p = 1.15
    mnerr = 0.0
    mxerr = 15000.0
    winstart = 25.0
    winlengt = 20.0
    mndiff = -3.0
    mxdiff = 3.0

    os.system("cp %s ./" % aqfile)
    bn = os.path.basename(aqfile)
    path = os.path.dirname(aqfile)
    print bn, path

    fp = open("tcas.cmd", "w")
    fp.write("%d\n" % (niter))
    fp.write("%.1f\n" % (index))
    fp.write("%.4f\n" % (p))
    fp.write("%.1f\t%.1f\n" % (mnerr, mxerr))
    fp.write("%.1f\t%.1f\n" % (winstart, winlengt))
    fp.write("%s\n" % (bn))
    fp.write("%.1f\n%.1f" % (mndiff, mxdiff))
    fp.close()

    subprocess.call("./tcas", shell=True)

    # move the output file into the original path
    file1 = bn.replace("rts", "asi")
    file2 = bn.replace("rts", "asf")
    file3 = bn.replace(".aq", ".ttr")
    for f in [file1, file2, file3]:
        os.system("mv %s %s/" % (f, path))

    os.system("rm %s" % bn)

def pgplot(aqfile, outpath="aspicking/figures", binname="./aqplot",
           device="/vps", width=5.0, scale=1.0,
           cheight=1.0, line_thickness=1,
           label=0):

    try:
        os.makedirs(outpath)
    except:
        pass

    bn = os.path.basename(aqfile)
    os.system("cp %s ./" % (aqfile))

    with open(aqfile, "r") as ff:
        lines = ff.readlines()
    nstation = float(lines[0].strip())
    height = nstation * 0.25 * abs(scale)

    fp = open("aqplot.in", "w")
    fp.write("%s\n" % bn)
    fp.write("%s\n" % device)
    fp.write("%.2f\t%.2f\n" % (width, height))
    fp.write("%.2f\n" % cheight)
    fp.write("%d\n" % line_thickness)
    fp.write("%d\n" % label)
    fp.close()

    subprocess.call(binname, shell=True)
    fn = outpath + "/" + bn
    fn = fn.replace(".aq", ".ps")
    os.system("mv pgplot.ps %s" % (fn))
    os.system("rm %s" % bn)

def aq2sac(aqfile, outpath="ts"):

    with open(aqfile, "r") as fp:
        lines = fp.readlines()

    nsta = int(lines[0].strip())
    row = lines[1].split()


    evla = float(row[0])
    evlo = float(row[1])
    evdp = float(row[2])

    row1 = lines[2].split()
    row2 = lines[3].split()
    origin = UTCDateTime(int(row1[0]), int(row1[1]), int(row1[2]),
                         int(row2[0]), int(row2[1]), float(row2[2]))
    starttime = float(row2[3])
    # print origin, starttime

    o = origin - starttime
    # print o

    ostr = str(o)
    ostr = ostr.replace(":", "")
    ostr = ostr.replace("T", "")
    ostr = ostr.replace("-", "")
    ostr = ostr[:14]
    # print ostr

    op = outpath+"/"+ostr
    try:
        os.makedirs(op)
    except:
        pass

    row = lines[4].split()
    delta = float(row[0])
    phase = row[1]

    # print delta, phase

    for i in range(nsta):
        i1 = 2*i + 5
        i2 = i1+1
        # print lines[i1].strip(), i1, i2

        row = lines[i1].split()
        # print row
        # print row[0]
        flat = int(row[0])
        npts = int(row[1])
        tshift = float(row[2])
        kstnm = row[3]

        row = lines[i2].split()
        data = []
        for i in row:
            a = float(i)
            data.append(a)

        data = np.array(data)
        # print data.size

        header = {"evla":evla, "evlo":evlo, "evdp": evdp, "kstnm":kstnm}
        sac = SACTrace(data=data, **header)

        fn = op + "/%s.sac" % (kstnm)
        # sac.reftime = origin
        sac.reftime = o
        sac.o = 0
        sac.b = starttime + tshift
        sac.delta = delta
        sac.write(dest=fn)
        print fn

        pass

    pass


def picks(ttrfile, wavepath, refpath):

    with open(ttrfile, "r") as fp:
        lines = fp.readlines()

    # find initial and final time shift of the reference station
    bn = os.path.basename(ttrfile)
    eventid = bn[3:17]
    print eventid
    print wavepath
    print refpath
    f = refpath + "/*." + eventid + ".sac"
    fs = glob.glob(f)
    if len(fs)>1:
        print "Multiple reference station found for ", eventid
    elif len(fs)<1:
        print "No reference station found for ", eventid
    reffile = fs[0]
    row = reffile.split(".")
    refstnm = row[1]
    for line in lines[8:]:
        # print line.strip()
        row = line.split()
        kstnm = row[1]
        final = float(row[2])
        initi = float(row[4])

        if kstnm==refstnm:
            refinit = initi
            refinal = final
            break

    reftr = read(reffile)[0]
    refpick = reftime_from_trace(reftr) + reftr.stats.sac.t2
    print refpick

    for line in lines[8:]:
        # print line.strip()
        row = line.split()
        kstnm = row[1]
        final = float(row[2])
        initi = float(row[4])


        fn = wavepath + "/" + eventid + "/*.%s.*" % kstnm
        files = glob.glob(fn)
        fn = files[0]
        if len(files)!=1:
            continue

        # read trace
        tr = read(fn)[0]
        reftime = reftime_from_trace(tr)

        tr.stats.sac.t2 = refpick + refinit + refinal - final + initi - reftime
        tr.write(filename=fn, format="SAC")



if __name__=="__main__":

    # nsta = len(files)
    # print nsta
    #
    # # aq2sac(aqfile="out.aq")
    # aq2sac(aqfile="asf1101559.aq", outpath="asf")


    # files = glob.glob("Event_2004_02_20_05_58_45/*BHZ*.sac")
    # files.sort()
    #
    # init_process(sacfiles=files)

    files = glob.glob("mydata/20040220055845/*.sac")
    files.sort()
    sac2aq(sacfiles=files)

    pass