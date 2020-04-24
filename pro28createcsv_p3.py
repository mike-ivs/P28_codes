# Project 28 database statistics
# /mike/project28/...
# Mike Laverick

#Unused modules seem 'necessary' to ignore sqrt error for file 397377_reduced
#due to invalid value of sqrt(usignal) and sqrt(bsignal)
###################################################################################################

import sys,glob,os,pylab
from pylab import *
import astropy.io.fits as pyfits
from scipy import interpolate
import math
import numpy as np
from ivs.catalogs import hermes, sesame
import time

###################################################################################################

#lil fix to change directory depending on machine in use (/lhome != /home)
banana = 0
outputdir = 'NULL'

while (banana != 1):
  currentworkplace = input('laptop or ivs?')
  if currentworkplace == 'laptop':
    outputdir = '/lhome/mike/project28/'
    print('Output directory set to /lhome/mike/project28/')
    banana = 1
  elif currentworkplace == 'ivs':
    outputdir = '/home/mike/python/'
    print('Output directory set to /home/mike/python/')
    banana = 1
  else:
    print('not a valid input, try again')

###################################################################################################

# obsOnly = true => just look at Hermes obs list, don't extract info from Simbad
obsOnly=0

#quickfix list of names that are not matched to the SIMBAD database, these are
#manually updated (try and automate the flagging of problem stars)
problematicNames = {}
problematicNames["siriusa"] = "Sirius A"
problematicNames["bd+062036b"] = "BD +06 2036"
problematicNames["v*ocet"] = "V* omi Cet"
problematicNames["ccdmj16054-1948ac"] = "CCDM J16054-1948AB"
problematicNames["ads513ab"] = "pi And"  #HD 3369
problematicNames["bd144369a"] = "bet Del"
problematicNames["uvcam"] = "V* UV Cam"
problematicNames["gpori"] = "V* GP Ori"
problematicNames["ccdmj11479+0816ac"] = "ccdmj11479+0816a"
problematicNames["bd+144369a"] = "bet Del"
problematicNames["ccdmj00219-2300ab"] = "HD 1766A"
problematicNames["ccdmj08246-0108b"] = "HD 70923"
problematicNames["ccdmj08112-1256ae"] = "* 19 Pup"
problematicNames["ads3824abc"] = "*14 Aur"
problematicNames["ccdmj21446+2539b"] = "*kap Peg"
problematicNames["bd+244463a"] = "*kap Peg"
problematicNames["ads4241abc"] = "*sig Ori"
problematicNames["ads4182ab"] = "HD 36960"
problematicNames["NAME REGULUS"] = "Regulus"
problematicNames["ads2888ab"] = "eps Per"
problematicNames["gj332b"] = "10 uma"
problematicNames["ccdmj15416+1940a"] = "iot Ser"
problematicNames["bd-074876a"] = "HD 179002"
problematicNames["bd+231170a"] = "1 Gem"
problematicNames["gj188a"] = "m tau"
problematicNames["hd37742j"] = "HD 37742"
problematicNames["hd224635j"] = "hd224635"



###################################################################################################

#def's name and directory of output csv
today = time.strftime("%Y%m%d",time.localtime())
overview = open(outputdir+"OverviewProgram28_"+today+".csv",'w')

#header changes if data checked with SIMBAD database
if obsOnly:
 line = "  unseq  \t  starname  \t  program  \t  date  \t airmass \t skyq \t exptime  \t ra \t dec \t observer \t obsmode \t bvcor \t bjd \t usignal \t usignal1sec  \t  bsignal \t bsignal1sec \t  vsignal \t vsignal1sec \t  rsignal \t rsignal1sec \t isignal \t isignal1sec \t pmtotal \t filename \t vrad \n"
else:
 line = "  unseq  \t  starname  \t  HD \t HIP \t Mv \t BV \t Sptype \t comments \t ra \t dec \t pmra \t pmdec \t pmq \t airmass \t exptime \t skyq \t usn \t usignal1sec  \t  bsn \t bsignal1sec \t  vsn \t vsignal1sec \t  rsn \t rsignal1sec \t isn \t isignal1sec \t pmtotal \t bvcor \t bjd \t program \t obsmode \t observer \t date \t filename \t vrad \n"

overview.write(line)

###################################################################################################
#functions to convert angles etc
def findNth(s,sub,N,replaceString="XXX"):
    """ findNth(s,sub,N,replaceString="XXX")
    Find position of nth occurance of substring sub in string s
    """
    from string import find, replace
    if s.find(replaceString)>=0:
        print("Input contains replaceString, please change parameter replaceString")
        return
    else:
        return s.replace(sub,replaceString,N-1).find(sub) - (len(replaceString)-len(sub))*(N-1)

def deg2alpha(deg,sep=":"):
    """
    alpha from deg to hh:mm:ss
    """
    hh = int(deg / 15)
    mdeg = deg - 15*hh
    mm = int(mdeg * 4)
    sdeg = (mdeg*4) - mm
    ss = sdeg *60.
    if ss<10.:
        ss = ("0"+"%4.2f" % (ss))
    else:
        ss = ("%5.2f" % (ss))
    if ss=="010.00": ss = "10.00"
    return str(hh).zfill(2)+sep+str(mm).zfill(2)+sep+ss

def deg2delta(deg,sep=":"):
    """
    delta from deg to hh:mm:ss
    """
    dd = int(deg)
    mdeg = abs(deg - dd)
    #print dd, mdeg
    mm = int(mdeg * 60.)
    sdeg = (mdeg*60.) - mm
    #print mm, sdeg
    ss = sdeg * 60.
    #print ss
    if ss < 10:
        ss = ("0"+"%4.2f" % (ss))
    else:
        ss = ("%5.2f" % (ss))
    #
    if (deg >= 0.):
        dd = " " + str(dd).zfill(2)
    elif -10.<deg<=-1.:
        dd = '-'+str(abs(dd)).zfill(2)
    elif -1<deg<0.:
        dd = '-00'
    else:
        dd = str(dd)
    if ss=="010.00": ss = "10.00"
    return dd+sep+str(mm).zfill(2)+sep+ss

###################################################################################################
#Written by Pierre Royer

fullProgram = hermes.search(prog_ID=28)
namefits = sorted(fullProgram['filename'])

#fullProgram = np.loadtxt('../programming/nonP28_obs_of_P28_objs.tsv',delimiter='\t',dtype=object)
#namefits = sorted(fullProgram[:,-1])

#print (namefits)
#sys.exit()
# Normalization factor for perfect sky conditions (arbitrary). Intervenes in skyq
norm = 30000.

for i in np.arange(0,len(namefits)):
    dirname    = namefits[i][:namefits[i].rfind('/')]
    shortname  = namefits[i][namefits[i].rfind('/'):]
    try:
        fitsfile2d = dirname + shortname[:9]+"_HRF_OBJ_ext_CosmicsRemoved.fits"
    #outputname = string.join(namefits[i].split('/')[-3:]).replace(' ','/')
    #
        ext = pyfits.open(fitsfile2d)
        header = ext[0].header
    except:
        print('file not found')
        continue
    try:
        starname = str(header['OBJECT'])
        lowered = starname.lower().replace(" ","")
        if lowered in list(problematicNames.keys()):
            print()
            print("Problematic Name Identified --> Replacing ", starname, " with ", problematicNames[lowered])
            starname = problematicNames[lowered]
    except:
        starname = "not_known"
    unseq = str(header['UNSEQ'])
    #print i, "\t", namefits[i], "\t", fitsfile2d, "\t", unseq, "\t --"+starname+"--"
    print(i, "\t", namefits[i], "\t", unseq, "\t --"+starname+"--")#, header["RA"], header["RA"].__class__
    if ("PROG_ID" in header) and (header['PROG_ID'] !="") :
        program = str(header['PROG_ID'])
    elif ("PROGRAM" in header) :
        program = str(header['PROGRAM'])
    else:
        program = str(0)
    program = program.replace(" ","")
    if "TELALT" in header:
        secansz = 1/math.cos((90-float(header['TELALT']))*math.pi/180.)
        airmass = secansz * (1-0.0012*(secansz**2-1))
    else:
        airmass = 1.2
    airmass = "%5.3f"%airmass
    if "DATE-AVG" in header:
        date = str(header['DATE-AVG'])
    elif "DATE-END" in header:
        date = str(header['DATE-END'])
    else:
        date = str("no date known")
    exptime = float(header['EXPTIME'])
    if "PMTOTAL" in header:
        pmtotal = float(header['PMTOTAL'])
    else:
        pmtotal = 0.0
    usignal = np.median(ext[0].data[2050:2250,51])
    bsignal = np.median(ext[0].data[2250:2350,41])
    vsignal = np.median(ext[0].data[2350:2450,24])
    rsignal = np.median(ext[0].data[2450:2650,15])
    isignal = np.median(ext[0].data[2300:3100,5])
    usignal1sec = "%9.2f"%(usignal/exptime)
    bsignal1sec = "%9.2f"%(bsignal/exptime)
    vsignal1sec = "%9.2f"%(vsignal/exptime)
    rsignal1sec = "%9.2f"%(rsignal/exptime)
    isignal1sec = "%9.2f"%(isignal/exptime)
    usn = "%7.1f"%sqrt(usignal)
    bsn = "%7.1f"%sqrt(bsignal)
    vsn = "%7.1f"%sqrt(vsignal)
    rsn = "%7.1f"%sqrt(rsignal)
    isn = "%7.1f"%sqrt(isignal)
    #
    try:
        obsmode = str(header['OBSMODE'])
        observer= str(header['OBSERVER'])
        ra      = deg2alpha(float(header["RA"]))
        dec     = deg2delta(float(header["DEC"]))
    except:
        print('KEYWORD missing')
        continue
    try:
        bvcor   = "%11.6f"%float(header['BVCOR'])
    except:
        bvcor  = " "
    try:
        bjd     = "%14.6f"%float(header['BJD'])
    except:
        bjd  = " "
    #
    if obsOnly:
      line = unseq + "\t" + starname + "\t" + program + "\t" + date + "\t" + str(airmass) + "\t" + str(exptime) + "\t" + \
        ra + "\t" + dec + "\t" + observer + "\t" + obsmode + "\t" + bvcor + "\t" + bjd + "\t" + \
        str(usn) + "\t" + str(usignal1sec) + "\t" + str(bsn) + "\t" + str(bsignal1sec) + "\t" + str(vsn) + "\t" + \
        str(vsignal1sec)+  "\t" + str(rsn) + "\t" + str(rsignal1sec)+ "\t"+  str(isn) + "\t" + str(isignal1sec) + "\t" + str(pmtotal)+ "\t" + namefits[i] + "\n"
      overview.write(line)
    else:
      target = sesame.search(starname)
      if target:
          alias = target["alias"]
          #print target["Vel"]["v"]
          try:
            hd = [n for n in alias if n[:2]=="HD"][0]
          except:
            hd = " "
          try:
            hip = [n for n in alias if n[:3]=="HIP"][0]
          except:
            hip = " "
          try:
            Mv = float(target['mag']['V']['v'])
            vref = 10.**(-0.4 * Mv)
            skyq = "%7.5f" % (float(vsignal1sec) / vref / norm)
          except:
            Mv = "90.0"
            skyq = "0.0"
          try:
            Mb = float(target['mag']['B']['v'])
            bv = Mb - Mv
          except:
            Mb = "90.0"
            bv = "90.0"
          try:
            pmra  = str(target['pm']['pmRA'])
            pmdec = str(target['pm']['pmDE'])
            pmq   = str(target['pm']['q'])
          except:
            pmra, pmdec, pmq = " "," "," "
          try:
            sptype = str(target['spType'])
          except:
            sptype = " "
          try:
            otype = str(target['otype'])
          except:
            otype = " "
          try:
            vrad = str(target['Vel']['v'])
          except:
            vrad = ' '

      else:
         print(i, "\t", namefits[i], "\t", fitsfile2d, "\t --"+starname+"-- \t", lowered, "-- NO SIMBAD MATCH")
         hd = hip = pmra = pmdec = pmq = sptype = otype = vrad = " "
         Mv = Mb = bv = "90.0"
      line = unseq + "\t" + starname + "\t" + hd + "\t" + hip + "\t" + str(Mv) + "\t" + str(bv) + "\t" + sptype  + "\t" + otype + "\t" + \
        ra + "\t" + dec + "\t" +  pmra + "\t" + pmdec + "\t" + pmq + "\t" + str(airmass) + "\t" + str(exptime) + "\t" + str(skyq) + "\t" + \
        str(usn) + "\t" + str(usignal1sec) + "\t" + str(bsn) + "\t" + str(bsignal1sec) + "\t" + \
        str(vsn) + "\t" + str(vsignal1sec) + "\t" + str(rsn) + "\t" + str(rsignal1sec) + "\t" + str(isn) + "\t" + str(isignal1sec) + "\t" + \
        str(pmtotal)+ "\t" + bvcor + "\t" + bjd + "\t" + program + "\t" + obsmode + "\t" + observer+ "\t" + date + "\t" + namefits[i]+ '\t' + vrad + "\n"
      overview.write(line)
    ext.close()

overview.close()
