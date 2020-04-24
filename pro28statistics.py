# Project 28 database statistics
# /mike/project28/...
# Mike Laverick

"""
Program28.csv info:
Column 0 : Unique index
Column 1 : Object name
Column 2 : HD catalog identifier
Column 3 : Hipparcos identifier
Column 4 : Visual magnitude M_V
Column 5 : Color B-V
Column 6 : Spectral Type
Column 7 : Comments = Simbad object type
Column 8 : Right Ascension hh:mm:ss.ss
Column 9 : Declination dd:
Column 10: Proper motion in RA
Column 11: Proper motion in dec
Column 12: Quality flag A = good
Column 13: Airmass
Column 14: Exposure Time
Column 15: Sky Quality 1 = very good 0 = dark clouds
Column 16: SN in U range
Column 17: U signal in 1 sec
Column 18: SN in B range
Column 19: B signal in 1 sec
Column 20: SN in V range
Column 21: V signal in 1 sec
Column 22: SN in R range
Column 23: R signal in 1 sec
Column 24: SN in I range
Column 25: I signal in 1 sec
Column 26: Total photomultiplier count
Column 27: Velocity correction
Column 28: Julian Date
Column 29: Program
Column 30: Observation Mode
Column 31: Observer
Column 32: Date
Column 33: Filename
Column 34: Vrad
NEW-column 35: reference name (ie HdXXXXX or starname)


        ###COMMENTS/TO DO/ETC###

duplicate issue between Sn and rejected fix this using the filter method

O and M type filtering does not work with the NONON search in SN plotting, BYPASSED by not allowing
spectype filter and spectype plot to run together (implement sub-spectype plot!)

"""
##########################################IMPORTS##################################################

import numpy as np
import matplotlib.pyplot as plt
import collections
import re
import pyfits
#import ivs.inout.fits as fits

###########################################DATA####################################################

# Defines path to project28 overview csv                                             ####################################
workdir = "./"                     #<----------------------------##  IF NOT MIKE THEN CHANGE THIS  ##
pro28datafile = workdir + "OverviewProgram28_20190702.csv" #MAKE SURE THIS EXISTS!!! ####################################

# imports data into np array
pro28data = np.loadtxt(pro28datafile,delimiter="\t",skiprows=1,dtype="object")

# creates column to contain 'reference name'
col35 = np.array([],dtype='str')

for row in range(pro28data.shape[0]):

    if pro28data[row,2] != ' ':
        col35 = np.append(col35,pro28data[row,2])

    else:
        col35 = np.append(col35,pro28data[row,1])

# creates data variable to be used throughout program. Splices column 35 to the master csv file
masterdata = np.column_stack((pro28data,col35))
data = masterdata
masterdata_rejected = np.array([])
data_rejected = masterdata_rejected

# creates list of unique names in csv file
uniquenamelist = []

for i in range(data.shape[0]):
    if data[i,35] not in uniquenamelist:
        uniquenamelist.append(data[i,35])
########################################FUNCTIONS##################################################

#--------------------------------------------------------------------------------------------------
# Defines function to filter data based upon minimum snr and minimum number of signal channels
def SNfilter(input_data):

    # default of 3 channels at 100sn
    nochannels = 3
    minsnr = 100
    # allows global variables data and data_rejected to be edited rather than creating local names
    global data
    global data_rejected
    #new lists to contain filtered data
    new_data = []
    new_data_rejected = []

    # probably a neater way to do this!! Checks each row for specified SN criteria, appends counter
    for row in range(input_data.shape[0]):

        nogreat=0
        if float(input_data[row,16]) >= minsnr:
            nogreat +=1
        if float(input_data[row,18]) >= minsnr:
            nogreat +=1
        if float(input_data[row,20]) >= minsnr:
            nogreat +=1
        if float(input_data[row,22]) >= minsnr:
            nogreat +=1
        if float(input_data[row,24]) >= minsnr:
            nogreat +=1

        # accepted = append.new_data, rejected = append.new_data_rejected
        if nogreat >= nochannels:
            new_data.append(input_data[row,:])
        else:
            new_data_rejected.append(input_data[row,:])
    # updates global data variables with SN filtered data
    data = np.asarray(new_data)
    data_rejected = np.asarray(new_data_rejected)

#--------------------------------------------------------------------------------------------------

def spectraltypefilter(input_data):
    # allows global variables data and data_rejected to be edited rather than creating local names
    global data
    global data_rejected
    #new lists to contain filtered data
    new_data = []
    new_data_rejected = []

    spectraltypechoice = 'NULL'
    letterlist = ['O','B','A','F','G','K','M']
    while spectraltypechoice == 'NULL':
        inputspectraltype = raw_input('Filter by spectral type: O,B,A,F,G,K or M?')
        for i in letterlist:
            if inputspectraltype == i in letterlist:
                for row in range(input_data.shape[0]):
                    if (input_data[row,6])[0] == i:
                        new_data.append(input_data[row,:])
                    else:
                        a = ' '
                spectraltypechoice = i
            else:
                a = ' '
    data = np.asarray(new_data)

#--------------------------------------------------------------------------------------------------

# defines the plotting function of the SN values of the program spectra
def plot_SNR(input_data):
    # defines the channels U,B,V,R,I
    U = input_data[:,16].astype('float')
    U = U[np.isfinite(U)]
    B = input_data[:,18].astype('float')
    B = B[np.isfinite(B)]
    V = input_data[:,20].astype('float')
    V = V[np.isfinite(V)]
    R = input_data[:,22].astype('float')
    R = R[np.isfinite(R)]
    I = input_data[:,24].astype('float')
    I = I[np.isfinite(I)]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # defines the bin range and steps for the histogram for each channel
    ax.hist(U,bins = np.arange(0,400,10),
            color='cyan',
            alpha=1,
            histtype='step',
            label='U')
    ax.hist(B,bins = np.arange(0,400,10),
            color='blue',
            alpha=1,
            histtype='step',
            label='B')
    ax.hist(V,bins = np.arange(0,400,10),
            color='green',
            alpha=1,
            histtype='step',
            label='V')
    ax.hist(R,bins = np.arange(0,400,10),
            color='red',
            alpha=1,
            histtype='step',
            label='R')
    ax.hist(I,bins = np.arange(0,400,10),
            color='black',
            alpha=1,
            histtype='step',
            label='I')
    # Creates the shaded area from 0-100 SN
    plt.axvspan(0,100,facecolor='red',alpha=0.1)
    plt.xlim([0,400])
    plt.xlabel("S/N")
    plt.ylabel("# of spectra")
    plt.legend()
    plt.title("Program 28: S/N distribution for observations")
    fig.savefig('SN_distrib.png',dpi=1200)

#--------------------------------------------------------------------------------------------------

# defines the plotting function for the number of observations per object
def plot_uniquestarobs(input_data):
    count1 = collections.Counter(input_data[:,35])
    count2 = np.asarray(count1.items())
    count3 = count2[:,1].astype('float')
    print len(count3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # defines the bin range and steps for the histogram for each channel
    ax.hist(count3,bins = np.arange(0,25,1),
            color='blue',
            alpha=1,
            histtype='step')
    plt.xlim([1,25])
    plt.xlabel("# of observations")
    #plt.yscale('log')
    plt.ylabel("# of unique objects")
    plt.title("Program 28: Number of observations per object")
    fig.savefig('no_obs_obj.png',dpi=1200)

#--------------------------------------------------------------------------------------------------

# Plot of spectral types etc
def plot_spectraltype(input_data):
    # creates lists of unique names and according spectral type
    list_uniquenames = []
    list_spectraltypes = []
    for i in range(input_data.shape[0]):
        if input_data[i,35] not in list_uniquenames:
            list_uniquenames.append(input_data[i,35])
            list_spectraltypes.append(input_data[i,6])

    np_spectraltypes = np.asarray(list_spectraltypes)

    expr_non = re.compile('^([OBAFGKM])([0-9.]{1,3})$')
    idx_non = np.array([expr_non.findall(i) for i in np_spectraltypes],dtype='bool')
    nonon = np_spectraltypes[idx_non]

    expr_I = re.compile('^([OBAFGKM])([0-9.]{1,3})([I]{1})([aA]|[bB])')
    idx_I = np.array([expr_I.findall(i) for i in np_spectraltypes],dtype='bool')
    I = np_spectraltypes[idx_I]

    expr_II = re.compile('^([OBAFGKM])([0-9.]{1,3})([I]{2}(?!I))')					# (\+|$)
    idx_II = np.array([expr_II.findall(i) for i in np_spectraltypes],dtype='bool')
    II = np_spectraltypes[idx_II]

    expr_III = re.compile('^([OBAFGKM])([0-9.]{1,3})([I]{3})')
    idx_III = np.array([expr_III.findall(i) for i in np_spectraltypes],dtype='bool')
    III = np_spectraltypes[idx_III]

    expr_IV = re.compile('^([OBAFGKM])([0-9.]{1,3})(IV)')
    idx_IV = np.array([expr_IV.findall(i) for i in np_spectraltypes],dtype='bool')
    IV = np_spectraltypes[idx_IV]

    expr_V = re.compile('^([OBAFGKM])([0-9.]{1,3})(\()?(V)')
    idx_V = np.array([expr_V.findall(i) for i in np_spectraltypes],dtype='bool')
    V = np_spectraltypes[idx_V]


    maglist = [I,II,III,IV,V]
    letterlist = ['O','B','A','F','G','K','M']

    # Find Spectral classes
    idx_O = np.array([i[0] == 'O' for i in np_spectraltypes])
    O = np_spectraltypes[idx_O]
    idx_B = np.array([i[0] == 'B' for i in np_spectraltypes])
    B = np_spectraltypes[idx_B]
    idx_A = np.array([i[0] == 'A' for i in np_spectraltypes])
    A = np_spectraltypes[idx_A]
    idx_F = np.array([i[0] == 'F' for i in np_spectraltypes])
    F = np_spectraltypes[idx_F]
    idx_G = np.array([i[0] == 'G' for i in np_spectraltypes])
    G = np_spectraltypes[idx_G]
    idx_K = np.array([i[0] == 'K' for i in np_spectraltypes])
    K = np_spectraltypes[idx_K]
    idx_M = np.array([i[0] == 'M' for i in np_spectraltypes])
    M = np_spectraltypes[idx_M]
    idx_undef = np.array([i[0] not in letterlist for i in np_spectraltypes])
    undef = np_spectraltypes[idx_undef]


    classes = np.array(['O','B','A','F','G','K','M','Other'])

    stats = {'spec': {'idx' : { 'O': idx_O, 'B': idx_B, 'A': idx_A, 'F': idx_F, 'G': idx_G, 'K': idx_K, 'M': idx_M, 'Other': idx_undef},
                      'O': O, 'B': B, 'A': A, 'F': F, 'G': G, 'K': K, 'M': M, 'Other': undef},
             'lum':  {'idx' : { 'I'  : idx_I, 'II' : idx_II, 'III': idx_III, 'IV' : idx_IV, 'V' : idx_V},
                      'I'  : I, 'II' : II, 'III': III, 'IV' : IV, 'V'  : V}}

    freq_I = np.array([(stats['spec']['idx'][i] & stats['lum']['idx']['I']).sum() for i in classes])
    freq_II = np.array([(stats['spec']['idx'][i] & stats['lum']['idx']['II']).sum() for i in classes])
    freq_III = np.array([(stats['spec']['idx'][i] & stats['lum']['idx']['III']).sum() for i in classes])
    freq_IV = np.array([(stats['spec']['idx'][i] & stats['lum']['idx']['IV']).sum() for i in classes])
    freq_V = np.array([(stats['spec']['idx'][i] & stats['lum']['idx']['V']).sum() for i in classes])
    freq_non = np.array([(stats['spec']['idx'][i] & ~(stats['lum']['idx']['I']|
                        stats['lum']['idx']['II']|
                        stats['lum']['idx']['III']|
                        stats['lum']['idx']['IV']|
                        stats['lum']['idx']['V'])).sum() for i in classes])
    freq_total = freq_I + freq_II + freq_III + freq_IV + freq_V + freq_non
    x = np.arange(len(classes))
    width = 1.0

    print 'Number of spectral counts in graph - ', sum(freq_total), ' / ', np_spectraltypes.shape

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(x,freq_non, width=width,color='grey',label='N/A')
    ax.bar(x,freq_I,bottom=freq_non, width=width,color='blue',label='I')
    ax.bar(x,freq_II,bottom=freq_non + freq_I,
           width=width,color='deepskyblue',label='II')
    ax.bar(x,freq_III,bottom=freq_non + freq_I+freq_II,
           width=width,color='palegreen',label='III')
    ax.bar(x,freq_IV,bottom=freq_non + freq_I+freq_II+freq_III,
           width=width,color='coral',label='IV')
    ax.bar(x,freq_V,bottom=freq_non + freq_I+freq_II+freq_III+freq_IV,
           width=width,color='darkred',label='V')
    ax.set_ylabel("# of unique stars")
    ax.set_xlabel("Spectral type of stars")
    ax.set_xticks(x)#+width/2)
    ax.set_xticklabels(classes)
    for i,element in enumerate(freq_total):
        plt.text(x[i],5+element,'%i'%element,ha='center',va='center')
    ax.set_title("Program 28: spectral type distribution")
    plt.legend()
    fig.savefig('spectral_distrib.svg')

#--------------------------------------------------------------------------------------------------

def make_type_csv(input_data,banana):

    typecsv = open(workdir+banana+"_type_stars.csv",'w')

    line = "  unseq  \t  starname  \t  HD \t HIP \t Mv \t BV \t Sptype \t comments \t ra \t dec \t pmra \t pmdec \t pmq \t airmass \t exptime \t skyq \t usn \t usignal1sec  \t  bsn \t bsignal1sec \t  vsn \t vsignal1sec \t  rsn \t rsignal1sec \t isn \t isignal1sec \t pmtotal \t bvcor \t bjd  \t program \t obsmode \t observer \t date \t filename \t vrad \t hom-starname \n"
    typecsv.write(line)
    for i in range(input_data.shape[0]):

      #if input_data[i,6][0] == banana:
      line = data[i,0] + '\t' + data[i,1] + '\t' + data[i,2] + '\t' + data[i,3] + '\t' + data[i,4] + '\t' + data[i,5] + '\t' + data[i,6] + '\t' + data[i,7] + '\t' + data[i,8] + '\t' + \
             data[i,9] + '\t' + data[i,10] + '\t' + data[i,11] + '\t' + data[i,12] + '\t' + data[i,13] + '\t' + data[i,14] + '\t' + data[i,15] + '\t' + data[i,16] + '\t' + data[i,17] + '\t' + \
             data[i,18] + '\t' + data[i,19] + '\t' + data[i,20] + '\t' + data[i,21] + '\t' + data[i,22] + '\t' + data[i,23] + '\t' + data[i,24] + '\t' + data[i,25] + '\t' + data[i,26] + '\t' + \
             data[i,27] + '\t' + data[i,28] + '\t' + data[i,29] + '\t' + data[i,30] + '\t' + data[i,31] + '\t' + data[i,32] + '\t' + data[i,33] + '\t' + data[i,34] + '\t' + data[i,35] + "\n"
      typecsv.write(line)
    typecsv.close()
#--------------------------------------------------------------------------------------------------

def mainfunction():
    global data
    global data_rejected
    SN_choice = raw_input('Do you wish to filter by SN? (Default currently 3 channels at 100+) enter  \'y\' or \'n\' ')
    if SN_choice == 'y':
        SNfilter(data)
        choice = raw_input('Do you want to work with rejected data? enter  \'y\' or \'n\' ')
        if choice == 'y':
            data = data_rejected

    csv_choice = raw_input('Do you wish to write out a spectral type filtered csv? enter \'y\' or \'n\' ')
    if csv_choice == 'y':
        sptype = raw_input('Which spectral type? \'O\',\'B\',\'A\',\'F\',\'G\',\'K\' or \'M\' ')
        make_type_csv(data, sptype)

    spectypechoice = raw_input('Do you wish to filter by spectral type? enter  \'y\' or \'n\' ')
    if spectypechoice == 'y':
        spectraltypefilter(data)
    SNRplot_choice = raw_input('Do you wish to plot SNRs? enter  \'y\' or \'n\' ')
    if SNRplot_choice == 'y':
        plot_SNR(data)
    uniqobschoice_choice = raw_input('Do you wish to plot no. of observations per unique star? enter  \'y\' or \'n\' ')
    if uniqobschoice_choice == 'y':
        plot_uniquestarobs(data)
    if spectypechoice != 'y':
        plotspectypechoice = raw_input('Do you wish to plot no. of stars by spectral type? enter  \'y\' or \'n\' ')
        if plotspectypechoice == 'y':
            plot_spectraltype(data)
    plt.show()

#--------------------------------------------------------------------------------------------------
mainfunction()
