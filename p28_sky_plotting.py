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
Column 34: vrad in km/s
NEW-column 35: reference name (ie HdXXXXX or starname)

"""
# #########################################IMPORTS##################################################

import numpy as np
import matplotlib.pyplot as plt
import re
import astropy.coordinates as coord
import astropy.units as unitss
# ##########################################DATA####################################################

# Defines path to project28 overview csv
pro28datafile = "./OverviewProgram28_20190702.csv"

# imports data into np array
pro28data = np.loadtxt(pro28datafile, delimiter="\t", skiprows=1, dtype="object")


# creates column to contain 'reference name'
col35 = np.array([], dtype='str')

for row in range(pro28data.shape[0]):
    print(pro28data[row, 2])
    if len(pro28data[row, 2]) != 0:
        col35 = np.append(col35, pro28data[row, 2])

    else:
        col35 = np.append(col35, pro28data[row, 1])

# creates data variable to be used throughout program. Splices column 34 to the master csv file
data = np.column_stack((pro28data, col35))

########################################FUNCTIONS##################################################

def SNfilter(input_data, N, min_SNR):

    new_data = []
    for row in input_data:

        if np.isnan(row[[16, 18, 20, 22, 24]].astype(float)).any():
            print('skipping NaN SNR: - ', row[0])
            continue

        if len(np.where(row[[16, 18, 20, 22, 24]].astype(float) > min_SNR)[0]) > N -1:
            new_data.append(row)
    return np.array(new_data)


# defines the plotting function of the SN values of the program spectra
def plot_SNR(input_data):
    # defines the channels U,B,V,R,I
    U = input_data[:,16].astype('float')
    B = input_data[:,18].astype('float')
    V = input_data[:,20].astype('float')
    R = input_data[:,22].astype('float')
    I = input_data[:,24].astype('float')
    mv = input_data[:,4].astype('float')

    list_uniquenames = []
    list_mv = []

    for i in range(input_data.shape[0]):
        if input_data[i,35] not in list_uniquenames:
            list_uniquenames.append(input_data[i,35])
            list_mv.append(float(input_data[i,4]))

    list_mv = np.clip(np.array(list_mv),min(list_mv),9.3)
    mv = np.clip(mv,min(mv),9.3)

    fig = plt.figure()
    ax = fig.add_subplot(231)
    ax1 = fig.add_subplot(232)
    ax2 = fig.add_subplot(233)
    ax3 = fig.add_subplot(234)
    ax4 = fig.add_subplot(235)
    ax5 = fig.add_subplot(236)

    #ax9 = ax.twinx()

    # defines the bin range and steps for the histogram for each channel
    ax.hist(list_mv,bins = np.arange(-2,10,0.3),color='k',alpha=1,histtype='step',label='objects')
    #ax9.hist(mv,bins = np.arange(-2,10,0.3),color='k',alpha=1,histtype='step',label='exposures')

    ax1.hist(U,bins = np.arange(0,400,10),color='cyan',alpha=1,histtype='step',label='U')
    ax2.hist(B,bins = np.arange(0,400,10),color='blue',alpha=1,histtype='step',label='B')
    ax3.hist(V,bins = np.arange(0,400,10),color='green',alpha=1,histtype='step',label='V')
    ax4.hist(R,bins = np.arange(0,400,10),color='red',alpha=1,histtype='step',label='R')
    ax5.hist(I,bins = np.arange(0,400,10),color='black',alpha=1,histtype='step',label='I')
    # Creates the shaded area from 0-100 SN

    ax.set_xlabel('Mv')
    ax.set_ylabel('Number of objects')

    xlabels = ['U','B','V','R','I']
    for idx,axis in enumerate([ax1,ax2,ax3,ax4,ax5]):
        axis.axvspan(0,100,facecolor='red',alpha=0.1)
        axis.axvline(x=100,color='gray',linestyle='--')
        axis.set_xlim([0,400])
        axis.set_xlabel("S/N in "+xlabels[idx])
        axis.set_ylabel("Number of spectra")
    #plt.legend()
    plt.show()

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

    normal_spt = '^([OBAFGKM])([0-9.]{1,3})'
    dash_spt = '^([OBAFGKM])([0-9.]{1,3})([-+:])'
    lum_spt = '^([OBAFGKM])([0-9.]{1,3})([-/])([0-9.]{1,3})'

    total_spt = '('+normal_spt+'|'+dash_spt +'|'+lum_spt+')'

    expr_non = re.compile(normal_spt+'$')
    idx_non = np.array([expr_non.findall(i) for i in np_spectraltypes],dtype='bool')
    nonon = set(np_spectraltypes[idx_non])

    expr_I = re.compile(total_spt+'(([I]{1})([aA]|[bB])|([I]{1})$)')
    idx_I = np.array([expr_I.findall(i) for i in np_spectraltypes],dtype='bool')
    I = set(np_spectraltypes[idx_I])

    expr_II = re.compile(total_spt+'([I]{2}(?!I))')					# (\+|$)
    idx_II = np.array([expr_II.findall(i) for i in np_spectraltypes],dtype='bool')
    II = set(np_spectraltypes[idx_II])

    expr_III = re.compile(total_spt+'([I]{3})')
    idx_III = np.array([expr_III.findall(i) for i in np_spectraltypes],dtype='bool')
    III = set(np_spectraltypes[idx_III])

    expr_IV = re.compile(total_spt+'(IV)')
    idx_IV = np.array([expr_IV.findall(i) for i in np_spectraltypes],dtype='bool')
    IV = set(np_spectraltypes[idx_IV])

    expr_V = re.compile(total_spt+'(\()?(V)')
    idx_V = np.array([expr_V.findall(i) for i in np_spectraltypes],dtype='bool')
    V = set(np_spectraltypes[idx_V])

    maglist = [I,II,III,IV,V]
    letterlist = ['O','B','A','F','G','K','M']

    # Find Spectral classes
    idx_O = np.array([i[0] == 'O' for i in np_spectraltypes])
    O = set(np_spectraltypes[idx_O])
    idx_B = np.array([i[0] == 'B' for i in np_spectraltypes])
    B = set(np_spectraltypes[idx_B])
    idx_A = np.array([i[0] == 'A' for i in np_spectraltypes])
    A = set(np_spectraltypes[idx_A])
    idx_F = np.array([i[0] == 'F' for i in np_spectraltypes])
    F = set(np_spectraltypes[idx_F])
    idx_G = np.array([i[0] == 'G' for i in np_spectraltypes])
    G = set(np_spectraltypes[idx_G])
    idx_K = np.array([i[0] == 'K' for i in np_spectraltypes])
    K = set(np_spectraltypes[idx_K])
    idx_M = np.array([i[0] == 'M' for i in np_spectraltypes])
    M = set(np_spectraltypes[idx_M])
    idx_undef = np.array([i[0] not in letterlist for i in np_spectraltypes])
    undef = set(np_spectraltypes[idx_undef])

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

    for i in classes:

        print((np_spectraltypes[stats['spec']['idx'][i] & ~(stats['lum']['idx']['I']|stats['lum']['idx']['II']|
                              stats['lum']['idx']['III']|stats['lum']['idx']['IV']|stats['lum']['idx']['V'])]))

        print(len(np_spectraltypes[stats['spec']['idx'][i] & ~(stats['lum']['idx']['I']|stats['lum']['idx']['II']|
                              stats['lum']['idx']['III']|stats['lum']['idx']['IV']|stats['lum']['idx']['V'])]))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.bar(x,freq_V,width=width,color='darkgreen',edgecolor='white',linewidth=0.5,label='V')
    ax.bar(x,freq_IV,bottom=freq_V,width=width,color='limegreen',edgecolor='white',linewidth=0.5,label='IV')
    ax.bar(x,freq_III,bottom=freq_V + freq_IV,
           width=width,color='skyblue',edgecolor='white',linewidth=0.5,label='III')
    ax.bar(x,freq_II,bottom=freq_V + freq_IV + freq_III,
           width=width,color='royalblue',edgecolor='white',linewidth=0.5,label='II')
    ax.bar(x,freq_I,bottom=freq_V + freq_IV + freq_III + freq_II,
           width=width,color='darkblue',edgecolor='white',linewidth=0.5,label='I')
    ax.bar(x,freq_non, bottom= freq_V + freq_IV + freq_III + freq_II + freq_I,
           width=width,color='grey',edgecolor='white',linewidth=0.5,label='N/A')

    ax.set_ylabel("Number of objects")
    ax.set_xlabel("Spectral type")
    ax.set_xticks(x)#+width/2)
    ax.set_xticklabels(classes)
    for i,element in enumerate(freq_total):
        plt.text(x[i],5+element,'%i'%element,ha='center',va='center')
    #ax.set_title("Program 28: spectral type distribution")
    plt.legend()
    plt.show()


def galactic_plot(input_data):

    list_uniquenames = []
    list_ra = []
    list_dec = []
    list_mv = []
    list_bv = []

    for i in range(input_data.shape[0]):
        if input_data[i,35] not in list_uniquenames:

            list_uniquenames.append(input_data[i,35])
            list_ra.append(coord.Angle(input_data[i,8],unit=unitss.hour).degree)
            list_dec.append(coord.Angle(input_data[i,9],unit=unitss.deg).degree)
            list_mv.append(float(input_data[i,4]))
            list_bv.append(float(input_data[i,5]))

    ra = coord.Angle(list_ra*unitss.degree)
    ra = ra.wrap_at(180*unitss.degree)
    dec = coord.Angle(list_dec*unitss.degree)


    list_mv = np.clip(np.array(list_mv),min(list_mv),9.3)
    list_mv -= min(list_mv)
    list_mv /= max(list_mv)
    list_mv = 1. - list_mv
    list_mv *= list_mv
    list_mv *= 80
    list_mv = np.clip(list_mv,5,70)


    list_bv = np.clip(np.array(list_bv),min(list_bv),2.28)

    cmap_min = min(list_bv)
    cmap_max = max(list_bv)

    list_bv -= min(list_bv)
    list_bv /= max(list_bv)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection="mollweide")
    scatt1 = ax.scatter(ra.radian, dec.radian,s=list_mv,c='royalblue',marker='o')

    #plt.colorbar()

    ax.grid(True)
    ax.grid(zorder=0)
    ax.set_axisbelow(True)

    print(min(list_bv),max(list_bv))

    plt.show()
########################################MAIN_ROUTINE##############################################

print('number of exposures:', data.shape[0],'\n','number of unique targets', len(set(data[:,35])))
print('------------')
data = SNfilter(data,2,100)
print('------------')

print('number of exposures:', data.shape[0],'\n','number of unique targets', len(set(data[:,35])))



#plot_SNR(data)
#plot_spectraltype(data)
galactic_plot(data)
