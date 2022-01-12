import os
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


# Step 1: combine all log data into one dataset and save dataset to file system
def combine_and_save_data(folder, outpath=r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/combined_datasets'):
    test_name = folder.split('/')[-1]
    logs = [l[2] for l in os.walk(folder)]
    if '.DS_Store' in logs[0]: logs[0].remove('.DS_Store')
    if 'combined_csv.csv' in logs[0]: logs[0].remove('combined_csv.csv')
    all_logs = [folder + '/' + l for l in logs[0]]
    combined_logs = pd.concat([pd.read_csv(l) for l in all_logs])
    datetime_string = combined_logs.date + " " + combined_logs.time # dd-mm-yy h:mm:ss
    # datetime_string.to_csv(outpath + "/" + test_name + "datetime.csv")
    combined_logs['Datetime'] = datetime_string
    combined_logs.dropna(subset = ["Datetime"], inplace=True)
    datetime_object = [datetime.strptime(d, '%d-%m-%Y %H:%M:%S') for d in combined_logs['Datetime']]
    combined_logs['Datetime'] = datetime_object
    datetime_object.sort()
    combined_logs.sort_values(by=['Datetime'],inplace=True)
    combined_logs.to_csv( outpath + "/" + test_name + ".csv", index=False, encoding='utf-8-sig')
    return combined_logs, datetime_object

# Step 2: get indices of StartV and EndV for characterization data
def get_StartEndIs(data, show_plot=False):
    v = data.BMU01_Cell_1_Voltage.to_numpy()
    # v = data.modbus_Voltage.to_numpy()    # for using module kernels
    # --------- convolution -----------
    # Do this once for each processing
    kernel_high = np.genfromtxt("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/kernel_high.csv",delimiter=',')
    length_high = len(kernel_high)
    kernel_low = np.genfromtxt("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/kernel_low.csv",delimiter=',')
    length_low = len(kernel_low)

    tmp_high = np.zeros(len(v) - length_high)
    for i in range(len(v) - length_high):
        tmp_high[i] = np.sum(np.abs(v[i:i+length_high] - kernel_high))
    tmp_low = np.zeros(len(v) - length_low)
    for i in range(len(v) - length_low):
        tmp_low[i] = np.sum(np.abs(v[i:i+length_low] - kernel_low))


    startpoints_high = signal.find_peaks(-tmp_high, height=-10,distance=100)
    # startpoints_high = signal.find_peaks(-tmp_high, height=-300,distance=100)
    startpoints_high = startpoints_high[0] + 399
    # startpoints_high = startpoints_high[0] + 103  # for using kernel_high_mod.csv
    startpoints_low = signal.find_peaks(-tmp_low, height=-10,distance=100)
    startpoints_low = startpoints_low[0] + 181
    # startpoints_low = startpoints_low[0] + 227    # for using kernel_low_mod.csv
    # ----------------------------------------

    # ------- get startIs and endIs
    def getStartEndVs(startpoints, highLow, z=20):
        V = np.zeros(len(startpoints))
        I = np.zeros(len(startpoints))
        z=20
        for i,s in enumerate(startpoints):
            v_tmp = v[s-z:s+z]
            if highLow=='high': spots = np.where(v_tmp==np.amax(v_tmp))[0]
            elif highLow=='low': spots = np.where(v_tmp==np.amin(v_tmp))[0]
            V[i] = v[s + (spots[-1]-z)]
            I[i] = s + (spots[-1]-z)
        return V, I
    [startVs, startIs] = getStartEndVs(startpoints_high, 'high')
    [endVs, endIs] = getStartEndVs(startpoints_low, 'low')

    # check that start and end points are in the right spot
    if show_plot:
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(v)
        ax.scatter(startIs,startVs,color='red')
        ax.scatter(endIs,endVs,color='red')
        plt.show()

    return startIs, endIs

# Step 3: get characterization dch ah
def getCharacDch(data, datetime_object, startI, endI):
    # ----------- get DCH1, CapDCH, DCH2
    # get Ah w/ trapezoidal integraion of current
    c = data.modbus_Current.to_numpy()
    p = data.modbus_Power.to_numpy()
    tmstmp = [t.timestamp()/3600 for t in datetime_object]
    ahDch = np.zeros(len(c)-1)
    whDch = np.zeros(len(p)-1)
    for i in range(len(c)-2):
        # calcualte ahDch for positive values of current only
        if c[i+1]>0: 
            ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(c[i+2]+c[i+1])
            whDch[i+1]= whDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(p[i+2]+p[i+1])
        else: 
            ahDch[i+1] = ahDch[i]
            whDch[i+1] = whDch[i]

    # fig, ax = plt.subplots(figsize=(8,6))
    # ax.plot(ahDch)
    # ax2=ax.twinx()
    # ax2.plot(c, color='red')
    # fig, ax = plt.subplots()
    # ax.plot(tmstmp)
    # plt.show()
    DCH1 = ahDch[startI]
    DCH2 = ahDch[-1] - ahDch[startI]
    CapDch = ahDch[endI]-ahDch[startI]

    WhDch1 = whDch[startI]
    WhDch2 = whDch[-1] - whDch[startI]
    return DCH1, DCH2, CapDch, WhDch1, WhDch2

# Step 4: get aging dch ah
def getAgingDch(data, datetime_object):
    c = data.modbus_Current.to_numpy()
    p = data.modbus_Power.to_numpy()
    tmstmp = [t.timestamp()/3600 for t in datetime_object]
    c[c<0] = 0      # we want to integrate discharge only, so disregard any positive current 
    p[p<0] = 0      # we want to integrate discharge only, so disregard any positive power 
    ahDch = np.trapz(c,tmstmp)
    whDch = np.trapz(p,tmstmp)
    return ahDch, whDch

# Step 5: get capacity
def getCap(v, CapDch, startI, endI, type='cell'):
    # convert start/end Vs to SOC
    soc_curve = pd.read_excel("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/SOCcurve.xlsx")
    soc = soc_curve.iloc[2:102,2].values.astype(float)
    if type=='cell': ocv = soc_curve.iloc[2:102,3].values.astype(float)
    else: ocv = soc_curve.iloc[2:102,4].values.astype(float)
    v_start = v[startI]
    v_end = v[endI]
    [startSOC, endSOC] = np.interp([v[startI], v[endI]], ocv, soc)

    # calculate cell capacity
    cap = 100/(startSOC-endSOC)*CapDch
    return cap, v_start, v_end


# User input: Char folder path and Aging folder path
char_folder = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/Char9'
aging_folder = ''
aging_folder = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/Aging12'


[data, datetime_object] = combine_and_save_data(char_folder)
[startIs, endIs] = get_StartEndIs(data, show_plot=True)
# Char1: n=1, Char2: n=0, Char3: no n
n=0 
if len(startIs)>n:
    startI = int(startIs[n])
else:
    # startI = 2781       # Char3
    # startI = 7860       # Char4
    # startI = 9542       # Char6
    # startI = 23926      # Char7
    startI = 6354       # Char8
    print('No startpoints found via convolution, startI was mannually assigned')
if len(endIs)>n: 
    endI = int(endIs[n])
else: 
    # endI = 4275     # Char2
    # endI = 5386     # Char3
    # endI = 10420    # Char4
    # endI = 4020     # Char5
    # endI = 11678    # Char6
    endI = 26057      # Char7
    endI = 8482       # Char8
    print('No endpoints found via convolution, endI was mannually assigned')

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(data['BMU01_Cell_1_Voltage'].to_numpy())
ax.scatter(startI,data['BMU01_Cell_1_Voltage'].iloc[startI],color='red')
ax.scatter(endI,data['BMU01_Cell_1_Voltage'].iloc[endI],color='red')
plt.show()

[DCH1, DCH2, CapDch, WhDch1, WhDch2] = getCharacDch(data, datetime_object, startI, endI)

if len(aging_folder)>0: 
    [aging_data, aging_datetime_object] = combine_and_save_data(aging_folder)
    [agingDch, agingWhDch] = getAgingDch(aging_data, aging_datetime_object)
    # add aging ahdch, whdch to dch1, WhDch1 from char 
    DCH1 += agingDch
    WhDch1 += agingWhDch

# loop through all modules
module_num = 8
cell_num = 12
module_caps = np.zeros(module_num)
cell_caps = np.zeros(module_num*cell_num)
module_Vs = np.zeros(module_num*2)
for m in range(module_num):      
    module = 'BMU0' + str(m+1)
    v_mod = data[module + '_CMA_Voltage'].to_numpy()
    [cap_mod, startv, endv] = getCap(v_mod, CapDch, startI, endI, type='module')
    module_caps[m] = cap_mod
    module_Vs[m] = startv
    module_Vs[m+8] = endv
    for c in range(cell_num):
        cell = module + '_Cell_' + str(c+1) + '_Voltage'
        v_cell = data[cell].to_numpy()
        [cap_cell, startv, endv] = getCap(v_cell, CapDch, startI, endI, type='cell')
        cell_caps[m*cell_num + c] = cap_cell

# concatenate into summary array
test_name = char_folder.split('/')
summary_updated = np.concatenate([np.array([test_name[-1], DCH1, DCH2, WhDch1, WhDch2]), module_caps, cell_caps])

# import summary csv, append new summary, save as csv
summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_summary.csv'
summary = pd.read_csv(summary_file)
summary.loc[len(summary)] = summary_updated
summary.to_csv(summary_file, index=False, encoding='utf-8-sig')

# concatenate CMA details into summary array
CMA_details_updated = np.concatenate([np.array([test_name[-1], CapDch]), module_Vs])

# import CMA details summary, append new summary, and save as csv
details_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_CMA_details.csv'
details = pd.read_csv(details_file)
details.loc[len(details)] = CMA_details_updated
details.to_csv(details_file, index=False, encoding='utf-8-sig')

    










# ------------------------- dev code --------------------------------------------------------------------------------
# file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/combined_datasets/Char_1.csv'
# data = pd.read_csv(file)

# v = data.BMU01_Cell_1_Voltage.to_numpy()
# # v = data.BMU05_Cell_1_Voltage.to_numpy()

# # --------- convolution -----------
# # Do this once for each processing
# kernel_high = np.genfromtxt("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/kernel_high.csv",delimiter=',')
# length_high = len(kernel_high)
# kernel_low = np.genfromtxt("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/kernel_low.csv",delimiter=',')
# length_low = len(kernel_low)

# tmp_high = np.zeros(len(v) - length_high)
# for i in range(len(v) - length_high):
#     tmp_high[i] = np.sum(np.abs(v[i:i+length_high] - kernel_high))
# tmp_low = np.zeros(len(v) - length_low)
# for i in range(len(v) - length_low):
#     tmp_low[i] = np.sum(np.abs(v[i:i+length_low] - kernel_low))

# startpoints_high = signal.find_peaks(-tmp_high, height=-10,distance=100)
# startpoints_high = startpoints_high[0] + 399
# startpoints_low = signal.find_peaks(-tmp_low, height=-10,distance=100)
# startpoints_low = startpoints_low[0] + 181
# # ----------------------------------------

# # ------- get startVs and endVs
# def getStartEndVs(startpoints, highLow, z=20):
#     V = np.zeros(len(startpoints))
#     I = np.zeros(len(startpoints))
#     z=20
#     for i,s in enumerate(startpoints):
#         v_tmp = v[s-z:s+z]
#         if highLow=='high': spots = np.where(v_tmp==np.amax(v_tmp))[0]
#         elif highLow=='low': spots = np.where(v_tmp==np.amin(v_tmp))[0]
#         V[i] = v[s + (spots[-1]-z)]
#         I[i] = s + (spots[-1]-z)
#     return V, I
# [startVs, startIs] = getStartEndVs(startpoints_high, 'high')
# [endVs, endIs] = getStartEndVs(startpoints_low, 'low')
# # check that start and end points are in the right spot
# fig, ax = plt.subplots(figsize=(8,6))
# ax.plot(v)
# ax.scatter(startIs,startVs,color='red')
# ax.scatter(endIs,endVs,color='red')
# # plt.show()
# # ----------------------------------------

# # ----------- get DCH1, CapDCH, DCH2
# n=1
# startV = startVs[n]
# endV = endVs[n]
# startI = int(startIs[n])
# endI = int(endIs[n])
# # get Ah w/ trapezoidal integraion of current
# v2 = data.BMU03_CMA_Voltage.to_numpy()
# c = data.modbus_Current.to_numpy()
# t = data.Datetime.to_numpy()
# datetime_object = [datetime.strptime(d, '%Y-%m-%d %H:%M:%S') for d in t]
# tmstmp = [t.timestamp()/3600 for t in datetime_object]
# ahDch = np.zeros(len(c)-1)
# for i in range(len(c)-2):
#     # calcualte ahDch for positive values of current only
#     if c[i+1]>0: 
#         ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(c[i+2]+c[i+1])
#     else: 
#         ahDch[i+1] = ahDch[i]

# DCH1 = ahDch[startI]
# DCH2 = ahDch[-1] - ahDch[startI]
# CapDch = ahDch[endI]-ahDch[startI]

# # convert start/end Vs to SOC
# soc_curve = pd.read_excel("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/SOCcurve.xlsx")
# soc = soc_curve.iloc[2:102,2].values.astype(float)
# ocv = soc_curve.iloc[2:102,3].values.astype(float)
# # test = np.interp(3.1, ocv, soc)
# [startSOC, endSOC] = np.interp([v2[startI], v2[endI]], ocv, soc)

# # calculate cell capacity
# cap = (startSOC-endSOC)/100*CapDch

# fig, ax = plt.subplots(figsize=(8,6))
# ax.plot(ah)
# ax.plot(c,color='red')
# ax2 = ax.twinx()
# ax2.plot(ahDch,color='black')
# plt.show()
