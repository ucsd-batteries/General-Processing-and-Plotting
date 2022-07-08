import numpy as np
import numpy
import pandas as pd
import datetime as dt
import os
import re
from openpyxl import load_workbook
import matplotlib.pyplot as plt
import scipy as sp
import pytz

def EOLparam(data_file_path, cell_num, summary_file):
    # import data
    data = pd.read_csv(data_file_path)

    # remove rows with duplicate timestamps
    data = data.drop_duplicates(subset=['Time'])
    # datatime = pd.to_datetime(data['Time'], format='%a %b %d %H:%M:%S %Z %Y')
    # seconds = datatime.total_seconds()
    # print(seconds)
    # get all cell voltages

    Ipack = data[' Pack Current'].to_numpy()

    # create Step array - each 'Step' corresponds to a sequential command in the Digatron test setup
    #    (eg "10 A for 5 seconds", "40 A for 10 minutes", etc)
    Step = []
    for i in range(len(Ipack)-1):
        if Ipack[i]!=0 and Ipack[i+1]==0:
            Step = np.append(Step,i+1)
        elif Ipack[i]==0 and Ipack[i+1]!=0:
            Step = np.append(Step,i+1)

    # determine test name
    test_name = data_file_path.split('_')
    test_name = test_name[-1].split('.')

    # first pack in series
    hightemp = data.iloc[:, 13].to_numpy()
    lowtemp = data.iloc[:, 14].to_numpy()
    highave = np.mean(hightemp)
    lowave = np.mean(lowtemp)
    aveT = (highave + lowave)/2
    maxT = np.max(hightemp)
    lowT = np.min(lowtemp)

    cell_voltages = data.iloc[:, 15 + n:15 + n + cell_num].to_numpy()
    # --------to remove range of values---------------
    # print(cell_voltages.shape)
    #
    cell_v1 = cell_voltages[0:int(Step[8]), :]
    length = len(cell_voltages)-1
    cell_v2 = cell_voltages[int(Step[10]):length, :]
    cell_voltages = np.concatenate((cell_v1, cell_v2), axis=0)
    numpy.savetxt('new_aging_27.csv', cell_voltages)
    #
    # print(cell_voltages.shape)
    #
    I1 = Ipack[0:int(Step[8])]
    length = len(Ipack)-1
    I2 = Ipack[int(Step[10]):length]
    Ipack = np.concatenate((I1,I2),axis=0)
    #
    # # must recreate steps
    Step = []
    for i in range(len(Ipack) - 1):
        if Ipack[i] != 0 and Ipack[i + 1] == 0:
            Step = np.append(Step, i + 1)
        elif Ipack[i] == 0 and Ipack[i + 1] != 0:
            Step = np.append(Step, i + 1)
    # ---------------------end of range elimination--------------------
    cell_res = data.iloc[:, 33:33 + cell_num].to_numpy()
    res_ave_ind = np.mean(cell_res,axis=1)
    res_std = np.std(res_ave_ind)
    res_ave = np.mean(res_ave_ind)


    #cell voltage imbalance
    CCCV_voltage = np.array(cell_voltages[int(Step[16]-1):int(Step[22]-1), :])
    std_indiv = np.std(CCCV_voltage, axis=1)
    std_final = np.mean(std_indiv, axis=0)
    vmean_indiv = np.mean(CCCV_voltage, axis=1)
    vmean = np.mean(vmean_indiv, axis=0)
    std_std = np.std(std_indiv, axis=0)

    summary = pd.read_csv(summary_file)
    summary_updated = np.concatenate([np.array([test_name[0], maxT, lowT, aveT, std_final, vmean, res_ave, res_std])])
    summary.loc[len(summary)] = summary_updated
    summary.to_csv(summary_file, index=False, encoding='utf-8-sig')
# ------------------------------------ General Tesla Info ------------------------------
# path to test summary file
summary_file = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/Historic_Summary.csv'

# number of total cells in test (3 modules each with 6 cells)
cell_num=18

# set n to signify which battery is being processed in series
# if cycle aging or first in calendar
n = 0
# second in calendar
#n = 6
# third in calendar
#n = 12


# ------------------------------------ Tesla Aging ------------------------------
# path to where raw test csv is stored
path = r'/Users/liamk/OneDrive/Desktop/Data to Process/'

# name of csv file
data_file = r'Tesla_Aging_28.csv'
data_file_path = path + data_file
# data_file2 = r'Tesla_Aging16_2.csv'
# data_file_path_2 = path + data_file2

EOLparam(data_file_path, cell_num, summary_file)

# # ------------------------- For Processing Multiple Files ---------------------
# folder = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Tesla/Tesla_Cycle/csvs/'
# # logs = [l[2] for l in os.walk(folder) if len(l[2])>len('.DS_Store')]
# logs = [l[2] for l in os.walk(folder)]
# if '.DS_Store' in logs[0]: logs[0].remove('.DS_Store')
# csvs = logs[0]
#
# for l in logs[0]:
#     data_file_path = folder + l
#     EOLparam(data_file_path, cell_num, summary_file)
#     print(l)