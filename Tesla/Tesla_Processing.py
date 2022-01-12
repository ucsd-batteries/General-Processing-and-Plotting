import numpy as np
import pandas as pd
import datetime as dt
import os
import re
from openpyxl import load_workbook

def getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=False):
    # This function imports raw test data, calculates test results (discharge amp hours and cell capacites), 
    # and then appends the test results to an ongoing summary file
    
    # import data
    data = pd.read_csv(data_file_path)

    # remove rows with duplicate timestamps
    data = data.drop_duplicates(subset=['Time'])

    # get all cell voltages
    cell_voltages = data.iloc[:,15:15+cell_num].to_numpy()

    # pack voltage and current
    Vpack = data[' Pack Voltage'].to_numpy()
    Ipack = data[' Pack Current'].to_numpy()

    # create Step array - each 'Step' corresponds to a sequential command in the Digatron test setup 
    #    (eg "10 A for 5 seconds", "40 A for 10 minutes", etc)
    Step = []
    for i in range(len(Ipack)-1):
        if Ipack[i]!=0 and Ipack[i+1]==0:
            Step = np.append(Step,i+1)
        elif Ipack[i]==0 and Ipack[i+1]!=0:
            Step = np.append(Step,i+1)

    # get StartVs, EndVs 
    StartVs = cell_voltages[int(Step[16]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
    EndVs = cell_voltages[int(Step[22]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge

    # interpolate OCV-SOC curve to get SOCstart and SOCend values corresponding to StartVs and EndVs
    soc_curve = pd.read_csv(soc_curve_file)
    ocv = soc_curve.iloc[:,0].values.astype(float)
    soc = soc_curve.iloc[:,1].values.astype(float)
    [startSOC, endSOC] = np.interp([StartVs, EndVs], ocv, soc)

    # integrate discharge current to get discharge amp hours
    datetime_object = [dt.datetime.strptime(d, '%a %b %d %H:%M:%S %Z %Y') for d in data['Time']]
    tmstmp = [t.timestamp()/3600 for t in datetime_object]
    ahDch = np.zeros(len(Ipack)-1)
    for i in range(len(Ipack)-2):
        if Ipack[i+1]>0:    # calcualte ahDch for positive values of current only (positive = discharge)
            ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(Ipack[i+2]+Ipack[i+1])
        else: 
            ahDch[i+1] = ahDch[i]
    DCH1 = ahDch[int(Step[16])]
    DCHTotal = ahDch[-1]
    DCH2 = DCHTotal - DCH1
    # based on CC (const curr.) as specified by the Digatron test 
    CapDchAh = cc*(tmstmp[int(Step[17])]-tmstmp[int(Step[16])])  

    # print test start and end dates for logging purposes
    print('start: ' + data['Time'].iloc[0])
    print('end: ' + data['Time'].iloc[-1]) 

    # calculate individual cell capacities
    Cap = np.zeros((len(startSOC,)))
    for i in range(len(Cap)):
        Cap[i] = (100/(startSOC[i]-endSOC[i]))*CapDchAh

    # determine test name
    test_name = data_file_path.split('_')
    test_name = test_name[-1].split('.')

    if isCalendarAging: 
        # get NP number
        x = re.findall("T1-\d+", data_file_path)
        TPname = x[-1]
        # import summary from correct sheet name
        summary = pd.read_excel(summary_file, sheet_name=TPname)
        test_name = 'char ' + str(summary.shape[0]+1)
        # get end_date
        end_date = datetime_object[-1]
        # get days_elapsed
        days_elapsed = (datetime_object[0] - summary['test_end_date'].iloc[-1]).days
        # import summary csv, append new summary, save as csv
        summary_updated = np.concatenate([np.array([test_name, end_date, days_elapsed, DCH1, DCH2]), Cap])
        summary.loc[len(summary)] = summary_updated
        # summary.to_excel(summary_file, sheet_name = NPname, index=False, encoding='utf-8-sig')

        # open excel sheet with ExcelWriter to avoid overwritting other sheets
        book = load_workbook(summary_file)
        writer = pd.ExcelWriter(summary_file, engine='openpyxl') 
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        summary.to_excel(writer, TPname, index=False)
        writer.save()

    else: 
        # import summary csv, append new summary, save as csv
        summary = pd.read_csv(summary_file)
        summary_updated = np.concatenate([np.array([test_name[0], DCH1, DCH2]), Cap])
        summary.loc[len(summary)] = summary_updated
        summary.to_csv(summary_file, index=False, encoding='utf-8-sig')

# ------------------------------------ General Tesla Info ------------------------------
# path to OCV-SOC csv file
soc_curve_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Tesla/Tesla_3s/TeslaOCVcurve.csv'

# Value of constant current during characteriations
cc = 40

# ------------------------------------ Tesla Aging ------------------------------
# path to test summary file
summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Tesla/Tesla_3s/Tesla_test_summary.csv'

# number of total cells in test (3 modules each with 6 cells)
cell_num=18     

# path to where raw test csv is stored
path = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Tesla/'

# name of csv file
data_file = r'cellvoltages_2021-05-03-14-19-59_T1-7-2-6_aging5.csv'
data_file_path = path + data_file

# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)


# ------------------------------------ Tesla Calendar Aging ------------------------------
# path to test summary file
summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Tesla/Tesla_Calendar/Calendar_Processed_Data.xlsx'
cell_num = 6

# path to where raw test csv is stored
path = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Tesla/Tesla_Calendar/'

# # T1-3
# data_file = r'cellvoltages_T1-3_50_char1.csv'
# data_file_path = path + 'T1-3/' + data_file
# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)
# print("T1-3 done")

# # T1-9
# data_file = r'cellvoltages_T1-9_75_char1.csv'
# data_file_path = path + 'T1-9/' + data_file
# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)
# print("T1-9 done")

# # T1-10
# data_file = r'cellvoltages_T1-10_90_char1.csv'
# data_file_path = path + 'T1-10/' + data_file
# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)
# print("T1-10 done")

# T1-12
data_file = r'cellvoltages_T1-12_100_char2-repeat.csv'
data_file_path = path + 'T1-12/' + data_file
getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)
print("T1-12 done")


# ------------------------- For Processing Multiple Files ---------------------
# folder = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Tesla/csvs/'
# # logs = [l[2] for l in os.walk(folder) if len(l[2])>len('.DS_Store')]
# logs = [l[2] for l in os.walk(folder)]
# if '.DS_Store' in logs[0]: logs[0].remove('.DS_Store')
# csvs = logs[0]

# for l in logs[0]:
#     data_file_path = folder + l
#     getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file)
#     print(l)

# os.system('say "ding ding ding ding Im done ding ding ding ding Im done"')