import numpy as np
import pandas as pd
import datetime as dt
import os
import re
from openpyxl import load_workbook
import matplotlib.pyplot as plt

def getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=False):
    # This function imports raw test data, calculates test results (discharge amp hours and cell capacites), 
    # and then appends the test results to an ongoing summary file
    
    # import data
    data = pd.read_csv(data_file_path)
    # data2= pd.read_csv(data_file_path_2)
    # data = data.append(data2)

    # remove rows with duplicate timestamps
    data = data.drop_duplicates(subset=['Time'])

    # get all cell voltages

    # first pack in series
    n = 0
    cell_voltages = data.iloc[:, 15 + n :15 + n + cell_num].to_numpy()

    # second pack in series
    # n = 6
    # cell_voltages = data.iloc[:, 15 + n:15 + n + cell_num].to_numpy()

    # third pack in series
    # n = 12
    # cell_voltages = data.iloc[:, 15 + n:15 + n + cell_num].to_numpy()


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


    # --------to remove range of values---------------
    # print(cell_voltages.shape)
    #
    # cell_v1 = cell_voltages[0:int(Step[8]), :]
    # length = len(cell_voltages)-1
    # cell_v2 = cell_voltages[int(Step[10]):length, :]
    # cell_voltages = np.concatenate((cell_v1, cell_v2), axis=0)
    #
    # print(cell_voltages.shape)
    #
    # I1 = Ipack[0:int(Step[8])]
    # length = len(Ipack)-1
    # I2 = Ipack[int(Step[10]):length]
    # Ipack = np.concatenate((I1,I2),axis=0)
    #
    # # must recreate steps
    # Step = []
    # for i in range(len(Ipack) - 1):
    #     if Ipack[i] != 0 and Ipack[i + 1] == 0:
    #         Step = np.append(Step, i + 1)
    #     elif Ipack[i] == 0 and Ipack[i + 1] != 0:
    #         Step = np.append(Step, i + 1)
    # ---------------------end of range elimination--------------------

    # get StartVs, EndVs
    start_idx = 16
    end_idx = start_idx + 6
    StartVs = cell_voltages[int(Step[16]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
    EndVs = cell_voltages[int(Step[22]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge
    # interpolate OCV-SOC curve to get SOCstart and SOCend values corresponding to StartVs and EndVs
    soc_curve = pd.read_csv(soc_curve_file)
    ocv = soc_curve.iloc[:,0].values.astype(float)
    soc = soc_curve.iloc[:,1].values.astype(float)
    ocv2=ocv[::-1]
    soc2=soc[::-1]

    endSOC= np.interp(EndVs, ocv2, soc2)
    startSOC= np.interp(StartVs, ocv2, soc2)



    # integrate discharge current to get discharge amp hours
    datetime_object = [dt.datetime.strptime(d[0:20]+d[24:28], '%a %b %d %H:%M:%S %Y') for d in data['Time']]
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
    print('st')
    print(startSOC)
    print('end')
    print(endSOC)
    print('dch')
    print(CapDchAh)
    if isCalendarAging: 
        # get NP number
        x = re.findall("T1-\d+", data_file_path)
        TPname = x[0]
        # import summary from correct sheet name
        summary = pd.read_excel(summary_file, sheet_name=TPname)
        test_name = 'char ' + str(summary.shape[0]+1)
        # get end_date
        end_date = datetime_object[-1]
        # get days_elapsed
        previous_test = dt.datetime.strptime(summary['test_end_date'].iloc[0], '%Y-%m-%d ')
        days_elapsed = (datetime_object[0] - previous_test).days
        # days_elapsed = 0
        # import summary csv, append new summary, save as csv
        summary_updated = np.concatenate([np.array([test_name, end_date.strftime('%Y-%m-%d '), days_elapsed, DCH1, DCH2]), Cap])
        # summary_updated = np.concatenate([np.array([test_name, end_date.strftime('%-m/%-d/%Y'), days_elapsed, DCH1, DCH2]), Cap])
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


        # ---- Plotting in case results are defective -------
    fig, ax = plt.subplots()
    ax.plot(Ipack, label='current')
    ax.plot(cell_voltages[:, 0], label='voltage')
    dots = [Ipack[int(s)] for s in Step[17:20]]
    ax.scatter(Step[17:20], dots, color='red')
    # ax.scatter(Step[start_idx], StartVs[0], color='black')
    # ax.scatter(Step[end_idx], EndVs[0], color='black')

    ax.legend()
    plt.show()
# ------------------------------------ General Tesla Info ------------------------------
# path to test summary file
summary_file = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/Tesla_test_summary.csv'

# path to OCV-SOC csv file
soc_curve_file = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/TeslaOCVcurve_matt.csv'

# number of total cells in test (3 modules each with 6 cells)
cell_num=18     

# Value of constant current during characteriations
cc = 40

# ------------------------------------ Tesla Aging ------------------------------
# path to where raw test csv is stored
path = r'/Users/liamk/OneDrive/Desktop/Data to Process/'

# name of csv file
data_file = r'T1-7-2-6_Aging21.csv'
data_file_path = path + data_file
# data_file2 = r'Tesla_Aging16_2.csv'
# data_file_path_2 = path + data_file2

getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)


# ------------------------------------ Tesla Calendar Aging ------------------------------
# path to test summary file
summary_file = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/Tesla_test_summary_calendar.xlsx'
cell_num = 6

# path to where raw test csv is stored
path = r'/Users/liamk/OneDrive/Desktop/Data to Process/'

# # T1-3
# data_file = r'December_T1-3_50.csv'
# data_file_path = path + 'T1-3/' + data_file
# data_file2 = r'December_T1-3_50_2.csv'
# data_file_path_2 = path + 'T1-3/' + data_file2

# # T1-9
# data_file = r'December_T1-9_75.csv'
# data_file_path = path + 'T1-9/' + data_file

# T1-10
data_file = r'December_T1-10_90.csv'
data_file_path = path + 'T1-10/' + data_file

# T1-12
# data_file = r'Decemeber_T1-12_100.csv'
# data_file_path = path + 'T1-12/' + data_file
#
# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)


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