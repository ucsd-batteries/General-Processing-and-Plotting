import numpy as np
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
import re
from openpyxl import load_workbook

def getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=False):
    # This function imports raw test data, calculates test results (discharge amp hours and cell capacites), 
    # and then appends the test results to an ongoing summary file
    
    # import data
    #data = pd.read_csv(data_file_path)
    data = pd.read_csv(data_file_path, error_bad_lines=False)
    #data2= pd.read_csv(data_file_path_2,error_bad_lines=False)     # use these lines to concatenate data if it was split in two parts  also have to make second file path at the bottom
    #data = data.append(data2)

    # check if csv has headers or not
    has_headers = True
    if data.columns[0] != 'Time':
        has_headers = False

    if has_headers:
        # remove rows with duplicate timestamps
        data = data.drop_duplicates(subset=['Time'])
        # pack voltage and current
        Vpack = data[' Pack Voltage'].to_numpy()
        Ipack = data[' Pack Current'].to_numpy()
        # get all cell voltages
        cell_voltages = data.iloc[:,15:15+cell_num].to_numpy()
    else: 
        # remove rows with duplicate timestamps
        data.drop_duplicates(subset=['0'])
        # pack voltage and current
        Vpack = data.iloc[:,2].to_numpy()
        Ipack = data.iloc[:,1].to_numpy()
        # get all cell voltages
        cell_voltages = data.iloc[:,3:3+cell_num].to_numpy()

    # create Step array - each 'Step' corresponds to a sequential command in the Digatron test setup 
    #    (eg "10 A for 5 seconds", "40 A for 10 minutes", etc)
    Step = []
    for i in range(len(Ipack)-1):
        if Ipack[i]!=0 and Ipack[i+1]==0:
            Step = np.append(Step,i+1)
        elif Ipack[i]==0 and Ipack[i+1]!=0:
            Step = np.append(Step,i+1)

    # get StartVs, EndVs 
    if has_headers:
        # change start_idx and end_idx as needed: originally 16 and 22 (16+6), respectively
        start_idx = 16
        end_idx = start_idx + 6
        StartVs = cell_voltages[int(Step[start_idx]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
        EndVs = cell_voltages[int(Step[end_idx]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge
    else:
        StartVs = cell_voltages[int(Step[2]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
        EndVs = cell_voltages[int(Step[3]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge

    # ---- Plotting in case results are defective ------- 
    # fig, ax = plt.subplots()
    # ax.plot(Ipack, label='current')
    # ax.plot(cell_voltages[:,0], label='voltage')
    # dots = [Ipack[int(s)] for s in Step]
    # ax.scatter(Step, dots, color='red')
    # ax.scatter(Step[start_idx], StartVs[0], color='black')
    # ax.scatter(Step[end_idx], EndVs[0], color='black')
    # ax.legend()
    # plt.show()


    # interpolate OCV-SOC curve to get SOCstart and SOCend values corresponding to StartVs and EndVs
    soc_curve = pd.read_csv(soc_curve_file)
    ocv = soc_curve.iloc[:,0].values.astype(float)
    soc = soc_curve.iloc[:,1].values.astype(float)
    [startSOC, endSOC] = np.interp([StartVs, EndVs], ocv, soc)   #multiple inputs gives multiple outputs

    # integrate discharge current to get discharge amp hours
    if has_headers:
        datetime_object = [dt.datetime.strptime(d[0:20]+d[24:28], '%a %b %d %H:%M:%S %Y') for d in data['Time']]
        print('start: ' + data['Time'].iloc[0])
        print('end: ' + data['Time'].iloc[-1])
        tmstmp = [t.timestamp()/3600 for t in datetime_object]
    else:
        tmstmp = data.iloc[:,0]
    ahDch = np.zeros(len(Ipack)-1)
    for i in range(len(Ipack)-2):
        if Ipack[i+1]>0:    # calcualte ahDch for positive values of current only (positive = discharge) #ahdch = amp hour discharge
            ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(Ipack[i+2]+Ipack[i+1]) #trapezoidal integration
        else: 
            ahDch[i+1] = ahDch[i]
    DCH1 = ahDch[int(Step[16])] #what is this?? --> how much was discharged up until the vs point
    #dch2 is then the discharge amount agter that
    
# ahdch is the total discharged over the whole graph
    
    #dch 1 and 2 are the ifferent discharger cylces, its used to distinguish capacity values for the state of health graph
    #(capacity drops after 1 discharge cycle, so on and so on)  soh is capacity over rated capacity
    #cycles is ah/rated capacity
    DCHTotal = ahDch[-1]  #total ah discharge
    DCH2 = DCHTotal - DCH1
    # C (const curr.) as specified by the associated Digatron test 
    CapDchAh = cc*(tmstmp[int(Step[start_idx+1])]-tmstmp[int(Step[start_idx])])   #cc is constant current, what current it was discharging at

    # calculate individual cell capacities
    Cap = np.zeros((len(startSOC,)))
    for i in range(len(Cap)):
        Cap[i] = (100/(startSOC[i]-endSOC[i]))*CapDchAh

    # determine test name
    test_name = data_file_path.split('_')
    test_name = test_name[-2] +'-' + test_name[-1].split('.')[0]

    if isCalendarAging: 
        # get NP number
        x = re.findall("NP\d", data_file_path)
        NPname = x[-1]
        NPname = "NP12"
        # import summary from correct sheet name
        summary = pd.read_excel(summary_file, sheet_name=NPname)
        test_name = 'char ' + str(summary.shape[0]+1)
        # get end_date
        end_date = datetime_object[-1]
        # get days_elapsed
        days_elapsed = (datetime_object[0] - summary['end_date'].iloc[-1]).days
        # import summary csv, append new summary, save as csv
        summary_updated = np.concatenate([np.array([test_name, end_date, days_elapsed, DCH1, DCH2]), Cap])
        summary.loc[len(summary)] = summary_updated
        #summary.to_excel(summary_file, sheet_name = NPname, index=False, encoding='utf-8-sig')

        # open excel sheet with ExcelWriter to avoid overwritting other sheets
        book = load_workbook(summary_file)
        writer = pd.ExcelWriter(summary_file, engine='openpyxl') 
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        summary.to_excel(writer, NPname, index=False)
        writer.save()
    else: 
        # import summary csv, append new summary, save as csv
        summary = pd.read_csv(summary_file)
        summary_updated = np.concatenate([np.array([test_name, DCH1, DCH2]), Cap])
        summary.loc[len(summary)] = summary_updated
        summary.to_csv(summary_file, index=False, encoding='utf-8-sig')

# ------------------------------------ General Nissan Info ------------------------------
# path to OCV-SOC csv file
soc_curve_file = r'C:/Users/amirs/Downloads/SOC_curve.csv'

# number of total cells in test (3 modules each with 6 cells)
cell_num = 16     

# Value of constant current during characteriations
cc = 20

# ------------------------------------ Nissan Pack 5 ------------------------------
# path to where raw test csv is stored
#r is like a backslash in other programs where it ignores what comes after the r' and treats it like text instead of a command
path = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/raw data/NP5/'

# name of csv file
data_file = r'cellvoltages_2022-01-03-12-31-01_NP5_char8.csv'
#data_file2 =r'cellvoltages_2021-09-27-16-08-23-NP5-Aging18_2.csv'
data_file_path = path + data_file
#data_file_path_2 = path + data_file2

# path to test summary file
summary_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/NP5_test_summary.csv'

getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)

# ------------------------------------ Nissan Pack 6 ------------------------------
# path to where raw test csv is stored
path = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/raw data/NP6/'

# name of csv file
data_file = r'cellvoltages_2022-01-03-12-29-36_NP6_char10.csv'
data_file_path = path + data_file

# path to test summary file
summary_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/NP6_test_summary.csv'

getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)

# ------------------------------------ Nissan Calendar Aging ------------------------------
# path to test summary file
summary_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/Calendar_Processed_Data.xlsx'

# path to where raw test csv is stored
path = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/raw data/Calendar/'

# NP7
data_file = r'cellvoltages_2021-12-09-09-34-49_NP7_100.csv'
data_file_path = path + 'NP7/' + data_file



# NP9
data_file = r'cellvoltages_2021-12-10-14-17-27_NP9_75 (1).csv'
data_file_path = path + 'NP9/' + data_file



# NP10
data_file = r'cellvoltages_2021-12-13-13-33-22_NP10_90.csv'
data_file_path = path + 'NP10/' + data_file


# NP12
data_file = r'cellvoltages_2021-12-11-16-39-30_NP12_50.csv'
data_file_path = path + 'NP12/' + data_file
data_file2 =r'cellvoltages_2021-12-12-16-27-08_NP12_50_2.csv'
data_file_path_2 = path +'NP12/' + data_file2

#getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)



# ------------------------- For Processing Multiple Files ---------------------
# folder = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Nissan/NP6/csvs/'
# # logs = [l[2] for l in os.walk(folder) if len(l[2])>len('.DS_Store')]
# logs = [l[2] for l in os.walk(folder)]
# if '.DS_Store' in logs[0]: logs[0].remove('.DS_Store')
# csvs = logs[0]

# for l in logs[0]:
#     data_file_path = folder + l
#     getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file,cc)
#     print(l)

# os.system('say "ding ding ding ding Im done ding ding ding ding Im done"')