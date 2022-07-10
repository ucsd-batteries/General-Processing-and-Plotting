from itertools import cycle
import numpy as np
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
import re
from openpyxl import load_workbook
# import matlab.engine #look at matlab documentation to install this module

def getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=False):
    # This function imports raw test data, calculates test results (discharge amp hours and cell capacites), 
    # and then appends the test results to an ongoing summary file
    
    # import data
    data = pd.read_csv(data_file_path, error_bad_lines=False)
    #data2= pd.read_csv(data_file_path_2,error_bad_lines=False)     # use these lines to concatenate data if it was split in two parts  also have to make second file path at the bottom
    #data = data.append(data2)
    rated_cap = 56.3
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
        offset =0 #if bms data is non standard
        offset2 =0
        # change start_idx and end_idx as needed: originally 16 and 22 (16+6), respectively
        start_idx = 16-offset
        end_idx = start_idx + 6 - offset2
        StartVs = cell_voltages[int(Step[start_idx]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
        EndVs = cell_voltages[int(Step[end_idx]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge
    else:
        StartVs = cell_voltages[int(Step[2]-1),:]  # get the voltage at the height of the CCCV curve after constant voltage
        EndVs = cell_voltages[int(Step[3]-1),:];  # get the voltage at the end of the CCCV curve after constant current discharge

    # ---- Plotting in case results are defective ------- 
    # fig, ax = plt.subplots()
    # # end_index = 300000
    # ax.plot(Ipack, label='current')
    # ax.plot(cell_voltages[:, 0], label='voltage')
    # dots = [Ipack[int(s)] for s in Step]
    # ax.scatter(Step, dots, color='red')
    # ax.scatter(Step[start_idx], StartVs[0], color='black')
    # ax.scatter(Step[end_idx], EndVs[0], color='black')
    # ax.legend()
    # plt.show()
    # return
    # ----------------------------------------------------

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
    DCH1 = ahDch[int(Step[start_idx])] #what is this?? --> how much was discharged up until the vs point
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
        NPname = re.findall("NP[0-9]+", data_file_path)[-1]
        print(f"Processing Calendar Aging {NPname}")

        # import summary from correct sheet name
        summary = pd.read_excel(summary_file, sheet_name=NPname)
        test_name = 'Characterization ' + str(summary.shape[0]+1)
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
        #Internal resistance
        drive_cycle_start = Step[25]
        drive_cycle_end = len(data)-1
        
        #remove first and last 10% of data
        t1 = int(0.1*(drive_cycle_end-drive_cycle_start))+int(drive_cycle_start)
        t2 = drive_cycle_end- int(0.1*(drive_cycle_end-drive_cycle_start))
        t = tmstmp[t1:t2]
        
        drive_cyc_cell_v = cell_voltages[t1:t2]
        drive_cyc_cell_i = Ipack[t1:t2]
        
        def ParamID(SOC_f,drive_cyc_cell_i,V,t):
            p = np.polyfit(soc/100,ocv,9) #fit OCV-SOC model with a 9-degree Power function polynomial. This degree is determined by a comprehensive comparison of the calculation and accuracy of polynomials with different degrees
            p = np.poly1d(p)
            N = len(t) # the data length
            '''The Thevenin battery model is applied to depict the battery's dynamic characteristics.
            Thus, Rp, Cp, R0 represents the polarization resistance, polarization capacitance, 
            the equivalent ohmic resistance, respectively.'''
            z = np.zeros((N,1)) #observation vector
            Phi = np.zeros((N-1,3)) #input and output sequences
            dt = 1 #sampling interval
            Up = np.zeros((1,N+1)) #Up means the polarization voltage. 
            error = np.zeros((1,N)) #voltage estimation error
            Lambda=1;# the forgetting factor, and its value usually from 0.95 to 1.The smaller the Â¦Ã‹, the stronger the tracking ability, but the greater the fluctuation.
            for i in range(N):
                if i== 0:
                    #observation vector=[OCV - measurement battery voltage]
                    z[i] = (p(SOC_f[0,i])-V[i])*Lambda**(N-i)
                else:
                    z[i] = (p(SOC_f[0,i])-V[i])*Lambda**(N-i)
                    #input and output sequences
                    Phi[i-1,:] = Lambda**(N-i)*np.array([z[i-1],drive_cyc_cell_i[i],drive_cyc_cell_i[i-1]])
            
                    
            #Parameter identification with RLS
            '''Theta indicates the parameters to be identified. These are the parameters 
             in the recursive least square (RLS) method instead of the battery parameters.  '''
            Theta = np.matmul(  np.matmul(  np.linalg.inv(np.matmul(np.transpose(Phi),Phi)),  np.transpose(Phi)),z[1:N])
            R0 = (Theta[1]-Theta[2])/(1+Theta[0])
            Rp = 2*(Theta[0]*Theta[1]+Theta[2])/((1-Theta[0])*(1+Theta[0]))
            Cp = (1+Theta[0])**2*dt/(4*(Theta[0]*Theta[1]+Theta[2]))

            
            
            #calculate voltage estimation error
            A = np.exp(-dt/(Rp*Cp))
            B = Rp * (1-np.exp(-dt/(Rp*Cp)))
            for j in range(N):
                Up[0,j+1] = A*Up[0,j]+B*drive_cyc_cell_i[j] #update Up
                error[0,j] = p(SOC_f[0,j])-Up[0,j+1]-int(R0)*drive_cyc_cell_i[j]-V[j] #calculate voltage estimation error according to the Thevenin battery model and OCV-SOCC model
                
            # mean absolute error of voltage estimation
            RMSE = np.sqrt(np.sum(error**2)/len(error))
            #root mean square error of voltage estimation
            MAE = np.mean(np.abs(error))
            
            
            return [R0,Rp,Cp,RMSE,MAE]
            
        
        R0_all = np.zeros((1,cell_num))
        Rp_all = np.zeros((1,cell_num))
        Cp_all = np.zeros((1,cell_num))
        RMSE_all = np.zeros((1,cell_num))
        MAE_all = np.zeros((1,cell_num))
        for cell in range(cell_num):
            # get current, voltage, and SOC
            V = drive_cyc_cell_v[:,cell]
            cap = Cap[cell]
            
            #integrate current
            ah = np.zeros((1,len(drive_cyc_cell_i)))
            for i in range(1,len(drive_cyc_cell_i)):
                ah[0,i] = ah[0,i-1]-.5*(t[i]-t[i-1])*(drive_cyc_cell_i[i]+drive_cyc_cell_i[i-1])
                
                
            SOC_f = ah/cap
            SOC_f = SOC_f + (abs(0.99-np.max(SOC_f)))
            [R0,Rp,Cp,RMSE,MAE] = ParamID(SOC_f,drive_cyc_cell_i,V,t)
            print(R0)
            

            
            
            
            R0_all[0,cell] =  R0
            Rp_all[0,cell] = Rp
            Cp_all[0,cell] = Cp
            RMSE_all[0,cell] = RMSE
            MAE_all[0,cell] = MAE
        R0_std = np.std(R0_all)
        
        
        '''
        R0_all = 16*[0]
        R0 = 0
        R0_std = 0
        #R0 = 0
        #R0_all = data.iloc[-1,31:31+16]
        '''
        # import summary csv, append new summary, save as csv

        summary = pd.read_csv(summary_file)
        summary_updated = np.concatenate([np.array([test_name, DCH1, DCH2]), Cap,np.array([np.mean(Cap/rated_cap)]),np.array([np.sum(R0_all),np.std(np.mean(data.iloc[int(Step[start_idx]):int(Step[start_idx+1]),15:31])),R0_std]),np.squeeze(R0_all)])
        summary.loc[len(summary)] = summary_updated
        summary.to_csv(summary_file, index=False, encoding='utf-8-sig')

# ------------------------------------ General Nissan Info ------------------------------
# path to OCV-SOC csv file
soc_curve_file = 'SOC_curve.csv'

# number of total cells in test (3 modules each with 6 cells)
cell_num = 16     

# Value of constant current during characteriations
cc = 20

# ------------------------------------ Nissan Pack 5 ------------------------------
# path to where raw test csv is stored
#r is like a backslash in other programs where it ignores what comes after the r' and treats it like text instead of a command
path = './'

# name of csv file
data_file = 'cellvoltages_2022-06-23-16-17-41_NP5_Aging45failed.csv'
data_file_path = path + data_file

# path to test summary file
summary_file = 'NP5_test_summary.csv'

# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)

# ------------------------------------ Nissan Pack 6 ------------------------------
# path to where raw test csv is stored
path = './'

# name of csv file
data_file = 'cellvoltages_2022-06-13-17-20-38_NP6_Char14.csv'
data_file_path = path + data_file

# path to test summary file
summary_file = 'NP6_test_summary.csv'

# getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)


# ------------------------------------ Nissan Pack 3 ------------------------------
# path to where raw test csv is stored
path = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/raw data/NP3/'

# name of csv file
data_file = r'cellvoltages_2021-03-03-12-37-22_NP3_Aging21.csv'
data_file_path = path + data_file

# path to test summary file
summary_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/NP3_test_summary.csv'

#getAhAndCaps(data_file_path, cell_num, soc_curve_file, summary_file, cc)

# ------------------------------------ Nissan Calendar Aging ------------------------------
# path to test summary file
summary_file = 'Calendar_Processed_Data.xlsx'

# path to where raw test csv is stored
path = './calendar_aging_data'

for file in os.listdir(path):
    calendar_data_path = os.path.join(path, file)
    if os.path.isfile(calendar_data_path):
        getAhAndCaps(calendar_data_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=True)
    # TODO: move the processed files from path to ./calendar_aging_data/processed


summary_file = './NP5_test_summary.csv'
path = './cycle_aging_data'
for file in os.listdir(path):
    cycle_data_path = os.path.join(path, file)
    print("Processing", cycle_data_path)


    if os.path.isfile(cycle_data_path):
        pack_num = re.findall("NP[0-9]+", cycle_data_path)[-1]
        print(pack_num)
        getAhAndCaps(cycle_data_path, cell_num, soc_curve_file, summary_file, cc, isCalendarAging=False) 
        print("Moving", path, "to", os.path.join(path, 'processed', file))
        os.rename(cycle_data_path, os.path.join(path, 'processed', file))
    else:
        print("Skipping", cycle_data_path)

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