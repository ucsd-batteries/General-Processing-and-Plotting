# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import datetime as dt



def getInternalResitance(filepath,cell_num,cc,soc_curve_file):
    R0 = np.zeros((1,cell_num))
    Rp_all =np.zeros((1,cell_num))
    Cp_all = np.zeros((1,cell_num))
    RMSE_all = np.zeros((1,cell_num))
    MAE_all = np.zeros((1,cell_num))
    data = pd.read_csv(data_file_path, error_bad_lines=False)
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
    #[steps, Caps, data] = getCellCap(data,has_headers,cell_num,StartVs,EndVs,Step)
    #cell_voltages = data.iloc[:,3:3+cell_num].to_numpy()
    if has_headers:
        datetime_object = [dt.datetime.strptime(d[0:20]+d[24:28], '%a %b %d %H:%M:%S %Y') for d in data['Time']]
        #print('start: ' + data['Time'].iloc[0])
        #print('end: ' + data['Time'].iloc[-1])
        tmstmp = [t.timestamp()/3600 for t in datetime_object]
    else:
        tmstmp = data.iloc[:,0]
    ahDch = np.zeros(len(Ipack)-1)
    for i in range(len(Ipack)-2):
        if Ipack[i+1]>0:    # calcualte ahDch for positive values of current only (positive = discharge) #ahdch = amp hour discharge
            ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(Ipack[i+2]+Ipack[i+1]) #trapezoidal integration
        else: 
            ahDch[i+1] = ahDch[i]
            
            
    soc_curve = pd.read_csv(soc_curve_file)
    ocv = soc_curve.iloc[:,0].values.astype(float)
    soc = soc_curve.iloc[:,1].values.astype(float)
    [startSOC, endSOC] = np.interp([StartVs, EndVs], ocv, soc) 
    
    CapDchAh = cc*(tmstmp[int(Step[start_idx+1])]-tmstmp[int(Step[start_idx])])   #cc is constant current, what current it was discharging at

    # calculate individual cell capacities
    Cap = np.zeros((len(startSOC,)))
    for i in range(len(Cap)):
        Cap[i] = (100/(startSOC[i]-endSOC[i]))*CapDchAh
        
        
        
    #New code
    #get data from drive cycle only
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
    
    return [R0_all,Rp_all,Cp_all,RMSE_all,MAE_all]
        
    
    
            
        
    

            
            
            

def getCellCap(data, has_headers,cell_num,StartVs,EndVs,Step):
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
    if has_headers:
        datetime_object = [dt.datetime.strptime(d[0:20]+d[24:28], '%a %b %d %H:%M:%S %Y') for d in data['Time']]
        #print('start: ' + data['Time'].iloc[0])
        #print('end: ' + data['Time'].iloc[-1])
        tmstmp = [t.timestamp()/3600 for t in datetime_object]
    else:
        tmstmp = data.iloc[:,0]
    ahDch = np.zeros(len(Ipack)-1)
    for i in range(len(Ipack)-2):
        if Ipack[i+1]>0:    # calcualte ahDch for positive values of current only (positive = discharge) #ahdch = amp hour discharge
            ahDch[i+1] = ahDch[i] + .5*(tmstmp[i+2]-tmstmp[i+1])*(Ipack[i+2]+Ipack[i+1]) #trapezoidal integration
        else: 
            ahDch[i+1] = ahDch[i]


if __name__ == '__main__':
#load data file
    soc_curve_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/SOC_curve.csv'
    
    # number of total cells in test (3 modules each with 6 cells)
    cell_num = 16     
    
    # Value of constant current during characteriations
    cc = 20


    path = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/raw data/NP5/'
    
    # name of csv file
    data_file = r'cellvoltages_2022-05-06-11-33-35_NP5_Aging41.csv'
    #data_file2 =r'cellvoltages_2021-09-27-16-08-23-NP5-Aging18_2.csv'
    data_file_path = path + data_file
    [R0_all,Rp_all,Cp_all,RMSE_all,MAE_all] = getInternalResitance(data_file_path,cell_num,cc,soc_curve_file)