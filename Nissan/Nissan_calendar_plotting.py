from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from numpy.lib.function_base import append
import pandas as pd
import datetime as dt

"""Plots the processed data in line graphs to show battery degradation"""

summary_file = './Calendar_Processed_Data.xlsx'

save_plot = True   # Change to True to save plot as .png

pack_names = ['NP7', 'NP9', 'NP10', 'NP12']
charge_levels = ['100%', '75%', '90%', '50%']
legend_loc = 'lower left'
rated_cap = 56.3

# used to make size of days and avg_soh matrixes
highestNumberOfChars = 0
for i, pack in enumerate(pack_names):
    tpnumber = pd.read_excel(summary_file, sheet_name=pack,index_col='Test').shape[0]
    if tpnumber > highestNumberOfChars:
        highestNumberOfChars = tpnumber
        

# makes a list of the packs that have fewer characterizations than the max number of characterizations        
fewerCharList = []
for i, pack in enumerate(pack_names):
    tpnumber = pd.read_excel(summary_file, sheet_name=pack,index_col='Test').shape[0]
    if tpnumber < highestNumberOfChars:
        fewerCharList.append(i)      


#initialize the days and avg_soh arrays 
days = np.zeros((4,highestNumberOfChars))
avg_soh = np.zeros((4,highestNumberOfChars))

# goes through each pack on the csv
for i, pack in enumerate(pack_names):
    summary_data = pd.read_excel(summary_file, sheet_name=pack,index_col='Test')
    for j in range(1,summary_data.shape[0]+1):
        avg_soh[i][j-1] = sum(summary_data.loc["Characterization "+str(j)][4:20])/16/rated_cap*100 #goes through ever characterization and sums up the soh values
        if j >= 1:
            days[i][j-1] = summary_data["days"][j-1]+days[i][j-2]   #adds the length of the gaps together to get cumulative day count for graphing
        else:
            days[i][j-1] = summary_data["days"][j-1]  

#  ---------------------------- Pack Comparison ----------------------------
fs = 15     # universal fontsize
outpath = './plots/'    # path for saving plot

fig, ax = plt.subplots(figsize=(12,6))
for i in range(4):
    if i in fewerCharList:
        ax.plot(days[i][0:-1], avg_soh[i][0:-1], label=pack_names[i] + ' ' + charge_levels[i], marker='.', markersize=12) #does not graph the last characterization if there are fewwer characterizations
    else:
        ax.plot(days[i], avg_soh[i], label=pack_names[i] + ' ' + charge_levels[i], marker='.', markersize=12)

ax.set_title('Nissan Calendar Aging Pack Comparison', fontsize=fs)
ax.set_ylabel('Pack Average SOH [%]', fontsize=fs)
ax.set_xlabel('Time Elapsed [days]', fontsize=fs)
ax.tick_params(direction='in')
plt.ylim([57,67])
ax.legend()
if save_plot: plt.savefig(outpath + 'Nissan Calendar Pack Comparison.jpg', bbox_inches='tight', dpi=1000)



plt.show()


# --------------------------- normalized pack comparison ----------------------------------

fig, ax = plt.subplots(figsize=(12,6))
for i in range(4):
    if i in fewerCharList:
        ax.plot(days[i][0:-1], avg_soh[i][0] - avg_soh[i][0:-1], label=pack_names[i] + ' ' + charge_levels[i], marker='.', markersize=12) #does not graph the last characterization if there are fewwer characterizations
    else:
        ax.plot(days[i], avg_soh[i][0] -avg_soh[i], label=pack_names[i] + ' ' + charge_levels[i], marker='.', markersize=12)

ax.set_title('Nissan Calendar Aging Pack Comparison', fontsize=fs)
ax.set_ylabel('Pack Average SOH Drop [%]', fontsize=fs)
ax.set_xlabel('Time Elapsed [days]', fontsize=fs)
ax.tick_params(direction='in')
ax.legend()
if save_plot: 
    plt.savefig(outpath + 'Nissan Calendar Pack Comparison Normalized.jpg', bbox_inches='tight', dpi=1000)

plt.show()
