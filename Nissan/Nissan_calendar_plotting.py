from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from numpy.lib.function_base import append
import pandas as pd
import datetime as dt



    # ---------------------------- import data from summary spreadsheet ---------------------------- 
    summary_data = pd.read_excel(summary_file, sheet_name=pack_name)
    rated_cap = 200
    cell_num = 6
    fs = 15     # universal fontsize
    outpath = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/plots/'    # path for saving plot

    days = summary_data['days'].to_numpy()

    #  ---------------------------- get cell capacities and convert to SOH (State of Health) ---------------------------- 
    idx1 = summary_data.columns.get_loc("cell_1")
    cell_caps = summary_data.iloc[:,idx1:idx1+cell_num].to_numpy()
    cell_soh = cell_caps/rated_cap*100      # convert to %
    soh_mean = np.mean(cell_soh,axis=1)
    soh_std = np.std(cell_soh, axis=1)

    #  ---------------------------- SOH summary plot ---------------------------- 
    fig, ax = plt.subplots(figsize=(12,6))
    # plot dashed vertical lines at each test index
    ax.vlines(days,50,100,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.6)

    # plot cells
    shapes = ['o']*cell_num
    cmap = list(plt.get_cmap("tab10").colors)
    cmap.pop(5)
    for i in range(cell_num):
        cell_name = 'Cell ' + str(i+1)   
        ax.scatter(days, cell_soh[:,i], label=cell_name, marker=shapes[i], color=cmap[i], alpha=.5)
    ax.set_ylabel("State of Health [%]", fontsize=fs)
    ax.set_title(pack_name + " State of Health Summary: Day " + str(round(max(days))), fontsize=fs)
    # legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fancybox=True, ncol=6, borderaxespad=0.)
    legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
    plt.gca().add_artist(legend1)

    # plot avg, lower bound, and upper bound
    line1 = ax.plot(days, soh_mean, '-k', label="Average SOH")
    line2 = ax.plot(days, soh_mean+soh_std, label="Upper Bound SOH", color=cmap[0])
    line3 = ax.plot(days, soh_mean-soh_std, label="Lower Bound SOH", color=cmap[2])

    # calculate and set left axis lims
    low = np.amin(cell_soh)
    high = np.amax(cell_soh)
    yrange = high-low
    ax.set_ylim([low-.1*yrange, high+.1*yrange])

    # plot standard deviation
    ax.set_xlabel("Time [days]", fontsize=fs)
    ax2 = ax.twinx()
    line4 = ax2.plot(days, soh_std, '--k', label="Standard Deviation", alpha=0.8)
    ax2.set_ylabel("Standard Deviation", fontsize=fs)

    # calculate and set right axis lims
    low = np.amin(soh_std)
    high = np.amax(soh_std)
    yrange = high-low
    ax2.set_ylim([0, high+.1*yrange])

    # put all lines in one legend on plot
    lns = line1+line2+line3+line4
    labls = [l.get_label() for l in lns]
    # plt.legend(lns, labls, loc='center', bbox_to_anchor=(1.22, .5), fancybox=True, ncol = 2, borderaxespad=0.)
    plt.legend(lns, labls, loc=legend_loc)
    ax.set_xlim([-2,days[-1]+4])

    # make things look nice
    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')
    # plt.subplots_adjust(left=None, bottom=.2, right=.77, top=None, wspace=None, hspace=None)
    plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

    test_names = summary_data['Test'].to_numpy()
    # plt.tight_layout()
    if save_plot: plt.savefig(outpath + 'Calendar ' + pack_name + ' SOH Summary.jpg', bbox_inches='tight', dpi=1000)

    return days, soh_mean
summary_file = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/Calendar_Processed_Data.xlsx'

save_plot = True   # Change to True to save plot as .png

pack_names = ['NP7', 'NP9', 'NP10', 'NP12']
#pack_names = ['NP9','NP12']
charge_levels = ['100%', '75%', '90%', '50%']
#charge_levels = ['75%','50%']
legend_loc = 'lower left'


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
        avg_soh[i][j-1] = sum(summary_data.loc["Characterization "+str(j)][4:20])/16 #goes through ever characterization and sums up the soh values
        if j >= 1:
            days[i][j-1] = summary_data["days"][j-1]+days[i][j-2]   #adds the length of the gaps together to get cumulative day count for graphing
        else:
            days[i][j-1] = summary_data["days"][j-1]  

#  ---------------------------- Pack Comparison ----------------------------
fs = 15     # universal fontsize
outpath = r'C:/Users/amirs/OneDrive - UC San Diego/college/research/Dr Tong ESS/Nissan cycle testing/plots/'    # path for saving plot
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
plt.ylim([32,37])
ax.legend()
if save_plot: plt.savefig(outpath + 'Nissan Calendar Pack Comparison.jpg', bbox_inches='tight', dpi=1000)



plt.show()
