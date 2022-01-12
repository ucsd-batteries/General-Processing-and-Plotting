from matplotlib import pyplot as plt
import matplotlib as mplt
import numpy as np
import pandas as pd
from datetime import datetime
import os
import re
import csaps
from scipy.interpolate import UnivariateSpline

from scipy.signal.ltisys import TransferFunctionDiscrete

# ------- Tool Functions ------------
def hours_from_raw(data, format='%d/%m/%y %H:%M:%S'):
    datetime_string = data.date + " " + data.time # dd-mm-yy h:mm:ss
    datetime_object = [datetime.strptime(d, format) for d in datetime_string]
    tmstmp = [t.timestamp()/3600 for t in datetime_object]
    return tmstmp

def get_SOC(v1):
    soc_curve = pd.read_excel("/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/SOCcurve.xlsx")
    soc = soc_curve.iloc[2:102,2].values.astype(float)
    ocv = soc_curve.iloc[2:102,3].values.astype(float)
    soc1 = np.interp(v1, ocv, soc)
    return soc1

def moving_average(a, n):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# ------- Plotting Functions ------------
def soh_summary(save_plot=False, show_plot=True):
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    fs = 12     # universal fontsize
    summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_summary.csv'
    summary = pd.read_csv(summary_file)
    dch1 = summary.DCH1
    dch2 = summary.DCH2

    # rated_cap = max(summary.iloc[0,3:-1])
    rated_cap = 113

    rows = summary.shape[0]
    cycles = np.zeros(rows)
    cycles[0] = dch1[0]
    for i in range(rows-1):
        cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])
    cycles = cycles/rated_cap

    # get colors for cells
    def rainbow():
        all_colors = [
        #reds
            (153,0,0), (165,0,0), (187,0,0), (204,0,0), (220,0,0),(232,0,0),(255,0,0), (255,25,25), (255,51,51), (255,70,70), (255,102,102), (255,127,127), #(255,153,153), (255,170,170), (255,204,204),
        #oranges 
        (153,76,0), (165,85,0), (187,93,0), (204,102,0), (215,110,0), (237,120,0), (255,128,0), (255,153,51), (255,168,70), (255,178,102), (255,190,127), (255,204,153), 
        #yellows 
        (153,153,0), (165,165,0), (187,187,0), (204,204,0), (218,218,0), (237,237,0), (255,255,0), (255,255,51), (255,255,70), (255,255,102), (255,255,127), (255,255,153),
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #pinks 
        (153,0,153), (174,0,174), (204,0,204), (220,0,220), (233,0,233), (255,0,255), (255,10,255), (255,25,255), (255,51,255), (255,102,255), (255,125,255), (255,153,255)] #(255,175,255), (255,187,255), (255,204,255)]
        return all_colors
    def cool_colors(): 
        all_colors = [
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125)] #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        return all_colors
    def rgb_to_dec(x):
        out = [i/255 for i in x]
        return out
    # cell_colors = rainbow()
    # shapes = ['o']*8   # for rainbow
    cell_colors = cool_colors()
    shapes = ['o', '^','s','*']*2 

    fig, ax = plt.subplots(figsize=(12,8))
    # add dotted vertical lines at each data point  
    ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)

    # plot cells
    cell_soh = summary.iloc[:,13:109]/rated_cap*100
    for i in range(8*12):
        bmu = int(np.floor(i/12)+1)
        cell_label = 'Cell ' + str(bmu) + '.' + str(np.mod(i,12)+1)
        ax.scatter(cycles, cell_soh.iloc[:,i], label=cell_label, color=rgb_to_dec(cell_colors[i]), marker=shapes[bmu-1])
    ax.set_ylabel('SOH [%]', size=fs)
    ax.set_xlabel('Cycles [-]', size=fs)
    ax.set_title('Cummins State of Health Summary: Cycle ' + str(round(max(cycles))), size=fs)
    legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), fancybox=True, ncol = 8, borderaxespad=0.)
    plt.gca().add_artist(legend1)

    # set y limit for cell SOH
    low = cell_soh.min().min()
    high = cell_soh.max().max()
    yrange = high - low
    ax.set_ylim([low-yrange*.1, min(100,high+yrange*.1)])

    # plot upper and lower SOH bounds, avg SOH
    avg_SOH = np.mean(summary.iloc[:,11:-1]/rated_cap*100, axis=1)
    std_SOH = np.std(summary.iloc[:,11:-1]/rated_cap*100,  axis=1)
    coefs = np.polyfit(cycles, avg_SOH, 3)
    x = np.linspace(min(cycles), max(cycles), 1000)
    y_poly = np.polyval(coefs, x)
    line1 = ax.plot(x, y_poly, '-k', label='Average SOH')
    # line1 = ax.plot(cycles, avg_SOH, '-k', label='Average SOH', marker='o')
    # line2 = ax.plot(cycles, avg_SOH+std_SOH, label='Upperbound SOH')
    # line3 = ax.plot(cycles, avg_SOH-std_SOH, '-g', label='Lowerbound SOH')

    # plot standard deviation on right axes
    ax2 = ax.twinx()
    ax2.set_ylabel('Standard Deviation [%]', size=fs)
    line4 = ax2.plot(cycles, std_SOH, '--k', label='Standard Deviation')

    # put all lines in one legend on plot
    # lns = line1+line2+line3+line4
    lns = line1+line4
    labls = [l.get_label() for l in lns]
    plt.legend(lns, labls, loc='lower left')
    # plt.legend(lns, labls, loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, ncol = 8, borderaxespad=0.)
    # ax.set_xlim([0,cycles[-1]+4])

    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')
    plt.subplots_adjust(left=None, bottom=.5, right=None, top=None, wspace=None, hspace=None)
    if save_plot: plt.savefig(outpath+' SOH Summary.jpg', dpi=1000, bbox_inches='tight')
    if show_plot: plt.show()

def time_series(save_plot=False, show_plot=True):
    def plot_p_soc_v_c(cell_v, cell_c, hours, save_plot, show_plot, x1, x2, title, show_legend=True):
        outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
        fs = 12     # universal font size

        hours = hours - hours[x1]
        cell_p = cell_v*cell_c
        cell_soc = get_SOC(cell_v)

        # current and voltage
        fig, ax = plt.subplots(figsize=(9,4))
        line1 = ax.plot(hours[x1:x2],cell_v[x1:x2], label='Voltage', color='blue', zorder=1)
        ax.spines['left'].set_color('blue')
        ax.tick_params(axis='y', colors='blue')
        ax.yaxis.label.set_color('blue')
        ax.set_ylabel('Cell Voltage [V]', size=fs)
        ax.set_ylim([2.1, 4.1])

        ax2 = ax.twinx()
        line2 = ax2.plot(hours[x1:x2],cell_c[x1:x2], label='Current', color='red', zorder=0)
        ax2.spines['right'].set_color('red')
        ax2.tick_params(axis='y', colors='red')
        ax2.yaxis.label.set_color('red')
        ax2.set_ylabel('Cell Current [A]', size=fs)
        ax2.set_ylim([-40, 120])

        # add dot at v1 and v2 & i1 and i2
        v1 = (5.39305, 3.97)
        v2 = (8.41655, 3.24)
        ax.scatter(v1[0], v1[1], color='k', zorder=10)
        ax.scatter(v2[0], v2[1], color='k', zorder=10)
        ax.annotate("V_1", xy=v1, xytext=(6,3.97), arrowprops=dict(arrowstyle="->"))
        ax.annotate("V_2", xy=v2, xytext=(9,3.24), arrowprops=dict(arrowstyle="->"))

        i1 = (5.39457, 35.040)
        i2 = (7.74329, 33.99)
        ax2.scatter(i1[0], i1[1], color='k', zorder=10)
        ax2.scatter(i2[0], i2[1], color='k', zorder=10)
        ax2.annotate("I_1", xy=i1, xytext=(4.5,15), arrowprops=dict(arrowstyle="->"))
        ax2.annotate("I_2", xy=i2, xytext=(8,15), arrowprops=dict(arrowstyle="->"))


        ax.set_title('Characterization: Cell-Level Voltage')

        ax.set_xlim(0,max(hours[x1:x2]))
        ax.set_xlabel('Time [hours]', size=fs)
        # ax.set_title('Cummins ' + title + ': Cell Voltage and Current')
        lns = line1+line2
        # lns = line1
        labls = [l.get_label() for l in lns]
        if show_legend: plt.legend(lns, labls, loc='upper left')
        if save_plot: plt.savefig(outpath+ title + ' C and V demo', dpi=1000)

        # # power and SOC
        # fig, ax = plt.subplots(figsize=(9,4))
        # line1 = ax.plot(hours[x1:x2],cell_soc[x1:x2], label='SOC', color='blue')
        # ax.spines['left'].set_color('blue')
        # ax.tick_params(axis='y', colors='blue')
        # ax.yaxis.label.set_color('blue')
        # ax.set_ylabel('Cell SOC [%]', size=fs)

        # ax2 = ax.twinx()
        # line2 = ax2.plot(hours[x1:x2],cell_p[x1:x2], label='Power', color='red')
        # ax2.spines['right'].set_color('red')
        # ax2.tick_params(axis='y', colors='red')
        # ax2.yaxis.label.set_color('red')
        # ax2.set_ylabel('Cell Power [W]', size=fs)

        # ax.set_xlim(0,max(hours[x1:x2]))
        # ax.set_xlabel('Time [hours]', size=fs)
        # ax.set_title('Cummins ' + title + ': Cell SOC and Power')
        # lns = line1+line2
        # labls = [l.get_label() for l in lns]
        # if show_legend: plt.legend(lns, labls, loc='upper left')
        # if save_plot: plt.savefig(outpath+title+' P and SOC', dpi=1000)
    
    # plot characterization
    char_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/Char_1/combined_csv.csv'
    char = pd.read_csv(char_file)
    cell_v = char.BMU01_Cell_10_Voltage.to_numpy()
    cell_c = char.modbus_Current.to_numpy()
    hours = np.array(hours_from_raw(char, format='%d-%m-%Y %H:%M:%S'))
    # plot_p_soc_v_c(cell_v, cell_c, hours, save_plot, show_plot, x1=4860, x2=10250, title='Characterization')
    plot_p_soc_v_c(cell_v, cell_c, hours, save_plot, show_plot, x1=4860, x2=12250, title='Characterization')


    # # plot aging
    # drive_cycle_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/Aging_1/combined_csv.csv'
    # drive_cycle = pd.read_csv(drive_cycle_file)
    # cell_v = drive_cycle.BMU01_Cell_10_Voltage.to_numpy()
    # cell_c = drive_cycle.modbus_Current.to_numpy()
    # hours = np.array(hours_from_raw(drive_cycle, format='%d/%m/%y %H:%M:%S'))
    # plot_p_soc_v_c(cell_v, cell_c, hours, save_plot, show_plot, x1=6211, x2=7900, title='Aging', show_legend=False)

    if show_plot: plt.show()

def monthly_meeting(save_plot=False, show_plot=True):
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    fs = 12     # universal fontsize
    summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_summary.csv'
    summary = pd.read_csv(summary_file)

    dch1 = summary.DCH1
    dch2 = summary.DCH2
    Wh1 = summary.WhDch1
    Wh2 = summary.WhDch2

    rated_cap = 113

    shapes = ['o', '^','s','*']*2 
    colors = []
    fills = ['green', 'green', 'green', 'green', 'none', 'none', 'none', 'none']
    all_colors = [(153,0,0), (153,76,0),(153,153,0),(0,153,0),(0,153,153),(0,0,153),(76,0,153),(153,0,153)]
    def rgb_to_dec(x):
        out = [i/255 for i in x]
        return out
        # #reds
        #     (153,0,0), (165,0,0), (187,0,0), (204,0,0), (220,0,0),(232,0,0),(255,0,0), (255,25,25), (255,51,51), (255,70,70), (255,102,102), (255,127,127), #(255,153,153), (255,170,170), (255,204,204),
        # #oranges 
        # (153,76,0), (165,85,0), (187,93,0), (204,102,0), (215,110,0), (237,120,0), (255,128,0), (255,153,51), (255,168,70), (255,178,102), (255,190,127), (255,204,153), 
        # #yellows 
        # (153,153,0), (165,165,0), (187,187,0), (204,204,0), (218,218,0), (237,237,0), (255,255,0), (255,255,51), (255,255,70), (255,255,102), (255,255,127), (255,255,153),
        # #greens 
        # (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        # #light_blues 
        # (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        # #blues 
        # (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        # #purples 
        # (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        # #pinks 
        # (153,0,153), (174,0,174), (204,0,204), (220,0,220), (233,0,233), (255,0,255), (255,10,255), (255,25,255), (255,51,255), (255,102,255), (255,125,255), (255,153,255)] #(255,175,255), (255,187,255), (255,204,255)]
        # return all_colors

    rows = summary.shape[0]
    cycles = np.zeros(rows)
    energy = np.zeros(rows)
    cycles[0] = dch1[0]/rated_cap
    energy[0] = Wh1[0]
    for i in range(rows-1):
        cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])/rated_cap
        energy[i+1] = energy[i] + Wh2[i] + Wh1[i+1]
    
    # plot cap vs cycle count 
    fig, ax = plt.subplots(figsize=(8,6))
    ax2=ax.twinx()
    module_caps = summary.iloc[:,5:13].to_numpy()
    cell_caps = summary.iloc[:,13:109].to_numpy()
    # ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)
    lns = []
    for i in range(8):
        module = "BMU0" + str(i+1)
        # ax.scatter(cycles, module_caps[:,i]/rated_cap*100, color=rgb_to_dec(all_colors[i]), label=module, marker='.', linewidth=3)
        coefs = np.polyfit(cycles, module_caps[:,i]/rated_cap*100, 4)
        y_poly = np.polyval(coefs, cycles)
        ax.plot(cycles, y_poly,color=rgb_to_dec(all_colors[i]), label=module + ' SOH')
        cell_std = np.std(cell_caps[:,i*8:i*8+12]/rated_cap*100,axis=1)
        coefs = np.polyfit(cycles, cell_std, 3)
        y_poly = np.polyval(coefs, cycles)
        lns += ax2.plot(cycles, y_poly,'--',linewidth=1.5, color=rgb_to_dec(all_colors[i]), label=module+ ' Std')

    ax.set_xlabel('Cycle Count [-]', fontsize=fs)
    ax.set_ylabel('State of Health [%]', fontsize=fs)
    ax.set_title('Module SOH vs Cycle Count', fontsize=fs)
    ax.set_ylim([85, 98])
    ax.tick_params(direction='in')
    ax.legend(loc='center', bbox_to_anchor=(1.3, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    ax2.set_ylabel('Standard Deviation [%]')
    ax2.set_ylim([0,3])
    plt.subplots_adjust(left=None, bottom=None, right=.83, top=None, wspace=None, hspace=None)
    plt.tight_layout()
    if save_plot: plt.savefig(outpath+'Monthly Cap vs Cycles no filter.jpg', dpi=1000, tight_layout=True)

    # # plot SOH vs energy discharge
    # fig, ax = plt.subplots(figsize=(8,6))
    # ax.vlines(energy, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)
    # for i in range(8):
    #     module = "BMU0" + str(i+1)
    #     ax.scatter(energy, module_caps[:,i]/rated_cap*100, label=module, marker='.', linewidth=3)
    #     coefs = np.polyfit(energy, module_caps[:,i], 4)
    #     y_poly = np.polyval(coefs, energy)
    #     plt.plot(energy, y_poly)
    # ax.set_xlabel('Energy Discharge [kWh]', fontsize=fs)
    # ax.set_ylabel('Module State of Health [%]', fontsize=fs)
    # ax.set_title('Module Capacity vs Energy Discharge', fontsize=fs)
    # ax.set_ylim([85, 100])
    # ax.tick_params(direction='in')
    # ax.legend(loc='center', bbox_to_anchor=(1.11, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    # plt.subplots_adjust(left=None, bottom=None, right=.83, top=None, wspace=None, hspace=None)
    # plt.tight_layout()
    # if save_plot: plt.savefig(outpath+'Monthly SOH vs kWh no filter.jpg', dpi=1000, tight_layout=True)

    # # Cell-Level Bar Plots
    # fig, axs = plt.subplots(2,4,figsize=(14,6))
    # cells = [str(i+1) for i in range(12)]
    # for i in range(8):
    #     row = int(np.floor(i/4))
    #     col = i%4
    #     cell_vals = cell_caps[-1, i*12:i*12+12]
    #     axs[row, col].bar(cells, cell_vals)
    #     axs[row, col].set_ylim([80, 115])
    #     axs[row, col].set_title('BMU0'+str(i+1))
    #     if row==1: axs[row, col].set_xlabel('Cell Number [-]')
    #     if col==0: axs[row, col].set_ylabel('Capacity [Ah]')
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.3)
    # plt.tight_layout()
    # if save_plot: plt.savefig(outpath+'Cell Caps.jpg', dpi=1000, tight_layout=True)
    
    # individual module plots
    fig, ax = plt.subplots(figsize=(12,6))
    mods = [5,8]
    colors = [(153,76,0),(0,153,153)]
    m = 7
    cell_vals = cell_caps[:, m*12:m*12+12]
    for i in range(12):
        ax.scatter(cycles, cell_vals[:,i]/rated_cap*100, label='BMU'+str(m+1)+'.Cell'+str(i+1), marker='.', linewidth=3)
        ax.plot(cycles, cell_vals[:,i]/rated_cap*100)
    plt.legend(loc='best')
    ax.set_xlabel('Cycles [-]')
    ax.set_ylabel('State of Health [%]')

    # system level soh vs cycles 
    fig, ax = plt.subplots(figsize=(12,6))
    module_SOH = module_caps/rated_cap*100
    cell_SOH = cell_caps/rated_cap*100
    ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)
    cell_std = np.std(cell_SOH, axis=1)
    mod_std = np.std(module_SOH, axis=1)
    for i in range(8):
        module = "BMU0" + str(i+1)
        ax.scatter(cycles, module_SOH[:,i], label=module, marker=shapes[i], facecolor=fills[i], edgecolor='k', linewidth=1)
    # legend1 = ax.legend(loc='center', bbox_to_anchor=(1.13, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    legend1 = ax.legend(loc='lower left', bbox_to_anchor=(.01, .02), fancybox=True, ncol = 1, borderaxespad=0.)
    plt.gca().add_artist(legend1)
    avg_SOH = np.mean(module_SOH, axis=1)
    # plot line of fit
    coefs = np.polyfit(cycles, avg_SOH, 4)
    x = np.linspace(min(cycles), max(cycles), 1000)
    y_poly = np.polyval(coefs, x)
    # moving average
    # y = moving_average(avg_SOH, 10)
    # x = cycles
    # smoothing spling
    # spl = UnivariateSpline(cycles, avg_SOH, k=4)
    # # sp = csaps.UnivariateCubicSmoothingSpline(cycles, avg_SOH, smooth=0.85)
    # y_poly = spl(x)
    line1 = plt.plot(x, y_poly, label="System SOH", color='black')
    # plot standard deviation
    ax2 = ax.twinx()
    ax2.spines['right'].set_color('blue')
    ax2.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_color('blue')
    line2 = ax2.plot(cycles, cell_std, linestyle='dotted', color='blue', label="Cell Std")
    line3 = ax2.plot(cycles, mod_std, linestyle='dashed', color='blue', label="Module Std")
    ax2.set_ylabel("Standard Deviation", fontsize=fs)
    ax2.set_ylim([1, 2.5])
    lns = line2+line3+line1
    labls = [l.get_label() for l in lns]
    # ax.legend(lns, labls, loc='lower left')
    # ax.legend(lns, labls, loc='lower left', bbox_to_anchor=(1.071, .15), fancybox=True, ncol = 1, borderaxespad=0.)
    ax.legend(lns, labls, loc='lower left', bbox_to_anchor=(.14, .02), fancybox=True, ncol = 1, borderaxespad=0.)
    # other plot formatting: 
    ax.set_xlabel('Cycle Count [-]', fontsize=fs)
    ax.set_ylabel('System SOH [%]', fontsize=fs)
    ax.set_title('System State of Health vs Cycle Count', fontsize=fs)
    ax.set_ylim([85, 98])
    ax.tick_params(direction='in')
    # legend
    # ax.legend(lns, labls, loc='center', bbox_to_anchor=(1.11, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    plt.subplots_adjust(left=None, bottom=None, right=.83, top=None, wspace=None, hspace=None)
    if save_plot: plt.savefig(outpath+'System SOH.jpg', dpi=1000, tight_layout=True)

    # # plot temp and cooling
    # fig, ax = plt.subplots(figsize=(10,5))
    # drive_cycle_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/combined_datasets/DutyCycle_4_5_6.csv'
    # drive_cycle = pd.read_csv(drive_cycle_file)
    # all_temps = np.zeros((drive_cycle.shape[0],8*7))
    # for i in range(8*7):
    #     bmu = int(np.floor(i/7)+1)
    #     col_name = 'BMU0' + str(bmu) + '_CMA_Temp' + str(np.mod(i,7)+1)
    #     all_temps[:,i] = drive_cycle[col_name]
    # avg_temp = np.mean(all_temps, axis=1)
    # fan = drive_cycle.isPumpFanOff.to_numpy()
    # fan = np.invert(fan)
    # hours = np.array(hours_from_raw(drive_cycle, format='%d-%m-%Y %H:%M:%S'))
    # # x1=7000
    # x1=0
    # x2=len(fan)
    # hours = hours - hours[x1]
    # line1 = ax.plot(hours[x1:x2],avg_temp[x1:x2], label='System Temperature', color='blue')
    # line2 = ax.plot(hours[x1:x2], 30*np.ones(fan[x1:x2].size), 'k--', label='Fan Temp Trigger')
    # line3 = ax.plot(hours[x1:x2], 50*np.ones(fan[x1:x2].size), 'k:', markersize=2,label='Danger Temp Trigger')
    # # ax.spines['left'].set_color('blue')
    # # ax.tick_params(axis='y', colors='blue')
    # # ax.yaxis.label.set_color('blue')
    # ax.set_ylabel('Temperature [°C]', size=fs)
    # ax.tick_params(direction='in')
    # ax.set_ylim([27, 52])
    # ax2 = ax.twinx()
    # line4 = ax2.plot(hours[x1:x2],fan[x1:x2], label='Fan Status', color='red')
    # ax2.spines['right'].set_color('red')
    # ax2.tick_params(axis='y', colors='red')
    # ax2.yaxis.label.set_color('red')
    # ax2.set_ylabel('Fan Status [-]', size=fs)
    # ax2.set_ylim([-.5, 8.5])
    # ax2.set_yticks([0,1])
    # ax2.set_yticklabels(['off', 'on'])
    # ax2.tick_params(direction='in')
    # ax.set_xlim(0,max(hours[x1:x2]))
    # ax.set_xlabel('Time [hours]', size=fs)
    # ax.set_title('Average Module Temperature and Cooling')
    # lns = line3+line1+line2+line4
    # labls = [l.get_label() for l in lns]
    # # plt.legend(lns, labls, loc='upper left')
    # ax.legend(lns, labls, loc='center', bbox_to_anchor=(1.24, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    # plt.subplots_adjust(left=None, bottom=None, right=.83, top=None, wspace=None, hspace=None)
    # plt.tight_layout()
    # if save_plot: plt.savefig(outpath+ 'temp and cooling.jpg', dpi=1000)

    # # plot characterization
    # char_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/combined_datasets/Char9.csv'
    # char = pd.read_csv(char_file)
    # mod_v = char.modbus_Voltage.to_numpy()
    # mod_c = char.modbus_Current.to_numpy()
    # hours = np.array(hours_from_raw(char, format='%d-%m-%Y %H:%M:%S'))
    # x1=4860
    # x2=10250
    # x1=0
    # x2=int(len(hours)-.4*len(hours))
    # hours = hours - hours[x1]
    # fig, ax = plt.subplots(figsize=(10,5))
    # # plot voltage
    # line1 = ax.plot(hours[x1:x2],mod_v[x1:x2], label='Voltage', color='blue')
    # ax.spines['left'].set_color('blue')
    # ax.tick_params(axis='y', colors='blue')
    # ax.yaxis.label.set_color('blue')
    # ax.set_ylabel('Voltage [V]', size=fs)
    # ax.tick_params(direction='in')
    # # plot current
    # ax2 = ax.twinx()
    # line2 = ax2.plot(hours[x1:x2],mod_c[x1:x2], label='Current', color='red')
    # ax2.spines['right'].set_color('red')
    # ax2.tick_params(axis='y', colors='red')
    # ax2.yaxis.label.set_color('red')
    # ax2.set_ylabel('Current [A]', size=fs)
    # ax2.tick_params(direction='in')
    # ax.set_xlim(0,max(hours[x1:x2]))
    # ax.set_xlim(0,15)
    # ax.set_xlabel('Time [hours]', size=fs)
    # ax.set_title('Cummins System Characterization: Voltage and Current')
    # lns = line1+line2
    # labls = [l.get_label() for l in lns]
    # plt.legend(lns, labls, loc='upper left')

    if save_plot: plt.savefig(outpath+ 'System Characterization C and V.jpg', dpi=1000)

    if show_plot: plt.show()

def all_cell_voltage(save_plot=False):
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/Aging15/combined_csv.csv'
    data = pd.read_csv(file)
    tmstmp = hours_from_raw(data, format='%d-%m-%Y %H:%M:%S')
    hours = [h - tmstmp[0] for h in tmstmp]
    # hours = [h - 18.5 for h in hours]
    # limits: 
    cell_lim = [3.1, 3.9]
    mod_lim = [37.25, 46.87]
    sys_lim = [297, 375]

    # cell fig
    fig_cell, ax_cell = plt.subplots()
    # # module fig
    # fig_mod, ax_mod = plt.subplots()

    # mod_vs = np.zeros((data.shape[0],8))
    # cell_vs = np.zeros((data.shape[0],8*12))
    # for m in range(8):      
    #     module = 'BMU0' + str(m+1)
    #     v_mod = data[module + '_CMA_Voltage'].to_numpy()
    #     mod_vs[:,m] = v_mod
    #     ax_mod.plot(hours, v_mod, label=module)
    # #     for c in range(12):
    # #         cell = module + '_Cell_' + str(c+1) + '_Voltage'
    # #         v_cell = data[cell].to_numpy()
    # #         cell_vs[:,m*12 + c] = v_cell
    # #         ax_cell.plot(hours, v_cell, label=cell)
    # # # ax_cell.set_xlim([0, 7])
    # # ax_cell.hlines(cell_lim, 0, 24, colors='k')
    # # ax_cell.set_title('Cell Voltages')
    # # ax_cell.set_xlabel('Time [hours]')
    # # ax_cell.set_ylabel('Cell Voltage [V]')
    # ax_mod.set_xlabel('Time [hours]')
    # ax_mod.legend()
    # ax_mod.set_ylabel('Module Voltage [V]')
    # ax_mod.set_xlabel('Time [hours]')
    # ax_mod.set_title('Module Voltages')
    # # ax_mod.set_xlim([0, 7])
    # ax_mod.hlines(mod_lim, 0, 24, colors='k')

    # fig_sys, ax_sys = plt.subplots()
    # ax_sys.plot(hours, data['modbus_Voltage'])
    # ax_sys.set_ylabel('Modbus Voltage [V]')
    # ax_sys.set_xlabel('Time [hours]')
    # ax_sys.set_title('Modbus Voltage')
    # # ax_sys.set_xlim([0, 7])
    # ax_sys.hlines(sys_lim, 0, 24, colors='k')
    
    m=8
    module = 'BMU0' + str(m)
    for c in range(12):
        cell = module + '_Cell_' + str(c+1) + '_Voltage'
        v_cell = data[cell].to_numpy()
        ax_cell.plot(hours, v_cell, label=cell)
    ax_cell.hlines(cell_lim, 0, 24, colors='k')
    ax_cell.set_ylabel('Cell Voltage [V]')
    ax_cell.set_xlabel('Test Interval [-]')
    ax_cell.set_title(module + ' Cells')
    # ax_cell.set_xlim([0, 7])
    ax_cell.legend(loc='center', bbox_to_anchor=(1.3, .5), fancybox=True, ncol = 1, borderaxespad=0.)
    plt.subplots_adjust(left=None, bottom=None, right=.6, top=None, wspace=None, hspace=None)

    # # max temp
    # fig_temp, ax_temp = plt.subplots()
    # dates = mplt.dates.date2num(data['Datetime'])
    # ax_temp.hlines(30,min(dates), max(dates), colors='k')
    # for m in range(8):      
    #     module = 'BMU0' + str(m+1)
    #     t_mod = data[module + '_CMA_Max_Temp'].to_numpy()
    #     ax_temp.plot_date(dates, t_mod, '-')
    #     # ax_temp.plot(tmstmp, t_mod, label=module)
    # ax_temp.legend(loc='best')
    # ax_temp.set_xlabel('Datetime [mm-dd hh]')
    # ax_temp.set_ylabel('Temp [°C]')
    # ax_temp.set_title('Module Max Temp')

    plt.show()
    # if save_plot: 
    #     ax_mod.savefig(outpath+ 'Module troubleshooting.jpg', dpi=1000)
    #     ax_cell.savefig(outpath + cell + 'troubleshooting.jpg', dpi=1000)

def trendline():
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    fs = 12     # universal fontsize
    summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_summary.csv'
    summary = pd.read_csv(summary_file)
    dch1 = summary.DCH1
    dch2 = summary.DCH2

    # rated_cap = max(summary.iloc[0,3:-1])
    rated_cap = 113

    rows = summary.shape[0]
    cycles = np.zeros(rows)
    cycles[0] = dch1[0]
    for i in range(rows-1):
        cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])
    cycles = cycles/rated_cap

    # get colors for cells
    def rainbow():
        all_colors = [
        #reds
            (153,0,0), (165,0,0), (187,0,0), (204,0,0), (220,0,0),(232,0,0),(255,0,0), (255,25,25), (255,51,51), (255,70,70), (255,102,102), (255,127,127), #(255,153,153), (255,170,170), (255,204,204),
        #oranges 
        (153,76,0), (165,85,0), (187,93,0), (204,102,0), (215,110,0), (237,120,0), (255,128,0), (255,153,51), (255,168,70), (255,178,102), (255,190,127), (255,204,153), 
        #yellows 
        (153,153,0), (165,165,0), (187,187,0), (204,204,0), (218,218,0), (237,237,0), (255,255,0), (255,255,51), (255,255,70), (255,255,102), (255,255,127), (255,255,153),
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #pinks 
        (153,0,153), (174,0,174), (204,0,204), (220,0,220), (233,0,233), (255,0,255), (255,10,255), (255,25,255), (255,51,255), (255,102,255), (255,125,255), (255,153,255)] #(255,175,255), (255,187,255), (255,204,255)]
        return all_colors
    def cool_colors(): 
        all_colors = [
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125)] #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        return all_colors
    def rgb_to_dec(x):
        out = [i/255 for i in x]
        return out
    # cell_colors = rainbow()
    # shapes = ['o']*8   # for rainbow
    cell_colors = cool_colors()
    shapes = ['o', '^','s','*']*2 

    fig, ax = plt.subplots(figsize=(12,8))
    # add dotted vertical lines at each data point  
    ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)

    # get upper and lower SOH bounds, avg SOH
    avg_SOH = np.mean(summary.iloc[:,11:-1]/rated_cap*100, axis=1)
    std_SOH = np.std(summary.iloc[:,11:-1]/rated_cap*100,  axis=1)
    upper = avg_SOH+std_SOH
    lower = avg_SOH+std_SOH

    # # plot cells
    # cell_soh = summary.iloc[:,13:109]/rated_cap*100
    # for i in range(8*12):
    #     bmu = int(np.floor(i/12)+1)
    #     cell_label = 'Cell ' + str(bmu) + '.' + str(np.mod(i,12)+1)
    #     ax.scatter(cycles, cell_soh.iloc[:,i], label=cell_label, color=rgb_to_dec(cell_colors[i]), marker=shapes[bmu-1])
    
    # plot module averages
    cell_soh = summary.iloc[:,13:109]/rated_cap*100
    for i in range(8*12):
        bmu = int(np.floor(i/12)+1)
        cell_label = 'Cell ' + str(bmu) + '.' + str(np.mod(i,12)+1)
        ax.scatter(cycles, cell_soh.iloc[:,i], label=cell_label, color=rgb_to_dec(cell_colors[i]), marker=shapes[bmu-1])

    ax.set_ylabel('SOH [%]', size=fs)
    ax.set_xlabel('Cycles [-]', size=fs)
    ax.set_title('Cummins State of Health Summary: Cycle ' + str(round(max(cycles))), size=fs)
    legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), fancybox=True, ncol = 8, borderaxespad=0.)
    plt.gca().add_artist(legend1)

    # set y limit for cell SOH
    low = cell_soh.min().min()
    high = cell_soh.max().max()
    yrange = high - low
    ax.set_ylim([low-yrange*.1, min(100,high+yrange*.1)])

    # line1 = ax.plot(cycles, avg_SOH, '-k', label='Average SOH', marker='o')
    # line2 = ax.plot(cycles, avg_SOH+std_SOH, label='Upperbound SOH')
    # line3 = ax.plot(cycles, avg_SOH-std_SOH, '-g', label='Lowerbound SOH')
    
    # # deg = 2
    # m1, m2, b = np.polyfit(cycles, avg_SOH, 2)
    # line1 = plt.plot(cycles, m1*cycles**2 + m2*cycles + b, label="Linear Regression")

    # deg = 1
    m, b = np.polyfit(cycles, avg_SOH, 1)
    line1 = plt.plot(cycles, m*cycles + b, label="SOH Regression")
    
    # plot standard deviation on right axes
    ax2 = ax.twinx()
    ax2.set_ylabel('Cell Standard Deviation [%]', size=fs)
    # line4 = ax2.plot(cycles, std_SOH, '--k', label='Standard Deviation', linewidth=1)
    m, b = np.polyfit(cycles, std_SOH, 1)
    line4 = plt.plot(cycles, m*cycles + b, '--k', label='Std Regression', linewidth=1)

    # put all lines in one legend on plot
    lns = line1+line4
    labls = [l.get_label() for l in lns]
    plt.legend(lns, labls, loc='upper right')
    # plt.legend(lns, labls, loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, ncol = 8, borderaxespad=0.)
    # ax.set_xlim([0,cycles[-1]+4])

    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')
    plt.subplots_adjust(left=None, bottom=.5, right=None, top=None, wspace=None, hspace=None)
    plt.show()
    # if save_plot: plt.savefig(outpath+' SOH Summary.jpg', dpi=1000, bbox_inches='tight')

def inLab_data(save_plot=False):
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    fs = 12     # universal fontsize

    summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/inLab_cell_summary.xlsx'
    BMUs = ['BMU1_BMU3', 'BMU2_BMU4', 'BMU5_BMU6']
    strength = ['Agressive', 'Normal', 'Mild']

    for bmu in range(3):
        mods = BMUs[bmu].split('_')
        summary = pd.read_excel(summary_file, sheet_name=BMUs[bmu])
        dch1 = summary['DCH1'].to_numpy()
        dch2 = summary['DCH2'].to_numpy()
        days = summary['added'].to_numpy()
        days_e = summary['eliminated'].to_numpy()
        # datetimes = [datetime.strptime(d, format) for d in summary['timestamps']]
        tmstmp = [t.timestamp()/3600/24 for t in summary['timestamps']]
        days = [d - tmstmp[0] for d in tmstmp]
        rated_cap = 108.820559325236    #max cell DchCapAh for first test

        rows = summary.shape[0]
        cycles = np.zeros(rows)
        days_cum = np.zeros((rows,2))
        cycles[0] = dch1[0]/rated_cap
        days_cum[0] = [days[0], days_e[0]]
        for i in range(rows-1):
            cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])/rated_cap
            days_cum[i+1] = days_cum[i] + [days[i+1],days_e[i+1]]

        bats = [summary['bat1'].to_numpy()/rated_cap*100, summary['bat2'].to_numpy()/rated_cap*100]

        # plot SOH vs cycle count 
        fig, ax = plt.subplots(figsize=(5,6))
        # ax.vlines(days_cum[:,0], 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)
        for i in range(2):
            module = mods[i]
            d = days_cum[:,0]
            d = np.array(days)
            ax.scatter(d, bats[i], label=module, marker='.', linewidth=3)
            # m, b = np.polyfit(d, bats[i], 1)
            coefs = np.polyfit(d, bats[i], 2)
            # plt.plot(d, m*d + b)
            x = np.linspace(0,70,1000)
            y_poly = np.polyval(coefs, x)
            plt.plot(x, y_poly)
            # plt.plot(x, m*x + b)
            # fit_string = 'SOH = ' + "{:.2f}".format(m*100) + 'e-2 x days + ' + "{:.2f}".format(b)
            print()
            print(y_poly[0]-y_poly[-1])
            print(y_poly[-1])
        ax.set_xlabel('Time [days]', fontsize=fs)
        ax.set_ylabel('State of Health [%]', fontsize=fs)
        ax.set_title(strength[bmu] + ' Duty Cycle', fontsize=fs)
        ax.set_ylim([np.amax(bats)-5, np.amax(bats)+.5])
        # ax.set_ylim([94, 100.5])
        # ax.set_xlim([-5, 290])
        # ax.set_xlim([0, 50])
        ax.tick_params(direction='in')
        # ax.legend(loc='center', bbox_to_anchor=(1.11, .5), fancybox=True, ncol = 1, borderaxespad=0.)
        ax.legend(loc='upper right')
        plt.subplots_adjust(left=None, bottom=None, right=.83, top=None, wspace=None, hspace=None)
        plt.tight_layout()
        if save_plot: plt.savefig(outpath+'In Lab ' + strength[bmu] + '.jpg', dpi=1000, tight_layout=True)
        plt.show()

def temps():
    # path to combined datasets
    combined_datasets = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Data/Cummins/combined_datasets'

    # get all files in combined datasests folder
    logs = [l[2] for l in os.walk(combined_datasets)]
    logs = logs[0]

    # # remove '.DS_Store' file from logs
    # if '.DS_Store' in logs[0]: logs[0].remove('.DS_Store')

    # put all Aging and Char tests together
    Agings = [l for l in logs if l[0]=='A']
    Chars = [l for l in logs if l[0]=='C']

    # sort Aging and Char arrays
    # sort step 1: use regular expressions to get number value in each file name
    Aging_nums = np.zeros(len(Agings))
    for i,f in enumerate(Agings):
        x = re.search("\d+", f)
        Aging_nums[i] = x.group(0)

    Char_nums = np.zeros(len(Chars))
    for i,f in enumerate(Chars):
        x = re.search("\d+", f)
        Char_nums[i] = x.group(0)

    # sort step 2: sort the Agings and Chars arrays based on their respective num arrays
    Agings = [x for _, x in sorted(zip(Aging_nums, Agings))]
    Chars = [x for _, x in sorted(zip(Char_nums, Chars))]

    # get average temp for each Aging + Char group over all modules
    avg_temp = np.zeros(len(Chars)) # initiate output avg temps array
    num_temp_sensors = 8*7      # 7 temp sensors for each of the 8 modules

    # since Char1 had no preceeding Aging test, the first avg temp is based only on Char1.csv
    char1 = pd.read_csv(combined_datasets+'/Char1.csv')
    char1_temps = np.zeros((char1.shape[0],num_temp_sensors))  # initiate array for all temps for this test
    for i in range(num_temp_sensors):
        bmu = int(np.floor(i/7)+1)
        col_name = 'BMU0' + str(bmu) + '_CMA_Temp' + str(np.mod(i,7)+1)
        char1_temps[:,i] = char1[col_name]
    avg_temp[0] = np.mean(np.mean(char1_temps, axis=1))

    # loop through other files 
    for i in range(len(Agings)):
        # get aging and char datasets: Aging1 goes with Char2
        aging_set = pd.read_csv(combined_datasets+'/'+Agings[i])
        char_set = pd.read_csv(combined_datasets+'/'+Chars[i+1])

        # initiate arrays for temps from both datasets
        aging_temps = np.zeros((aging_set.shape[0],num_temp_sensors))  
        char_temps = np.zeros((char_set.shape[0],num_temp_sensors))  

        # get all temp values for each dataset
        for j in range(num_temp_sensors):
            bmu = int(np.floor(j/7)+1)
            col_name = 'BMU0' + str(bmu) + '_CMA_Temp' + str(np.mod(j,7)+1)
            aging_temps[:,j] = aging_set[col_name]
            char_temps[:,j] = char_set[col_name]
        
        # combine temp arrays into one
        all_temps = np.mean(aging_temps, axis=1).tolist() + np.mean(char_temps, axis=1).tolist()

        # add to output temp array
        avg_temp[i+1] = sum(all_temps)/len(all_temps)

    print()
    # summary = pd.read_csv(summary_file)

def cell_level(save_plot=False, show_plot=True):
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Cummins/'
    fs = 12     # universal fontsize
    outdoor_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/testing_summary.csv'
    inLab_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Cummins/in-lab_cell_summary.xlsx'
    outdoor = pd.read_csv(outdoor_file)

    # get colors for cells
    def rainbow():
        all_colors = [
        #reds
            (153,0,0), (165,0,0), (187,0,0), (204,0,0), (220,0,0),(232,0,0),(255,0,0), (255,25,25), (255,51,51), (255,70,70), (255,102,102), (255,127,127), #(255,153,153), (255,170,170), (255,204,204),
        #oranges 
        (153,76,0), (165,85,0), (187,93,0), (204,102,0), (215,110,0), (237,120,0), (255,128,0), (255,153,51), (255,168,70), (255,178,102), (255,190,127), (255,204,153), 
        #yellows 
        (153,153,0), (165,165,0), (187,187,0), (204,204,0), (218,218,0), (237,237,0), (255,255,0), (255,255,51), (255,255,70), (255,255,102), (255,255,127), (255,255,153),
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #pinks 
        (153,0,153), (174,0,174), (204,0,204), (220,0,220), (233,0,233), (255,0,255), (255,10,255), (255,25,255), (255,51,255), (255,102,255), (255,125,255), (255,153,255)] #(255,175,255), (255,187,255), (255,204,255)]
        return all_colors
    def cool_colors(): 
        all_colors = [
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204), (110,0,232), (127,0,255), (130,13,255), (138,25,255), (153,51,255), (165,80,255), (178,102,255), (187,125,255), (204,153,255), #(211,175,255), (215,187,255), (229,204,255), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), (0,186,186), (0,204,204), (0,218,218),(0,237,237),(0,255,255), (25,255,255), (51,255,255), (102,255,255), (125,255,255), #(153,255,255), (175,255,255), (187,255,255), (204,255,255), 
        #greens 
        (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125)] #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        return all_colors
    def inLab_cool_colors(): 
        all_colors = [
        #greens 
        (0,153,0), (0,164,0), (0,175,0), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #light_blues 
        (0,153,153), (0,164,164),(0,175,175), 
        #blues 
        (0,0,153), (0,0,174), (0,0,204), 
        #purples 
        (76,0,153), (85,0,174), (102,0,204)]
        return all_colors

    def contrast_colors():
        all_colors = [
        # #oranges 
        # (153,76,0), (165,85,0), (187,93,0), (204,102,0), (215,110,0), (237,120,0), (255,128,0), (255,153,51), (255,168,70), (255,178,102), (255,190,127), (255,204,153), 
        # #greens 
        # (0,153,0), (0,164,0), (0,175,0), (0,186,0), (0,198,0), (0,215,0), (0,230,0), (0,255,0), (25,255,25), (51,255,51), (102,255,102), (125,255,125), #(153,255,153), (175,255,175), (187,255,187), (204,255,204), 
        #reds
            (153,0,0), (165,0,0), (187,0,0), (204,0,0), (220,0,0),(232,0,0),(255,0,0), (255,25,25), (255,51,51), (255,70,70), (255,102,102), (255,127,127), #(255,153,153), (255,170,170), (255,204,204),
        #blues 
        (0,0,153), (0,0,174), (0,0,204), (0,0,255), (25,25,255), (51,51,255), (102,102,255), (110,110,255), (123,123,255), (135,135,255), (153,153,255), (175,175,255), #(187,187,255), (204,204,255), 
        ]
        return all_colors
    def rgb_to_dec(x):
        out = [i/255 for i in x]
        return out
    # cell_colors = rainbow()
    # shapes = ['o']*8   # for rainbow
    # cell_colors = cool_colors()
    cell_colors = contrast_colors()
    inLab_colors = inLab_cool_colors()
    shapes = ['o', '^','s','*']*2 

    # # ------- in lab processing -------
    # BMUs = ['BMU1_BMU3', 'BMU2_BMU4', 'BMU5_BMU6']
    # bNums = [[1,3], [2,4], [5,6]]
    # strength = ['Agressive', 'Normal', 'Mild']
    # titles = ['(A)', '(B)', '(C)']
    # for bmu in range(3):
    #     mods = BMUs[bmu].split('_')
    #     mods = bNums[bmu]
    #     summary = pd.read_excel(inLab_file, sheet_name=BMUs[bmu])
    #     dch1 = summary['DCH1'].to_numpy()
    #     dch2 = summary['DCH2'].to_numpy()
    #     rated_cap = 108.820559325236    #max cell DchCapAh for first test

    #     rows = summary.shape[0]
    #     cycles = np.zeros(rows)
    #     cycles[0] = dch1[0]/rated_cap
    #     for i in range(rows-1):
    #         cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])/rated_cap

    #     cell_soh = summary.iloc[:,3:27]/rated_cap*100

        # plot SOH vs cycle count 
        # fig, ax = plt.subplots(figsize=(9,6))
        # # fig, ax = plt.subplots(figsize=(3,6))
        # ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)
        # for i in range(2*12):
        #     bmuNum = int(np.floor(i/12)+1)
        #     cell_label = 'Cell ' + str(mods[bmuNum-1]) + '.' + str(np.mod(i,12)+1)
        #     ax.scatter(cycles, cell_soh.iloc[:,i], label=cell_label, color=rgb_to_dec(cell_colors[i]), marker=shapes[bmuNum-1])
        # ax.set_ylabel('SOH [%]', size=fs)
        # ax.set_xlabel('Cycles [-]', size=fs)
        # # ax.set_title('In-Lab: ' + strength[bmu-1] + ' Cell-Level SOH', size=fs)
        # ax.set_title(titles[bmu])
        # legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.68, -0.25), fancybox=True, ncol = 4, borderaxespad=0.)
        # plt.gca().add_artist(legend1)
        # low = cell_soh.min().min()
        # high = cell_soh.max().max() 
        # yrange = high - low
        # ax.set_ylim([low-yrange*.1, min(100,high+yrange*.1)])
        # ax.tick_params(direction='in')

        # ## plot upper and lower SOH bounds, avg SOH
        # # avg_SOH = np.mean(cell_soh, axis=1)
        # # std_SOH = np.std(cell_soh,  axis=1)
        # # coefs = np.polyfit(cycles, avg_SOH, 3)
        # # x = np.linspace(min(cycles), max(cycles), 1000)
        # # y_poly = np.polyval(coefs, x)
        # # line1 = ax.plot(x, y_poly, '-k', label='Average SOH')

        # ## plot standard deviation on right axes
        # # ax2 = ax.twinx()
        # # ax2.set_ylabel('Standard Deviation [%]', size=fs)
        # # line4 = ax2.plot(cycles, std_SOH, '--k', label='Standard Deviation')
        # # low = min(std_SOH)
        # # high = max(std_SOH)
        # # yrange = high - low
        # # ax2.set_ylim([low-yrange*.4, min(100,high+yrange*.1)])
        # # ax2.set_ylim([.1, .6])

        # # plot bmus separately 
        # ax2 = ax.twinx()
        # lns = []
        # for i in range(2):
        #     cell_sohi = cell_soh.iloc[:,i*12:i*12+12]
        #     avg_SOHi = np.mean(cell_sohi, axis=1)
        #     std_SOHi = np.std(cell_sohi,  axis=1)
        #     coefs = np.polyfit(cycles, avg_SOHi, 3)
        #     x = np.linspace(min(cycles), max(cycles), 1000)
        #     y_poly = np.polyval(coefs, x)
        #     lns += ax.plot(x, y_poly, '-', color=rgb_to_dec(cell_colors[i*12]), label='BMU0' + str(mods[i]) +  ' SOH')
        #     # line1[i] = ax.plot(x, y_poly, '-k', label='Average SOH')
 
        #     ax2.set_ylabel('Standard Deviation [%]', size=fs)
        #     lns +=  ax2.plot(cycles, std_SOHi, '--', color=rgb_to_dec(cell_colors[i*12]), label='BMU0' + str(mods[i])+ ' Standard Deviation')
        #     low = min(std_SOHi)
        #     high = max(std_SOHi)
        #     yrange = high - low
        #     # ax2.set_ylim([low-yrange*.4, min(100,high+yrange*.1)])
        #     # ax2.set_ylim([.1, .6])

        # # low = min(std_SOH)
        # # high = max(std_SOH)
        # # yrange = high - low
        # # ax2.set_ylim([low-yrange*.4, min(100,high+yrange*.1)])
        # ax2.set_ylim([.01, .55])
        
        # # put all lines in one legend on plot
        # # lns = line1 + line4
        # labls = [l.get_label() for l in lns]
        # # plt.legend(lns, labls, loc='lower left')
        # ax.legend(lns, labls, loc='upper center', bbox_to_anchor=(.15, -.25), fancybox=True, ncol = 1, borderaxespad=0.)

        # # legend2 = ax.legend(lns, labls,loc='upper center', bbox_to_anchor=(0.2, -0.25), fancybox=True, ncol = 1, borderaxespad=0.)
        # # plt.gca().add_artist(legend2)

        # plt.subplots_adjust(left=None, bottom=.5, right=None, top=None, wspace=None, hspace=None)
        # if save_plot: plt.savefig(outpath+' Cell-Level ' + strength[bmu] + '.jpg', dpi=1000, bbox_inches='tight')
        # plt.show()

    # ------- outdoor -------
    dch1 = outdoor.DCH1
    dch2 = outdoor.DCH2

    # rated_cap = max(outdoor.iloc[0,3:-1])
    rated_cap = 113

    rows = outdoor.shape[0]
    cycles = np.zeros(rows)
    cycles[0] = dch1[0]
    for i in range(rows-1):
        cycles[i+1] = cycles[i] + (dch2[i] + dch1[i+1])
    cycles = cycles/rated_cap


    fig, ax = plt.subplots(figsize=(12,8))
    # add dotted vertical lines at each data point  
    ax.vlines(cycles, 0, 120, linestyles='dashed', linewidth=1, colors='grey', zorder=0)

    # plot cells
    cell_soh = outdoor.iloc[:,13:109]/rated_cap*100
    for i in range(8*12):
        bmu = int(np.floor(i/12)+1)
        cell_label = 'Cell ' + str(bmu) + '.' + str(np.mod(i,12)+1)
        ax.scatter(cycles, cell_soh.iloc[:,i], label=cell_label, color=rgb_to_dec(cell_colors[i]), marker=shapes[bmu-1])
    ax.set_ylabel('SOH [%]', size=fs)
    ax.set_xlabel('Cycles [-]', size=fs)
    ax.set_title('Cummins State of Health outdoor: Cycle ' + str(round(max(cycles))), size=fs)
    legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), fancybox=True, ncol = 8, borderaxespad=0.)
    plt.gca().add_artist(legend1)

    # set y limit for cell SOH
    low = cell_soh.min().min()
    high = cell_soh.max().max()
    yrange = high - low
    ax.set_ylim([low-yrange*.1, min(100,high+yrange*.1)])

    # plot upper and lower SOH bounds, avg SOH
    avg_SOH = np.mean(outdoor.iloc[:,11:-1]/rated_cap*100, axis=1)
    std_SOH = np.std(outdoor.iloc[:,11:-1]/rated_cap*100,  axis=1)
    coefs = np.polyfit(cycles, avg_SOH, 4)
    x = np.linspace(min(cycles), max(cycles), 1000)
    y_poly = np.polyval(coefs, x)
    line1 = ax.plot(x, y_poly, '-k', label='Average SOH')
    # line1 = ax.plot(cycles, avg_SOH, '-k', label='Average SOH', marker='o')
    # line2 = ax.plot(cycles, avg_SOH+std_SOH, label='Upperbound SOH')
    # line3 = ax.plot(cycles, avg_SOH-std_SOH, '-g', label='Lowerbound SOH')

    # plot standard deviation on right axes
    ax2 = ax.twinx()
    ax2.set_ylabel('Standard Deviation [%]', size=fs)
    # ax2.set_ylim([.1, .6])
    line4 = ax2.plot(cycles, std_SOH, '--k', label='Standard Deviation')

    # put all lines in one legend on plot
    # lns = line1+line2+line3+line4
    lns = line1+line4
    labls = [l.get_label() for l in lns]
    plt.legend(lns, labls, loc='lower left')
    # plt.legend(lns, labls, loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, ncol = 8, borderaxespad=0.)
    # ax.set_xlim([0,cycles[-1]+4])

    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')
    plt.subplots_adjust(left=None, bottom=.5, right=None, top=None, wspace=None, hspace=None)
    if save_plot: plt.savefig(outpath+' SOH Summary.jpg', dpi=1000, bbox_inches='tight')
    if show_plot: plt.show()
    

# soh_summary(save_plot=True, show_plot=False)
# time_series(save_plot=True, show_plot=False)
monthly_meeting(save_plot=False, show_plot=True)
# all_cell_voltage(save_plot=False)
# trendline()
# inLab_data(save_plot=False)
# temps()
# cell_level(save_plot=True, show_plot=False)

