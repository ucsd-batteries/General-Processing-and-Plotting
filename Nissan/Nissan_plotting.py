from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import datetime as dt

def SOH_Summary(summary_file, pack_name, legend_loc='lower right', save_plot=False):
    # ---------------------------- import data from summary spreadsheet ---------------------------- 
    summary_data = pd.read_csv(summary_file)
    rated_cap = 56.3
    cell_num = 16
    fs = 15     # universal fontsize
    outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Nissan/'    # path for saving plot

    #  ---------------------------- calculate cycles ---------------------------- 
    rows = summary_data.shape[0]
    dch1 = summary_data['DCH1'].to_numpy()
    dch2 = summary_data['DCH2'].to_numpy()
    cycles = np.zeros((rows,))
    cycles[0] = dch1[0]
    for i in range(rows-1):
        cycles[i+1] = cycles[i] + dch1[i+1] + dch2[i]
    cycles = cycles/rated_cap

    #  ---------------------------- get cell capacities and convert to SOH (State of Health) ---------------------------- 
    idx1 = summary_data.columns.get_loc("cell_1")
    cell_caps = summary_data.iloc[:,idx1:idx1+cell_num].to_numpy()
    cell_soh = cell_caps/rated_cap*100      # convert to %
    soh_mean = np.mean(cell_soh,axis=1)
    soh_std = np.std(cell_soh, axis=1)

    #  ---------------------------- SOH summary plot ---------------------------- 
    fig, ax = plt.subplots(figsize=(12,6))
    # plot dashed vertical lines at each test index
    ax.vlines(cycles,50,100,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.6)

    # plot cells
    half_cell_num = int(cell_num/2)
    shapes = ['o']*half_cell_num + ['s']*half_cell_num
    cmap = list(plt.get_cmap("tab10").colors)
    cmap.pop(5)
    for i in range(cell_num):
        cell_name = 'Cell ' + str(i+1)   
        ax.scatter(cycles, cell_soh[:,i], label=cell_name, marker=shapes[i], color=cmap[i%half_cell_num], alpha=.5)
    ax.set_ylabel("State of Health [%]", fontsize=fs)
    ax.set_title(pack_name + " State of Health Summary: Cycle " + str(round(max(cycles))) + " of 2190", fontsize=fs)
    # legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fancybox=True, ncol=6, borderaxespad=0.)
    legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
    plt.gca().add_artist(legend1)

    # plot avg, lower bound, and upper bound
    line1 = ax.plot(cycles, soh_mean, '-k', label="Average SOH")
    line2 = ax.plot(cycles, soh_mean+soh_std, label="Upper Bound SOH", color=cmap[0])
    line3 = ax.plot(cycles, soh_mean-soh_std, label="Lower Bound SOH", color=cmap[2])

    # calculate and set left axis lims
    low = np.amin(cell_soh)
    high = np.amax(cell_soh)
    yrange = high-low
    ax.set_ylim([low-.1*yrange, high+.1*yrange])

    # plot standard deviation
    ax.set_xlabel("Cycles [-]", fontsize=fs)
    ax2 = ax.twinx()
    line4 = ax2.plot(cycles, soh_std, '--k', label="Standard Deviation", alpha=0.8)
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
    ax.set_xlim([0,cycles[-1]+4])

    # make things look nice
    ax.tick_params(direction='in')
    ax2.tick_params(direction='in')
    # plt.subplots_adjust(left=None, bottom=.2, right=.77, top=None, wspace=None, hspace=None)
    plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

    test_names = summary_data['test'].to_numpy()
    # plt.tight_layout()
    if save_plot: plt.savefig(outpath + test_names[-1] + ' SOH Summary.jpg', bbox_inches='tight', dpi=1000)

    return cycles

save_plot = True   # Change to True to save plot as .png
#  ---------------------------- NP5 Summary Plot ---------------------------- 
summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/NP5/NP5_test_summary.csv'
pack_name = 'NP5'  
legend_loc = 'lower right'
NP5_cycles = SOH_Summary(summary_file, pack_name, legend_loc, save_plot)

#  ---------------------------- NP6 Summary Plot ---------------------------- 
summary_file = r'/Users/quiana/Documents/UCSD/CER/Data_Processing/Processing/Nissan/NP6/NP6_test_summary.csv'
pack_name = 'NP6'  
legend_loc = 'upper right'
NP6_cycles = SOH_Summary(summary_file, pack_name, legend_loc, save_plot)

#  ---------------------------- gantt chart ---------------------------- 
outpath = r'/Users/quiana/Documents/UCSD/CER/Plots/Nissan/'    # path for saving plot
task = ['Nissan Pack 3', 'Nissan Pack 6', 'Nissan Pack 5']
status = ['Terminated', 'In Progress', 'In Progress']   # testing status for NP3, NP6, NP5; respectively
# print(NP5_cycles[-1])
# print(NP6_cycles[-1])
progress = [463/2190, NP6_cycles[-1]/2190, NP5_cycles[-1]/2190] # NP3, NP5, NP6

# estimate end date
def est_end_date(cycles, days_per_normal_test_cycle, delay):
    # days_per_normal_test_cycle: estimate this from testing progress spreadsheet on drive
    # delay: accumulate delays (from testing progress spreadsheet)
    cycles_per_normal_test_cycle = cycles[-1] - cycles[-1-5]
    cycles_per_day = cycles_per_normal_test_cycle/days_per_normal_test_cycle
    cycles_remaining = np.ceil((6*365)-cycles[-1])
    days_remaining = cycles_remaining/cycles_per_day+delay
    end_date = dt.datetime.now() + dt.timedelta(days=days_remaining)
    return end_date
NP3_end_date = dt.datetime(2023, 6, 1)
NP5_end_date = est_end_date(NP5_cycles, days_per_normal_test_cycle=34.5, delay=23)
NP6_end_date = est_end_date(NP6_cycles, days_per_normal_test_cycle=35.8, delay=108)
end_date = [NP3_end_date, NP6_end_date, NP5_end_date]
end_date_str = [str(d.month)+'-'+str(d.day)+'-'+str(d.year) for d in end_date]

#Convert dates to datetime format
start_date = [dt.datetime(2020, 7, 7), dt.datetime(2020, 8, 19), dt.datetime(2021, 3, 18)]
start_date_str = [str(d.month)+'-'+str(d.day)+'-'+str(d.year) for d in start_date]
start_num = [(d - dt.datetime(2020, 7, 7)).days for d in start_date]
start_number = [d - min(start_num) for d in start_num]
total_days = [(end_date[i]-start_date[i]).days + 1 for i in range(len(start_date))]
completed = [progress[i]*total_days[i] for i in range(len(progress))]
x_ticks = [round(i) for i in np.linspace(0, max(total_days)+1, 20, endpoint=True)]
x_labels=[(dt.datetime(2020, 7, 7)+dt.timedelta(days=i)).strftime('%d-%b') for i in x_ticks]

# plt.figure(figsize=(8,1.5))
fig, ax = plt.subplots(figsize=(12,2.2))
ax.set_title('Gantt Chart', size=10)
#Darker bar for completed part
ax.barh(y=task, left=start_number, width=completed, alpha=1, color='green', label='Completed')
#Light bar for entire task
ax.barh(y=task, left=start_number, width=total_days, alpha=0.4, color='green', label='Estimated Duration')
plt.gca().invert_yaxis()
plt.xticks(ticks=x_ticks[::3], labels=x_labels[::3], fontsize=9, rotation=-30, ha='left')
ax.grid(axis='x')
ax.legend(loc='center', bbox_to_anchor=(0.5, -.8), fancybox=True, ncol=2, borderaxespad=0., fontsize=10)
plt.subplots_adjust(left=.37, bottom=.35, right=None, top=.7, wspace=None, hspace=None)
ax.axes.yaxis.set_visible(False)

# add tabulated vlaues
col_labels=['Task','Start Date','End Date', 'Progress', 'Status']
table_vals=[[task[i], start_date_str[i], end_date_str[i], str(round(progress[i]*100,2))+'%', status[i]] for i in range(len(task))]
# the rectangle is where I want to place the table
the_table = plt.table(cellText=table_vals, colLabels=col_labels, cellLoc='center', colWidths=[.25, .20, .20, .16, .21],
                loc='center',bbox=[-.64, 0, .64, 1.33])
the_table.set_fontsize(10)
the_table.auto_set_font_size(False)
# add rectangle around plot so that title appears to be a table column
xlim_left, xlim_right = plt.xlim()
ax.add_patch(Rectangle((0,-1.55), xlim_right, 1, edgecolor='black', fill=False, clip_on=False))

if save_plot: plt.savefig(outpath + ' Nissan Gantt Chart.jpg', bbox_inches='tight', dpi=1000)

plt.show()
