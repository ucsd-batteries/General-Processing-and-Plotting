from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import datetime as dt
from matplotlib.lines import Line2D
import statistics
from math import pi
from matplotlib.table import table
from matplotlib.patches import Patch
# ---------------------------- import data from summary spreadsheet ---------------------------- 
summary_file = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/Tesla_test_summary.xlsx'
summary_data = pd.read_excel(summary_file)
rated_cap = 200
fs = 15     # universal fontsize
save_plot = False   # Change to True to save plot as .png
outpath = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/'    # path for saving plot

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
idx2 = summary_data.columns.get_loc("cell_18")
cell_caps = summary_data.iloc[:,idx1:idx2+1].to_numpy()
cell_soh = cell_caps/rated_cap*100      # convert to %
soh_mean = np.mean(cell_soh,axis=1)
soh_std = np.std(cell_soh, axis=1)


#  ---------------------------- SOH summary plot ---------------------------- 
fig, ax = plt.subplots(figsize=(12,6))
# plot dashed vertical lines at each test index
ax.vlines(cycles,50,100,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.8)

# plot cells
shapes = ['o']*9 + ['s']*9
cmap = list(plt.get_cmap("tab10").colors)
cmap.pop(5)
for i in range(18):
    cell_name = 'Cell ' + str(int(np.floor(i/6)+1)) + '.' + str(int(i%6 + 1))   # Module.Cell
    ax.scatter(cycles, cell_soh[:,i], label=cell_name, marker=shapes[i], color=cmap[i%9], alpha=.6)
ax.set_ylabel("State of Health [%]", fontsize=fs)
ax.set_title("Tesla State of Health Summary: Cycle " + str(round(max(cycles))) + " of 2200", fontsize=fs)
# legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fancybox=True, ncol=6, borderaxespad=0.)
legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
plt.gca().add_artist(legend1)

# plot avg, lower bound, and upper bound
line1 = ax.plot(cycles, soh_mean, '-ok', label="Average SOH")
line2 = ax.plot(cycles, soh_mean+soh_std, label="Upper Bound SOH", color=cmap[0])
line3 = ax.plot(cycles, soh_mean-soh_std, label="Lower Bound SOH", color=cmap[1])

# calculate and set left axis lims
low = np.amin(cell_soh)
high = np.amax(cell_soh)
yrange = high-low
ax.set_ylim([low-.1*yrange, high+.1*yrange])

# plot standard deviation
ax.set_xlabel("Cycles [-]", fontsize=fs)
ax2 = ax.twinx()
line4 = ax2.plot(cycles, soh_std, '--k', label="Standard Deviation")
ax2.set_ylabel("Standard Deviation", fontsize=fs)

# calculate and set right axis lims
low = np.amin(soh_std)
high = np.amax(soh_std)
yrange = high-low
ax2.set_ylim([low-.1*yrange, high+yrange])

# put all lines in one legend on plot
lns = line1+line2+line3+line4
labls = [l.get_label() for l in lns]
# plt.legend(lns, labls, loc='center', bbox_to_anchor=(1.22, .5), fancybox=True, ncol = 2, borderaxespad=0.)
plt.legend(lns, labls, loc='upper right')
ax.set_xlim([0,cycles[-1]+4])

# make things look nice
ax.tick_params(direction='in')
ax2.tick_params(direction='in')
# plt.subplots_adjust(left=None, bottom=.2, right=.77, top=None, wspace=None, hspace=None)
plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

test_names = summary_data['test'].to_numpy()
# plt.tight_layout()
if save_plot: plt.savefig(outpath + test_names[-1] + ' SOH Summary.jpg', dpi=1000)

#  ---------------------------- gantt chart ---------------------------- 
task = 'Tesla 3s Aging'
progress = cycles[-1]/2190

# estimate end date
days_per_normal_test_cycle = 43      # estimate this from testing progress spreadsheet on drive
delay = 16+97                        # accumulate delays (from testing progress spreadsheet)
cycles_per_normal_test_cycle = cycles[-1] - cycles[-1-5]
cycles_per_day = cycles_per_normal_test_cycle/days_per_normal_test_cycle
cycles_remaining = np.ceil((6*365)-cycles[-1])
days_remaining = cycles_remaining/cycles_per_day+delay
end_date = dt.datetime.now() + dt.timedelta(days=days_remaining)

#Convert dates to datetime format
start_date = dt.datetime(2020, 11, 16)
start_number = 0
total_days = (end_date-start_date).days + 1
completed = progress*total_days
x_ticks = [round(i) for i in np.linspace(0, total_days+1, 20, endpoint=True)]
x_labels=[(start_date+dt.timedelta(days=i)).strftime('%d-%b') for i in x_ticks]

# plt.figure(figsize=(8,1.5))
fig, ax = plt.subplots(figsize=(12,1.5))
ax.set_title('Gantt Chart', size=10)
#Darker bar for completed part
ax.barh(y=task, left=start_number, width=completed, alpha=1, color='green', label='Completed')
#Light bar for entire task
ax.barh(y=task, left=start_number, width=total_days, alpha=0.4, color='green', label='Estimated Duration')
plt.gca().invert_yaxis()
plt.xticks(ticks=x_ticks[::3], labels=x_labels[::3], fontsize=9, rotation=-30, ha='left')
ax.grid(axis='x')
ax.legend(loc='center', bbox_to_anchor=(0.5, -2.6), fancybox=True, ncol=2, borderaxespad=0., fontsize=10)
plt.subplots_adjust(left=.35, bottom=.55, right=None, top=.7, wspace=None, hspace=None)
ax.axes.yaxis.set_visible(False)

# add tabulated vlaues
col_labels=['Task','Start Date','End Date', 'Progress']
table_vals=[['Tesla 3s Aging',str(start_date.month)+'-'+str(start_date.day)+'-'+str(start_date.year),
                str(end_date.month)+'-'+str(end_date.day)+'-'+str(end_date.year), str(round(progress*100,2))+'%']]
# the rectangle is where I want to place the table
the_table = plt.table(cellText=table_vals, colLabels=col_labels, cellLoc='center', colWidths=[.22, .20, .20, .15],
                fontsize=25, loc='center',bbox=[-.57, 0, .57, 2])
# add rectangle around plot so that title appears to be a table column
xlim_left, xlim_right = plt.xlim()
ax.add_patch(Rectangle((0,-1.32), xlim_right, 1.76, edgecolor='black', fill=False, clip_on=False))

test_names = summary_data['test'].to_numpy()
if save_plot: plt.savefig(outpath + test_names[-1] + ' Gantt Chart.jpg', dpi=1000)

# ------------------ characterization chart -------------------------
cell_std = np.std(cell_caps, axis=1)
end_limits = [200, 60, .8]
resistance = 120
current_val = [resistance, soh_mean[-1], cell_std[-1]]
polar_data = [resistance/end_limits[0]*100, end_limits[1]/soh_mean[-1]*100, cell_std[-1]/end_limits[2]*100]

fig, ax = plt.subplots(figsize=(6, 6))
ax = plt.subplot(projection='polar')
# data = [63, 85, 43]
startangle = 90
colors = ['#4393E5', '#43BAE5', '#7AE6EA']
xs = [(i * pi * 2) / 100 for i in polar_data]
ys = [-0.2, 1, 2.2]
left = (startangle * pi * 2) / 360  # this is to control where the bar starts
# plot bars and points at the end to make them round
for i, x in enumerate(xs):
    ax.barh(ys[i], x, left=left, height=1, color=colors[i])
    ax.scatter(x + left, ys[i], s=350, color=colors[i], zorder=2)
    ax.scatter(left, ys[i], s=350, color=colors[i], zorder=2)

plt.ylim(-4, 4)
# legend
lbl1 = ''
legend_elements = [Line2D([0], [0], marker='o', color='w', label='Resistance: {val} %'.format(val=round(polar_data[0],2)), markerfacecolor='#4393E5', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='SOH : {val} %'.format(val=round(polar_data[1],2)), markerfacecolor='#43BAE5', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='Cell Imbalance: {val} %'.format(val=(round(polar_data[2],2))), markerfacecolor='#7AE6EA', markersize=10)]
ax.legend(handles=legend_elements, loc='best', frameon=True, prop={'size': 7})
ax.set_title('Battery Chararcteristics', size=12)
# clear ticks, grids, spines
plt.xticks([])
plt.yticks([])
ax.spines.clear()
a = np.array([current_val])
b = np.array([end_limits])
c = np.array([polar_data])
data2 = np.round(np.concatenate((a.transpose(), b.transpose(), c.transpose()), axis=1),2)
rows = ['Resistance [% increase]', 'SOH', 'Cell Imbalance']
columns = ['Current Value', 'End Limit', 'Completion %']
the_table = plt.table(cellText=data2,
                      rowLabels=rows,
                      colLabels=columns,
                      loc='bottom',)
the_table.auto_set_font_size(False)
the_table.set_fontsize(11)
plt.subplots_adjust(bottom=.2, hspace=.6)
plt.show()
