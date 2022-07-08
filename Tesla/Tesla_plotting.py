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
summary_file2 = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/Tesla_Historic.xlsx'
summary_data = pd.read_excel(summary_file)
summary_data2 = pd.read_excel(summary_file2)
rated_cap = 200
fs = 15     # universal fontsize
save_plot = False   # Change to True to save plot as .png
outpath = r'/Users/liamk/OneDrive/Desktop/ESS LAB/Data/Tesla/'    # path for saving plot

#  ---------------------------- calculate cycles ----------------------------
# cycle for SOH
rows = summary_data.shape[0]
dch1 = summary_data['DCH1'].to_numpy()
dch2 = summary_data['DCH2'].to_numpy()
cycles1 = np.zeros((rows,))
cycles1[0] = dch1[0]
for i in range(rows-1):
    cycles1[i+1] = cycles1[i] + dch1[i+1] + dch2[i]
cycles1 = cycles1/rated_cap

# cycle for other plots
rows = summary_data2.shape[0]
dch1 = summary_data2['DCH1'].to_numpy()
dch2 = summary_data2['DCH2'].to_numpy()
cycles2 = np.zeros((rows,))
cycles2[0] = dch1[0]
for i in range(rows-1):
    cycles2[i+1] = cycles2[i] + dch1[i+1] + dch2[i]
cycles2 = cycles2/rated_cap

#  ---------------------------- get cell capacities and convert to SOH (State of Health) ----------------------------
idx1 = summary_data.columns.get_loc("cell_1")
idx2 = summary_data.columns.get_loc("cell_18")
idx3 = summary_data.columns.get_loc("cell_6")
idx4 = summary_data.columns.get_loc("cell_12")

cell_caps = summary_data.iloc[:,idx1:idx2+1].to_numpy()

cell_soh = cell_caps/rated_cap*100      # convert to %
soh_mean = np.mean(cell_soh,axis=1)
soh_std = np.std(cell_soh, axis=1)
p1_soh = cell_soh[:, 0:6]
p2_soh = cell_soh[:, 6:12]
p3_soh = cell_soh[:, 12:18]
p1_mean = np.mean(p1_soh, axis=1)
p2_mean = np.mean(p1_soh, axis=1)
p3_mean = np.mean(p1_soh, axis=1)

V_std = summary_data2['ave voltage std'].to_numpy()
V_ave = summary_data2['ave voltage'].to_numpy()

r_idx1 = summary_data2.columns.get_loc("Cell 1")
r_idx2 = summary_data2.columns.get_loc("Cell 18")
r_idx3 = summary_data2.columns.get_loc("Cell 6")
r_idx4 = summary_data2.columns.get_loc("Cell 12")

rcell_1 = summary_data2.iloc[:, r_idx1:r_idx3+1].to_numpy()
rcell_2 = summary_data2.iloc[:, r_idx3+1:r_idx4+1].to_numpy()
rcell_3 = summary_data2.iloc[:, r_idx4+1:r_idx2+1].to_numpy()

rpack_1 = np.sum(rcell_1, axis=1)
rpack_2 = np.sum(rcell_2, axis=1)
rpack_3 = np.sum(rcell_3, axis=1)
# print(rpack_1)
print(np.size(rpack_1))
rstd_1 = np.std(rcell_1, axis=1)
rstd_2 = np.std(rcell_2, axis=1)
rstd_3 = np.std(rcell_3, axis=1)
rstd_all = [rstd_1, rstd_2, rstd_3]
rstd_tot = np.mean(rstd_all, axis=0)
rave_1 = np.mean(rcell_1, axis=1)
rave_2 = np.mean(rcell_2, axis=1)
rave_3 = np.mean(rcell_3, axis=1)

rave_all = [rave_1, rave_2, rave_3]
rave_tot = np.mean(rave_all, axis=0)
# print(rave_tot)
print(np.size(rave_tot))
print(np.size(rstd_tot))
# print(rstd_tot)
cell_res = summary_data2.iloc[:, r_idx1:r_idx2+1].to_numpy()
total_res = summary_data2['Total'].to_numpy()
res_std = summary_data2['Res STD'].to_numpy()
res_ave = np.mean(cell_res, axis=1)


# -------------------------State of Health------------------------------
fig, ax = plt.subplots(figsize=(12,6))
# plot dashed vertical lines at each test index
ax.vlines(cycles1,50,100,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.8)

# plot cells
shapes = ['o']*9 + ['s']*9
cmap = list(plt.get_cmap("tab10").colors)
cmap.pop(5)
for i in range(18):
    cell_name = 'Cell ' + str(int(np.floor(i/6)+1)) + '.' + str(int(i%6 + 1))   # Module.Cell
    ax.scatter(cycles1, cell_soh[:,i], label=cell_name, marker=shapes[i], color=cmap[i%9], alpha=.6)
ax.set_ylabel("State of Health [%]", fontsize=fs)
ax.set_title("Tesla State of Health Summary: Cycle " + str(round(max(cycles1))) + " of 1500", fontsize=fs)
# legend1 = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fancybox=True, ncol=6, borderaxespad=0.)
legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
plt.gca().add_artist(legend1)
print(np.size(soh_mean))
print(np.size(soh_std))
# plot avg, lower bound, and upper bound
# line1 = ax.plot(cycles1, soh_mean, '-ok', label="Average SOH")
line2 = ax.plot(cycles1, soh_mean+soh_std, label="Upper Bound SOH", color=cmap[0])
line3 = ax.plot(cycles1, soh_mean-soh_std, label="Lower Bound SOH", color=cmap[1])
line5 = ax.plot(cycles1, p1_mean, '-', label="Pack 7 SOH", color=cmap[2])
line6 = ax.plot(cycles1, p2_mean, '-', label="Pack 2 SOH", color=cmap[3])
line7 = ax.plot(cycles1, p3_mean, '-', label="Pack 6 SOH", color=cmap[4])

# line5 = plt.plot(x, soh_interp, '-', label="polyfit")
# calculate and set left axis lims
low = np.amin(cell_soh)
high = np.amax(cell_soh)
yrange = high-low
ax.set_ylim([low-.1*yrange, high+.1*yrange])

# plot standard deviation
ax.set_xlabel("Cycles [-]", fontsize=fs)
ax2 = ax.twinx()
line4 = ax2.plot(cycles1, soh_std, '--k', label="Standard Deviation")
ax2.set_ylabel("Standard Deviation", fontsize=fs)

# calculate and set right axis lims
low = np.amin(soh_std)
high = np.amax(soh_std)
yrange = high-low
ax2.set_ylim([low-.1*yrange, high+yrange])

# put all lines in one legend on plot
lns = line2+line3+line4+line5+line6+line7
labls = [l.get_label() for l in lns]
# plt.legend(lns, labls, loc='center', bbox_to_anchor=(1.22, .5), fancybox=True, ncol = 2, borderaxespad=0.)
plt.legend(lns, labls, loc='upper right')
ax.set_xlim([0,cycles1[-1]+4])

# make things look nice
ax.tick_params(direction='in')
ax2.tick_params(direction='in')
# plt.subplots_adjust(left=None, bottom=.2, right=.77, top=None, wspace=None, hspace=None)
plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

test_names = summary_data['test'].to_numpy()
# plt.tight_layout()
if save_plot: plt.savefig(outpath + test_names[-1] + ' SOH Summary.jpg', dpi=1000)

# ----------------------------Pack Res ----------------------------------

upper = np.multiply(np.ones((len(cycles2),1)), .02)
fig, ax = plt.subplots(figsize=(12, 6))
ax.vlines(cycles2,0,.025,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.8)
shapes = ['o']*9 + ['s']*9
cmap = list(plt.get_cmap("tab10").colors)
cmap.pop(5)

line1 = ax.plot(cycles2, upper, '-', label="Upper Bound Resistance",color=cmap[0])
line2 = ax.plot(cycles2, rpack_1, '-', label='Pack 7 Resistance', marker='o', alpha=.6,color=cmap[1])
line3 = ax.plot(cycles2, rpack_2, '-', label='Pack 2 Resistance', marker='o', alpha=.6, color=cmap[2])
line4 = ax.plot(cycles2, rpack_3, '-', label='Pack 6 Resistance', marker='o', alpha=.6, color=cmap[3])
# line5 = ax.plot(cycles2, rave_tot-rstd_tot, '-', label="Lower Bound Res", color=cmap[5])
# line6 = ax.plot(cycles2, rave_tot+rstd_tot, '-', label="Upper Bound Res", color=cmap[6])
# line7 = ax.plot(cycles2, rave_tot, '-', label="Upper Bound Res", color=cmap[7])
lns = line1+line2+line3+line4
labls = [l.get_label() for l in lns]

ax.set_ylabel("Resistance [Ohms]", fontsize='15')
ax.set_xlabel('Cycles [-]', fontsize='15')
ax.set_title("Tesla Cell Resistance: Cycle " + str(round(max(cycles1))) + " of 1500", fontsize='19')
plt.legend(lns, labls, loc='upper right')
plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

for i in range(18):
    cell_name = 'Cell ' + str(int(np.floor(i/6)+1)) + '.' + str(int(i%6 + 1))   # Module.Cell
    ax.scatter(cycles2, cell_res[:,i], label=cell_name, marker=shapes[i], color=cmap[i%9], alpha=.6)


legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
plt.gca().add_artist(legend1)

ax.tick_params(direction='in')

# -----------------------Res STD ----------------------------------

fig, ax = plt.subplots(figsize=(12, 6))
upper = np.multiply(np.ones((len(cycles2),1)), .0005)
ax.vlines(cycles2,0,.0006,colors='grey',linestyles='dashed', linewidth=1, zorder=0, alpha=.8)
shapes = ['o']*9 + ['s']*9
cmap = list(plt.get_cmap("tab10").colors)
cmap.pop(5)
line1 = ax.plot(cycles2, rstd_1, '--', label="Pack 7 STD", color=cmap[0])
line2 = ax.plot(cycles2, rstd_2, '--', label="Pack 2 STD", color=cmap[1])
line3 = ax.plot(cycles2, rstd_3, '--', label="Pack 6 STD", color=cmap[2])
line4 = ax.plot(cycles2, upper, '-', label="Upper Bound Resistance STD")
lns = line1 + line2 + line3 + line4
labls = [l.get_label() for l in lns]
ax.set_ylabel("Resistance STD [Ohms]", fontsize='15')
ax.set_xlabel('Cycles [-]', fontsize='15')
ax.set_title("Tesla Cell Resistance STD: Cycle " + str(round(max(cycles1))) + " of 1500", fontsize='19')
plt.legend(lns, labls, loc='upper right')
plt.subplots_adjust(left=None, bottom=None, right=.8, top=None, wspace=None, hspace=None)

ax.tick_params(direction='in')

# for i in range(18):
#     cell_name = 'Cell ' + str(int(np.floor(i/6)+1)) + '.' + str(int(i%6 + 1))   # Module.Cell
#     ax.scatter(cycles2, cell_res[:,i], label=cell_name, marker=shapes[i], color=cmap[i%9], alpha=.6)
#
#
# legend1 = ax.legend(loc='center', bbox_to_anchor=(1.17, .5), fancybox=True, ncol=1, borderaxespad=0.)
# plt.gca().add_artist(legend1)

plt.show()