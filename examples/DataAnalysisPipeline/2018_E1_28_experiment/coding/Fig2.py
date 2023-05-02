
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale
from DataAnalysis.Statistics.PlotStatistics import amplitud_histograms

# Figure 02


sns.despine();sns.set(context='paper', style='ticks');plt.grid(0);colors = sns.color_palette('muted')
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 


cells=[10,38,6]
conc_labels = {'an_SpryCNP_ESL': 'WT-ESL'}


fig = plt.figure(constrained_layout=False, figsize=(8.27/2, 11.69/4))
gs_main = gridspec.GridSpec(nrows=1, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.95, hspace=0.0)
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs_main[0], wspace=0.4, hspace=0.1)


inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells), ncols=1, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)


conc = filtmin_conc.data
traze_value = filtmin_conc.traze_value
th = 1; order = 6
ylim = filtmin_conc.ylim
M = filtmin_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)


df = conc['an_SpryCNP_ESL'];df.sort_index(inplace=True);

for k,c in enumerate(cells):
    ax = plt.subplot(inner_gs0[k]); color=colors[k]
    
    time = len(df.loc[c][traze_value].values)
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors[k], linewidth=0.5);
    ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=0.8, alpha=1)
    ax.set_xlim([0, 171]);
    ax.set_ylim(  [-6,21])
    #ax.text(0.8, 0.1, conc_labels[c], ha='center', va='center', transform=ax.transAxes, fontsize=5)


    if k==len(cells)-1:
        set_scale(ax,[0, 170], [-6,21]); ax.set_ylabel('Intensity \n (a.u.)',fontsize=10,rotation=270); ax.set_xlabel('Frames',fontsize=10)
        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(1.01,0.5)
    else:
        silent_ax(ax)


k = 0
ax = plt.subplot(inner_gs1[k]); color=colors[k]
range_values = [0, 170]; num_bins = 50

n, bins, patches = ax.hist(conc[C][plot_value], num_bins, range=range_values[j], density=1,
                           facecolor=color[i], alpha=0.8, cumulative=cumulative)
ax.set_xlim(range_values)
ax.grid(False)

set_scale(ax,[0, 170], [0,0.001]); ax.set_ylabel('Ocurrences',fontsize=10); ax.set_xlabel('Intensity (a.u.)',fontsize=10)
ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(1.01,0.5)




fig.savefig('/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Paper_Figs/Figure2.pdf', format='pdf')
