
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

# Figure 01


sns.despine();sns.set(context='paper', style='ticks');plt.grid(0);colors = sns.color_palette('muted')
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 



#
cells_conc = {'an_WT_ESL': 16,'an_WT_N2chHep':4}
#conc_labels = {'an_0ng': 'KO-0ng of FGF','an_2_5ng':'KO-2.5ng of FGF','an_5ng':'KO-5ng of FGF',
#               'an_10ng':'KO-10ng of FGF','an_20ng': 'KO-20ng of FGF','an_WT_ESL':
#                   'WT 1','an_WT_N2ChHep':'WT 2' }

#cells_conc = {'an_0ng': 19,'an_2_5ng':11,'an_5ng':48,'an_10ng':54,'an_20ng': 51}

conc_labels = {'an_0ng': 'KO-0ng of FGF','an_2_5ng':'KO-2.5ng of FGF','an_5ng':'KO-5ng of FGF',
               'an_10ng':'KO-10ng of FGF','an_20ng': 'KO-20ng of FGF','an_WT_ESL':
                  'WT-ESL','an_WT_N2chHep':'WT-N2chHep' }



fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69/6))
gs_main = gridspec.GridSpec(nrows=1, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.95, hspace=0.0)
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[0], wspace=0.4, hspace=0.1)


inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc), ncols=1, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc), ncols=1, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)
inner_gs2 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc), ncols=1, subplot_spec=gs0[2], wspace=0.0, hspace=0.0)

#conc_aux = TimeSerieDC.TimeSerieDC(data_folder,dataset_name,save_folder)
conc = conc_aux.data
traze_value = conc_aux.traze_value
#th = 2500; order = 2
ylim = conc_aux.ylim
M = conc_aux.traze_max_values

for k,c in enumerate(cells_conc):
    df = conc[c];df.sort_index(inplace=True); ax = plt.subplot(inner_gs0[k]); color=colors[k]
    
    time = len(df.loc[cells_conc[c]][traze_value].values)
    ax.plot(np.arange(time),df.loc[cells_conc[c]][traze_value].values, color=colors[k], linewidth=0.5);
    ax.set_xlim([0, None]);
    ax.set_ylim( [35000,65000])
    ax.text(0.8, 0.1, conc_labels[c], ha='center', va='center', transform=ax.transAxes, fontsize=5)


    if k==len(cells_conc)-1:
        set_scale(ax,[0, 400],[35000,65000]); ax.set_ylabel('Intensity \n (a.u.)',fontsize=10); ax.set_xlabel('Frames',fontsize=10)
        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
    else:
        silent_ax(ax)



conc = filtmin_conc.data
peaks_conc = filtmin_conc.peaks
traze_value =filtmin_conc.traze_value
th = 2500; order = 2
ylim = filtmin_conc.ylim
M = filtmin_conc.traze_max_values


for k,c in enumerate(cells_conc):
    df = conc[c];df.sort_index(inplace=True); ax = plt.subplot(inner_gs1[k]); color=colors[k]
    
    time = len(df.loc[cells_conc[c]][traze_value].values)
    ax.axhline(y = 0, linewidth=0.5, color='gray',alpha = 0.8)
    ax.plot(np.arange(time),df.loc[cells_conc[c]][traze_value].values, color=colors[k], linewidth=0.5);
    ax.set_xlim([0, None]);
    ax.set_ylim( [-2000, 6000])

    if k==len(cells_conc)-1:
        set_scale(ax,[0, 400],[-2000,0,6000]); ax.set_ylabel('Intensity \n (a.u.)',fontsize=10); ax.set_xlabel('Frames',fontsize=10)
        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
    else:
        silent_ax(ax)


    ax1 = plt.subplot(inner_gs2[k]);
    ax1.hist(df[traze_value], 50, density=1, facecolor=color, alpha=0.8); 
    ax1.set_xlim([-2000,6000]);ax1.set_ylim([0,0.0005])


    if k==len(cells_conc)-1:
        set_scale(ax1,[-2000,0,6000],[0,0.0005]); ax1.set_xlabel('Intensity \n (a.u.)',fontsize=10); ax1.set_ylabel('Ocurrences',fontsize=10,rotation = 270)
        ax1.xaxis.set_label_coords(0.5, -0.35);ax1.yaxis.set_label_coords(1.01,0.5)
    else:
        silent_ax(ax1)



fig.savefig('/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Paper_Figs/Figure1.pdf', format='pdf')

