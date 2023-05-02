#New Fig. 2 Supp 4: Autocorrelation plots of serum + LIF data. 

import sys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from DataAnalysis.Preprocesing.PlotData import silent_ax
import matplotlib
import numpy as np
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.5

matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 0.25
matplotlib.rcParams['xtick.minor.size'] = 10
matplotlib.rcParams['xtick.minor.width'] = 0.25

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')
from main_oscillations import autocorr_normed
from exponential_dm_main import download_data

conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')

#%%
ylim = [-1,1]
xlim = [-1,330] #110 min
long=1;cols=5
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'
save_name ='NEW_Fig2Supp4'

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))

fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=4, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True);
handles = [None] * len(conc_data)

Cols = 5
for k,C in enumerate(conc_labels):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[C];
    df.sort_index(inplace=True);
    cells = df.index.get_level_values(0).unique()

    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        data = df.query('cell ==' + str(cell))
        AC = autocorr_normed(data.sm_MEAN_INTENSITY.values - np.mean(data.sm_MEAN_INTENSITY.values))#data.detrended.values)
    
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
    
        if cells[-1] >= cell:
            handles[k],= ax.plot(AC, color=color, linewidth=0.4,label=label)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            #set_scale(ax, [], ylim);
            ax.set_ylabel('autocorrelation', fontsize=8);
            ax.set_xlabel('time lag (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels([-1,1], fontsize=6)

        else:
            silent_ax(ax)

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + save_name + '.pdf', format='pdf')
