#New Fig. 2 Supp 3: Details of stochastic pulsing hypothesis. 
#In the first row, the time between pulses (ex data), waiting time (simulation) and pulse duration (ex data) histograms. 
#Both waiting time and time between pulses histograms include the exponential fit. 

#On the second row, the simulated traces (population hypothesis). Third row, shuffle hypothesis traces. 
#Bottom row, experimental data.  All traces with pulses (maxima).
#%%

import sys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale
import matplotlib
import numpy as np
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.75

matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.major.width'] = 0.75
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['xtick.minor.width'] = 0.75
plt.rcdefaults()


sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')

from exponential_dm_main import download_data,get_linear_fit,get_population_time_series,get_single_cell_time_series,find_dm,get_single_cell_time_series_moving_boundaries,get_single_cells_pulses_shuffle
#%%

# =============================================================================
# calculating the input variables for population simulation
# =============================================================================
conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')

find_dm(conc_data['an_WT_ESL'])

bins1 = plt.hist(conc_data['an_WT_ESL']['dm_'],bins=np.arange(2,200,6), density=1);
popt,pcov = get_linear_fit(bins1)
lambda_ = popt[0]   

exponential = np.random.exponential(1/lambda_,len(conc_data['an_WT_ESL']['dt_peaks']))
orange =  sns.color_palette()[1]
colors =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))
colors = [colors[i] for i in [15,25]]


 #%%
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=4, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.5,wspace= 0.0)

gs_main_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[0], wspace=0.5, hspace=0.0)
# =============================================================================
# dm
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'dm_'
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_main_row_1[0]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,200,6), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    #ax1.axvspan(0,dm_max*60/20, color='gray',alpha = 0.4,linewidth=0.0)
    ax1.plot(lambda_*np.exp(-lambda_*np.arange(200)),color = orange,alpha = 1,linewidth=0.8)

    ax1.set_xlim([0,200]);ax1.set_ylim([0,0.045])

    mode_dm =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_dm = 'mode (m): '+str(np.round(mode_dm*20/60,2))+' \n'
    ax1.text(0.65, 0.90, mode_dm, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    mode_dm =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_dm = 'mode (frames): '+str(np.round(mode_dm,2))+' \n'
    ax1.text(0.65, 0.70, mode_dm, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(df[col3].quantile(0.25),2))+' ; '+str(np.round(df[col3].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75),2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    ax1.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(df[col3]),2)) +'\n std : '+str(np.round(np.std(df[col3]),2)), ha='center', va='center', transform=ax1.transAxes, fontsize=6)

set_scale(ax1,np.arange(2,200,60),[0,0.045]); ax1.set_xlabel('time between pulses (min)',fontsize=10); ax1.set_ylabel('probability density (1/min)',fontsize=10,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=8)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.045]],fontsize=8)


ax1.tick_params(labelsize=8,direction='out', pad=1,length=2)
# =============================================================================
# exponential distrobution
# 
# =============================================================================
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_main_row_1[1]); 
    

    bins3 = ax3.hist(exponential,bins=np.arange(0,200,6), density=1, facecolor=orange, alpha=0.8,linewidth=0.0 )
    ax3.plot(lambda_*np.exp(-lambda_*np.arange(400)),color = orange,alpha = 1,linewidth=0.8)
    ax3.set_xlim([0,200]);ax3.set_ylim([0,0.045]);

    mode_exp =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_exp = 'mode (m): '+str(np.round(mode_exp*20/60,2))+' \n'
    ax3.text(0.65, 0.90, mode_exp, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(np.quantile(exponential,0.25)*20/60,2))+' ; '+str(np.round(np.quantile(exponential,0.50)*20/60,2))+' ; '
                                  +str(np.round(np.quantile(exponential,0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
 

    mode_exp =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_exp = 'mode (f): '+str(np.round(mode_exp,2))+' \n'
    ax3.text(0.65, 0.7, mode_exp, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(np.quantile(exponential,0.25),2))+' ; '+str(np.round(np.quantile(exponential,0.50),2))+' ; '
                                  +str(np.round(np.quantile(exponential,0.75),2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
 

    ax3.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(exponential),2)) +'\n std : '+str(np.round(np.std(exponential),2)), ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    

ax3.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=8)
ax3.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.045]],fontsize=8)

set_scale(ax3,np.arange(2,200,60),[0,0.045]); ax3.set_xlabel(r'waiting time (min)',fontsize=10); ax3.set_ylabel('probability density (1/min)',fontsize=10,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=8)
ax3.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.045]],fontsize=8)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.12,0.5)
ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_main_row_1[2]); 
    
    bins3 = ax3.hist(df[col2],bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.1]);

    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode (m): '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(0.8, 1.07, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.8, 1.05,r'$Q$ (m): '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    
    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode (f): '+str(np.round(mode_dt,2))+' \n'
    ax3.text(0.8, 0.92, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.8, 0.9,r'$Q$ (f): '+str(np.round(df[col2].quantile(0.25),2))+' ; '+str(np.round(df[col2].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75),2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
  

set_scale(ax3,np.arange(0,100,15),[0,0.1]); ax3.set_xlabel('duration (min)',fontsize=10); ax3.set_ylabel('probability density (1/min)',fontsize=10,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=8)
ax3.set_yticklabels([np.round(y/(20/60),2)  for y in [0,0.1]],fontsize=8)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.0,0.5)


ax3.tick_params(labelsize=8,direction='out', pad=1,length=2)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
#
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/population_data/'
save_name = 'population0'

save_folder_save_name = save_folder + save_name
df_forecast = download_data(save_folder_save_name + '.pkl')

ylim = [-0.2,1.1]
xlim = [-15,330]
long=1;Cols=5
 
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rc('axes.spines', top=True, bottom=True, left=True, right=True);

k = 0
label = 'an_WT_ESL'
Rows = len_cells[k] // Cols;
if len_cells[k] % Cols != 0:
    Rows = Rows + 1
inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_main[2], wspace=0, hspace=0.0)

df = df_forecast
cells = df.index.get_level_values(0).unique()


for zo,cell in enumerate(cells):
    ax = plt.subplot(inner_gs[zo])
    ax.set_ylim(ylim);
    ax.set_xlim(xlim)
    data = df.query('cell ==' + str(cell))

    if zo==0: ax.text(0.0, 1.4, 'HoPM', ha='left', va='center', transform=ax.transAxes, fontsize=10,color = color)

    if cells[-1] >= cell:
        ax.plot(data.exponential_IPI.values ,color=orange, linewidth=0.4,label=label)
        ax.plot(data.amp_peaks.values,'o',color = orange,markersize=1)
        #ax.text(0.05, 0.2, str(data.amp_peaks.count()), ha='center', va='center', transform=ax.transAxes, fontsize=6,color = orange)
        ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
    ax.grid(False);
    ax.tick_params(labelsize=8, direction='out', pad=1, length=2)

    
    if zo==((Rows-1)*Cols):
        ax.set_ylabel('amplitude (a.u)', fontsize=10);
        ax.set_xlabel('time (min)', fontsize=10)
        ax.xaxis.set_label_coords(0.5, -1);
        ax.yaxis.set_label_coords(-0.2, 1)
        ax.set_xticks([0,xlim[-1]])
        ax.set_yticks(ylim)
        ax.set_xticklabels([int(i/3) for i in [0,xlim[-1]]], fontsize=6)
        #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
        ax.set_yticklabels(ylim, fontsize=8)


    else:
        silent_ax(ax)

            
# =============================================================================
# plotting time series for shuffle
# =============================================================================
# =============================================================================
# =============================================================================

save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/data/'
save_name = 'shuffle_0'

ylim = [-0.2,1.1]
xlim = [-15,330]
long=1;Cols=5


len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))



k = 0
label = 'single cell hyp.'

Rows = len_cells[k] // Cols;
if len_cells[k] % Cols != 0:
    Rows = Rows + 1
    inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_main[3], wspace=0, hspace=0.0)

color = sns.color_palette()[3]
df_forecast = download_data(save_folder+save_name+'.pkl')
df = df_forecast
cells = df.index.get_level_values(0).unique()


for zo,cell in enumerate(cells):
    ax = plt.subplot(inner_gs[zo])
    ax.set_ylim(ylim);
    ax.set_xlim(xlim)
    data = df.query('cell ==' + str(cell))

    if zo==0: ax.text(0.0, 1.4, 'HePM', ha='left', va='center', transform=ax.transAxes, fontsize=10,color = color)

    if cells[-1] >= cell:
        ax.plot(data.exponential_IPI.values ,color=color, linewidth=0.4,label=label)
        ax.plot(data.amp_peaks.values,'o',color = color,markersize=1)
        #ax.text(0.05, 0.2, str(data.amp_peaks.count()), ha='center', va='center', transform=ax.transAxes, fontsize=6,color = color)
        ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
    ax.grid(False);
    ax.tick_params(labelsize=8, direction='out', pad=1, length=2)

    
    if zo==((Rows-1)*Cols):
        ax.set_ylabel('amplitude (a.u)', fontsize=8);
        ax.set_xlabel('time (min)', fontsize=8)
        ax.xaxis.set_label_coords(0.5, -1);
        ax.yaxis.set_label_coords(-0.2, 1)
        ax.set_xticks([0,xlim[-1]])
        ax.set_yticks(ylim)
        ax.set_xticklabels([int(i/3) for i in [0,xlim[-1]]], fontsize=8)
        #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
        ax.set_yticklabels(ylim, fontsize=8)


    else:
        silent_ax(ax)


#aca van las celulas posta
ylim = [41500,63000]
k = 0
green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
Rows = len_cells[k] // Cols;
if len_cells[k] % Cols != 0:
    Rows = Rows + 1
inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_main[1], wspace=0, hspace=0.0)

df = conc_data[c]
cells = df.index.get_level_values(0).unique()
    

for zo,cell in enumerate(cells):
    ax = plt.subplot(inner_gs[zo])
    ax.set_ylim(ylim);
    ax.set_xlim(xlim)
    data = df.query('cell ==' + str(cell))

    if zo==0: ax.text(0.0, 1.4, 'serum + LIF', ha='left', va='center', transform=ax.transAxes, fontsize=10,color = green)

    if cells[-1] >= cell:
        ax.plot(data.sm_MEAN_INTENSITY.values ,color=green, linewidth=0.4,label=label)
        ax.plot(data.max_.values,'o',color = green,markersize=1)
        ax.text(0.05, 0.2, str(data.max_.count()), ha='center', va='center', transform=ax.transAxes, fontsize=6,color = green)
        ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
    ax.grid(False);
    ax.tick_params(labelsize=8, direction='out', pad=1, length=2)

    
    if zo==((Rows-1)*Cols):
        ax.set_ylabel('KTR signal (a.u)', fontsize=10);
        ax.set_xlabel('time (min)', fontsize=10)
        ax.xaxis.set_label_coords(0.5, -1);
        ax.yaxis.set_label_coords(-0.3, 1)
        ax.set_xticks([0,xlim[-1]])
        ax.set_yticks(ylim)
        ax.set_xticklabels([int(i/3) for i in [0,xlim[-1]]], fontsize=8)
        ax.set_yticklabels([int(y/100) for y in ylim], fontsize=8)


    else:   
        silent_ax(ax)
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/NEW_Fig2Supp3.pdf', format='pdf')

