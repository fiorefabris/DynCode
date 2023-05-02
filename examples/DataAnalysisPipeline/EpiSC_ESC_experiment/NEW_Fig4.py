import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 
import sys
import matplotlib
import os

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

plt.rcdefaults()
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.75

matplotlib.rcParams['xtick.major.size'] = 2
matplotlib.rcParams['xtick.major.width'] = 0.75
matplotlib.rcParams['xtick.minor.size'] = 2
matplotlib.rcParams['xtick.minor.width'] = 0.75
matplotlib.rcParams['ytick.major.width'] = 0.75
matplotlib.rcParams['ytick.major.width'] = 0.75


sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')
from plotting_main import load_file,load_file5,set_scale_bars, joint_duration, is_consecutive
sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')
from plotting_main import  save_file

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from consecutive_main import consecutive_non_cumulative,makeColours,consecutive_cumulative,mean_consecutive_value,filter_activity_cells, mask_low_variance


def replace_zeros(cell_population_list):
    outcome = []
    for i in cell_population_list:
        if i == 0:
            outcome.append(None)
        else:
            outcome.append(i)
    return(outcome)


#%%

fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=4, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.5,wspace= 0.0)

#%% EPISC

epi_main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/EpiSC/data/'
epi_conc_data = {}
epi_color = sns.color_palette()[5]
epi_key =[ 'an_EpiSC_FAX_PD03', 'an_EpiSC_FAX']
keywords = ['CN','WT']
epi_TH = ('33.1315210671255', '305.08474576271186')
epi_slope_th = 33
epi_amp_th = 305

for j,k in enumerate(keywords):
    filename = epi_main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    epi_conc_data[epi_key[j]] = aux[epi_TH]
    
if False: save_file(epi_conc_data,'/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/epi_conc_data.pkl')

#%%ESC

esc_main_folder ='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/'
esc_conc_data = {}
esc_color = sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
esc_key =[ 'an_ESC_FAX_PD03', 'an_ESC_FAX']
keywords = ['CN','WT']
esc_TH =('57.25643730066835', '406.77966101694915')
esc_slope_th = 57
esc_amp_th = 406

for j,k in enumerate(keywords):
    filename = esc_main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    esc_conc_data[esc_key[j]] = aux[esc_TH]

if False: save_file(esc_conc_data,'/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/esc_conc_data_.pkl')

#%%
# =============================================================================
# EpiSC - row 1; 
# =============================================================================
gs_main_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[0], wspace=0.5, hspace=0.0, width_ratios = [1,1.5,0.5]  )

# ==================================== traces =========================================
ylim = [61000,66000]
M = max([epi_conc_data[key].max().FRAME for key in epi_conc_data.keys()])
epi_cells_conc = {'an_EpiSC_FAX':[27,46,50],'an_EpiSC_FAX_PD03':[5,14,25]} 
epi_conc_labels = {'an_EpiSC_FAX_PD03':'MEKi','an_EpiSC_FAX': 'FAX'}
box=[];handles = [None]*2
c = 'an_EpiSC_FAX'
df = epi_conc_data[c];df.sort_index(inplace=True);
inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(epi_cells_conc[c]), ncols=1, subplot_spec=gs_main_row_1[0], wspace=0.1, hspace=0.0)

for w,z in enumerate(epi_cells_conc[c]):
    ax = plt.subplot(inner_gs0[w])

    ax.plot(df.loc[epi_cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=epi_color, linewidth=0.75);
    handles[0], = ax.plot(df.loc[epi_cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=1, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
    handles[1], = ax.plot(df.loc[epi_cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=1, alpha=1,label = 'Peak minima')
    ax.set_xlim([0, M]);
    ax.set_ylim(ylim)

    if  (w == len(epi_cells_conc[c])-1) :
       
       set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
       ax.set_ylabel('KTR signal \n 50 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
       ax.xaxis.set_label_coords(0.2, -0.15);
       ax.yaxis.set_label_coords(-0.05,0.6);
    else:
        set_scale_bars(ax)

# ==================================== activity  =========================================

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_main_row_1[1], wspace=0.8, hspace=0.75)

for k,c in enumerate(epi_cells_conc):
    ax1 = plt.subplot(gs0_inner[k]);
    df = epi_conc_data[c]
    activity = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=epi_color,alpha=0.8,linewidth=0.0)
    if k == 0: plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)
    ax1.set_xlim([-1,df.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xlabel( epi_conc_labels[c],fontsize=8); 
    ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
    ax1.set_yticks([0,50,100])
    ax1.tick_params(labelsize=8,direction='out', length=2)
    ax1.xaxis.set_label_coords(0.5,-0.06)
    
    only_pulsing_cells = '// # activity cells = ' + str( sum((activity == 100)*1))
    pulsing_cells = '// # pulsing cells = ' + str(len(activity) - sum((activity == 0)*1))
    only_silent_cells = '// # silent cells = ' + str( sum((silent == 100)*1))
    print(c,only_pulsing_cells + pulsing_cells + only_silent_cells)


ax1.set_ylabel('fraction of\ncell track' ,fontsize=8); 
ax1.yaxis.set_label_coords(-0.17,1.3)

# ==================================== mean activity  =========================================

plt.rc('axes.spines', top=False, bottom=True, left=False, right=False); 

ax1 = plt.subplot(gs_main_row_1[2]);
activity = []; silent = []
activity_err = []; silent_err = []

for c in epi_cells_conc:
    df = epi_conc_data[c]
    activity_ = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))

    
p1 = ax1.barh(np.arange(len(epi_cells_conc)),width = silent,xerr=silent_err,left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
p2 = ax1.barh(np.arange(len(epi_cells_conc)),width = activity,left=silent,xerr = activity_err,color=epi_color,alpha=0.8,linewidth=0.0,height=0.6)

labels = ['FAX','MEKi']
plt.xlabel('fraction of\ncell track' ,fontsize=8); 
plt.xticks([0,50,100]);plt.xlim([0,100]);
plt.yticks(np.arange(len(epi_cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=8,direction='out',length=2)
ax1.xaxis.set_label_coords(0.45, -0.25);ax1.yaxis.set_label_coords(-0.20,0.50)
ax1.invert_yaxis()

# =============================================================================
# ESC - row 2; 
# =============================================================================
gs_main_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[1], wspace=0.5, hspace=0.0, width_ratios = [1,1.5,0.5]  )
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 

# ==================================== traces =========================================
ylim = [58000,65000]

M = max([esc_conc_data[key].max().FRAME for key in esc_conc_data.keys()])
esc_cells_conc = {'an_ESC_FAX':[40,33,38],'an_ESC_FAX_PD03':[5,14,25]} 
esc_conc_labels = {'an_ESC_FAX':' FAX ','an_ESC_FAX_PD03': 'MEKi' }
cell_number = ['A','B','C','D','E']

box=[];handles = [None]*2
c = 'an_ESC_FAX'
df = esc_conc_data[c];df.sort_index(inplace=True)

inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(esc_cells_conc[c]), ncols=1, subplot_spec=gs_main_row_2[0], wspace=0.1, hspace=0.0)

for w,z in enumerate(esc_cells_conc[c]):
    ax = plt.subplot(inner_gs0[w])

    ax.plot(df.loc[esc_cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=esc_color, linewidth=0.75);
    handles[0], = ax.plot(df.loc[esc_cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=1, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
    handles[1], = ax.plot(df.loc[esc_cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=1, alpha=1,label = 'Peak minima')


    ax.set_xlim([0, M]);
    ax.set_ylim(ylim)

    if  (w == len(esc_cells_conc[c])-1) :
       
       set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
       ax.set_ylabel('KTR signal \n 70 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
       ax.xaxis.set_label_coords(0.2, -0.15);
       ax.yaxis.set_label_coords(-0.05,0.6);
    else:
        set_scale_bars(ax)

# ==================================== activity  =========================================

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_main_row_2[1], wspace=0.8, hspace=0.75)

for k,c in enumerate(esc_cells_conc):
    ax1 = plt.subplot(gs0_inner[k]);
    df = esc_conc_data[c]
    activity = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=esc_color,alpha=0.8,linewidth=0.0)
    if k == 0: plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)
    ax1.set_xlim([-1,df.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xlabel( esc_conc_labels[c],fontsize=8); 
    ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
    ax1.set_yticks([0,50,100])
    ax1.tick_params(labelsize=8,direction='out',length=2)

    ax1.xaxis.set_label_coords(0.5,-0.06)
    
    only_pulsing_cells = '// # activity cells = ' + str( sum((activity == 100)*1))
    pulsing_cells = '// # pulsing cells = ' + str(len(activity) - sum((activity == 0)*1))
    only_silent_cells = '// # silent cells = ' + str( sum((silent == 100)*1))
    print(c,only_pulsing_cells + pulsing_cells + only_silent_cells)


ax1.set_ylabel('fraction of\ncell track' ,fontsize=8); 
ax1.yaxis.set_label_coords(-0.17,1.3)


# ==================================== mean activity  =========================================

plt.rc('axes.spines', top=False, bottom=True, left=False, right=False); 

ax1 = plt.subplot(gs_main_row_2[2]);
activity = []; silent = []
activity_err = []; silent_err = []


for c in esc_cells_conc:
    df = esc_conc_data[c]
    activity_ = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))

    
p1 = ax1.barh(np.arange(len(epi_cells_conc)),width = silent,xerr=silent_err,left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
p2 = ax1.barh(np.arange(len(epi_cells_conc)),width = activity,left=silent,xerr = activity_err,color=esc_color,alpha=0.8,linewidth=0.0,height=0.6)

labels = ['FAX','MEKi']
plt.xlabel('fraction of\ncell track' ,fontsize=8); 
plt.xticks([0,50,100]);plt.xlim([0,100]);
plt.yticks(np.arange(len(epi_cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=8,direction='out',length=2)
ax1.xaxis.set_label_coords(0.45, -0.25);ax1.yaxis.set_label_coords(-0.20,0.50)
ax1.invert_yaxis()

# =============================================================================
# EpiSC - row 3; 
# =============================================================================

gs_main_row_3_4 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=4, subplot_spec=gs_main[2:4], wspace=0.5, hspace=0.5) #, width_ratios = [1,1.5,0.5]
gs_main_row_3 = [gs_main_row_3_4[0,0], gs_main_row_3_4[0,1], gs_main_row_3_4[0,2]]

#gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[2], wspace=0.5, hspace=0.0 )
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

# ==================================== histograms  =========================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'IPI' #ax1
col2='dt_peaks'#ax3
k = 0
c = 'an_EpiSC_FAX'

df = epi_conc_data[c]
ax1 = plt.subplot(gs_main_row_3[1]); ax3 = plt.subplot(gs_main_row_3[0]); 
    
bins1 = ax1.hist(df[col3],bins=np.arange(2,140,4),range=(2,None), density=1, facecolor=epi_color, alpha=0.8,linewidth=0.0 )
bins3 = ax3.hist(df[col2],bins=np.arange(2,24,2),range=(2,None), density=1, facecolor=epi_color, alpha=0.8,linewidth=0.0 )

ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

ax1.set_xlim([0,140]);ax1.set_ylim([0,0.1])
ax3.set_xlim([0,30]);ax3.set_ylim([0,0.20]);

    
mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
ax1.text(0.75, 0.8, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
ax1.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                              +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
ax3.text(0.75, 0.8, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
ax3.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
         +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)

set_scale(ax1,np.arange(0,140,30),[0,0.1]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.2);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,140,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(x,1) for x in [0,0.1*3]],fontsize=6)

set_scale(ax3,np.arange(0,36,6),[0,0.20]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,36,6)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.2);ax3.yaxis.set_label_coords(-0.01,0.5)
ax3.set_yticklabels([np.round(x,1) for x in [0,0.20*3]],fontsize=6)


for ax_ in [ax1,ax3]: ax_.tick_params(labelsize=8,direction='out', length=2)

# ==================================== correlation  =========================================



consecutive_ = is_consecutive(df)
total_consecutive     = consecutive_.sum()
total_non_consecutive = len(consecutive_) - total_consecutive


ax3 = plt.subplot(gs_main_row_3[2]); 
ax3.plot(df['IPI'].dropna().values-joint_duration(df),joint_duration(df),'o',markersize=1.5,alpha=0.6,color = epi_color)
ax3.set_xlim([-5, 60]);
ax3.set_ylim( [-5, 60])
outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-joint_duration(df),joint_duration(df))]
outliers_sum = np.sum([1 for i,j in outliers if i > 60 or j > 60])
print('epi outliers : ' + str(outliers_sum))
print('epi total_data: ' + str(len(outliers)))
    
set_scale(ax3,np.arange(0,220,30),np.arange(0,360,30)); 
        
X = np.arange(0,150)
ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.75,linestyle='--') 

x_label = 'Interpulse intervals - ' + 'silence interval between pulses ' + '\n (min)'
x_label = 'joint duration (min)'

ax3.set_xlabel('silence (min)',fontsize=8,rotation=0); #between pulses 
ax3.set_ylabel(x_label,fontsize=8); 

ax3.xaxis.set_label_coords(0.50, -0.2);ax3.yaxis.set_label_coords(-0.30,0.32)
ax3.set_xticklabels([int(x) for x in np.arange(0,220,30)*20/60],fontsize=6)
ax3.set_yticklabels([int(x) for x in np.arange(0,360,30)*20/60],fontsize=6)
ax3.tick_params(labelsize=8,direction='out',length=2)

ax3.set_aspect(1)
text_ = 'total: ' + str(int(len(consecutive_))) + '\n '+ r'$\leq 0.5 x :$'+ str(np.round(100*total_consecutive/len(consecutive_),1)) + ' % \n '+ r'$> 0.5 x:$' + str(np.round(100*total_non_consecutive/len(consecutive_),1)) + ' % \n '+ r'$ \frac{\leq 0.5 x}{> 0.5 x} :$ ' + str(np.round(total_consecutive/total_non_consecutive,2))    
text_ = 'consecutive : '+ str(np.round(100*total_consecutive/len(consecutive_),1))+ ' %'

ax3.text(0.70, 1.07, text_, ha='center', va='center', transform=ax3.transAxes, fontsize=6,color = 'black')

# ==================================== seq of consecutive pulses  =========================================
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 


# =============================================================================
# ESC - row 4; 
# =============================================================================
gs_main_row_4 = [gs_main_row_3_4[1,0], gs_main_row_3_4[1,1], gs_main_row_3_4[1,2]]

#gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[3], wspace=0.5, hspace=0.0 )
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
# ==================================== histograms  =========================================
col3 = 'IPI' #ax1
col2='dt_peaks'#ax3
c = 'an_ESC_FAX'
k = 0


df = esc_conc_data[c]
ax1 = plt.subplot(gs_main_row_4[1]);ax3 = plt.subplot(gs_main_row_4[0]); 
bins1 = ax1.hist(df[col3],bins=np.arange(2,140,4),range=(2,None), density=1, facecolor=esc_color, alpha=0.8,linewidth=0.0 )
bins3 = ax3.hist(df[col2],bins=np.arange(2,30,2),range=(2,None), density=1, facecolor=esc_color, alpha=0.8,linewidth=0.0 )

ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

ax1.set_xlim([0,140]);ax1.set_ylim([0,0.1])
ax3.set_xlim([0,30]);ax3.set_ylim([0,0.20]);
    
mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
ax1.text(0.75, 0.8, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
ax1.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                              +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)


mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
ax3.text(0.75, 0.8, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
ax3.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                              +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)

set_scale(ax1,np.arange(0,140,30),[0,0.1]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.2);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,140,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(x,1) for x in [0,0.1*3]],fontsize=6)

set_scale(ax3,np.arange(0,36,6),[0,0.20]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,36,6)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.2);ax3.yaxis.set_label_coords(-0.01,0.5)
ax3.set_yticklabels([np.round(x,1) for x in [0,0.20*3]],fontsize=6)


for ax_ in [ax1,ax3]: ax_.tick_params(labelsize=8,direction='out',length=2)

# ==================================== correlation  =========================================

#
consecutive_ = is_consecutive(df)
total_consecutive     = consecutive_.sum()
total_non_consecutive = len(consecutive_) - total_consecutive


ax3 = plt.subplot(gs_main_row_4[2]); 
ax3.plot(df['IPI'].dropna().values-joint_duration(df),joint_duration(df),'o',markersize=1.5,alpha=0.6,color = esc_color)
ax3.set_xlim([-5, 60]);
ax3.set_ylim( [-5, 60])
outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-joint_duration(df),joint_duration(df))]
outliers_sum = np.sum([1 for i,j in outliers if i > 60 or j > 60])
print('esc outliers : ' + str(outliers_sum))
print('esc total_data: ' + str(len(outliers)))
    
set_scale(ax3,np.arange(0,220,30),np.arange(0,360,30)); 
        
X = np.arange(0,150)
ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.75,linestyle='--') 

x_label = 'Interpulse intervals - ' + 'silence interval between pulses ' + '\n (min)'
x_label = 'joint duration (min)'

ax3.set_xlabel('silence (min)',fontsize=8,rotation=0); #between pulses 
ax3.set_ylabel(x_label,fontsize=8); 

ax3.xaxis.set_label_coords(0.50, -0.2);ax3.yaxis.set_label_coords(-0.30,0.32)
ax3.set_xticklabels([int(x) for x in np.arange(0,220,30)*20/60],fontsize=6)
ax3.set_yticklabels([int(x) for x in np.arange(0,360,30)*20/60],fontsize=6)
ax3.tick_params(labelsize=8,direction='out',length=2)

ax3.set_aspect(1)
text_ = 'total: ' + str(int(len(consecutive_))) + '\n '+ r'$\leq 0.5 x :$'+ str(np.round(100*total_consecutive/len(consecutive_),1)) + ' % \n '+ r'$> 0.5 x:$' + str(np.round(100*total_non_consecutive/len(consecutive_),1)) + ' % \n '+ r'$ \frac{\leq 0.5 x}{> 0.5 x} :$ ' + str(np.round(total_consecutive/total_non_consecutive,2))    
text_ = 'consecutive : '+ str(np.round(100*total_consecutive/len(consecutive_),1))+ ' %'
ax3.text(0.70, 1.07, text_, ha='center', va='center', transform=ax3.transAxes, fontsize=6,color = 'black')
    #%%

# =============================================================================
# consecutiveness
# =============================================================================
save_folder_data = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/'
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);


ax1 = plt.subplot(gs_main_row_3_4[0,3])

hypothesis_list_df = [load_file(save_folder_data + 'epi_conc_data.pkl')['an_EpiSC_FAX'],load_file(save_folder_data + 'esc_conc_data.pkl')['an_ESC_FAX']]
labels_hypothesis_list_df =[ 'EpiSC in FAX','ESC in FAX']
colors = [epi_color,esc_color]
shuffle_data_folder = {'EpiSC in FAX':save_folder_data+'epi_shuffle_data/' ,'ESC in FAX':save_folder_data+'esc_shuffle_data/'}

for i,df_consecutive in enumerate(hypothesis_list_df):
    ### experimental data
    label = labels_hypothesis_list_df[i]
    if label == 'EpiSC in FAX':
         alpha = 1
    else:
        alpha = 0.3
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    cell_population_list = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    norm = len(df_consecutive.index.get_level_values(0).unique())
    cell_population_list = replace_zeros([i/norm for i in cell_population_list])
    
    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.75, marker = "." , markersize=4, alpha=alpha,color =colors[i],label=labels_hypothesis_list_df[i])
    cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []
    
    if label == 'EpiSC in FAX':
        for df_file_name in os.listdir(shuffle_data_folder[label]) :
            
                df_consecutive = load_file5(shuffle_data_folder[label]+df_file_name)    
                consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
                box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
                box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
               
        norm = len(df_consecutive.index.get_level_values(0).unique()) # la cantidad de células en un trial
        print('norm: ',norm)
        box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
        ax1.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean/norm, color=colors[i], linewidth=0.75, marker = "^", markersize=4, alpha=1,label= 'HePM '+label)
        ax1.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),mask_low_variance(box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std)/norm,(box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std)/norm,color =colors[i],alpha = 0.2)        

ax1.set_xlim([0, 10.5]);
ax1.set_yscale('log')
ax1.set_ylim([10e-3,15])
ax1.xaxis.set_label_coords(0.5, -0.09);ax1.yaxis.set_label_coords(-0.1,0.5);
ax1.set_xticks([1,4,7,10])
ax1.set_ylabel('counts x trace',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
ax1.legend(fontsize=5, ncol=1,framealpha=0, fancybox=True)


# =============================================================================
# =============================================================================

save_folder_data = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/'
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);


ax1 = plt.subplot(gs_main_row_3_4[1,3])

hypothesis_list_df = [load_file(save_folder_data + 'epi_conc_data.pkl')['an_EpiSC_FAX'],load_file(save_folder_data + 'esc_conc_data.pkl')['an_ESC_FAX']]
labels_hypothesis_list_df =[ 'EpiSC in FAX','ESC in FAX']
colors = [epi_color,esc_color]
shuffle_data_folder = {'EpiSC in FAX':save_folder_data+'epi_shuffle_data/' ,'ESC in FAX':save_folder_data+'esc_shuffle_data/'}

for i,df_consecutive in enumerate(hypothesis_list_df):
    ### experimental data
    label = labels_hypothesis_list_df[i]
    if label == 'ESC in FAX':
        alpha = 1
        consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
        cell_population_list = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
        norm = len(df_consecutive.index.get_level_values(0).unique())
        cell_population_list = replace_zeros([i/norm for i in cell_population_list])
        
        ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.75, marker = "." , markersize=4, alpha=alpha,color =colors[i],label=labels_hypothesis_list_df[i])
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []
        
        if True:
            for df_file_name in os.listdir(shuffle_data_folder[label]) :
                
                    df_consecutive = load_file5(shuffle_data_folder[label]+df_file_name)    
                    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
                    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
                    box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
                   
            norm = len(df_consecutive.index.get_level_values(0).unique()) # la cantidad de células en un trial
            print('norm: ',norm)
            box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
            ax1.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean/norm, color=colors[i], linewidth=0.75, marker = "^", markersize=3, alpha=1,label= 'HePM '+label)
            ax1.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),mask_low_variance(box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std)/norm,(box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std)/norm,color =colors[i],alpha = 0.2)        
            

ax1.set_xlim([0, 10.5]);
ax1.set_yscale('log')
ax1.set_ylim([10e-3,15])
ax1.xaxis.set_label_coords(0.5, -0.09);ax1.yaxis.set_label_coords(-0.1,0.5);
ax1.set_xticks([1,4,7,10])
ax1.set_ylabel('counts x trace',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
ax1.legend(fontsize=5, ncol=1,framealpha=0, fancybox=True)
plt.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/Fig. 4.pdf', format='pdf')

#%%
total_pulses_epi = epi_conc_data['an_EpiSC_FAX'].max_.count()
total_pulses_epi_meki = epi_conc_data['an_EpiSC_FAX_PD03'].max_.count()
total_pulses_esc =  esc_conc_data['an_ESC_FAX'].max_.count()

epi_pairs_pulses = epi_conc_data['an_EpiSC_FAX'].IPI.count()
esc_pairs_pulses = esc_conc_data['an_ESC_FAX'].IPI.count()
#%%
save_folder_data = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/'
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 

ax1 = plt.subplot(gs_main_row_3_4[0:2,3])

hypothesis_list_df = [load_file(save_folder_data + 'epi_conc_data.pkl')['an_EpiSC_FAX'],load_file(save_folder_data + 'esc_conc_data.pkl')['an_ESC_FAX']]
labels_hypothesis_list_df =[ 'EpiSC in FAX','ESC in FAX']
colors = [epi_color,esc_color]
shuffle_data_folder = {'EpiSC in FAX':save_folder_data+'epi_shuffle_data/' ,'ESC in FAX':save_folder_data+'esc_shuffle_data/'}

for i,df_consecutive in enumerate(hypothesis_list_df):
    ### experimental data
    label = labels_hypothesis_list_df[i]

    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    cell_population_list = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    norm = len(df_consecutive.index.get_level_values(0).unique())
    cell_population_list = replace_zeros([i/norm for i in cell_population_list])
    
    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.75, marker = "." , markersize=2, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])

    
    cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []

    for df_file_name in os.listdir(shuffle_data_folder[label]) :
        
            df_consecutive = load_file5(shuffle_data_folder[label]+df_file_name)    
            consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
           
    norm = len(df_consecutive.index.get_level_values(0).unique()) # la cantidad de células en un trial
    print('norm: ',norm)
    box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
    ax1.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean/norm, color=colors[i], linewidth=0.75, marker = "^", markersize=2, alpha=1,label= 'HePM '+label)
    ax1.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),(box_plot_consecutive_cumulative_trials_mean-box_plot_consecutive_cumulative_trials_std)/norm,(box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std)/norm,color =colors[i],alpha = 0.2)        

ax1.set_xlim([0, 7.5]);
#ax1.set_yscale('log')
#ax1.set_ylim([0.01,10])
ax1.xaxis.set_label_coords(0.5, -0.09);ax1.yaxis.set_label_coords(-0.1,0.5);
ax1.set_xticks([1,3,5,7])
ax1.set_ylabel('counts x trace',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
        
ax1.legend(fontsize=5, ncol=1,framealpha=0, fancybox=True)


ax3.text(0.70, 1.07, text_, ha='center', va='center', transform=ax3.transAxes, fontsize=6,color = 'black')
