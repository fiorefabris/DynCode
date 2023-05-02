import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 
import sys
import os

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot



import matplotlib
plt.rcdefaults()

matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.75

matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.major.width'] = 0.75
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['xtick.minor.width'] = 0.75


sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from consecutive_main import consecutive_non_cumulative,consecutive_cumulative, mean_consecutive_value,mask_low_variance#,first_element_box,second_element_box

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from exponential_dm_main import download_data

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')
from plotting_main import load_file,save_file,get_linear_detrend,set_scale_bars,joint_duration,is_consecutive


        
        #%%    DOWNLOAD WILD TYPE DATA
filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/2019_E8_2_WTonly.pkl'
conc = load_file(filename)

main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
conc_data = {}
key = ['an_WT_ESL_PD03','an_WT_ESL']
keywords = ['NC','WT']
TH=('187.17241379310343', '2068.9655172413795')
amp_th = 2068
slope_th = 187
for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)


#for key in list(conc_data.keys()):
#    dm =find_dm(conc_data[key])

def ax_fontsize(fontsize,ax,rotation=None):
    plt.setp(ax.get_xticklabels(), rotation=rotation, fontsize=fontsize)



#%% EL GRIDSPECT
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'

#Figura principal
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.5,wspace= 0.0)

#3 filas en toda la hoja superior
gs_main_rows = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_main[0], wspace=0.0, hspace=0.7)

#%primer fila, tres columnas
gs_row0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main_rows[0], wspace=0.5, hspace=0.5, width_ratios= [1.30,1.3,0.45, 0.7])

#%segunda fila, cuatro columnas
gs_row1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_main_rows[1], wspace=0.0, hspace=0.0)




#%%
# =============================================================================
# Selected traces - ~ 10 per condition; 
# =============================================================================
ylim=[43000,63000]

M = max([conc_data[key].max().FRAME for key in conc_data.keys()])
l = ['*','**','***','****','*****','*','**','***','****','*****']
traze_label = [0,1,2]
#cells_conc = {'an_WT_ESL':[19,3,31,18,23,35,48,65],'an_WT_ESL_PD03':[5,14,25,31,39,45,61,65]} 
cells_conc = {'an_WT_ESL':[40,24,52],'an_WT_ESL_PD03':[5,14,25,31,39,45,61,65]} 
conc_labels = {'an_WT_ESL_PD03':'MEKi','an_WT_ESL': 'serum + LIF'}
cell_number = ['A','B','C','D','E']

box=[];handles = [None]*2
colors =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))
colors = [colors[i] for i in [15,25]]


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
k = 0; c = 'an_WT_ESL'
df = conc_data[c];df.sort_index(inplace=True); color=colors[k]

inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc[c]), ncols=1, subplot_spec=gs_row0[0], wspace=0.1, hspace=0.0)

for w,z in enumerate(cells_conc[c]):
    ax = plt.subplot(inner_gs0[w])

    ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
    handles[0], = ax.plot(df.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
    handles[1], = ax.plot(df.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')


    ax.set_xlim([0, M]);
    ax.set_ylim(ylim)

    if  (w == len(cells_conc[c])-1) :
       
       set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
       ax.set_ylabel('KTR signal \n 200 a.u.',fontsize=10); ax.set_xlabel('20 minutes',fontsize=10)
       ax.xaxis.set_label_coords(0.2, -0.15);
       ax.yaxis.set_label_coords(-0.05,0.6);
    else:
        set_scale_bars(ax)

#    if w == 0:  
        #ax.text(0.2, 1.3, conc_labels[c], ha='center', va='center', transform=ax.transAxes, fontsize=10,color = colors[k])
    if ((w in traze_label) and (c == 'an_WT_ESL') ):
        txt = l.pop(0)
        ax.text(0.01, 0.7, txt, ha='left', va='center', transform=ax.transAxes, fontsize=10,color = colors[k])


#################
#histogramas

inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_row1[0], wspace=0.6, hspace=0.0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col1='amp_peaks'#ax2
col3 = 'IPI' #ax1
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(inner_gs1[k,2]);ax2 = plt.subplot(inner_gs1[k,0]); ax3 = plt.subplot(inner_gs1[k,1]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    bins3 = ax3.hist(df[col2],bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    bins2 = ax2.hist(df[col1],bins=20,range =(amp_th,13000), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )

    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax2.axvspan(0,amp_th, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.035])
    ax2.set_xlim([0,13000]);ax2.set_ylim([0,0.00045]);
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.1]);
    
    #mode_IPI = mode(df[col3])[0][0]*20/60; mode_IPI = 'mode: '+str(np.round(mode_IPI,2))+' \n'
    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
    #ax1.text(0.75, 0.75, mode_IPI + ' Q25: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+'\n Q50: '+str(np.round(df[col3].quantile(0.50)*20/60,2)) +'\n Q75: '
#                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=10)
    ax1.text(0.75, 0.8, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)


    #mode_amp = mode(df[col1])[0][0]; mode_amp = 'mode: '+str(np.round(mode_amp/100,1))+' \n'
    mode_amp =  (bins2[1][np.argmax(bins2[0])] + bins2[1][np.argmax(bins2[0])+1])/2; mode_amp = 'mode: '+str(np.round(mode_amp/100,1))+' \n'
#    ax2.text(0.75, 0.70, mode_amp + ' Q25: '+str(np.round(df[col1].quantile(0.25)/100,1))+'\n Q50: '+str(np.round(df[col1].quantile(0.50)/100,1)) +'\n Q75: '
#                                  +str(np.round(df[col1].quantile(0.75)/100,1)) , ha='center', va='center', transform=ax2.transAxes, fontsize=9)

    ax2.text(0.75, 0.8, mode_amp, ha='center', va='center', transform=ax2.transAxes, fontsize=6)
    ax2.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col1].quantile(0.25)/100,1))+' ; '+str(np.round(df[col1].quantile(0.50)/100,1))+' ; '
                                  +str(np.round(df[col1].quantile(0.75)/100,1)) , ha='center', va='center', transform=ax2.transAxes, fontsize=6)



    #mode_dt = mode(df[col2])[0][0]*20/60 ; mode_dt = 'mode: '+str(np.round(mode_dt,2))+' \n'
    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(0.75, 0.8, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.6, 1.07,r'$Q$: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)

#    ax3.text(0.75, 0.75,mode_dt+' Q25: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+'\n Q50: '+str(np.round(df[col2].quantile(0.50)*20/60,2)) +'\n Q75: '
#                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=10)
    

    if k == 0: silent_ax(ax1),silent_ax(ax2);silent_ax(ax3);



set_scale(ax1,np.arange(0,152,30),[0,0.035]); ax1.set_xlabel('IPI (min)',fontsize=10); ax1.set_ylabel('probability\ndensity (1/min)',fontsize=10,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.4);ax1.yaxis.set_label_coords(-0.1,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,152,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(x*60/20,2) for x in [0,0.035]],fontsize=6)


set_scale(ax2,[amp_th,13000],[0,0.00045]); ax2.set_xlabel('amplitude (a.u.)',fontsize=10); ax2.set_ylabel('probability\ndensity (1/a.u.)',fontsize=10,rotation = 90)
ax2.xaxis.set_label_coords(0.5, -0.4);ax2.yaxis.set_label_coords(-0.1,0.5)
ax2.set_xticklabels([amp_th/100,13000/100],fontsize=6)
ax2.set_yticklabels([np.round(x*100,3) for x in [0,0.00045]],fontsize=6)

set_scale(ax3,np.arange(0,100,15),[0,0.1]); ax3.set_xlabel('duration (min)',fontsize=10); ax3.set_ylabel('probability\ndensity (1/min)',fontsize=10,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.4);ax3.yaxis.set_label_coords(-0.1,0.5)
ax3.set_yticklabels([np.round(x*60/20,2) for x in [0,0.1]],fontsize=6)


for ax_ in [ax1,ax2,ax3]: ax_.tick_params(labelsize=6,direction='out', pad=1,length=2)



# =============================================================================
# correlation // consecutiveness
# =============================================================================


for k,c in enumerate(ipi_conc):
    df = conc_data[c]

    consecutive_ = is_consecutive(df)
    total_consecutive     = consecutive_.sum()
    total_non_consecutive = len(consecutive_) - total_consecutive


    ax3 = plt.subplot(inner_gs1[3]); 
    ax3.plot(df['IPI'].dropna().values-joint_duration(df),joint_duration(df),'o',markersize=1.5,alpha=0.6,color = colors[k+1])
    ax3.set_xlim([-5, 60]);
    ax3.set_ylim( [-5, 60])
    outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-joint_duration(df),joint_duration(df))]
    outliers_sum = np.sum([1 for i,j in outliers if i > 60 or j > 60])
    print('outliers : ' + str(outliers_sum))
    print('total_data: ' + str(len(outliers)))
    
    set_scale(ax3,np.arange(0,220,30),np.arange(0,360,30)); 
        
    X = np.arange(0,150)
    #ax3.fill_between(X,X*0.5,0,color='darkgray',alpha=0.2,linewidth=0.0) 
    ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.8,linestyle='--') 

    x_label = 'Interpulse intervals - ' + 'silence interval between pulses ' + '\n (min)'
    x_label = 'joint duration (min)'

    ax3.set_xlabel('silence (min)',fontsize=10,rotation=0); #between pulses 
    ax3.set_ylabel(x_label,fontsize=10); 

    ax3.xaxis.set_label_coords(0.50, -0.2);ax3.yaxis.set_label_coords(-0.30,0.32)
    ax3.set_xticklabels([int(x) for x in np.arange(0,220,30)*20/60],fontsize=6)
    ax3.set_yticklabels([int(x) for x in np.arange(0,360,30)*20/60],fontsize=6)
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

    ax3.set_aspect(1)
    text_ = 'total: ' + str(int(len(consecutive_))) + '\n '+ r'$\leq 0.5 x :$'+ str(np.round(100*total_consecutive/len(consecutive_),1)) + ' % \n '+ r'$> 0.5 x:$' + str(np.round(100*total_non_consecutive/len(consecutive_),1)) + ' % \n '+ r'$ \frac{\leq 0.5 x}{> 0.5 x} :$ ' + str(np.round(total_consecutive/total_non_consecutive,2))    
    text_ = 'consecutive : '+ str(np.round(100*total_consecutive/len(consecutive_),1))+ ' %'

    ax3.text(0.70, 1.07, text_, ha='center', va='center', transform=ax3.transAxes, fontsize=6,color = 'black')

# =============================================================================
# # activity
# =============================================================================
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_row0[1], wspace=0.8, hspace=0.75)

for k,c in enumerate(cells_conc):
    ax1 = plt.subplot(gs0_inner[k]);
    df = conc_data[c]
    activity = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=colors[k],alpha=0.8,linewidth=0.0)
    if k == 0: plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)
    ax1.set_xlim([-1,df.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xlabel( conc_labels[c],fontsize=10); 
    ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
    ax1.set_yticks([0,50,100])
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

    ax1.xaxis.set_label_coords(0.5,-0.06)
    
    only_pulsing_cells = '// # activity cells = ' + str( sum((activity == 100)*1))
    pulsing_cells = '// # pulsing cells = ' + str(len(activity) - sum((activity == 0)*1))
    only_silent_cells = '// # silent cells = ' + str( sum((silent == 100)*1))
    print(c,only_pulsing_cells + pulsing_cells + only_silent_cells)


ax1.set_ylabel('% of cell track' ,fontsize=10); 
ax1.yaxis.set_label_coords(-0.17,1.3)


# =============================================================================
# mean activity
# =============================================================================
plt.rc('axes.spines', top=False, bottom=True, left=False, right=False); 
ax1 = plt.subplot(gs_row0[2]);

activity = []; silent = []
activity_err = []; silent_err = []

for c in cells_conc:
    df = conc_data[c]
    activity_ = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))

    
p1 = ax1.barh(np.arange(len(cells_conc)),width = silent,xerr=silent_err,left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
p2 = ax1.barh(np.arange(len(cells_conc)),width = activity,left=silent,xerr = activity_err,color=colors,alpha=0.8,linewidth=0.0,height=0.6)
#plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.95, 1), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)
#plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 1.4), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)

labels = ['serum \n + LIF','MEKi']
plt.xlabel('% of cell track' ,fontsize=10); 
#plt.xlabel('[FGF4] (ng/ml)',fontsize=10); 
plt.xticks([0,50,100]);plt.xlim([0,100]);
plt.yticks(np.arange(len(cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.xaxis.set_label_coords(0.45, -0.25);ax1.yaxis.set_label_coords(-0.20,0.50)
ax1.invert_yaxis()

# =============================================================================
# =============================================================================
# ================================ el otro dataset ============================
# =============================================================================
# =============================================================================
 


#Segunda hoja principal
gs_main_rows_2 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[1], wspace=0.5, hspace=0.5)

#primera fila
#gs_row_2_0 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main_rows_2[0], wspace=0.5, hspace=0.5)#, width_ratios= [1.30,1.3,0.45, 0.7]

#segunda fila
#gs_row_2_1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main_rows_2[1], wspace=0.5, hspace=0.0)


c = 'an_WT_ESL'
save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/population_data/'
save_folder_4='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/'
hypothesis_list_df_folder = [False,save_folder_1,save_folder_4+'data/']
labels_hypothesis_list_df =[ 'Serum + LIF','HoPM','HePM']
hypothesis_list_df = [conc_data[c],download_data(save_folder_1 + 'population0.pkl'),download_data(save_folder_4+'data/shuffle_1.pkl')]

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
colors =  [green,sns.color_palette()[1],sns.color_palette()[3]]

#activity plot
gs_0_0 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=2, subplot_spec=gs_main_rows_2[0,:2], wspace=0.5, hspace=0.5)#, width_ratios= [1.30,1.3,0.45, 0.7]

ax1 = plt.subplot(gs_0_0[1,0])
ax2 = plt.subplot(gs_0_0[1,1])

ax_list = [None,ax1,ax2]
leg = []
for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]
    color = colors[i]
    ax = ax_list[i]
    
    activity = df_consecutive['dt_peaks'].groupby(level='cell').sum() / df_consecutive['FRAME'].groupby(level='cell').count() *  100   
    activity_index = np.argsort(activity.values)[::-1]
    activity = [activity[j] for j in activity_index]
    
    silent = np.ones(len(activity)) * 100 - activity
    
    if i == 0: 
        activity_experiment = activity
        silent_experiment = silent
    
    else:    
        p1 = ax.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=color,alpha=0.8,linewidth=0.0)
        ax.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),activity_experiment,bottom=silent_experiment,width=0.8,color=colors[0],alpha=0.3,linewidth=0.0)
    
        leg.append(plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1))
        
        ax.set_xlim([-1,df_consecutive.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
        #ax1.hlines(50,y,len(activity) - y,color = 'black')
    
        ax.set_xlabel( labels_hypothesis_list_df[i],fontsize=10); 
        ax.set_xticks([0,df_consecutive.index.get_level_values(0)[-1]])
        ax.set_yticks([0,50,100])
        ax.set_yticklabels([0,0.5,1],fontsize=8)
 
        ax.xaxis.set_label_coords(0.5,-0.06)
    
        ax.set_ylabel('% of cell track' ,fontsize=10); 
        ax.yaxis.set_label_coords(-0.13,0.05)

#ax1.add_artist(leg[0])
#%
#boxplots
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 


ax1 = plt.subplot(gs_main_rows_2[1,0])
ax2 = plt.subplot(gs_main_rows_2[1,1])
ax_list = [ax1,ax2]
colors =  [sns.color_palette()[1],sns.color_palette()[3]]

data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons/'

iso_con_population = download_data(data_folder + 'iso_con_population_2.pkl')
iso_con_shuffle = download_data(data_folder + 'iso_con_shuffle_2.pkl')



#population boxplot
for n,j in enumerate([iso_con_population,iso_con_shuffle]):
    ax = ax_list[n]
    X = [np.ones(len(j[i]))*(i+1) for i in range(0,len(j))]
    bp = ax.boxplot(j,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )
    ax.axhline(y = 1,color = green,linewidth=0.5,linestyle = 'dashed')

    for i,box_ in enumerate(bp['boxes']):
         box_.set( color=colors[n], linewidth=0.0,facecolor=colors[n],alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp['whiskers']):
        whisker.set(color=colors[n//2],linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp['caps']):
        cap.set(color=colors[n//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp['medians']):
        median.set(color=colors[n],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians

    for i in range(len(X)):
        xA = np.random.normal(0, 0.1, len(j[i])), 
        ax.scatter(xA+X[i],j[i], alpha=1,s = 0.2,color='black',edgecolors='black',linewidths=0.0)
           
    #cambiar las labels por las qze soi    
    ax.tick_params(axis='x', labelsize=8,length=2); 
    ax.tick_params(axis='y', labelsize=8,length=2)
    ax.set_xticklabels(['isolated','consecutive'],rotation = 0)
    #ax.set_xlabel('isolated pulses',fontsize=8)
    ax.set_ylabel('normalized \n counts',fontsize=10)
    ax.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
    #xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
    #ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
    ax.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax.set_ylim([0.0,2.0])

    ax.legend(fontsize=10, ncol=1, framealpha=0, fancybox=True)
    ax.tick_params(labelsize=6,direction='out', pad=1,length=2)

#consecutivenes plot

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
colors =  [green,sns.color_palette()[1],sns.color_palette()[3]]
ax2 = plt.subplot(gs_main_rows_2[:,2])


for i,df_consecutive_folder in enumerate(hypothesis_list_df_folder):
    label = labels_hypothesis_list_df[i]
    
    if df_consecutive_folder is False:
        df_consecutive = conc_data[c]
        

        consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
        box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
        
        norm = len(df_consecutive.index.get_level_values(0).unique())#box_plot_consecutive_cumulative[0]
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),[i/norm for i in box_plot_consecutive_cumulative], color=colors[i], linewidth=0.75, marker = "." , markersize=5, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.set_yscale('log')
        ax2.set_xlim([0, 13]); ax2.set_ylim([10e-3,15])
        ax2.set_ylabel('counts x trace',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses (100%) ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.2,0.5);
        ax2.set_xticks([1,4,7,10,13])
        

    else:
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []
        for df_file_name in os.listdir(df_consecutive_folder) :
        
            df_consecutive = download_data(df_consecutive_folder+df_file_name)
    
            consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
           
        norm = len(df_consecutive.index.get_level_values(0).unique()) # la cantidad de células en un trial
        box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
        #norm = box_plot_consecutive_cumulative_trials_mean[0]
        print(box_plot_consecutive_cumulative_trials_mean,norm)
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean/norm, color=colors[i], linewidth=0.75, marker = "." , markersize=5, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),mask_low_variance(box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std)/norm,(box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std)/norm,color =colors[i],alpha = 0.2)        
        ax2.set_yscale('log')
        ax2.set_xlim([0,10.5]); ax2.set_ylim([10e-3,15])
        ax2.set_ylabel('counts x trace',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.2,0.5);
        ax2.set_xticks([1,4,7,10])
        
ax2.legend(fontsize=8, ncol=1,loc=1, framealpha=0, fancybox=True)
ax2.tick_params(labelsize=8,direction='out', pad=1,length=2)


plt.savefig(save_folder+'figure2.pdf', format='pdf')

#%%
# =============================================================================
# #%% calcular consecutive boxplot
# 
# =============================================================================
c = 'an_WT_ESL'
save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/population_data/'
save_folder_4='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/'
hypothesis_list_df_folder = [False,save_folder_1,save_folder_4+'data/']
labels_hypothesis_list_df =[ 'Serum + LIF','HoPM','HePM']

#labels_hypothesis_list_df =[ 'Serum + LIF','HoPM']
#hypothesis_list_df_folder = [False,save_folder_1]


box_total = []
box_isolated = []
box_consecutive = []
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons/'

for i,df_consecutive_folder in enumerate(hypothesis_list_df_folder):
    label = labels_hypothesis_list_df[i]
    
    if df_consecutive_folder is False: #este es el dataset
        df_consecutive = conc_data[c]
        consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
        
        isolated_N = sum(consecutive_non_cumulative_obj.is_isolated_box())
        box_isolated.append(1)
        
        consecutive_N = sum(consecutive_non_cumulative_obj.count_consecutive_pulses_number())
        box_consecutive.append(1)
        
        total_N = sum(consecutive_non_cumulative_obj.get_number_of_pulses())
        box_total.append(1)
        
    
    else:
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []; box_plot_consecutive_cumulative_trials
        for df_file_name in os.listdir(df_consecutive_folder) :
            
            df_consecutive = download_data(df_consecutive_folder + df_file_name)
            consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
        
            isolated = sum(consecutive_non_cumulative_obj.is_isolated_box())
            box_isolated.append(isolated/isolated_N)
        
            consecutive = sum(consecutive_non_cumulative_obj.count_consecutive_pulses_number())
            box_consecutive.append(consecutive/consecutive_N)
        
            total = sum(consecutive_non_cumulative_obj.get_number_of_pulses())
            box_total.append(total/total_N)
 
                 
            
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
                    
        if label == 'HoPM':
            iso_con_population = [box_isolated,box_consecutive]
            save_file(iso_con_population,save_folder + 'iso_con_population_2.pkl')
        if label == 'HePM':
            iso_con_shuffle = [box_isolated,box_consecutive]
            save_file(iso_con_shuffle,save_folder + 'iso_con_shuffle_2.pkl')