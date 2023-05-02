import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 


from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)


import pickle
import sys


sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
colors = sns.color_palette('muted')
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
%matplotlib inline
matplotlib.rcParams['savefig.transparent'] = True



def set_scale_bars(ax, x_bar=None, y_bar=None, xunits='', yunits='', x_scale=1, y_scale=1, round_x=False,round_y=False):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if x_bar is not None:
        x0, xf = x_bar[0], x_bar[1]
        ax.spines['bottom'].set_bounds(x0, xf)
        ax.set_xticks([(xf + x0) / 2])
        scale_x = (xf - x0) // x_scale
        #if round_x:       ax.set_xticklabels([str(int(scale_x)) + xunits])#, rotation='vertical', va='center')
        #else:         ax.set_xticklabels([str(scale_x) + xunits])#, rotation='vertical', va='center')
        ax.set_xticklabels([])
        ax.xaxis.set_tick_params(length=0)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
    if y_bar is not None:
        y0, yf = y_bar[0], y_bar[1]
        ax.spines['left'].set_bounds(y0, yf)
        ax.set_yticks([(yf + y0) / 2])
        scale = (yf - y0) // y_scale
        if round_y:
            scale = int(scale)
        #ax.set_yticklabels([str(scale) + yunits], rotation='vertical', va='center')
        ax.set_yticklabels([])
        ax.yaxis.set_tick_params(length=0)
    else:
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
        
        
#%%
#Time correlation 3: con el dt miti miti

def mixed_dt(df):
    values_dt_mixed = []
    for cells, data in df.groupby(level='cell'):
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
    return(values_dt_mixed)


def is_consecutive(df):
    # discrimina los intervalos de consecutividad vs no consecutividad
    y = df['IPI'].dropna().values-mixed_dt(df)
    x = mixed_dt(df)
    consecutive = np.zeros(len(x))
    for i in range(len(x)):
        if x[i]*0.5 >= y[i]:
            consecutive[i] = 1
    return(consecutive)


###########################################################################################################################################
########################aca cargamos los datos ##########################################################################################
###########################################################################################################################################
###########################################################################################################################################

save_folder
df_forecast


#%% EL GRIDSPECT

#Figura principal
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.0,wspace=0)

#3 filas en toda la hoja
gs_main_rows = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_main[0], wspace=0.0, hspace=0.7)

#%primer fila, tres columnas
gs_row0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main_rows[0], wspace=0.5, hspace=0.5, width_ratios= [1.30,1.3,0.45, 0.7])

#%segunda fila, cuatro columnas
gs_row1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_main_rows[1], wspace=0.0, hspace=0.0)



#%%
# =============================================================================
# Selected traces - ~ 10 per condition; 
# =============================================================================
ylim = [-0.2,1.1]
xlim = [0,330]

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
color=colors[k]

inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc[c]), ncols=1, subplot_spec=gs_row0[0], wspace=0.1, hspace=0.0)

for w,z in enumerate(cells_conc[c]):
    ax = plt.subplot(inner_gs0[w])

    ax.plot(df_forecast.loc[cells_conc[c][w]]['exponential_IPI'].values, color=colors[k], linewidth=0.7);
    handles[0], = ax.plot(df_forecast.loc[cells_conc[c][w]]['amp_peaks'].values,linestyle="None", marker = "." , color='blue', markersize=4, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
    handles[1], = ax.plot(df_forecast.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=3, alpha=1,label = 'Peak minima')


    ax.set_xlim([0, M]);
    ax.set_ylim(ylim)

    if  (w == len(cells_conc[c])-1) :
       
       set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
       ax.set_ylabel('amplitude \n 1 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
       ax.xaxis.set_label_coords(0.2, -0.15);
       ax.yaxis.set_label_coords(-0.05,0.6);
    else:
        set_scale_bars(ax)

#    if w == 0:  
        #ax.text(0.2, 1.3, conc_labels[c], ha='center', va='center', transform=ax.transAxes, fontsize=8,color = colors[k])
    if ((w in traze_label) and (c == 'an_WT_ESL') ):
        txt = l.pop(0)
        ax.text(0.01, 0.7, txt, ha='left', va='center', transform=ax.transAxes, fontsize=8,color = colors[k])


#################
#histogramas

inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_row1[0], wspace=0.6, hspace=0.0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col1='amp_peaks'#ax2
col3 = 'IPI' #ax1
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    ax1 = plt.subplot(inner_gs1[k,2]); ax3 = plt.subplot(inner_gs1[k,1]); 
    
    bins1 = ax1.hist(df_forecast[col3].dropna(),bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    bins3 = ax3.hist(df_forecast[col2].dropna(),bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )

#    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.035])
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.1]);
    
    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
    ax1.text(0.75, 0.8, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.6, 1.07,r'$Q$: '+str(np.round(df_forecast[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df_forecast[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df_forecast[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)


    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(0.75, 0.8, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.6, 1.07,r'$Q$: '+str(np.round(df_forecast[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df_forecast[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df_forecast[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)

    if k == 0: silent_ax(ax1);silent_ax(ax3);



set_scale(ax1,np.arange(0,152,30),[0,0.035]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.2);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,152,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.035]],fontsize=6)


set_scale(ax3,np.arange(0,100,15),[0,0.1]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
ax3.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.1]],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.2);ax3.yaxis.set_label_coords(-0.01,0.5)


for ax_ in [ax1,ax3]: ax_.tick_params(labelsize=6,direction='out', pad=1,length=2)



# =============================================================================
# correlation // consecutiveness
# =============================================================================


for k,c in enumerate(ipi_conc):

    consecutive_ = is_consecutive(df_forecast)
    total_consecutive     = consecutive_.sum()
    total_non_consecutive = len(consecutive_) - total_consecutive


    ax3 = plt.subplot(inner_gs1[3]); 
    ax3.plot(df_forecast['IPI'].dropna().values-mixed_dt(df_forecast),mixed_dt(df_forecast),'o',markersize=1.5,alpha=0.6,color = colors[k+1])
    ax3.set_xlim([-5, 60]);
    ax3.set_ylim( [-5, 60])
    outliers = [[i,j] for i,j in zip(df_forecast['IPI'].dropna().values-mixed_dt(df_forecast),mixed_dt(df_forecast))]
    outliers_sum = np.sum([1 for i,j in outliers if i > 60 or j > 60])
    print('outliers : ' + str(outliers_sum))
    print('total_data: ' + str(len(outliers)))
    
    set_scale(ax3,np.arange(0,220,30),np.arange(0,360,30)); 
        
    X = np.arange(0,150)
    #ax3.fill_between(X,X*0.5,0,color='darkgray',alpha=0.2,linewidth=0.0) 
    ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.8,linestyle='--') 

    x_label = 'Interpulse intervals - ' + 'silence interval between pulses ' + '\n (min)'
    x_label = 'joint duration (min)'

    ax3.set_xlabel('silence (min)',fontsize=8,rotation=0); #between pulses 
    ax3.set_ylabel(x_label,fontsize=8); 

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
    
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_row0[1], wspace=0.8, hspace=0.75)

for k,c in enumerate(ipi_conc):
    ax1 = plt.subplot(gs0_inner[k]);
    activity = df_forecast['dt_peaks'].groupby(level='cell').sum() / df_forecast['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df_forecast.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df_forecast.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=colors[k],alpha=0.8,linewidth=0.0)
    if k == 0: plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)
    ax1.set_xlim([-1,df_forecast.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xlabel( conc_labels[c],fontsize=8); 
    ax1.set_xticks([0,df_forecast.index.get_level_values(0)[-1]])
    ax1.set_yticks([0,50,100])
    ax1.set_yticklabels([0,0.5,1],fontsize=6)
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

    ax1.xaxis.set_label_coords(0.5,-0.06)
    
    only_pulsing_cells = '// # activity cells = ' + str( sum((activity == 100)*1))
    pulsing_cells = '// # pulsing cells = ' + str(len(activity) - sum((activity == 0)*1))
    only_silent_cells = '// # silent cells = ' + str( sum((silent == 100)*1))
    print(c,only_pulsing_cells + pulsing_cells + only_silent_cells)


ax1.set_ylabel('fraction of cell track' ,fontsize=8); 
ax1.yaxis.set_label_coords(-0.13,0.05)


# =============================================================================
# mean activity
# =============================================================================

ax1 = plt.subplot(gs_row0[2]);

activity = []; silent = []
activity_err = []; silent_err = []

for c in ipi_conc:
    activity_ = df_forecast['dt_peaks'].groupby(level='cell').sum() / df_forecast['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df_forecast.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df_forecast.index.get_level_values(0).unique())))

    
p1 = ax1.barh(np.arange(len(ipi_conc)),width = silent,xerr=silent_err,left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
p2 = ax1.barh(np.arange(len(ipi_conc)),width = activity,left=silent,xerr = activity_err,color=colors,alpha=0.8,linewidth=0.0,height=0.6)
#plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.95, 1), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)

labels = ['serum \n + LIF',' ']
plt.xlabel('fraction of cell track' ,fontsize=8); 
plt.xticks([0,50,100]);plt.xlim([0,100]); ax1.set_xticklabels([0,0.5,1],fontsize=6)
plt.yticks(np.arange(len(cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.xaxis.set_label_coords(0.45, -0.25);ax1.yaxis.set_label_coords(-0.20,0.50)
ax1.invert_yaxis()

# =============================================================================
# 
# =============================================================================
plt.savefig(save_folder+save_name+'Figure2_st.pdf', format='pdf')