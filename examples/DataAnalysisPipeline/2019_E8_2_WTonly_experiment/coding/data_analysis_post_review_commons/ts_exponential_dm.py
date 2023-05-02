import sys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from DataAnalysis.Preprocesing.PlotData import silent_ax

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')

from exponential_dm_main import download_data,get_linear_fit,get_population_time_series,get_single_cell_time_series,find_dm,get_single_cell_time_series_moving_boundaries,get_single_cells_pulses_shuffle
#%%

# =============================================================================
# calculating the input variables for population simulation
# =============================================================================

find_dm(conc_data['an_WT_ESL'])

bins1 = plt.hist(conc_data['an_WT_ESL']['dm_'],bins=np.arange(2,200,6), density=1);

popt,pcov = get_linear_fit(bins1)
dm_max = 0

#popt,pcov = get_linear_fit(bins1,20)
#dm_max = 20 #minutes

lambda_ = popt[0]   
exponential = np.random.exponential(1/lambda_,len(conc_data['an_WT_ESL']['dt_peaks']))


save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_04_2022/'
save_name = 'population'

# =============================================================================
# calculating the time series population simulation
# =============================================================================
ipi_conc = ['an_WT_ESL']
save_folder_save_name = save_folder + save_name
get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)

#%%

# =============================================================================
# plotting funtions for population simulation
# =============================================================================

# =============================================================================
# =============================================================================
# HYPOTHESIS PLOT


plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}
colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_0 = gridspec.GridSpec(nrows=4, ncols=3, figure=fig, wspace=0.3, hspace=0.3)

# =============================================================================
# Pulse rate
# =============================================================================

box = []
for k,c in enumerate(conc_labels):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    box.append((df['amp_peaks'].groupby(level='cell').count() / df['amp_peaks'].groupby(level='cell').size()).to_list()) 
    
X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
ax2 = plt.subplot(gs_0[0,0]);

for i in range(len(X)):
    xA = np.random.normal(0, 0.1, len(box[i])), 
    ax2.scatter(xA+X[i],box[i],  alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
bp = ax2.boxplot(box,whis=[5, 95],vert=True,patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp['fliers']):
    flyer.set(markeredgecolor=colors[i])## change color and linewidth of the medians


ax2.tick_params(axis='x', labelsize=6,length=2); ax2.tick_params(axis='y', labelsize=6,length=2)
ax2.set_xticklabels([conc_labels[c] for c in conc_labels],rotation = 0,fontsize=8)
label_peaks = r'$\rm{\left(min^{-1}\right)}$'
ax2.set_ylabel('pulse rate ' + label_peaks,fontsize=8)
ax2.yaxis.set_label_coords(-0.01,0.55)

ax2.text(0.015, 0.85, 'mean (m): ' + str(np.round(np.mean([i/(20/60) for i in box[0]]),2)) + '\n std : ' + str(np.round(np.std([i/(20/60) for i in box[0]]),2)), ha='left', va='baseline', transform=ax2.transAxes, fontsize=6)
ax2.text(0.015, 0.70, 'mean (f): ' + str(np.round(np.mean([i/1 for i in box[0]]),2)) + '\n std : ' + str(np.round(np.std([i/1 for i in box[0]]),2)), ha='left', va='baseline', transform=ax2.transAxes, fontsize=6)

ax2.text(0.55, 0.85, 'mean (m): ' + str(np.round(np.mean([i/(20/60) for i in box[1]]),2))+ '\n std : ' + str(np.round(np.std([i/(20/60) for i in box[1]]),2)), ha='left', va='baseline', transform=ax2.transAxes, fontsize=6)
ax2.text(0.55, 0.70, 'mean (f): ' + str(np.round(np.mean([i/1 for i in box[1]]),2))+ '\n std : ' + str(np.round(np.std([i/1 for i in box[1]]),2)), ha='left', va='baseline', transform=ax2.transAxes, fontsize=6)

yticks = ax2.get_yticks(); yticks[-1] = 0.06 
minutes_ylim = yticks[-1]/(20/60)
ax2.set_yticklabels([0,minutes_ylim],fontsize=6)
ax2.set_ylim([-0.005,yticks[-1]]); ax2.set_yticks([0,yticks[-1]])
ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
#  dt duration histogram
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_0[0,1]); 
    
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
  

set_scale(ax3,np.arange(0,100,15),[0,0.1]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
ax3.set_yticklabels([np.round(y/(20/60),2)  for y in [0,0.1]],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.0,0.5)


ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# IPI
# =============================================================================
        


plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'IPI'
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[1,0]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.035])

    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode (m): '+str(np.round(mode_IPI*20/60,2))+' \n'
    ax1.text(0.65, 0.90, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode (f): '+str(np.round(mode_IPI,2))+' \n'
    ax1.text(0.65, 0.70, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(df[col3].quantile(0.25),2))+' ; '+str(np.round(df[col3].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75),2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    ax1.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(df[col3]),2)) +'\n std : '+str(np.round(np.std(df[col3]),2)), ha='center', va='center', transform=ax1.transAxes, fontsize=6)

set_scale(ax1,np.arange(0,152,30),[0,0.035]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,152,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.035]],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# dm
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'dm_'
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[1,1]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,200,6), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.axvspan(0,dm_max*60/20, color='gray',alpha = 0.4,linewidth=0.0)
    ax1.plot(lambda_*np.exp(-lambda_*np.arange(200)),color = 'blue',alpha = 1,linewidth=0.5)

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

set_scale(ax1,np.arange(2,200,60),[0,0.045]); ax1.set_xlabel('time between pulses (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.045]],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# exponential distrobution
# 
# =============================================================================
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_0[1,2]); 
    

    bins3 = ax3.hist(exponential,bins=np.arange(0,200,6), density=1, facecolor='blue', alpha=0.5,linewidth=0.0 )
    ax3.axvspan(0,np.round(df[col2].quantile(0.50),2), color='red',alpha = 0.2,linewidth=0.0)
    ax3.plot(lambda_*np.exp(-lambda_*np.arange(400)),color = 'blue',alpha = 1,linewidth=0.5)
    ax3.set_xlim([0,200]);ax3.set_ylim([0,0.03]);

    mode_exp =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_exp = 'mode (m): '+str(np.round(mode_exp*20/60,2))+' \n'
    ax3.text(0.65, 0.90, mode_exp, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(np.quantile(exponential,0.25)*20/60,2))+' ; '+str(np.round(np.quantile(exponential,0.50)*20/60,2))+' ; '
                                  +str(np.round(np.quantile(exponential,0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
 

    mode_exp =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_exp = 'mode (f): '+str(np.round(mode_exp,2))+' \n'
    ax3.text(0.65, 0.7, mode_exp, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(np.quantile(exponential,0.25),2))+' ; '+str(np.round(np.quantile(exponential,0.50),2))+' ; '
                                  +str(np.round(np.quantile(exponential,0.75),2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
 

    ax3.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(exponential),2)) +'\n std : '+str(np.round(np.std(exponential),2)), ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    


set_scale(ax3,np.arange(0,200,50),[0,0.01,0.02,0.03]); ax3.set_xlabel(r'$T_i$ (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,200,50)*20/60],fontsize=6)
ax3.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.01,0.02,0.03]],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.01,0.5)
ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)


# =============================================================================
# #log plots
#  IPI log plot
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'IPI'
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[2,0]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_yscale('log')
    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.035])


set_scale(ax1,np.arange(0,152,30),ax1.get_yticks()); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.15,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,152,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),4) for y in ax1.get_yticks()],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# dm
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'dm_'
ipi_conc = ['an_WT_ESL']


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[2,1]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,200,6), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.plot(lambda_*np.exp(-lambda_*np.arange(200)),color = 'blue',alpha = 1,linewidth=0.5)
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.axvspan(0,dm_max*60/20, color='gray',alpha = 0.4,linewidth=0.0)
    ax1.set_yscale('log')
    ax1.set_xlim([0,200]);ax1.set_ylim([0,0.045])


set_scale(ax1,np.arange(2,200,60),ax1.get_yticks()); ax1.set_xlabel('time between pulses (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.15,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),4) for y in ax1.get_yticks()],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# exponential distrobution
# 
# =============================================================================
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
col2='dt_peaks'#ax3
ipi_conc = ['an_WT_ESL']

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_0[2,2]); 
    
    bins3 = ax3.hist(exponential,bins=np.arange(0,200,6), density=1, facecolor='blue', alpha=0.5,linewidth=0.0 )
    ax3.plot(lambda_*np.exp(-lambda_*np.arange(200)),color = 'blue',alpha = 1,linewidth=0.5)
    ax3.axvspan(0,np.round(df[col2].quantile(0.50),2), color='red',alpha = 0.2,linewidth=0.0)
    ax3.set_yscale('log')
    ax3.set_xlim([0,200]);ax3.set_ylim([0,0.03]);

set_scale(ax3,np.arange(0,200,50),ax3.get_yticks()); ax3.set_xlabel(r'$T_i$ (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,200,50)*20/60],fontsize=6)
ax3.set_yticklabels([np.round(y/(20/60),4) for y in ax3.get_yticks()],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.15,0.5)
ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder+save_name + '_hypothesis'+ '.pdf', format='pdf')



#%%
# =============================================================================
# =============================================================================
# TIME SERIES PLOT

df_forecast = download_data(save_folder_save_name + '.pkl')

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5
 
colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);


for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = df_forecast
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.exponential_IPI.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.amp_peaks.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)

            
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()
#%%

# =============================================================================
# calculating the time series single_cell simulation
# =============================================================================

''' Time series forecasting with individual cell parametrization module. 

We want to generate single cells from single cells statistics by generating an exponential
dristibution of silent intervals from single cells mean value of silent intervals. 
'''
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/'
save_name= 'single_cell_silences_new_cdc'
ipi_conc = ['an_WT_ESL']
get_single_cell_time_series(conc_data,ipi_conc,save_folder+save_name,True)
    
    
conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')
#%%
# =============================================================================
# plotting funtions for single cell simulation
# =============================================================================

# =============================================================================
# =============================================================================
# TIME SERIES PLOT

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);



for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = download_data(save_folder+save_name+'.pkl')
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.exponential_IPI.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.amp_peaks.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)
     
#aca van las celulas posta
ylim = [45000,75000]

for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k+1], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[c]
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.sm_MEAN_INTENSITY.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.max_.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.3, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)



#aca van los histogramas
ylim = [0,0.05]
xlim = [0,200]
for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[2], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[c]
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            lambda_ = np.mean(data.dm_.dropna().values)[0]
            if len(data.dt_peaks.dropna().values)> 0: 
                exponential_cell = np.random.exponential(lambda_,1000)
            else: 
                exponential_cell = np.array([])

            if exponential_cell.size > 0: 
                ax.hist(exponential_cell ,bins=30,density = True,color=color, linewidth=0.4,label=label)
                ax.plot(1/lambda_*np.exp(-1/lambda_*np.arange(600)),color = 'blue',alpha = 1,linewidth=0.5)

            ax.text(0.95, 0.2, str(cell) , ha='center', va='center', transform=ax.transAxes, fontsize=6)
            ax.text(0.7, 0.6, r'$T_i = $'+str(np.round(lambda_/3,2)) + 'min' , ha='center', va='center', transform=ax.transAxes, fontsize=4)

        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('counts ', fontsize=8);
            ax.set_xlabel(r'$T_i (min)$', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i*20/60 for i in [0,xlim[-1]]], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)


#aca va el histograma de lambda
ylim = [0,0.02]
xlim = [0,300]

inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_row[3:-1], wspace=0, hspace=0.0)
color = 'blue';
plt.rcdefaults()


ax = plt.subplot(inner_gs[0])
ax.set_ylim(ylim);
ax.set_xlim(xlim)
  


lambda_box = []
for c in ipi_conc:    
    for cell,data in df.groupby(level= 'cell'):
        lambda_ = np.mean(data.dm_.dropna().values)
        lambda_box.append(lambda_)      
    
ax.hist(lambda_box  ,bins=20,density = True,color=color, linewidth=0, alpha =0.8, label=label)
ax.grid(False);
ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
ax.set_ylabel('probability (1/min)', fontsize=8);
ax.set_xlabel(r'$T_i (min)$', fontsize=8)
ax.xaxis.set_label_coords(0.5, -0.1);
ax.yaxis.set_label_coords(-0.05, 0.5)
ax.set_xticks([xlim[0],25,40,xlim[1]])
ax.set_yticks(ylim)
ax.set_xticklabels([np.round(i*20/60,2) for i in [xlim[0],25,40,xlim[1]]], fontsize=6)
ax.set_yticklabels([i * 3 for i in ylim], fontsize=6)


         
            
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()

#%%

# =============================================================================
# calculating the time series single_cell simulation with moving boundaries
# =============================================================================

''' Time series forecasting with individual cell parametrization module. 

We want to generate single cells from single cells statistics by generating an exponential
dristibution of silent intervals from single cells mean value of silent intervals. 
'''
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/moving_bc/'
save_name= 'single_cell_silences_moving_bc'
ipi_conc = ['an_WT_ESL']
dm_cells,lambda_box = get_single_cell_time_series_moving_boundaries(conc_data,ipi_conc,save_folder+save_name)
    
    
conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')


#%%
# =============================================================================
# plotting funtions for single cell simulation
# =============================================================================

# =============================================================================
# =============================================================================
# TIME SERIES PLOT

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);



for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = download_data(save_folder+save_name+'.pkl')
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.exponential_IPI.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.amp_peaks.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)
     
#aca van las celulas posta
ylim = [45000,75000]

for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k+1], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[c]
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.sm_MEAN_INTENSITY.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.max_.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.3, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)



#aca van los histogramas
ylim = [0,0.05]
xlim = [0,200]
for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[2], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[c]
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            lambda_ = lambda_box[zo]
            if len(data.dt_peaks.dropna().values)> 0: 
                exponential_cell = np.random.exponential(lambda_,1000)
            else: 
                exponential_cell = np.array([])

            if exponential_cell.size > 0: 
                ax.hist(exponential_cell[np.isfinite(exponential_cell)] ,bins=30,density = True,color=color, linewidth=0.4,label=label)
                ax.plot(1/lambda_*np.exp(-1/lambda_*np.arange(600)),color = 'blue',alpha = 1,linewidth=0.5)

            ax.text(0.95, 0.2, str(cell) , ha='center', va='center', transform=ax.transAxes, fontsize=6)
            ax.text(0.7, 0.6, r'$T_i = $'+str(np.round(lambda_/3,2)) + 'min' , ha='center', va='center', transform=ax.transAxes, fontsize=4)

        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('counts ', fontsize=8);
            ax.set_xlabel(r'$T_i (min)$', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i*20/60 for i in [0,xlim[-1]]], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)


#aca va el histograma de lambda
ylim = [0,0.02]
xlim = [0,300]

inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_row[3:-1], wspace=0, hspace=0.0)
color = 'blue';
plt.rcdefaults()


ax = plt.subplot(inner_gs[0])
ax.set_ylim(ylim);
ax.set_xlim(xlim)
  
    
ax.hist(lambda_box  ,bins=20,density = True,color=color, linewidth=0, alpha =0.8, label=label)
ax.grid(False);
ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
ax.set_ylabel('probability (1/min)', fontsize=8);
ax.set_xlabel(r'$T_i (min)$', fontsize=8)
ax.xaxis.set_label_coords(0.5, -0.1);
ax.yaxis.set_label_coords(-0.05, 0.5)
ax.set_xticks([xlim[0],25,40,xlim[1]])
ax.set_yticks(ylim)
ax.set_xticklabels([np.round(i*20/60,2) for i in [xlim[0],25,40,xlim[1]]], fontsize=6)
ax.set_yticklabels([i * 3 for i in ylim], fontsize=6)


         
            
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()

#%%
# =============================================================================
# calculating the time series shuffle
# =============================================================================

''' Time series forecasting with individual cell parametrization module. 

We want to generate single cells from single cells statistics by trying to fit the number of pulses of each single cell. 
'''
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/data/'
save_name= 'shuffle_1'
ipi_conc = ['an_WT_ESL']
get_single_cells_pulses_shuffle(conc_data,ipi_conc,save_folder+save_name)
    
    
conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')



#%%
# =============================================================================
# plotting funtions for shuffle
# =============================================================================

# =============================================================================
# =============================================================================
# TIME SERIES PLOT

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]
conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);



for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df_forecast = download_data(save_folder+save_name+'.pkl')
    df = df_forecast
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.exponential_IPI.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.amp_peaks.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)
     
#aca van las celulas posta
ylim = [45000,75000]

for k,C in enumerate(ipi_conc):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k+1], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[c]
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        data = df.query('cell ==' + str(cell))

    
        if cells[-1] >= cell:
            ax.plot(data.sm_MEAN_INTENSITY.values ,color=color, linewidth=0.4,label=label)
            ax.plot(data.max_.values,'o',color = color,markersize=1)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('amplitude (a.u)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.3, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()

#%%

# =============================================================================
# plotting funtions for proobe
# =============================================================================
''' La idea es probar si esto es poisson o no. '''

#aca vamos a generar las series temporales
#%%
from exponential_dm_main import generar_poisson, shuffle_poison,get_silent_poisson
import numpy as np

conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/neg_control/'

save_folder_save_name_1 = save_folder + 'neg_control'
df = conc_data['an_WT_ESL']
L = 10000000
generar_poisson(df,L,save_folder_save_name_1)

#aca generamos lo que viene despues
save_folder_save_name_2 = save_folder + 'neg_control_shuffled'
dm_1 = shuffle_poison(save_folder_save_name_1,save_folder_save_name_2)

#%%
#aca vamos a plotearlas
# =============================================================================
# =============================================================================

ylim = [-0.2,1.1]
l = 1000
xlim = [0,l]

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [25, 5]]
plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.0)
gs_row = gridspec.GridSpec(nrows=3, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);

color = colors[0];
df_forecast = download_data(save_folder_save_name_1+'.pkl')
ax = plt.subplot(gs_row[0])
ax.set_ylim(ylim);
ax.set_xlim(xlim)    
ax.plot(df_forecast.exponential_IPI.values[0:l] ,color=color, linewidth=0.4)
ax.plot(df_forecast.amp_peaks.values[0:l],'o',color = color,markersize=1)
ax.grid(False);
ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
ax.set_ylabel('amplitude (a.u)', fontsize=8);
ax.set_xlabel('time (min)', fontsize=8)
#ax.xaxis.set_label_coords(0.5, -1);
#ax.yaxis.set_label_coords(-0.2, 1)
ax.set_xticks([0,xlim[-1]])
ax.set_yticks(ylim)
#ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
#ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
#ax.set_yticklabels(ylim, fontsize=6)
 #       else:
  #          silent_ax(ax)


inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs_row[2], wspace=0.5, hspace=0.0)
ax3 = plt.subplot(inner_gs[0]); 
    
bins3 = ax3.hist(get_silent_poisson(df_forecast).dropna().values,bins=30, density=1, facecolor=color, alpha=0.8,linewidth=0.0 )
popt,pcov = get_linear_fit(bins3)
lambda_ = popt[0] 
ax3.set_ylim([0,0.02]);ax3.set_xlim([0,400])
ax3.plot(lambda_*np.exp(-lambda_*np.arange(400)),color = 'blue',alpha = 1,linewidth=0.5)
ax3.set_xlabel('amplitude (a.u)', fontsize=8);
ax3.set_ylabel('frequency (a.u.)', fontsize=8)
ax3.text(0.8, 1.07, 'lambda: '+str(lambda_), ha='center', va='center', transform=ax3.transAxes, fontsize=6)

#set_scale(ax3,np.arange(0,100,15),[0,0.1]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
#ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
#ax3.set_yticklabels([np.round(y/(20/60),2)  for y in [0,0.1]],fontsize=6)
#ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.0,0.5)


#ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# =============================================================================   
#van las series tenporaes post ahora

color = colors[1];
df_forecast = download_data(save_folder_save_name_2+'.pkl')


ax = plt.subplot(gs_row[1])
ax.set_ylim(ylim);
ax.set_xlim(xlim)
ax.plot(df_forecast.exponential_IPI.values[0:l] ,color=color, linewidth=0.4)
ax.plot(df_forecast.amp_peaks.values[0:l],'o',color = color,markersize=1)
ax.grid(False);
ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
ax.set_ylabel('amplitude (a.u)', fontsize=8);
ax.set_xlabel('time (min)', fontsize=8)
#ax.xaxis.set_label_coords(0.5, -1);
#ax.yaxis.set_label_coords(-0.2, 1)
ax.set_xticks([0,xlim[-1]])
ax.set_yticks(ylim)

ax3 = plt.subplot(inner_gs[1]); 
bins3 = ax3.hist(get_silent_poisson(df_forecast).dropna().values,bins=30, density=1, facecolor= color, alpha=0.8,linewidth=0.0 )
popt,pcov = get_linear_fit(bins3)
lambda_ = popt[0] 
ax3.set_ylim([0,0.02]);ax3.set_xlim([0,400])
ax3.plot(lambda_*np.exp(-lambda_*np.arange(400)),color = 'blue',alpha = 1,linewidth=0.5)
ax3.text(0.8, 1.07, 'lambda: '+str(lambda_), ha='center', va='center', transform=ax3.transAxes, fontsize=6)
ax3.set_xlabel('amplitude (a.u)', fontsize=8);
ax3.set_ylabel('frequency (a.u.)', fontsize=8)


fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()


#%%
###########################################################################################################################################
###########################################################################################################################################111111111111111111
################ ESTO ES PARA GENERAR DATOS INDIVUDUALES .PKL QUE CONTRIBUYEN AL PLOT DE CONSECUTIVIDAD ###################################
###########################################################################################################################################
###########################################################################################################################################

save_folder  = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/population_data/'
save_name= 'population'
bins1 = plt.hist(conc_data['an_WT_ESL']['dm_'],bins=np.arange(2,200,6), density=1);

dm_max = 0
lambda_ = popt[0] 

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    exponential = np.random.exponential(1/lambda_,len(conc_data['an_WT_ESL']['dt_peaks']))
    get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)

#%%
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/population_20_min_data/'
save_name = 'population_20_min'
popt,pcov = get_linear_fit(bins1,20)
dm_max = 20 #minutes
lambda_ = popt[0] 

for i in range(100):
    save_folder_save_name = save_folder + save_name + str(i)
    exponential = np.random.exponential(1/lambda_,len(conc_data['an_WT_ESL']['dt_peaks']))
    get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)
    
    
    
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/single_cell_silences_data/'
save_name = 'single_cell_silences_data'

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    get_single_cell_time_series(conc_data,ipi_conc,save_folder_save_name)
  #%%
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/new_cdc_data/'
save_name = 'new_cdc_'

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    get_single_cell_time_series(conc_data,ipi_conc,save_folder_save_name,True)

#%%

save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/moving_bc/moving_bc_data/'
save_name = 'moving_bc'

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    dm_cells,lambda_box = get_single_cell_time_series_moving_boundaries(conc_data,ipi_conc,save_folder_save_name)


#%%
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/data/'
save_name = 'shuffle_'

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    get_single_cells_pulses_shuffle(conc_data,ipi_conc,save_folder_save_name)
#%%
