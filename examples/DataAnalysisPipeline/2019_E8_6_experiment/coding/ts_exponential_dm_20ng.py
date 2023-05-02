import sys

sys.path.append('/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from exponential_dm_main import download_data,get_linear_fit,get_population_time_series,get_single_cell_time_series,find_dm
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/figures_05_2022/'


def mixed_dt(df): #Esto en realidad hay que exportarlo de la figura 2 para ser consistentes!
    values_dt_mixed = []
    for cells, data in df.groupby(level='cell'):
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
            #print(dt_raise,dt_fall)
    return(values_dt_mixed)
    
    
    
ipi_conc = ['an_20ng']

find_dm(conc_data['an_20ng'])
bins1 = plt.hist(conc_data['an_20ng']['dm_'],bins=np.arange(2,200,6), density=1);

popt,pcov = get_linear_fit(bins1)
dm_max = 0
popt,pcov = get_linear_fit(bins1,20)
dm_max = 20 #minutes
lambda_ = popt[0]   
exponential = np.random.exponential(1/lambda_,len(conc_data['an_20ng']['dt_peaks']))


save_name = 'population_20_min'
#save_name = 'population'

# =============================================================================
# calculating the time series population simulation
# =============================================================================
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
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}
colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]

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


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax3 = plt.subplot(gs_0[0,1]); 
    
    bins3 = ax3.hist(df[col2],bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.15]);

    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode (m): '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(0.8, 1.07, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.8, 1.05,r'$Q$ (m): '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    
    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode (f): '+str(np.round(mode_dt,2))+' \n'
    ax3.text(0.8, 0.92, mode_dt, ha='center', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(0.8, 0.9,r'$Q$ (f): '+str(np.round(df[col2].quantile(0.25),2))+' ; '+str(np.round(df[col2].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75),2)) , ha='center', va='center', transform=ax3.transAxes, fontsize=6)
  

set_scale(ax3,np.arange(0,100,15),[0,0.15]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
ax3.set_yticklabels([np.round(y/(20/60),2)  for y in [0,0.15]],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.07);ax3.yaxis.set_label_coords(-0.0,0.5)


ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# IPI
# =============================================================================
        


plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'IPI'


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[1,0]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.065])

    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode (m): '+str(np.round(mode_IPI*20/60,2))+' \n'
    ax1.text(0.65, 0.90, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode (f): '+str(np.round(mode_IPI,2))+' \n'
    ax1.text(0.65, 0.70, mode_IPI, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(df[col3].quantile(0.25),2))+' ; '+str(np.round(df[col3].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75),2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    ax1.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(df[col3]),2)) +'\n std : '+str(np.round(np.std(df[col3]),2)), ha='center', va='center', transform=ax1.transAxes, fontsize=6)

set_scale(ax1,np.arange(0,152,30),[0,0.065]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,152,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.065]],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
# =============================================================================
# dm
# =============================================================================

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

col3 = 'dm_'


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[1,1]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,200,6), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.axvspan(0,dm_max*60/20, color='gray',alpha = 0.4,linewidth=0.0)
    ax1.plot(lambda_*np.exp(-lambda_*np.arange(200)),color = 'blue',alpha = 1,linewidth=0.5)

    ax1.set_xlim([0,200]);ax1.set_ylim([0,0.065])

    mode_dm =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_dm = 'mode (m): '+str(np.round(mode_dm*20/60,2))+' \n'
    ax1.text(0.65, 0.90, mode_dm, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.87,r'$Q$ (m): '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    mode_dm =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_dm = 'mode (frames): '+str(np.round(mode_dm,2))+' \n'
    ax1.text(0.65, 0.70, mode_dm, ha='center', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(0.65, 0.67,r'$Q$ (f): '+str(np.round(df[col3].quantile(0.25),2))+' ; '+str(np.round(df[col3].quantile(0.50),2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75),2)) , ha='center', va='center', transform=ax1.transAxes, fontsize=6)

    ax1.text(0.65, 0.5, 'mean (f) : ' + str(np.round(np.mean(df[col3]),2)) +'\n std : '+str(np.round(np.std(df[col3]),2)), ha='center', va='center', transform=ax1.transAxes, fontsize=6)

set_scale(ax1,np.arange(2,200,60),[0,0.065]); ax1.set_xlabel('time between pulses (min)',fontsize=8); ax1.set_ylabel('probability',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.55, -0.07);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(2,200,60)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(y/(20/60),2) for y in [0,0.065]],fontsize=6)


ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
# =============================================================================
# exponential distrobution
# 
# =============================================================================
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
col2='dt_peaks'#ax3

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


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs_0[2,0]); 
    
    bins1 = ax1.hist(df[col3],bins=np.arange(2,152,6),range=(2,None), density=1, facecolor=colors[k], alpha=0.8,linewidth=0.0 )
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    ax1.set_yscale('log')
    ax1.set_xlim([0,152]);ax1.set_ylim([0,0.15])


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


''' Plotting the time series that we generated prevously'''

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5

colors =  sns.color_palette(sns.color_palette("Blues_d",len(conc_data)))

len_cells = []
for c in ['an_20ng']:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=5, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);



for k,C in enumerate(['an_20ng']):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[0];
    df = download_data(save_folder_save_name+'.pkl')
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
fig.savefig(save_folder + str(save_name) + '_time_series_.pdf', format='pdf')
plt.show()
plt.close()
#%%
df_forecast = download_data(save_folder_save_name+'.pkl')

fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
#Main 6 filas

gs_main0 = gridspec.GridSpec(nrows=2, ncols=1, figure=fig)#,width_ratios= [1,4],height_ratios=[1,4]); 
gs_main0.update(left=0.1, right=0.95, bottom=0.1, top=1, hspace=0.0,wspace = 0.0)

gs_main = gridspec.GridSpecFromSubplotSpec(nrows=6, ncols=1, subplot_spec=gs_main0[0], wspace=0.0, hspace=0.5)

#gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); 
#gs_main.update(left=0.1, right=0.95, bottom=0.1, top=1, hspace=0.5,wspace = 0.0)
gs_main_ = [gs for gs in gs_main]


# =============================================================================
# 
# =============================================================================
#Lugar 0: esquema 
#lugar 1 traza2
#gs_0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main_[1], wspace=0.5, hspace=0.5)





ylim = [-0.1, 0.1]
M = max([conc_data[key].max().FRAME for key in conc_data.keys()])

cells_conc = {'an_0ng':[1,32,53,55,59],'an_2-5ng':[0,8,19,33,15],'an_5ng':[1,4,49,53,55],'an_20ng':[2,6,17,26,64]} 
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

box=[];handles = [None]*len(conc_labels)
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]



# =============================================================================
# trazas
# =============================================================================
inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=len(cells_conc), subplot_spec=gs_main_[1], wspace=0.2, hspace=0.4)

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 

for k,c in enumerate(['an_20ng']):
    df = df_forecast;df.sort_index(inplace=True); color=colors[k]

    inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc[c]), ncols=1, subplot_spec=inner_gs0[k], wspace=0.0, hspace=0.0)

    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(inner_inner_gs0[w]);       
        ax.patch.set_alpha(0.0)
        ax.patch.set_linewidth(0.0)
    
        ax.plot(df.loc[cells_conc[c][w]]['exponential_IPI'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ylim_=[-0.1,1.1]
        ax.set_ylim(ylim_)

        #if (k==len(cells_conc)-1) and (w == len(cells_conc[c])-1) :
        if (k==0) and (w == len(cells_conc[c])-1):   
           set_scale_bars(ax, x_bar=(0,60), y_bar=ylim_, xunits='', yunits='', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
           ax.set_ylabel('     KTR signal \n 200 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
           ax.xaxis.set_label_coords(0.3, -0.5);
           ax.yaxis.set_label_coords(-0.08,2.2);
           
        elif (k!=0) and (w == len(cells_conc[c])-1): 
            set_scale_bars(ax, x_bar=(0,60), y_bar=ylim_, xunits='', yunits='', x_scale=60/20, y_scale=100, round_x = True,round_y=True)

        else:
            set_scale_bars(ax)

        if w == 0: 
            ax.text(0.05, 1.5,  conc_labels[c], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = colors[k])


# =============================================================================
# statitics
# =============================================================================

ipi_conc = ['an_20ng']
gs3_inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[3:6], wspace=0.4, hspace=0.0)
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=len(ipi_conc), ncols=3, subplot_spec=gs3_inner[0:3], wspace=0.4, hspace=0.4)
col3 = 'IPI'#ax1
col2='dt_peaks'#ax3
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 



for k,c in enumerate(ipi_conc):
    df = df_forecast
    ax1 = plt.subplot(gs0_inner[k,1]); ax3 = plt.subplot(gs0_inner[k,0]); #ax4 = plt.subplot(gs0[k,2]);
    bins1= ax1.hist(df[col3].dropna(),bins=np.arange(2,80,4),range=(2,None), density=1, facecolor=colors[k+1], alpha=0.8,linewidth=0); 
    bins3= ax3.hist(df[col2].dropna(),bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k+1], alpha=0.8,linewidth=0); 


    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.set_xlim([0,80]);ax1.set_ylim([0,0.07])
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.12]);


    #mode_IPI = mode(df[col3])[0][0]*20/60; mode_IPI = 'mode: '+str(np.round(mode_IPI,2))+' \n'
    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
    ax1.text(1, 0.80, mode_IPI, ha='right', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(1, 1.05,r'$Q$: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='right', va='center', transform=ax1.transAxes, fontsize=6)

    #mode_dt = mode(df[col2])[0][0]*20/60 ; mode_dt = 'mode: '+str(np.round(mode_dt,2))+' \n'
    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(1, 0.80, mode_dt, ha='right', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(1, 1.05,r'$Q$: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='right', va='center', transform=ax3.transAxes, fontsize=6)

    if k != len(ipi_conc)-1:
        set_scale(ax1,np.arange(0,80,30),[0,0.07]) ;ax1.set_xticklabels([None,None]);ax1.set_yticklabels([None,None]);
        set_scale(ax3,np.arange(0,46,15),[0,0.12]);ax3.set_xticklabels([None,None]);ax3.set_yticklabels([None,None]);
        
    
    ax3.text(0.05, 0.9,  conc_labels[c], ha='left', va='center', transform=ax3.transAxes, fontsize=6,color = colors[k+1])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
set_scale(ax1,np.arange(0,80,30),[0,0.07]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('counts',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.5, -0.3);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,80,30)*20/60],fontsize=6)


set_scale(ax3,np.arange(0,46,15),[0,0.12]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('counts',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,46,15)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.3);ax3.yaxis.set_label_coords(-0.01,0.5)

for ax_ in [ax1,ax3]: ax_.tick_params(labelsize=6,direction='out', pad=1,length=2)



# =============================================================================
# silence vs mixed dt
# =============================================================================


for k,c in enumerate(ipi_conc):
    df = df_forecast
    


    ax3 = plt.subplot(gs0_inner[k,2]); 
    ax3.plot(df['IPI'].dropna().values-mixed_dt(df),mixed_dt(df),'o',markersize=1.5,alpha=0.5,color = colors[k+1])
    ax3.set_xlim([-5, 60]);
    ax3.set_ylim( [-5, 60])
    outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-mixed_dt(df),mixed_dt(df))]
    outliers_sum = np.sum([1 for i,j in outliers if i>60 or j >60])
    
    #print(c,sum(1* (df['IPI'].dropna().values-mixed_dt(df) < 0))) 
    
    print('outliers : ' + str(outliers_sum))
    print('total_data: ' + str(len(outliers)))
    
    X = np.arange(0,200)
    ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.8,linestyle='--') 
   
    
    x_label = 'joint duration (min)'
    if k ==2: ax3.set_xlabel('silence (min)',fontsize=8); 
    if k ==2: ax3.set_ylabel(x_label,fontsize=8); 
    if k != len(ipi_conc)-1:
            set_scale(ax3,np.arange(0,240,30),np.arange(0,240,30));ax3.set_xticklabels([None,None]);ax3.set_yticklabels([None,None]);


    ax3.set_aspect(1)
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
set_scale(ax3,np.arange(0,240,30),np.arange(0,240,30))
ax3.set_xticklabels([int(x) for x in np.arange(0,240,30)*20/60],fontsize=6)
ax3.set_yticklabels([int(x) for x in np.arange(0,240,30)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.3); ax3.yaxis.set_label_coords(-0.30,0.5)

ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)



# =============================================================================
# pulse density
# =============================================================================
#####

#gs3_inner_inner[0].remove()
gs1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=4, subplot_spec=gs_main[2:4], wspace=0.6, hspace=0.5,width_ratios=[1.5,0.2,0.6,0.9])


ylim = ylim

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
box=[];handles = [None]*2


for k,c in enumerate(ipi_conc):
    df = df_forecast;df.sort_index(inplace=True); color=colors[k]
    box.append((df['amp_peaks'].groupby(level='cell').count() / df['amp_peaks'].groupby(level='cell').size()).to_list()) 

#label = [conc_labels[c] for c in cells_conc]
label = [ 20]
X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
gs2_ = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs1[0:2,3], wspace=0.0, hspace=0.0,height_ratios=[1.9,0.1])

ax2 = plt.subplot(gs2_[0]);

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
ax2.set_xticklabels(label,rotation = 0)
label_peaks = r'$\rm{\left(min^{-1}\right)}$'
ax2.set_ylabel('pulse density ' + label_peaks,fontsize=8)
ax2.set_xlabel('[FGF4] (ng/ml)',fontsize=8)
ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.01,0.55)


yticks = ax2.get_yticks(); yticks[-1] = 0.06 
minutes_ylim = yticks[-1]/(20/60)
ax2.set_yticklabels([0,minutes_ylim],fontsize=6)
ax2.set_ylim([-0.005,yticks[-1]]); ax2.set_yticks([0,yticks[-1]])


ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
#ax2.invert_yaxis()

# =============================================================================
# activity
# =============================================================================

gs0_inner_ = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc), ncols=1, subplot_spec=gs1[0,0:2], wspace=0.0, hspace=0.0)


for k,c in enumerate(ipi_conc):
    ax1 = plt.subplot(gs0_inner_[k]);
    df = df_forecast
    activity = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=colors[k],alpha=0.8,linewidth=0.0)
    #if k == 0: plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.6), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)
    ax1.set_xlim([-1,df.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
    
    if k != len(ipi_conc):
        silent_ax(ax1);
    else:
        ax1.set_yticks([0,100]); ax1.set_xticks([]);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2);


    #ax1.set_xlabel( conc_labels[c] + ' concentration' ,fontsize=10); 
    #ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
ax1.set_ylabel('% of cell track' ,fontsize=8); 
ax1.yaxis.set_label_coords(-0.05,1.5)
ax1.set_xlabel('individual cells per condition' ,fontsize=8); 
ax1.xaxis.set_label_coords(0.5,-0.15)


# =============================================================================
# mean activity
# =============================================================================

ax1 = plt.subplot(gs1[0,2]);

activity = []; silent = []
activity_err = []; silent_err = []

for c in ipi_conc:

    activity_ = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))



p1 = ax1.barh(y = np.arange(len(ipi_conc)),width = silent,xerr=silent_err,left = 0, color='darkgray',alpha=0.5,linewidth=0.0 )
p2 = ax1.barh(y = np.arange(len(ipi_conc)),width = activity,xerr=activity_err,left = silent, color=colors,alpha=0.8,linewidth=0.0 )
plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(-0.1, 1.43), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)



#p1 = ax1.bar(np.arange(len(cells_conc)),silent,yerr=silent_err,width=0.8,color='darkgray',alpha=0.5,linewidth=0.4)
#p2 = ax1.bar(np.arange(len(cells_conc)),activity,bottom=silent,yerr = activity_err,width=0.8,color=colors,alpha=0.8,linewidth=0.4)
#plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 1.4), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)

labels = [20]
plt.xlabel('% of cell track' ,fontsize=8); 
plt.ylabel('[FGF4] (ng/ml)',fontsize=8); 
plt.xticks([0,50,100]);plt.xlim([0,100]);plt.ylim([-0.55,3.5])
plt.yticks(np.arange(len(cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.xaxis.set_label_coords(0.5, -0.2);ax1.yaxis.set_label_coords(-0.24,0.5)
ax1.invert_yaxis()


plt.savefig(save_folder+ save_name +  '_figure3.pdf', format='pdf',transparent=True)

#%%

# =============================================================================
# calculating the time series single_cell simulation
# =============================================================================

''' Time series forecasting with individual cell parametrization module. 

We want to generate single cells from single cells statistics by generating an exponential
dristibution of silent intervals from single cells mean value of silent intervals. 
'''
save_name= 'single_cell_silences'
get_single_cell_time_series(conc_data,ipi_conc,save_folder+save_name)
    
    

#%%
# =============================================================================
# plotting funtions for single cell simulation
# =============================================================================

# =============================================================================
# =============================================================================
# TIME SERIES PLOT
df_forecast = download_data(save_folder+save_name+'.pkl')

ylim = [-0.2,1.1]
xlim = [0,330]
long=1;Cols=5

cells_conc = {'an_0ng':[1,32,53,55,59],'an_2-5ng':[0,8,19,33,15],'an_5ng':[1,4,49,53,55],'an_20ng':[2,6,17,26,64]} 
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

box=[];handles = [None]*len(conc_labels)
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]


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
            ax.yaxis.set_label_coords(-0.2, 1)
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
    Rows = len_cells[-1] // Cols;
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
            lambda_ = np.mean(data.dm_.dropna().values)
            if len(data.dt_peaks.dropna().values)> 0: 
                exponential_cell = np.random.exponential(lambda_,1000)
            else: 
                exponential_cell = np.array([])

            if exponential_cell.size > 0: 
                ax.hist(exponential_cell ,bins=30,density = True,color=color, linewidth=0.4,label=label)
                ax.plot(1/lambda_*np.exp(-1/lambda_*np.arange(600)),color = 'blue',alpha = 1,linewidth=0.5)

            ax.text(0.95, 0.2, str(cell) , ha='center', va='center', transform=ax.transAxes, fontsize=6)
            ax.text(0.7, 0.8, r'$\lambda = $'+str(np.round(lambda_/3,2)) + 'min' , ha='center', va='center', transform=ax.transAxes, fontsize=6)
            #CAMBIAR A FRAMES
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('counts ', fontsize=8);
            ax.set_xlabel(r'$T_i (min)$', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            #ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            ax.set_xticklabels([i*20/60 for i in [0,xlim[-1]]], fontsize=6)
            ax.set_yticklabels(ylim, fontsize=6)


        else:
            silent_ax(ax)


#aca va el histograma de lambda
ylim = [0,0.04]
xlim = [0,210]

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

    
ax.hist(lambda_box  ,bins=30,density = True,color=color, linewidth=0, alpha =0.8, label=label)
ax.grid(False);
ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
ax.set_ylabel('probability (1/min)', fontsize=8);
ax.set_xlabel(r'$T_i (min)$', fontsize=8)
ax.xaxis.set_label_coords(0.5, -0.1);
ax.yaxis.set_label_coords(-0.05, 0.5)
ax.set_xticks([xlim[0],15,40,xlim[1]])
ax.set_yticks(ylim)
ax.set_xticklabels([np.round(i*20/60,2) for i in [xlim[0],15,40,xlim[1]]], fontsize=6)
ax.set_yticklabels([i * 3 for i in ylim], fontsize=6)


         
            
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.show()
plt.close()
#%%
###########################################################################################################################################
###########################################################################################################################################
################ ESTO ES PARA GENERAR DATOS INDIVUDUALES .PKL QUE CONTRIBUYEN AL PLOT DE CONSECUTIVIDAD ###################################
###########################################################################################################################################
###########################################################################################################################################

save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/figures_05_2022/data/population/'
save_name = 'population'
bins1 = plt.hist(conc_data['an_20ng']['dm_'],bins=np.arange(2,200,6), density=1);

dm_max = 0
lambda_ = popt[0] 

for i in range(100):
    save_folder_save_name = save_folder + save_name + str(i)
    exponential = np.random.exponential(1/lambda_,len(conc_data['an_20ng']['dt_peaks']))
    get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)

save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/figures_05_2022/data/population_20_min/'
save_name = 'population_20_min'
popt,pcov = get_linear_fit(bins1,20)
dm_max = 20 #minutes
lambda_ = popt[0] 

for i in range(100):
    save_folder_save_name = save_folder + save_name + str(i)
    exponential = np.random.exponential(1/lambda_,len(conc_data['an_WT_ESL']['dt_peaks']))
    get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)
    
    
    
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/figures_05_2022/data/single_cell_silences_data/'
save_name = 'single_cell_silences_data'

for i in range(200):
    save_folder_save_name = save_folder + save_name + str(i)
    get_single_cell_time_series(conc_data,ipi_conc,save_folder_save_name)
