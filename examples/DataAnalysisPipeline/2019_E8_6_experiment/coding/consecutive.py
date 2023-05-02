'''Para que esto funcione, tengo que primero correr la figura 3bis, y después el modulo WT de consecutive'''
''' Me da el plot de consecutividad para distintos valores de FGF'''

import sys

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')

#Esto es para plotear los datos de consecutividad
from consecutive_main import consecutive_non_cumulative,makeColours,consecutive_cumulative        
from matplotlib import cm

c = 'an_20ng'
#%%
hypothesis_list_df = [conc_data[c],download_data(save_folder + 'population_20_min.pkl'),download_data(save_folder + 'population.pkl'),download_data(save_folder + 'single_cell_silences.pkl')]
labels_hypothesis_list_df =[ 'experiment_'+c,'population_20_min','population','single_cell_silences']
colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs_main[0,0])
ax2 = plt.subplot(gs_main[0,1])


for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]

    consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
    cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()

    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
    ax1.set_xlim([0, 18]);
    ax1.set_yscale('log')
    ax1.set_ylim([-1,150])
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.1,0.5);
    
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 18]); ax2.set_ylim([-1,150])
    ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax2.xaxis.set_label_coords(0.5, -0.08);
    ax2.yaxis.set_label_coords(-0.1,0.5);



ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
plt.savefig(save_folder + '_FGF20ng_consecutive_pulses_.pdf', format='pdf')

#%% aca lo que hacemos es plotear todo junto, solamente el 80 porciento de las cosas! (seleccionadas como las que tienen la 
#actividad del medio)


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)


gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=len(hypothesis_list_df)-2, subplot_spec=gs_main[0], wspace=0.5, hspace=0.3)
ax1_list = [gs0_inner[0,0],gs0_inner[0,1],gs0_inner[1,0],gs0_inner[1,1]]


activity_filter_cells_list = []
for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]
    color = colors[i]
    ax1 = plt.subplot(ax1_list[i])
    
    activity = df_consecutive['dt_peaks'].groupby(level='cell').sum() / df_consecutive['FRAME'].groupby(level='cell').count() *  100   
    activity_index = np.argsort(activity.values)[::-1]
    activity = [activity[j] for j in activity_index]

    y = len(activity)*2/10 // 2 #lo que le quiero sacar de cada lado
    activity_filter_cells = []
    
    for n,(ix,cell) in enumerate(zip(activity_index,activity)):
        if (n > y) and (n<= len(activity) - y):
            activity_filter_cells.append(ix)
    
    silent = np.ones(len(activity)) * 100 - activity
    if i == 0: 
        activity_experiment = activity
        silent_experiment = silent
        
    p1 = ax1.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=color,alpha=0.8,linewidth=0.0)
    ax1.bar(np.arange(0,len(df_consecutive.index.get_level_values(0).unique())),activity_experiment,bottom=silent_experiment,width=0.8,color=colors[0],alpha=0.2,linewidth=0.0)
    if i == 0: plt.legend((p1[0], p2[0]), ('non pulsing', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.3), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)
    ax1.set_xlim([-1,df_consecutive.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.hlines(50,y,len(activity) - y,color = 'black')
    
    ax1.set_xlabel( labels_hypothesis_list_df[i],fontsize=8); 
    ax1.set_xticks([0,df_consecutive.index.get_level_values(0)[-1]])
    ax1.set_yticks([0,50,100])
    ax1.set_yticklabels([0,0.5,1],fontsize=6)
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

    ax1.xaxis.set_label_coords(0.5,-0.06)
    
    if i == 0:
        ax1.set_ylabel('fraction of cell track' ,fontsize=8); 
        ax1.yaxis.set_label_coords(-0.13,0.05)


    activity_filter_cells_list.append(activity_filter_cells)
#activity_filter_cells_list = [[1,68],[1,68],[1,68],[1,68]]
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs_main[1], wspace=0.5, hspace=0.0)
  
ax1 = plt.subplot(gs0_inner[0])
ax2 = plt.subplot(gs0_inner[1])

for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]
    
    consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive[df_consecutive.index.get_level_values(0).isin(activity_filter_cells_list[i])])
    df.query('cell == ' + str(cell))
    cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()

    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
    ax1.set_xlim([0, 16]);
    ax1.set_yscale('log')
    ax1.set_ylim([-1,150])
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.1,0.5);
    
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive[df_consecutive.index.get_level_values(0).isin(activity_filter_cells_list[i])])
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 16]); ax2.set_ylim([-1,150])
    ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax2.xaxis.set_label_coords(0.5, -0.08);
    ax2.yaxis.set_label_coords(-0.1,0.5);

ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)

plt.savefig(save_folder + 'FGF20ng_activity_consecutive_pulses_.pdf', format='pdf')


############################################################################################################################
###############OLD. REVISAR!
############################################################################################################################

    #%% Esto es para los datos reales! (experimentales)
   
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 


hypothesis_list_df = [conc_data[c] for c in conc_labels]
labels_hypothesis_list_df = [conc_labels[c] for c in conc_labels]

save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/'

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs_main[0,0])
ax2 = plt.subplot(gs_main[0,1])


for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]

    consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
    cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()

    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
    ax1.set_xlim([0, 10]);
    ax1.set_yscale('log')
    ax1.set_ylim([-1,350])
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.09);ax1.yaxis.set_label_coords(-0.1,0.5);
    
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 20]); ax2.set_ylim([-1,600])
    ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax2.xaxis.set_label_coords(0.5, -0.09);
    ax2.yaxis.set_label_coords(-0.1,0.5);



ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
plt.savefig(save_folder + 'consecutive_pulses_experimental_FGF.pdf', format='pdf')
