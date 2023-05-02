#me lleve todo a consecutive_main
import sys
import os
sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from consecutive_main import consecutive_non_cumulative,makeColours,consecutive_cumulative,mean_consecutive_value,filter_activity_cells,first_element_box,second_element_box
from exponential_dm_main import download_data
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


save_folder ='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/'
conc_data = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')

#%%
#Esto es para plotear los datos de consecutividad
        
hypothesis_list_df = [conc_data[c],download_data(save_folder + 'population.pkl'),download_data(save_folder + 'population_20_min.pkl'),download_data(save_folder + 'single_cell_silences.pkl'),download_data(save_folder + 'single_cell_silences_new_cdc.pkl'),download_data(save_folder + 'single_cell_silences_moving_bc.pkl')]
labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','population_20_min','single_cell_silences','single_cell_silences_new_bc','single_cell_silences_moving_bc']
#%%
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
    ax1.set_xlim([0, 15]);
    ax1.set_yscale('log')
    ax1.set_ylim([-1,350])
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.09);ax1.yaxis.set_label_coords(-0.1,0.5);
    ax1.set_xticks([0,3,6,9,12,15])

    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 15]); ax2.set_ylim([-1,500])
    ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax2.xaxis.set_label_coords(0.5, -0.09);
    ax2.yaxis.set_label_coords(-0.1,0.5);
    ax2.set_xticks([0,3,6,9,12,15])




ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
plt.savefig(save_folder_4 + 'consecutive_pulses_.pdf', format='pdf')

#%% aca lo que hacemos es plotear todo junto, solamente el 80 porciento de las cosas! (seleccionadas como las que tienen la 
#actividad del medio)

#hypothesis_list_df = [conc_data[c],download_data(save_folder + 'population.pkl'),download_data(save_folder + 'population_20_min.pkl'),download_data(save_folder + 'single_cell_silences.pkl'),download_data(save_folder + 'single_cell_silences_new_cdc.pkl')]
#labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','population_20_min','single_cell_silences','single_cell_silences_new_bc']
#colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))


c = 'an_WT_ESL'
save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/'
#save_folder_2 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/'
#save_folder_3 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/moving_bc/'
save_folder_4='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/'

#hypothesis_list_df_folder = [False,save_folder_1+'population_data/',save_folder_1+'population_20_min_data/',save_folder_1+'single_cell_silences_data/',save_folder_2+'new_cdc_data/',save_folder_3 + 'moving_bc_data/']
#labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','population_20_min','single_cell_silences','single_cell_silences_new_bc','single_cell_silences_moving_bc']
hypothesis_list_df = [conc_data[c],download_data(save_folder + 'population.pkl'),download_data(save_folder_4+'data/shuffle_1.pkl')]
labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','shuffle']


colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)


gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=2, subplot_spec=gs_main[0], wspace=0.5, hspace=0.3)
ax1_list = [gs0_inner[0,0],gs0_inner[1,0],gs0_inner[1,1],gs0_inner[2,0],gs0_inner[2,1]]


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
    #df.query('cell == ' + str(cell))
    cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()

    ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
    ax1.set_xlim([0, 10]);
    ax1.set_yscale('log')
    ax1.set_ylim([-1,350])
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
    ax1.set_xticks([0,2,4,6,8,10])

    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive[df_consecutive.index.get_level_values(0).isin(activity_filter_cells_list[i])])
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 10]); ax2.set_ylim([-1,500])
    ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=10)
    ax2.xaxis.set_label_coords(0.5, -0.08);
    ax2.yaxis.set_label_coords(-0.2,0.5);
    ax2.set_xticks([0,2,4,6,8,10])


ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)

plt.savefig(save_folder_4 + 'consecutive_pulses_.pdf', format='pdf')
#%%
###########################################################################################################################################
###########################################################################################################################################
################esto es para plotear el grafico de consecutividad con promedio! ######################################
###########################################################################################################################################
###########################################################################################################################################
#save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/'
#save_folder_2 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/'
#save_folder_3 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/moving_bc/'
#hypothesis_list_df_folder = [False,save_folder_1+'population_data/',save_folder_1+'population_20_min_data/',save_folder_1+'single_cell_silences_data/',save_folder_2+'new_cdc_data/',save_folder_3 + 'moving_bc_data/']
#labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','population_20_min','single_cell_silences','single_cell_silences_new_bc','single_cell_silences_moving_bc']


c = 'an_WT_ESL'
save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/'
save_folder_4='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/'

hypothesis_list_df_folder = [False,save_folder_1+'population_data/',save_folder_4+'data/']
labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','shuffle']
colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs_main[0,0])
ax2 = plt.subplot(gs_main[0,1])


for i,df_consecutive_folder in enumerate(hypothesis_list_df_folder):
    label = labels_hypothesis_list_df[i]
    
    if df_consecutive_folder is False:
        df_consecutive = conc_data[c]
        
        consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
        cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
    
        ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
        ax1.set_xlim([0, 15]);
        ax1.set_yscale('log')
        ax1.set_ylim([-1,350])
        ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses (100%) ',fontsize=10)
        ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
        ax1.set_xticks([0,3,6,9,12,15])

        consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
        box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
        
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.set_yscale('log')
        ax2.set_xlim([0, 15]); ax2.set_ylim([-1,500])
        ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses (100%) ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.2,0.5);
        ax2.set_xticks([0,3,6,9,12,15])

    else:
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []
        for df_file_name in os.listdir(df_consecutive_folder) :
        
            df_consecutive = download_data(df_consecutive_folder+df_file_name)
            
            consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
            cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
            cell_population_list_trials.append(cell_population_list)
            
            consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
            
           
        cell_population_list_trials_mean,cell_population_list_trials_std = mean_consecutive_value(cell_population_list_trials)
        ax1.plot(np.arange(1,len(cell_population_list_trials_mean)+1),cell_population_list_trials_mean, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
        ax1.fill_between(np.arange(1,len(cell_population_list_trials_mean)+1),cell_population_list_trials_mean-cell_population_list_trials_std,cell_population_list_trials_mean+cell_population_list_trials_std,color =colors[i],alpha = 0.2)
        ax1.set_xlim([0, 15]);
        ax1.set_yscale('log')
        ax1.set_ylim([-1,350])
        ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses (100%) ',fontsize=10)
        ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
        ax1.set_xticks([0,3,6,9,12,15])
        

        box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean-box_plot_consecutive_cumulative_trials_std,box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std,color =colors[i],alpha = 0.2)        
        ax2.set_yscale('log')
        ax2.set_xlim([0, 15]); ax2.set_ylim([-1,500])
        ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses (100%) ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.2,0.5);
        ax2.set_xticks([0,3,6,9,12,15])
        
ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
ax1 = plt.subplot(gs_main[1,0])
ax2 = plt.subplot(gs_main[1,1])


for i,df_consecutive_folder in enumerate(hypothesis_list_df_folder):
    label = labels_hypothesis_list_df[i]
    
    if df_consecutive_folder is False:
        df_consecutive = conc_data[c]
        activity_filter_cells = filter_activity_cells(df_consecutive)
        df_consecutive = df_consecutive[df_consecutive.index.get_level_values(0).isin(activity_filter_cells)]
        
        consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
        cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
    
        ax1.plot(np.arange(1,len(cell_population_list)+1),cell_population_list, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
        ax1.set_xlim([0, 15]);
        ax1.set_yscale('log')
        ax1.set_ylim([-1,350])
        ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses (80%) ',fontsize=10)
        ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
        ax1.set_xticks([0,3,6,9,12,15])
        
        consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
        box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
        
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),box_plot_consecutive_cumulative, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.set_yscale('log')
        ax2.set_xlim([0, 15]); ax2.set_ylim([-1,500])
        ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses (80%) ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.2,0.5);
        ax2.set_xticks([0,3,6,9,12,15])


    else:
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []
        for df_file_name in os.listdir(df_consecutive_folder) :
        
            df_consecutive = download_data(df_consecutive_folder + df_file_name)
            activity_filter_cells = filter_activity_cells(df_consecutive)
            
            df_consecutive = df_consecutive[df_consecutive.index.get_level_values(0).isin(activity_filter_cells)]
            consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
            cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
            cell_population_list_trials.append(cell_population_list)
            
            consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
            
           
        cell_population_list_trials_mean,cell_population_list_trials_std = mean_consecutive_value(cell_population_list_trials)
        ax1.plot(np.arange(1,len(cell_population_list_trials_mean)+1),cell_population_list_trials_mean, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[i],label=labels_hypothesis_list_df[i])
        ax1.fill_between(np.arange(1,len(cell_population_list_trials_mean)+1),cell_population_list_trials_mean-cell_population_list_trials_std,cell_population_list_trials_mean+cell_population_list_trials_std,color =colors[i],alpha = 0.2)
        ax1.set_xlim([0, 15]);
        ax1.set_yscale('log')
        ax1.set_ylim([-1,350])
        ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses (80%) ',fontsize=10)
        ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
        ax1.set_xticks([0,3,6,9,12,15])
        

        box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std = mean_consecutive_value(box_plot_consecutive_cumulative_trials)        
        ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean, color=colors[i], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=labels_hypothesis_list_df[i])
        ax2.fill_between(np.arange(1,len(box_plot_consecutive_cumulative_trials_mean)+1),box_plot_consecutive_cumulative_trials_mean-box_plot_consecutive_cumulative_trials_std,box_plot_consecutive_cumulative_trials_mean+box_plot_consecutive_cumulative_trials_std,color =colors[i],alpha = 0.2)        
        ax2.set_yscale('log')
        ax2.set_xlim([0, 15]); ax2.set_ylim([-1,500])
        ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses (80%) ',fontsize=10)
        ax2.xaxis.set_label_coords(0.5, -0.08);
        ax2.yaxis.set_label_coords(-0.1,0.5);
        ax2.set_xticks([0,3,6,9,12,15])

        
ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
plt.savefig(save_folder + 'consecutive_pulses_WT_mean_all.pdf', format='pdf')
#%%
#4) Hacer plots de 1) Pulsos aislados y 2) pulsos totales y 3) total pulsos consecutivos 
# =============================================================================
# boxplot consecutiveness
# =============================================================================

#si sumas box1 y box2, te puede dar uno (si hay pico) o cero (si no hay picos)
    # Filtramos la celulas que tienen cero picos en orden de poder hacer estadistica solo de picos

c = 'an_WT_ESL'
save_folder_1 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/figures_05_2022/'
#save_folder_2 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/'
#save_folder_3 = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/moving_bc/'
save_folder_4='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/shuffle/'

#hypothesis_list_df_folder = [False,save_folder_1+'population_data/',save_folder_1+'population_20_min_data/',save_folder_1+'single_cell_silences_data/',save_folder_2+'new_cdc_data/',save_folder_3 + 'moving_bc_data/']
#labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','population_20_min','single_cell_silences','single_cell_silences_new_bc','single_cell_silences_moving_bc']
hypothesis_list_df_folder = [False,save_folder_1+'population_data/',save_folder_4+'data/']
labels_hypothesis_list_df =[ 'experiment_an_WT_ESL','population','shuffle']


colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))

box_total = []
box_isolated = []
box_consecutive = []
for i,df_consecutive_folder in enumerate(hypothesis_list_df_folder):
    label = labels_hypothesis_list_df[i]
    
    if df_consecutive_folder is False:
        df_consecutive = conc_data[c]
        
        consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
        cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
        isolated_N = cell_population_list[0]
        
        box_isolated.append(cell_population_list[0]/isolated_N)

        consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
        box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
        total_N = box_plot_consecutive_cumulative[0]
        box_total.append([box_plot_consecutive_cumulative[0]/total_N])
        
        total_consecutive = box_plot_consecutive_cumulative[1]
        box_consecutive.append([box_plot_consecutive_cumulative[1]]/total_consecutive)
        
    
    else:
        cell_population_list_trials = [];box_plot_consecutive_cumulative_trials = []; box_plot_consecutive_cumulative_trials
        for df_file_name in os.listdir(df_consecutive_folder) :
        
            df_consecutive = download_data(df_consecutive_folder + df_file_name)
            
            consecutive_non_cumulative_obj = consecutive_non_cumulative(df_consecutive)
            cell_population_list = consecutive_non_cumulative_obj.get_consecutive_trains_of_pulses()
            cell_population_list_trials.append(cell_population_list)
            
            consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
            box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
            box_plot_consecutive_cumulative_trials.append(box_plot_consecutive_cumulative)
            
        box_isolated.append([i/isolated_N for i in first_element_box(cell_population_list_trials)])
        box_total.append([i/total_N for i in first_element_box(box_plot_consecutive_cumulative_trials)])
        box_consecutive.append([i/total_consecutive for i in second_element_box(box_plot_consecutive_cumulative_trials)])
        
 #%%
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=3, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs_main[0,0])
ax2 = plt.subplot(gs_main[1,0])
ax3 = plt.subplot(gs_main[2,0])


X1 = [np.ones(len(box_total[i]))*(i+1) for i in range(0,len(box_total))]
bp1 = ax1.boxplot(box_total,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp1['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp1['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp1['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp1['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp1['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X1)):
    xA = np.random.normal(0, 0.1, len(box_total[i])), 
    ax1.scatter(xA+X1[i],box_total[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

ax1.tick_params(axis='x', labelsize=8,length=2); 
ax1.tick_params(axis='y', labelsize=8,length=2)
ax1.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax1.set_xlabel('total pulses',fontsize=8)
ax1.set_ylabel('counts / experimental value',fontsize=8)
ax1.xaxis.set_label_coords(0.5, -0.12);ax1.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax1.get_yticks(); ax1.set_ylim([-0.05,1.05]); ax1.set_yticks([0,1])
#ax1.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)

box_isolated[0] = [1]
X2 = [np.ones(len(box_isolated[i]))*(i+1) for i in range(0,len(box_isolated))]
bp2 = ax2.boxplot(box_isolated,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp2['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp2['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp2['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp2['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp2['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X2)):
    xA = np.random.normal(0, 0.1, len(box_isolated[i])), 
    ax2.scatter(xA+X2[i],box_isolated[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
           
        
ax2.tick_params(axis='x', labelsize=8,length=2); 
ax2.tick_params(axis='y', labelsize=8,length=2)
ax2.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax2.set_xlabel('isolated pulses',fontsize=8)
ax2.set_ylabel('counts / experimental value',fontsize=8)
ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
#ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)


box_consecutive[0] = [1]
X3 = [np.ones(len(box_consecutive[i]))*(i+1) for i in range(0,len(box_consecutive))]
bp3 = ax3.boxplot(box_consecutive,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp3['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp3['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp3['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp3['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp3['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X3)):
    xA = np.random.normal(0, 0.1, len(box_consecutive[i])), 
    ax3.scatter(xA+X2[i],box_consecutive[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
           
        
ax3.tick_params(axis='x', labelsize=8,length=2); 
ax3.tick_params(axis='y', labelsize=8,length=2)
ax3.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax3.set_xlabel('consecutive pulses',fontsize=8)
ax3.set_ylabel('counts / experimental value',fontsize=8)
ax3.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
#ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)



plt.savefig(save_folder_4 + 'boxplot.pdf', format='pdf')

#%% est plotea en horizontal y algunas condiciones!


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
gs_main = gridspec.GridSpec(nrows=3, ncols=3, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs_main[0,0])
ax2 = plt.subplot(gs_main[0,1])
ax3 = plt.subplot(gs_main[0,2])

box_total_ = [box_total[0],box_total[-1]]

X1 = [np.ones(len(box_total_[i]))*(i+1) for i in range(0,len(box_total_))]
bp1 = ax1.boxplot(box_total_,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp1['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp1['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp1['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp1['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp1['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X1)):
    xA = np.random.normal(0, 0.1, len(box_total_[i])), 
    ax1.scatter(xA+X1[i],box_total_[i], alpha=1,s = 0.6,color='black',edgecolors='black',linewidths=0.0)

ax1.tick_params(axis='x', labelsize=8,length=2); 
ax1.tick_params(axis='y', labelsize=8,length=2)
ax1.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax1.set_xlabel('total pulses',fontsize=8)
ax1.set_ylabel('counts / experimental value',fontsize=8)
ax1.xaxis.set_label_coords(0.5, -0.12);ax1.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax1.get_yticks(); ax1.set_ylim([-0.05,1.05]); ax1.set_yticks([0,1])
#ax1.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.set_ylim([0.6,1.8])

box_isolated_ = [box_isolated[0],box_isolated[-1]]
X2 = [np.ones(len(box_isolated_[i]))*(i+1) for i in range(0,len(box_isolated_))]
bp2 = ax2.boxplot(box_isolated_,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp2['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp2['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp2['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp2['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp2['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X2)):
    xA = np.random.normal(0, 0.1, len(box_isolated_[i])), 
    ax2.scatter(xA+X2[i],box_isolated_[i], alpha=1,s = 0.6,color='black',edgecolors='black',linewidths=0.0)
           
        
ax2.tick_params(axis='x', labelsize=8,length=2); 
ax2.tick_params(axis='y', labelsize=8,length=2)
ax2.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax2.set_xlabel('isolated pulses',fontsize=8)
ax2.set_ylabel('counts / experimental value',fontsize=8)
ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
#ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax2.set_ylim([0.6,1.8])


box_consecutive_ = [box_consecutive[0],box_consecutive[-1]]
X3 = [np.ones(len(box_consecutive_[i]))*(i+1) for i in range(0,len(box_consecutive_))]
bp3 = ax3.boxplot(box_consecutive_,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp3['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp3['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp3['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp3['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp3['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians

for i in range(len(X3)):
    xA = np.random.normal(0, 0.1, len(box_consecutive_[i])), 
    ax3.scatter(xA+X2[i],box_consecutive_[i], alpha=1,s = 0.6,color='black',edgecolors='black',linewidths=0.0)
           
        
ax3.tick_params(axis='x', labelsize=8,length=2); 
ax3.tick_params(axis='y', labelsize=8,length=2)
ax3.set_xticklabels(labels_hypothesis_list_df,rotation = 0)
ax3.set_xlabel('consecutive pulses',fontsize=8)
ax3.set_ylabel('counts / experimental value',fontsize=8)
ax3.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
#xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
#ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje
ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax3.set_ylim([0.6,1.8])



plt.savefig(save_folder_4 + 'boxplot_aux.pdf', format='pdf')