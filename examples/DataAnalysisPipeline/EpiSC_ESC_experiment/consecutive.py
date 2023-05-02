import sys
import os
sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from consecutive_main import consecutive_non_cumulative,makeColours,consecutive_cumulative,mean_consecutive_value,filter_activity_cells#,first_element_box,second_element_box
from exponential_dm_main import download_data
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


save_folder ='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/EpiSC/figures/'
#%%

#Esto es para plotear los datos de consecutividad
save_folder1 =  '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/'
save_folder2 =  '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/EpiSC/data/'
hypothesis_list_df = [download_data(save_folder1 + 'conc_data.pkl')['an_ESC_FAX'],download_data(save_folder2 + 'conc_data.pkl')['an_EpiSC_FAX']]
        
labels_hypothesis_list_df =[ 'ESL','EpiSL']
colors =  cm.get_cmap('viridis', len(labels_hypothesis_list_df))(np.linspace(0,1,len(labels_hypothesis_list_df)))
colors = ['g','r']


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
plt.savefig(save_folder + 'consecutive_pulses_.pdf', format='pdf')