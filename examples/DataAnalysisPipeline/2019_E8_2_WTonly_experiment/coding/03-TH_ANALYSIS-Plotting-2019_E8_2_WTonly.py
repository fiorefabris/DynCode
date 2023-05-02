import pickle5 as pickle
import pandas as pd

#%%
##########################################
### CHARGING THE FILES ###################
##########################################

data_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'
dataset_name='2019_E8_2_WTonly'

conc_df = pd.read_pickle(data_folder + dataset_name + '.pkl')

filename = data_folder + 'df_results_WT_Norm_.pkl'
infile =  open(filename,'rb')
df_result_WT =pickle.load(infile)
infile.close()

filename = data_folder + 'df_results_CN_Norm_.pkl'
infile =  open(filename,'rb')
df_result_CN = pickle.load(infile)
infile.close()

#%%

condition = 'NC';df_result = df_result_CN
get_heatmap(df_result,condition,save=True,save_folder=save_folder)

condition = 'WT';df_result = df_result_WT
get_heatmap(df_result,condition,save=True,save_folder=save_folder)

#%%
heatmap_values = df_result_WT#heatmap
heatmap_condition = 'WT'

lcurves_values = df_result_CN.values #curvas de nivel
lcurves_condition = 'CN'

get_level_curves(heatmap_values,heatmap_condition,lcurves_values,
                 lcurves_condition,save=True,
                 save_folder=save_folder,
                 plot_lcurves=False,save_lcurves=True)

#%%
heatmap_values = df_result_CN #heatmap
heatmap_condition = 'CN'

gradient_values = df_result_CN #heatmap
gradient_condition = 'CN'

get_level_curves(heatmap_values,heatmap_condition,lcurves_values,
                 lcurves_condition,save=True,
                 save_folder=save_folder,
                 plot_lcurves=False,save_lcurves=True)

#%%
