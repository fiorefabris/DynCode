'''Threshold analysis pipeline for the FGF newest dataset :)'''


import pickle
import sys
#sys.path.insert(0, '/home/fiore/Documents/DO2019/analysis_2019_E8_2_WTonly/')
#from __Th_Analysis import  get_heatmap

def download_data(filename):
    infile =  open(filename,'rb')
    results = pickle.load(infile)
    infile.close()
    return(results)

#%%
##########################################
### CHARGING THE FILES ###################
##########################################

data_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/files/'
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/'
dataset_name='2019_E8_6'

filename = data_folder + dataset_name + '.pkl'
conc_df = download_data(filename)

filename = data_folder + 'df_results_20ng_Norm_.pkl'
df_result_20ng = download_data(filename)

filename = data_folder + 'df_results_5ng_Norm_.pkl'
df_result_5ng = download_data(filename)

filename = data_folder + 'df_results_2_5ng_Norm_.pkl'
df_result_2_5ng = download_data(filename)

filename = data_folder + 'df_results_0ng_Norm_.pkl'
df_result_NC = download_data(filename)

results = [df_result_NC,df_result_2_5ng,df_result_5ng,df_result_20ng]
conditions = ['NC','2_5ng','5ng','20ng']

#%%
results = [df_result_NC,df_result_2_5ng,df_result_5ng,df_result_20ng]
conditions = ['NC','2_5ng','5ng','20ng']
for (condition,df_result) in zip(conditions,results):
    get_heatmap(df_result,condition,save=True,save_folder=save_folder)

#%%

lcurves_values = df_result_NC.values #curvas de nivel
lcurves_condition = 'NC'


for (heatmap_condition,heatmap_values) in zip(conditions,results):
    get_level_curves(heatmap_values,heatmap_condition,lcurves_values,
                 lcurves_condition,save=True,
                 save_folder=save_folder,
                 plot_lcurves=False,save_lcurves=True)


    