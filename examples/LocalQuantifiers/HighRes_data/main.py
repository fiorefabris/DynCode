from LocalAnalysis import Quantifiers,Plot_Quantifiers


#conc = conc_data.data
traze_value = 'MEAN_INTENSITY'
local_time_window_len = 30
delay = 1

for local_time_window_len in [5,10,20,50,100]: # 30
    Quantifiers.Calculate_Quantifiers(conc,traze_value,local_time_window_len,delay)
#%%

def get_cell_values(conc):
    len_cells,M = [],[]
    
    for c in conc:
        df = conc[c]
        df.sort_index(inplace=True)
        cells = df.index.get_level_values(0).unique()
        len_cells.append(len(cells))
        M.append(df.index.get_level_values(1).max())
    return(len_cells,M)

#%%
len_cells,M=get_cell_values(conc)
ylim = [58000,65000]
list_local_quantifiers = ['N30_PE','N30_F','N30_C','N30_LOCAL_VAR']
local_quantifier = 'N30_PE'
list_ylim = [[-0.15,0.15],[-0.02,0.6],[-3,0.2],[-22000,70000]]

list_ylim = [[-0.10,0.20],[-0.02,0.6],[-20000,100000],[0,30]]
list_ylim = [[-0.10,0.4],[-0.02,1.2]]
list_ylim = [[-0.10,1],[-0.02,2.4]]

values_list = ['_PE','_F','_LOCAL_VAR','_CV2']
values_list = ['_PE']
local_time_window_len_list = [5,10,20,30,50,100]
local_time_window_len_list = [5,10]
local_time_window_len_list = [5]

#%%
for local_time_window_len in local_time_window_len_list: # 30
    for z,local_quantifier in enumerate(values_list):
        local_quantifier_ ='N'+str(local_time_window_len) +local_quantifier
        print(local_quantifier_,'------------------------------------------------------------------------')
        save_folder_name = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/local_quantifiers/ESC_'+local_quantifier_
        Plot_Quantifiers.plot_local_value(conc,traze_value,local_quantifier_,len_cells, np.max(M),ylim ,list_ylim[z],save=True,save_folder_name=save_folder_name)
        
#%%

plot_value_list = ['_N50_PE','_N50_F','_N50_LOCAL_VAR','_N50_CV2']
list_ylim = [[0,0.20],[0,0.6],[0,100000],[0,3000]]
save_folder_name = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/local_quantifiers/ESC_N50NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
Plot_Quantifiers.plot_local_value_all(conc,traze_value,plot_value_list, len_cells, np.max(M), ylim ,list_ylim,save=True,save_folder_name=save_folder_name)
