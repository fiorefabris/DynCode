
from DataAnalysis import TimeSerieDC

data_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Files/'
data_folder ='/home/fabris/Documents/Dyncode/low_res_datasets/2018_E1_28/Data_Files/'
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Sup_Figs/'
save_folder ='/home/fabris/Documents/Dyncode/low_res_datasets/2018_E1_28/Sup_Figs/'

dataset_name='2018_E1_28'

#plotea raw data 
conc = TimeSerieDC.TimeSerieDC(data_folder,dataset_name,save_folder,clean_last_points = False)
conc_order = ['an_0ng','an_2_5ng','an_5ng','an_10ng','an_20ng','an_WT_ESL','an_WT_N2chHep']
conc_order = ['an_WT_N2chHep','an_WT_ESL','an_0ng']

conc.def_conc_order(conc_order)
conc.get_values();

#ylim = [24962,63899] # conc.min_int, conc.max_int
ylim = [25000,70000]
var_th = 700
var_ylim = [649,5904]  # conc.min_var, conc.max_var

save = True
conc.set_plot_values(ylim = ylim, var_th = var_th, var_ylim = var_ylim,long_plot = 0.5,cols_plot=5)
conc.plot_data(save=save)



