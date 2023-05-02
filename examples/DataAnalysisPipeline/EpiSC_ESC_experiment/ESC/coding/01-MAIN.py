from DataAnalysis import TimeSerieDC


from DataAnalysis.Preprocesing.PlotData import silent_ax,set_scale
from DataAnalysis.Preprocesing.DownloadData import drop_cells, clean_cells
from DataAnalysis.Detrending import Smoothing 
from DataAnalysis.FindPeaks import Threshold_analysis
from DataAnalysis.Detrending.MinPolDT import Search_min , Fit_pol_min
from DataAnalysis.FindPeaks.Find_peaks_v2 import Full_Find_Peaks

import pickle
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
#%%
data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/raw_data/'
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/'
dataset_name='ESC'


conc = TimeSerieDC.TimeSerieDC(data_folder,dataset_name,save_folder,clean_last_points=False)

#%%
conc.traze_value = 'MEAN_INTENSITY'
ylim = [58000,65000]; var_th = 4 ;var_ylim = [1155,5686];save = False
conc.set_plot_values(ylim = ylim, var_th = var_th, var_ylim = var_ylim,long_plot = 0.5,cols_plot=5)
conc.plot_data(save=True,pointplot = 0,description_label = False)


#%% CLEANING
clean = {'an_ESC_FAX_PD03':[],'an_ESC_FAX':[1]} #[0,1,5] ,[1,5]
conc.data = clean_cells(conc.data,clean)

#drop = {'an_WT_ESL_PD03':[],'an_WT_ESL':[]}
#conc.data = drop_cells(conc.data,drop)

conc.get_values()
conc.plot_data(save=True,pointplot = False,description_label = False)

#%% SMOOTHING (FALTA FRECUENCIA)

conc.traze_value = 'MEAN_INTENSITY'
Smoothing.smooth_signal(conc.data,conc.traze_value,3)
conc.traze_value = 'sm_MEAN_INTENSITY'
conc.plot_data(save=True,pointplot = False,description_label = False)


#%% PEAK FINDING

conc.filter_name = 'sm'
conc.traze_value = 'sm_MEAN_INTENSITY'
conc.dist_bet_peaks_order = 1
conc.order_mins=1

conc.data, conc.peaks = Full_Find_Peaks(conc.data, conc.dist_bet_peaks_order,conc.traze_value ) 
conc.data = Search_min(conc.data, [conc.order_mins],traze_value = 'sm_MEAN_INTENSITY' )
conc.plot_peaks_v2(save=True,plot_mins = True)

#%%


filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/ESC.pkl'
outfile = open(filename,'wb')
pickle.dump(conc.data, outfile)
outfile.close()
