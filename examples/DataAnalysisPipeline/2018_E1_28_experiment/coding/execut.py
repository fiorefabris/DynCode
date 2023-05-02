# =============================================================================
# Pipeline Analysis for Dataset 2018_E1_28
# =============================================================================

from DataAnalysis import TimeSerieDC

data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/low_res_datasets/2018_E1_28/Data_Files/'
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/low_res_datasets/2018_E1_28/Sup_Figs/'
dataset_name='2018_E1_28'

conc = TimeSerieDC.TimeSerieDC(data_folder,dataset_name,save_folder)
conc_order = ['an_0ng','an_2_5ng','an_5ng','an_10ng','an_20ng','an_WT_ESL','an_WT_N2chHep']
conc_order = ['an_WT_ESL','an_WT_N2chHep','an_0ng']

conc.def_conc_order(conc_order)

ylim = [24962,63899] # conc.min_int, conc.max_int
var_th = 700
var_ylim = [649,5904]  # conc.min_var, conc.max_var

save = False
conc.set_plot_values(ylim = ylim, var_th = var_th, var_ylim = var_ylim,long_plot = 0.5,cols_plot=5)
conc.plot_data(save=save)
conc.plot_quality_tracking_test(save)

# =============================================================================
# Min Polynomial Detrend
# =============================================================================


filtmin_conc = TimeSerieDC.MinDT_TimeSerieDC(data_folder,dataset_name,save_folder,clean_last_points=True)
filtmin_conc.def_conc_order(conc_order)
filtmin_conc.set_plot_values(ylim = ylim, var_th = var_th, var_ylim = var_ylim,long_plot = 1,cols_plot=5)

orders_min = [2] #[2,3,4] 	 
delta = 3
filtmin_conc.min_detrend(orders_min,delta,save)
filter_value = 2
filtmin_conc.set_filter_value(filter_value)
#filtmin_conc.plot_pol_order() #mejorar si es necesario
filtmin_conc.plot_min_detrend(save=save)
 

##

ylim = [-4000,4000] #[-3100,8300]
filtmin_conc.set_plot_values(ylim=ylim)
filtmin_conc.plot_min_detrend_data(save=save)

# =============================================================================
# #Butter_Filter detrend
# =============================================================================

buterfilt_conc = TimeSerieDC.BFiltDT_TimeSerieDC(data_folder,dataset_name,save_folder)

buterfilt_conc.def_conc_order(conc_order)
buterfilt_conc.bfilt_detrend(0.025,0.6)
ylim = [24962,63899] 

buterfilt_conc.set_plot_values(ylim=ylim, long_plot = 2,cols_plot=7)
buterfilt_conc.plot_bpfilt_detrend(save=save)

##

ylim = [-4000,4000] 
buterfilt_conc.set_plot_values(ylim=ylim,long_plot = 2,cols_plot=7)
buterfilt_conc.plot_bpfilt_detrended_data(save=save)

# =============================================================================
# #Statistics
# =============================================================================
from DataAnalysis import statistics

stats=statistics.statistics(filtmin_conc,buterfilt_conc)
range_values = [(40000,63899),(-4000,4000),(-4000,4000)]
stats.set_plot_values(range_values=range_values)
stats.histograms(False,True)
stats.histograms(True,True)

stats.set_statistics(True)
stats.plot_statistics(True)

# =============================================================================
# #Find_Peaks
# =============================================================================

# #Find_Peaks Threshold set up
import numpy as np
from DataAnalysis.FindPeaks import Threshold_analysis
Thresholds = np.arange(0,3000,100) 
Thresholds = np.array([2500,2400])


orders = [2]#orders = [2, 3, 4]
filtmin_conc.find_peaks(Thresholds,orders)
buterfilt_conc.find_peaks(Thresholds,orders)

th_parameter = 0.01
test_conditions = ["an_0ng","an_WT_ESL"]
save = True
save_name = filtmin_conc.dataset_name + '_threshold_analysis_'+str(th_parameter)
save_folder = filtmin_conc.save_folder
description = '(up) Density of peaks as a function of the Threshold \n (down) The 0ng/WT ratio  & WT-0ng difference of peaks densities as a function of the Threshold' 

dic_filtmin = {'an_WT_ESL':filtmin_conc.peaks['an_WT_ESL'],'an_0ng':filtmin_conc.peaks['an_0ng']}
dic_buterfilt = {'an_WT_ESL':buterfilt_conc.peaks['an_WT_ESL'],'an_0ng':buterfilt_conc.peaks['an_0ng']}



Threshold_analysis.plot_threshold_analysis([dic_filtmin,dic_buterfilt],
                                           [filtmin_conc.traze_value,buterfilt_conc.traze_value]
               ,Thresholds,th_parameter,2,test_conditions,save,save_name = save_name,save_folder = save_folder,description = description)
   
    
##end
filtmin_conc.save_folder = '/home/Fiore/Documents/Documents/Trabajo/Doctorado/Datasets/2018_E1_28/Data_Analysis/Th_0.01/' + str(th_parameter) + '_'
buterfilt_conc.save_folder = '/home/Fiore/Documents/Documents/Trabajo/Doctorado/Datasets/2018_E1_28/Data_Analysis/Th_0.01/' + str(th_parameter) + '_'

filtmin_conc.distance_between_peaks(2500,2)
filtmin_conc.plot_peaks(True)
filtmin_conc.plot_time_bet_peaks_hist((5.25,29.75),num_bins = np.arange(5.25,29.75,1.75),cumulative = False,save=False)
filtmin_conc.plot_time_bet_peaks_hist_all((5.25,29.75),num_bins = np.arange(5.25,29.75,1.75),cumulative = False,save=False)
filtmin_conc.plot_peak_statistics(4000,(5.25,29.75),num_bins = np.arange(5.25,29.75,1.75),save=False)

buterfilt_conc.distance_between_peaks(2400,2)
buterfilt_conc.plot_peaks(True)
buterfilt_conc.plot_time_bet_peaks_hist((4.375,30.625),num_bins = np.arange(4.375,30.625,1.75),cumulative = False,save=True)
buterfilt_conc.plot_time_bet_peaks_hist_all((4.375,30.625),num_bins = np.arange(4.375,30.625,1.75),cumulative = False,save=True)
buterfilt_conc.plot_peak_statistics(4000,(4.375,30.625),num_bins = np.arange(4.375,30.625,1.75),save=True)


# =============================================================================
# #Quantifiers (running finding peaks is not necesary)
# =============================================================================


values = [conc,filtmin_conc,buterfilt_conc]
conc.traze_value = 'MEAN_INTENSITY'
for value in values:
    value.set_local_quantifiers(30,1)
    value.set_local_quantifiers(15,1)
    value.save_folder= '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Quantifiers/'
    for plot_value in ['PE','F','C','LOCAL_VAR']:
        value.plot_local_value(plot_value,save=True)


