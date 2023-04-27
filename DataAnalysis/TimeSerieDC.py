import matplotlib.pyplot as plt
import numpy as np

from DataAnalysis.Preprocesing.DownloadData import download_data
from DataAnalysis.Preprocesing.QualityTrackingTest import quality_tracking_test
from DataAnalysis.Preprocesing.PlotData import grid_plot
import warnings

from DataAnalysis.Detrending.MinPolDT import Search_min,Fit_pol_min
from DataAnalysis.Detrending.PlotDetrend import plot_detrend,heatmap_plot

from DataAnalysis.Detrending.ButterFilterDetrend import butter_filter

from DataAnalysis.FindPeaks.FindPeaks import Full_Find_Peaks
from DataAnalysis.FindPeaks.FindPeaks import Full_Dist_Between_Peaks
from DataAnalysis.FindPeaks.Plot_FindPeaks import peaks_grid_plot, plot_time_between_peaks_st,plot_time_between_peaks_st2,plot_time_between_peaks_hist, plot_dist_bet_peaks_hist_all


from LocalAnalysis.Quantifiers import Calculate_Quantifiers
from LocalAnalysis.Plot_Quantifiers import plot_local_value
from DataAnalysis.FindPeaks.Find_peaks_v2 import peaks_grid_plot_v2

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("talk")


class TimeSerieDC:

    def get_values(self):
        '''Calculates the amount of cells for each concentration & the maximum value of the intensity for each condition
        '''
        len_cells,M = [],[]

        for c in self.data:
            df = self.data[c]
            df.sort_index(inplace=True)
            cells = df.index.get_level_values(0).unique()
            len_cells.append(len(cells))
            M.append(df.index.get_level_values(1).max())

        self.len_cells = len_cells
        self.traze_max_values = M

        # Calculates the minima and maxima of MEAN_INTENSITY column
        aux_min,aux_max = [],[]
        for c in self.data:
            df = self.data[c]
            df.sort_index(inplace=True)
            aux_min.append(df.MEAN_INTENSITY.min())
            aux_max.append(df.MEAN_INTENSITY.max())
        self.min_int = np.min(aux_min)
        self.max_int = np.max(aux_max)

        # Calculates the minima and maxima of STANDARD_DEVIATION column
        aux_min,aux_max = [],[]
        for c in self.data:
            df = self.data[c]
            df.sort_index(inplace=True)
            aux_min.append(df.STANDARD_DEVIATION.min())
            aux_max.append(df.STANDARD_DEVIATION.max())
        self.min_var = np.min(aux_min)
        self.max_var = np.max(aux_max)

        #Calculates the minima and maxima of the filter value column
        if self.filter_name == 'min_detrend':
            aux_min, aux_max = [], []
            for c in self.data:
                df = self.data[c]
                df.sort_index(inplace=True)
                aux_min.append(df.DT_MEAN_INTENSITY.min())
                aux_max.append(df.DT_MEAN_INTENSITY.max())
            self.min_int_filt = np.min(aux_min)
            self.max_int_filt = np.max(aux_max)
        elif self.filter_name == 'butterfilt_detrend':
            aux_min, aux_max = [], []
            for c in self.data:
                df = self.data[c]
                df.sort_index(inplace=True)
                aux_min.append(df.MEAN_INTENSITY_BP.min())
                aux_max.append(df.MEAN_INTENSITY_BP.max())
            self.min_int_filt = np.min(aux_min)
            self.max_int_filt = np.max(aux_max)
        else:
            self.min_int_filt = self.min_int_filt
            self.max_int_filt = self.max_int_filt

    def __init__(self, data_folder,dataset_name,save_folder = None,clean_last_points = True):
        self.data = download_data(data_folder,True, save_folder=save_folder,save_name=dataset_name,clean_last_points = clean_last_points)
        self.conc_order = [c for c in self.data]
        self.len_cells = []
        self.traze_max_values = []
        self.min_int = None
        self.max_int = None
        self.min_var = None
        self.max_var = None
        self.min_int_filt = None
        self.max_int_filt = None
        self.dataset_name= dataset_name

        self.long_plot = 1
        self.cols_plot = 5
        self.filter_name = None


        self.save_folder = save_folder

        self.get_values()


        self.ylim = []
        self.var_ylim = []
        self.var_th = []

        self.peaks = None
        self.thresholds = None
        self.traze_value = 'MEAN_INTENSITY'
        self.dist_bet_peaks = None
        self.total_peaks = None
        self.dist_bet_peaks_th = None
        self.dist_bet_peaks_order = None

        self.permutations = None
        self.PE_ymax = None
        self.F_ymax = None
        self.C_ymax = None
        self.LocVAr_ymax = None

    def def_conc_order(self, conc_order):
        '''Order the dictionary conc with the order stablished on conc_order(list of strings)
        '''
        aux = {}
        for c in conc_order:
            aux[c] = self.data[c]
        self.data = aux
        self.conc_order = conc_order

    def set_plot_values(self,**kwargs):
        """Sets the limits of the plots
        args:
            ylim: [ylow,yhigh]
            var_th: y_th
            var_ylim: [var_ylow,var_yhigh] """

        self.ylim = kwargs.get('ylim',[])
        self.var_th = kwargs.get('var_th', [])
        self.var_ylim = kwargs.get('var_ylim', [])
        self.long_plot = kwargs.get('long_plot', self.long_plot)
        self.cols_plot = kwargs.get('cols_plot', self.cols_plot)

    def set_save_folder(self,save_folder):
        self.save_folder = save_folder

    def plot_data(self,save = False,pointplot = False,description_label=True):
        """Plot the raw data

            save_name: dataset_name + _raw_data"""
        if self.ylim == []:
            warnings.warn('ylim is not specified');return
        description = 'Time Series for different FGF-KO conditions. \n  y_lim = ' + str(self.ylim) + \
                      '\n  Concentration order = ' \
                      + str(self.conc_order)


        save_name = self.dataset_name + '_raw_data'


        grid_plot(self.data, self.len_cells, self.traze_max_values, self.ylim, self.traze_value,
                  description, pointplot = pointplot, save=save,long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name,
                  description_label = description_label)

    def plot_quality_tracking_test(self,save=False,description_label=False):
        """Plot the quality tracking test

            save_name: dataset_name + _quality_tracking_test"""

        if self.ylim == []:
            warnings.warn('ylim is not specified'); return
        if self.var_ylim == []:
            warnings.warn('var_ylim is not specified') ; return
        if self.var_th == []:
            warnings.warn('var_th is not specified') ; return

        description = 'Quality Tracking test for different FGF-KO conditions. \n y_lim = ' + str(
            self.ylim) + \
        '\n Concentration order = ' + str(self.conc_order)

        save_name= self.dataset_name + '_quality_tracking_test'

        quality_tracking_test(self.data, self.len_cells, self.traze_max_values, self.ylim,self.var_th,self.var_ylim, description,
                              save=save,long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name
                              ,description_label=description_label)


    def find_mins(self,orders):
        if self.traze_value != None:
            self.data= Search_min(self.data,orders, self.traze_value)
        else:
            print('ERROR: traze value is missing')

    ########

    def find_peaks(self,Thresholds,orders):
        if self.traze_value != None:
            self.data, self.peaks = Full_Find_Peaks(self.data, Thresholds,orders, self.traze_value)
            self.thresholds = Thresholds
        else:
            print('ERROR: traze value is missing')

    def plot_peaks(self,save=False,**kwargs):
        plot_mins = kwargs.get('plot_mins', False)
        order = kwargs.get('order_mins', 2)


        description = 'Peak finding for minima detrending \n ylim:' + str(self.ylim) + str(self.conc_order) + ' \n th: ' + str(
            self.dist_bet_peaks_th) + ' - time_window: ' + str(int(self.dist_bet_peaks_order) * 2 + 1)+ 'frames'

        save_name= self.dataset_name + self.filter_name + '_peaks_TS'

        peaks_grid_plot(self.data, self.peaks, self.len_cells, self.traze_max_values, self.ylim,
                        self.dist_bet_peaks_th, self.dist_bet_peaks_order, self.traze_value, save=save,
                        long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name, description=description,plot_mins = plot_mins,order_mins=order)

    def plot_peaks_v2(self,save=False,**kwargs):

        plot_mins = kwargs.get('plot_mins', True)
        #order = kwargs.get('order_mins', 2)
        save_name= self.dataset_name + self.filter_name + '_peaks_TS'
        peaks_grid_plot_v2(self.data, self.peaks, self.len_cells, self.traze_max_values, self.ylim
                        , self.dist_bet_peaks_order, self.traze_value,self.order_mins, save=save,
                        long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name,plot_mins = plot_mins)


    def distance_between_peaks(self,th,order):
        if self.peaks != None:
            self.dist_bet_peaks, self.total_peaks = Full_Dist_Between_Peaks(self.data, self.peaks, th, order, self.traze_value)
            self.dist_bet_peaks_th = th
            self.dist_bet_peaks_order = order
        else:
            print('ERROR: peaks are missing')

    def plot_peak_statistics(self,vmax,x_range,num_bins,save=False):
        save_name = self.dataset_name + self.filter_name + '_time_bet_peaks_st'
        plot_time_between_peaks_st(self.data, self.peaks, self.dist_bet_peaks, self.total_peaks, self.dist_bet_peaks_th,
                                   self.traze_value, self.dist_bet_peaks_order,
                                       vmax,x_range,num_bins, save=save, save_folder = self.save_folder, save_name=save_name)

        save_name = self.dataset_name + self.filter_name + '_time_bet_peaks_st2'
        plot_time_between_peaks_st2(self.data, self.peaks, self.dist_bet_peaks, self.dist_bet_peaks_th
                                    ,self.traze_value, self.dist_bet_peaks_order,x_range,num_bins,save=save
                                    ,save_folder=self.save_folder,save_name=save_name)

    def plot_time_bet_peaks_hist(self,x_range,num_bins=7,cumulative=False,save=False):
        save_name = self.dataset_name + self.filter_name + '_time_bet_peaks_hist'

        plot_time_between_peaks_hist(self.data, self.peaks, self.dist_bet_peaks, self.total_peaks, self.dist_bet_peaks_th,
                                   self.traze_value, self.dist_bet_peaks_order,x_range,num_bins,cumulative=cumulative,save=save,save_folder = self.save_folder, save_name=save_name)

    def plot_time_bet_peaks_hist_all(self,x_range,num_bins=7,cumulative=False,save=False):
        save_name = self.dataset_name + self.filter_name + '_time_bet_peaks_hist_all'
        plot_dist_bet_peaks_hist_all(self.data,  self.dist_bet_peaks, self.total_peaks, self.dist_bet_peaks_th,
                                   self.traze_value, self.dist_bet_peaks_order,x_range,num_bins,cumulative = cumulative,save=save,save_folder = self.save_folder, save_name=save_name)

    #########
    # Quantifiers
    #########
    def get_quantifiers_values(self):
        PE_ymax, F_ymax, C_ymax, LocVAr_ymax = [], [],[], []

        for c in self.data:
            df = self.data[c]
            df.sort_index(inplace=True)
            PE_ymax.append(df[self.traze_value + '_PE'].max())
            F_ymax.append(df[self.traze_value + '_F'].max())
            C_ymax.append(df[self.traze_value + '_C'].max())
            LocVAr_ymax.append(df[self.traze_value + '_LOCAL_VAR'].max())

        self.PE_ymax = PE_ymax
        self.F_ymax = F_ymax
        self.C_ymax = C_ymax
        self.LocVAr_ymax = LocVAr_ymax

    def set_local_quantifiers(self,local_time_window_len,delay):
        if self.traze_value == None:
            self.traze_value = 'MEAN_INTENSITY'
        self.permutations,self.data = Calculate_Quantifiers(self.data,self.traze_value,local_time_window_len,delay)
        #self.get_quantifiers_values()


    def plot_local_value(self,plot_value,save=False):
        """
        Plots the local quantifiers.
        :param plot_value: PE, F, C or LOCAL_VAR
        :return:
        """
        #plot_value = str(plot_value)
        description = 'Quality Tracking test for different FGF-KO conditions. \n y_lim = ' + str(
            self.ylim) + \
                      '\n Concentration order = ' + str(self.conc_order)

        save_name = self.dataset_name + '_'+ self.traze_value + '_'+  plot_value+ '_TS'

        if plot_value == 'N30_PE':            plot_value_ylim = np.mean(self.PE_ymax)
        elif plot_value == 'N30_F':           plot_value_ylim = np.mean(self.F_ymax)
        elif plot_value == 'N30_C':           plot_value_ylim = np.mean(self.C_ymax)
        elif plot_value == 'N30_LOCAL_VAR':     plot_value_ylim = np.mean(self.LocVAr_ymax)

        plot_local_value(self.data, self.traze_value, plot_value,self.len_cells, self.traze_max_values, self.ylim,plot_value_ylim,
                         description, save=save, long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder,
                         save_name=save_name)



class MinDT_TimeSerieDC(TimeSerieDC):
    def __init__(self,data_folder,dataset_name, save_folder = None,clean_last_points=False):
        super().__init__(data_folder,dataset_name, save_folder,clean_last_points)
        self.delta = None
        self.orders_mins = None
        self.filter_value = None

        self.table = None
        self.filter_flag = False
        self.filter_name = 'min_detrend'

        self.traze_value = 'DT_MEAN_INTENSITY'

    def min_detrend(self,orders_min,delta,save=False,version=None):
        """Set the minima polynomial detrending.

        args:
            delta: int
            orders_mins: list

        return:
            self.table contains relevant information about the detrend"""

        self.delta = delta
        self.orders_mins = orders_min
        self.data = Search_min(self.data, orders_min)
        self.data, self.table = Fit_pol_min(self.data, self.delta, self.orders_mins,save, save_folder=self.save_folder,
                                            save_name= self.dataset_name, version = version)

        self.filter_flag = True

    def set_filter_value(self, filter_value):
        if self.filter_flag:
            self.filter_value = str(filter_value)
            plot_value_filter = 'FITM_O' + self.filter_value

            for c in self.data:
                self.data[c]['DT_MEAN_INTENSITY'] = (self.data[c]['MEAN_INTENSITY'] - self.data[c][plot_value_filter])
                self.data[c].sort_index(inplace=True)
        else:
            print('ERROR: No detrend information. Please Run min_detrend function.')

    def plot_pol_order(self):
        """Plot the order of each polynomial fit"""
        colors = sns.color_palette('colorblind')
        fig, axs = plt.subplots(len(self.data), sharex=True, sharey=True, figsize=(8.27, 11.69))
        k = 0

        for conc, data in self.table.groupby(level='conc'):
            ax = axs[k]
            ax.set_title(conc)
            ax.scatter(list(data.loc[conc][str(self.filter_value) +'_Pol_Fit_Order'].index),
            data[str(self.filter_value) +'_Pol_Fit_Order'], color=colors[k])
            k = k + 1

    def plot_min_detrend(self , save= False,description_label = True):
        """Plot the raw data & the detrending

            save_name: dataset_name + _min_detrend"""

        if self.ylim == []:
            warnings.warn('ylim is not specified');return
        if self.filter_value == None:
            warnings.warn('filter_value is not specified');return

        save_name = self.dataset_name + '_min_detrend'
        plot_value_filter = 'FITM_O' + self.filter_value

        description = 'Time Series for different FGF-KO conditions & the Polynomial Detrend. \n  y_lim = '+\
                      str(self.ylim) + '\n Concentration order = ' + str(self.conc_order) + ' \n delta = ' + str(self.delta)+\
                      ' Time window = ' + str(int(self.filter_value)*2+1) + ' frames'

        plot_detrend(self.data, self.len_cells, self.traze_max_values, 'MEAN_INTENSITY', plot_value_filter,
                     description, self.ylim, save=save,long=self.long_plot,cols=self.cols_plot,
                     save_folder= self.save_folder, save_name=save_name,description_label=description_label)

    def plot_min_detrend_data(self,traze_value = 'DT_MEAN_INTENSITY',save=False,pointplot=False):
        """Plot the detrended data with the polynomial local minima as detrend strategy

            save_name: dataset_name + _min_detrended_TS"""
        if self.ylim == []:
            warnings.warn('ylim is not specified');return
        if self.filter_value == None:
            warnings.warn('filter_value is not specified');return

        save_name = self.dataset_name + '_min_detrended_TS'


        description = 'Detrended time series using Polynomial local minima as the detrending strategy . \n  y_lim = ' +\
                      str(self.ylim) + ' \n Concentration order = ' + str(self.conc_order) + ' \n delta = ' + str(self.delta) +\
                      ' Time window = ' + str(int(self.filter_value)*2+1) + ' frames'

        grid_plot(self.data, self.len_cells, self.traze_max_values, self.ylim,traze_value, description,pointplot=pointplot,
                  save=save,long=self.long_plot,cols=self.cols_plot, save_folder= self.save_folder, save_name=save_name)

    def min_detrend_heatmap(self,save=False):
        """Plot the detrended data with the polynomial local minima as detrend strategy

            save_name: dataset_name + _min_detrended_TS_heatmap"""
        if self.ylim == []:
            warnings.warn('ylim is not specified');return
        if self.filter_value == None:
            warnings.warn('filter_value is not specified');return


        save_name = self.dataset_name + '_min_detrended_TS_heatmap'

        description = 'Heatmap of the detrended time series using Polynomial local minima as the detrending strategy. \n y_lim = ' \
                      + str(self.ylim) + ' Concentration order = ' + str(self.conc_order)  + ' \n delta = ' + str(self.delta) +\
                      ' Time window = '  + str(int(self.filter_value)*2+1) + ' frames'

        heatmap_plot(self.data, 'DT_MEAN_INTENSITY', self.ylim, description, save, save_folder=self.save_folder,
                     save_name=save_name)


class BFiltDT_TimeSerieDC(TimeSerieDC):

    def __init__(self,data_folder,dataset_name, save_folder = None,clean_last_points=False):
        super().__init__(data_folder,dataset_name, save_folder)
        self.n=4
        self.wcut = [0.1,0.6]

        self.traze_value = 'MEAN_INTENSITY_BP'
        self.filter_name = 'butterfilt_detrend'



    def bfilt_detrend(self,low = 0.1, high=0.6):
        self.data = butter_filter(self.data,low,high) #    print('Band Pass Filter  n = 4 \t wcut = [0.1,0.6]')

        for c in self.data:  # Estaria bueno agregarla en la funcion filtrado
            self.data[c]['BP_Filt'] = (self.data[c]['MEAN_INTENSITY'] - self.data[c]['MEAN_INTENSITY_BP'])
            self.data[c].sort_index(inplace=True)

    def plot_bpfilt_detrend(self, save=False):
        """Plot the raw data & the detrending

            save_name: dataset_name + _BPFilt_detrend"""

        if self.ylim == []:
            warnings.warn('ylim is not specified')
            return


        save_name = self.dataset_name + '_BPFilt_detrend'

        description = 'Time Series for different FGF-KO conditions & the Butter Filter Detrend. \n y_lim = ' + str(
            self.ylim) + ' \n Concentration order = ' + str(self.conc_order)

        plot_detrend(self.data, self.len_cells, self.traze_max_values, 'MEAN_INTENSITY', 'BP_Filt', description, self.ylim,
                     save=save,long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name)

    def plot_bpfilt_detrended_data(self,traze_value = 'MEAN_INTENSITY_BP',save=False):
        """Plot the detrended data with the Butter Band Pass Filter detrending strategy

            save_name: dataset_name + _BPFilt_detrended_TS"""
        if self.ylim == []:
            warnings.warn('ylim is not specified')
            return

        save_name = self.dataset_name +'_BPFilt_detrended_TS'

        description = 'Detrended time series using Butter Band Pass Filter as the detrending strategy. \n y_lim = ' + str(
            self.ylim) + ' \n Concentration order = ' + str(self.conc_order)

        grid_plot(self.data, self.len_cells, self.traze_max_values,self.ylim, traze_value, description,
                  save=save,long=self.long_plot,cols=self.cols_plot, save_folder=self.save_folder, save_name=save_name)

    def bpfilt_detrend_heatmap(self,save=False):
        """Plot the detrended data with the polynomial local minima as detrend strategy

                    save_name: dataset_name + _BPFilt_detrended_TS_heatmap"""

        if self.ylim == []:
            warnings.warn('ylim is not specified');return

        save_name = self.dataset_name + '_BPFilt_detrended_TS_heatmap'

        description = 'Heatmap of the detrended time series using Band Pass Filter as the detrending strategy. \n y_lim = ' \
                      + str(self.ylim) + ' Concentration order = ' + str(self.conc_order)

        heatmap_plot(self.data, 'MEAN_INTENSITY_BP', self.ylim, description, save, save_folder=self.save_folder,
                     save_name=save_name)



