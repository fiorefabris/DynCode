from DataAnalysis.Statistics.PlotStatistics import amplitud_histograms,plot_statistics
from DataAnalysis.Statistics.Statistics import percentiles,make_statistics
from DataAnalysis.Panels.panel_detrending import plot_detrend_panel
from DataAnalysis.Panels.panel_find_peaks import plot_find_peaks_panel
import pandas as pd
import numpy as np


def check_compatibility(MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
    aux = []
    for c in MinDT_TimeSerieDC.data:
        aux.append(
            (MinDT_TimeSerieDC.data[c].MEAN_INTENSITY == BFiltDT_TimeSerieDC.data[c].MEAN_INTENSITY).all())
    return (bool(np.prod(aux)))



class statistics():

    def set_detrend_data(self, MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
        conc = {}
        for c in MinDT_TimeSerieDC.data:
            df = MinDT_TimeSerieDC.data[c]
            df['MEAN_INTENSITY_BP'] = BFiltDT_TimeSerieDC.data[c].MEAN_INTENSITY_BP
            conc[c] = df
        self.data = conc

    def __init__(self, MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
        if check_compatibility(MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
            self.set_detrend_data(MinDT_TimeSerieDC,BFiltDT_TimeSerieDC)
            self.range_values = None
            self.doc = None
            self.percentiles = None

            self.save_folder = MinDT_TimeSerieDC.save_folder
            self.dataset_name = MinDT_TimeSerieDC.dataset_name
            self.conc_order = MinDT_TimeSerieDC.conc_order
            self.len_cells = MinDT_TimeSerieDC.len_cells
            self.traze_max_values = MinDT_TimeSerieDC.traze_max_values

        else:
            print('ERROR: no compatibility between variables')


    def set_plot_values(self, **kwargs):
        """Sets the limits of the plots
               args:
                   range_values: [(ylow_raw,yhigh_raw),(ylow_mindt,yhigh_mindt),(ylow_butterdt,yhigh_butterdt)] """
        self.range_values = kwargs.get('range_values',self.range_values)

    def histograms(self,cumulative=False,save=False):
        """Plot the amplitude histograms of the raw & detrended data

            save_name: dataset_name + _raw_data"""

        if cumulative:
            save_name = self.dataset_name + 'cum_amp_histograms'
        else:
            save_name = self.dataset_name + '_amp_histograms'

        description = 'Histograms of amplitudes. Range values = ' + str(self.range_values)
        amplitud_histograms(self.data ,description ,cumulative,save ,'MEAN_INTENSITY' ,'DT_MEAN_INTENSITY' ,'MEAN_INTENSITY_BP',
                            range_values = self.range_values,save_folder= self.save_folder ,save_name = save_name)

    def set_statistics(self,save=False):
        """ Calculates the statistics of the amplitud values.
       if save=True, save the results on the data_chart.
        """

        self.doc = make_statistics(self.data ,save ,'MEAN_INTENSITY' ,'DT_MEAN_INTENSITY' ,'MEAN_INTENSITY_BP',
                                   save_folder = self.save_folder,save_name=self.dataset_name)

        self.percentiles = percentiles(self.data,save ,'MEAN_INTENSITY' ,'DT_MEAN_INTENSITY','MEAN_INTENSITY_BP',
                                       save_folder=self.save_folder,save_name=self.dataset_name)

    def plot_statistics(self,save=False):

        save_name = self.dataset_name + '_statistics'
        description = ' '
        plot_statistics(self.doc,self.percentiles ,description ,save, save_folder = self.save_folder, save_name = save_name)

class panel_plots():
    def set_detrend_data(self, MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
        conc = {}
        keep_cols = ['ID', 'TRACK_ID', 'FRAME', 'CELL_DIVISION_TIME', 'MEAN_INTENSITY',
       'STANDARD_DEVIATION','FITM_O' + MinDT_TimeSerieDC.filter_value,'DT_MEAN_INTENSITY','BP_Filt',
       'MEAN_INTENSITY_BP']
        drop_cols=['ID', 'TRACK_ID', 'FRAME', 'CELL_DIVISION_TIME', 'MEAN_INTENSITY',
       'STANDARD_DEVIATION']
        for c in MinDT_TimeSerieDC.data:
            df = pd.concat([MinDT_TimeSerieDC.data[c],BFiltDT_TimeSerieDC.data[c].drop(drop_cols,axis=1,inplace=False)],axis=1)
            conc[c] = df
        self.data = conc

    def __init__(self,cells_conc,conc_labels,filter_labels, MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
        if check_compatibility(MinDT_TimeSerieDC, BFiltDT_TimeSerieDC):
            self.set_detrend_data(MinDT_TimeSerieDC,BFiltDT_TimeSerieDC)

            self.save_folder = MinDT_TimeSerieDC.save_folder
            self.dataset_name = MinDT_TimeSerieDC.dataset_name
            self.conc_order = MinDT_TimeSerieDC.conc_order
            self.len_cells = MinDT_TimeSerieDC.len_cells
            self.traze_max_values = MinDT_TimeSerieDC.traze_max_values

            self.cells_conc = cells_conc
            self.conc_labels = conc_labels
            self.filter_labels = filter_labels
            self.plot_value_filter = ['FITM_O' + MinDT_TimeSerieDC.filter_value,'BP_Filt']
            self.traze_values = [MinDT_TimeSerieDC.traze_value,BFiltDT_TimeSerieDC.traze_value]
            self.ylim = None

            self.peaks = [MinDT_TimeSerieDC.peaks, BFiltDT_TimeSerieDC.peaks]

    def set_ylim(self,ylim_conc,ylim_DT,ylim_BP):
        #they must be three lists
        self.ylim = [ylim_conc,ylim_DT,ylim_BP]

    def detrend(self,save=False):
        #poner un error si no tenes ylim
        description = None
        save_name = self.dataset_name + 'detrend_panel'
        plot_detrend_panel(self.data,self.cells_conc,self.conc_labels,self.plot_value_filter,self.traze_values , self.traze_max_values, self.ylim,
                           'MEAN_INTENSITY', save=save, save_folder= self.save_folder, save_name=save_name)


    def find_peaks(self,th_values_dt,Thresholds,order,save=False):
        """

        :param th_values_dt: the thresholds we want to plot
        :param Thresholds: The test thresholds
        :param order: the order of the find peaks
        :return: find_peaks panel
        """
        test_conditions = ["an_0ng", "an_WT_ESL"]
        test_conc_labels = [self.conc_labels[c] for c in test_conditions]
        save_name = self.dataset_name + 'find_peaks_panel'

        plot_find_peaks_panel(self.data, self.peaks, test_conditions, self.traze_values, th_values_dt, Thresholds,
     order, self.cells_conc, test_conc_labels,self.filter_labels ,self.traze_values , self.traze_max_values, self.ylim, save=save,
                              save_folder=self.save_folder, save_name=save_name)