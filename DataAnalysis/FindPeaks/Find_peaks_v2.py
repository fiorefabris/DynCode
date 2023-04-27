import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
from DataAnalysis.Preprocesing.PlotData import silent_ax


def Find_Peaks(Data_Peaks, order):
    ''' Searchs for local maxima over an Abs_Thr of Data of an given time window (Order).
    Returns
    - Data_Peaks : a Pandas Series with the values and positions of the local maxima
    - Number of peaks , time-lenght of Data & density of peaks
    '''

    #### antes habia otro!


    M = Data_Peaks.min()
    for Shift in range(1, order + 1):  # Llega hasta order
        Data_Peaks = Data_Peaks.where(
            (Data_Peaks > (Data_Peaks.shift(Shift)).fillna(M)) & (Data_Peaks > (Data_Peaks.shift(- Shift)).fillna(M)),M)

    Data_Peaks = Data_Peaks.mask(Data_Peaks == M)  # Data_Peaks = Data_Peaks.mask(Data_Peaks == m,None)
    return (Data_Peaks, Data_Peaks.count(), len(Data_Peaks), Data_Peaks.count() / len(Data_Peaks))


def Full_Find_Peaks(conc,order, traze_value):
    '''
    -Traze value : El nombre de la columna del Df que quiera analizar
    -order : orden del algoritmo de búsqueda

    RETURNS:
        - conc with peaks information
        - peaks_df: dictionary of df. Each df has the information of Number_Peaks, Frames, Dens_Peaks for each concentration.
    '''

    peaks_conc = {}

    for c in conc:
        df = conc[c]
        peaks_df = pd.DataFrame()

        auxNP, auxF, auxDP, auxCell = [], [], [], []
        aux_peaks = []

        for cells, data in df.groupby(level='cell'):
            Peaks, Number_Peaks, Frames, Dens_Peaks = Find_Peaks(data[traze_value], order)

            aux_peaks.extend(Peaks)
            auxNP.append(Number_Peaks)
            auxF.append(Frames)
            auxDP.append(Dens_Peaks)
            auxCell.append(cells)

        df[str(traze_value) + '_PEAKS_O' + str(order)] = np.array(aux_peaks)
        peaks_df[str(traze_value) + '_NP_O' + str(order)] = auxNP
        peaks_df[str(traze_value) + '_F_O' + str(order)] = auxF
        peaks_df[str(traze_value) + '_DP_O' + str(order)] = auxDP

        peaks_df['Cell'] = auxCell
        peaks_df = peaks_df.set_index(['Cell'], drop=True)
        peaks_conc[c] = peaks_df

    return (conc, peaks_conc)



def peaks_grid_plot_v2(conc, peaks_conc, len_cells, M, ylim, order, traze_value,order_mins, save=False,cols = 5,long=1, **kwargs):

    col = str(traze_value) + '_PEAKS_O' + str(order)
    colors = sns.color_palette('muted')
    k = 0; z = 0
    Tot = np.sum(len_cells); Cols = cols;  Rows = Tot // Cols;
    if Tot % Cols != 0:
        Rows = Rows + 1
    M = np.max(M)

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*2, 11.69*2*long))
    axs = axs.ravel()
    plot_mins = kwargs.get('plot_mins', True)

    for c in conc:
        df = conc[c]
        df.sort_index(inplace=True)
        peaks_df = peaks_conc[c]
        peaks_df.sort_index(inplace=True)
        color = colors[k];  k = k + 1
        cells = df.index.get_level_values(0).unique()
        plt.rcdefaults()

        for cell in cells:
            ax = axs[z];z = z + 1
            ax.set_ylim(ylim)
            ax.set_xlim([0, M])

            ax.grid(False)
            silent_ax(ax)


            if cells[-1] >= cell:
                ax.plot(df.loc[cell][traze_value], color=color, linewidth=0.7,alpha=0.7)
                ax.text(0.9, 0.1, 'Cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
                ax.plot(df.loc[cell][col], 'o', color='r', markersize=0.3, alpha=1)
                if plot_mins:
                    col_mins = 'MIN_O' + str(order_mins)
                    ax.plot(df.loc[cell][col_mins], 'o', color='k', markersize=0.3, alpha=1)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)

    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')

