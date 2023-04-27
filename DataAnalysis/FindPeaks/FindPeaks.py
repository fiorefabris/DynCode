import numpy as np
import pandas as pd



def Find_Peaks(Data, Abs_Thr, Order):
    ''' Searchs for local maxima over an Abs_Thr of Data of an given time window (Order).
    Returns
    - Data_Peaks : a Pandas Series with the values and positions of the local maxima
    - Number of peaks , time-lenght of Data & density of peaks
    '''

    #### antes habia otro!

    if Abs_Thr != None:
        Data_Peaks = Data.where(Data > Abs_Thr, 0)
        M = 0
    else:
        Data_Peaks = Data
        M = Data.min()

    for Shift in range(1, Order + 1):  # Llega hasta order
        Data_Peaks = Data_Peaks.where(
            (Data_Peaks > (Data_Peaks.shift(Shift)).fillna(M)) & (Data_Peaks > (Data_Peaks.shift(- Shift)).fillna(M)),M)

    Data_Peaks = Data_Peaks.mask(Data_Peaks == M)  # Data_Peaks = Data_Peaks.mask(Data_Peaks == m,None)
    return (Data_Peaks, Data_Peaks.count(), len(Data_Peaks), Data_Peaks.count() / len(Data_Peaks))



def Distance_Between_Peaks(data, col):
    ''' Calculates the frames between peaks.
    Input:
    - Pandas Series with peaks [in column named "col"] (if there is no peak, Nan)
    '''
    dt_Peaks = []
    t_Peak = 0

    for Index, Row in data.iterrows():
        if not np.isnan(Row[col]):
            dt_Peaks.append(Index[1] - t_Peak)
            t_Peak = Index[1]

    if len(dt_Peaks) > 0:
        del dt_Peaks[0]
    return (dt_Peaks)


def Full_Find_Peaks(conc, Thresholds,orders, traze_value):
    '''
    -Thresholds  : un array de los valores que quiera del Threshold
    -Traze value : La columna del Df que quiera analizar
    -orders = lista de ordenes que quiero ver

    RETURNS:
        - conc with peaks information
        - peaks_df: dictionary of df. Each df has the information of Number_Peaks, Frames, Dens_Peaks for each concentration.
    '''

    peaks_conc = {}

    for c in conc:
        df = conc[c]
        peaks_df = pd.DataFrame()

        for order in orders:
            for th in Thresholds:
                auxNP, auxF, auxDP, auxCell = [], [], [], []
                aux_peaks = []

                for cells, data in df.groupby(level='cell'):
                    Peaks, Number_Peaks, Frames, Dens_Peaks = Find_Peaks(data[traze_value], th, order)

                    aux_peaks.extend(Peaks)
                    auxNP.append(Number_Peaks)
                    auxF.append(Frames)
                    auxDP.append(Dens_Peaks)
                    auxCell.append(cells)

                df[str(traze_value) + '_' + str(th) + 'PEAKS_O' + str(order)] = np.array(aux_peaks)
                peaks_df[str(traze_value) + '_' + str(th) + 'NP_O' + str(order)] = auxNP
                peaks_df[str(traze_value) + '_' + str(th) + 'F_O' + str(order)] = auxF
                peaks_df[str(traze_value) + '_' + str(th) + 'DP_O' + str(order)] = auxDP

        peaks_df['Cell'] = auxCell
        peaks_df = peaks_df.set_index(['Cell'], drop=True)
        peaks_conc[c] = peaks_df

    return (conc, peaks_conc)


def Full_Dist_Between_Peaks(conc, peaks_conc, th, order, traze_value):
    '''
    Returns the distance between peaks in time units (minutes) and in frames for the th threshold,
    the order  and the traze value indicated
    '''

    col = str(traze_value) + '_' + str(th) + 'PEAKS_O' + str(order)
    dist_bet_peaks_conc = {}
    total_peaks_conc  = {}

    for c in conc:
        df = conc[c]
        aux_dist_bet_peaks, aux_cell = [], []

        dist_bet_peaks_df = pd.DataFrame()
        total_peaks_conc[c] = np.sum(peaks_conc[c].loc[:, str(traze_value) + '_' + str(th) + 'NP_O' + str(order)])

        for cells, data in df.groupby(level='cell'):
            dt_Peaks = Distance_Between_Peaks(data, col)
            aux_dist_bet_peaks.extend(dt_Peaks)
            aux_cell.extend([cells] * len(dt_Peaks))

        #aux_dist_bet_peaks_time = [i * 105 / 60 for i in aux_dist_bet_peaks]
        aux_dist_bet_peaks_time = [i * 20 / 60 for i in aux_dist_bet_peaks]

        dist_bet_peaks_df['dist_bet_peaks_' + str(col)] = np.array(aux_dist_bet_peaks)
        dist_bet_peaks_df['time_dist_bet_peaks_' + str(col)] = np.array(aux_dist_bet_peaks_time)
        dist_bet_peaks_df['Cell'] = np.array(aux_cell)
        dist_bet_peaks_df = dist_bet_peaks_df.set_index(['Cell'], drop=True)

        dist_bet_peaks_conc[c] = dist_bet_peaks_df


    return (dist_bet_peaks_conc, total_peaks_conc)
