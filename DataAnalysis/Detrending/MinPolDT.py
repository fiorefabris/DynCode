import numpy as np
import pandas as pd
from openpyxl import load_workbook


def Find_Min(Data, Order):
    ''' Searchs for local minima of Data of an given time window (Order).
    Returns a Pandas Series with the values and positions of the local minima
    '''
    Data_Peaks = Data
    m = Data.max()
    for Shift in range(1, Order + 1):  # Llega hasta order
        Data_Peaks = Data_Peaks.where(
            (Data_Peaks < (Data.shift(Shift)).fillna(m)) & (Data_Peaks < (Data.shift(- Shift)).fillna(m)), m)
    Data_Peaks = Data_Peaks.mask(Data_Peaks == m, None)
    return (Data_Peaks.values)


def Search_min(conc, orders_min,traze_value = 'MEAN_INTENSITY' ):
    ''' Searchs for local minima of each Time Serie of conc of an given time window
    (it size is detailed in the list orders_min (list)).

    RETURN - conc with aditional columns where the minima values are detailed [named, for example, 'MIN_O3'
            for window size 3 + 1 + 3 frames]

    WARNING - the non-minima points are None points
    '''
    for C in conc:
        df = conc[C]
        for order in orders_min:
            aux = []
            for cells, data in df.groupby(level='cell'):
                Mins = Find_Min(data[traze_value], order)
                aux.extend(Mins)
            df['MIN_O' + str(order)] = np.array(aux)
    return (conc)


def Fit_pol_min(conc, delta, orders_min,save=False,**kwargs):
    ''' Makes a polinomial fit over the given minima & the first and last point of the Time Serie.

    degree of the polinomial function = (2 + #of minima)/delta

    RETURN - conc with aditional columns where the polinomia points are detailed [named, for example, 'FITM_O'
            for minima of window size 3 + 1 + 3 frames]
            - Specifications are write down on data_chart.xlsx file: the order of the polinomia, the number of minima
            and the polonomia coefficients for each time serie.

    '''
    # Fitea un polinomio de orden deg con los minimos.

    auxconc, auxcell = [], [];
    for C in conc:
        df = conc[C]
        for cells, data in df.groupby(level='cell'):
            auxconc.append(C)
            auxcell.append(cells)

    index = [auxconc, auxcell]
    index_tuples = list(zip(*index))
    index = pd.MultiIndex.from_tuples(index_tuples, names=['conc', 'cell'])
    table = pd.DataFrame(index=index)

    for C in conc:
        df = conc[C]
        for order in orders_min:
            aux, auxdeg, auxz, auxmins = [], [], [], [];
            for cells, data in df.groupby(level='cell'):
                points = [data['MEAN_INTENSITY'].iloc[0]] + data['MIN_O' + str(order)].dropna().tolist() + [
                    np.array(data['MEAN_INTENSITY'].iloc[-1])]
                time = [data['MEAN_INTENSITY'].index.get_level_values(1)[0]] + data[
                    'MIN_O' + str(order)].dropna().index.get_level_values(1).tolist() + [
                           data['MEAN_INTENSITY'].index.get_level_values(1)[-1]]

                deg =int(len(time) / delta) #poner 1 cuando queres una recta
                #deg = 1
                auxdeg.append(deg)
                # print('Concentración ' + str(C)+ ';Célula '+str(cells)+ ';Polinomio de orden '+str(deg))

                z = np.polyfit(time, points, deg=deg)
                p = np.poly1d(z)
                aux.extend(p(np.array(data['MIN_O' + str(order)].index.get_level_values(1))))
                auxz.append(z)
                auxmins.append(len(data['MIN_O' + str(order)].dropna().tolist()))

            df['FITM_O' + str(order)] = aux
            #table.loc[C, str(order) + '_Pol_Fit_Order'] = np.ravel(auxdeg)
            #table.loc[C, str(order) + '_Pol_Fit_Coef'] = np.ravel(auxz)
            #table.loc[C, str(order) + '_Number_of_Mins'] = np.ravel(auxmins)

    if False:
        if save:
            save_folder = kwargs.get('save_folder', None)
            save_name = kwargs.get('save_name', None)
            book = load_workbook(str(save_folder) + str(save_name)+ '.xlsx')
            writer = pd.ExcelWriter(str(save_folder) + str(save_name)+ '.xlsx', engine='openpyxl')
            writer.book = book
            version = kwargs.get('version', None)
            if version != None:
                table.to_excel(writer, sheet_name=str(version) + '_MinPolDT_Data')
            else:
                table.to_excel(writer, sheet_name='_MinPolDT_Data')
            writer.save()

    return (conc,table)
