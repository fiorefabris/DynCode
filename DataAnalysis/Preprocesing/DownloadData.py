import os
import pandas as pd


def download_data(data_folder, save=True, **kwargs):
    '''Returns a dictionary which has an entry per condition. Each entry has a DataFrame with (cells,frames) as index.
    Specifications are write down on data_chart.xlsx file if save=True

      INPUTS: - data_folder: a str with folder address

    WARNINGS: -  Time Series with more than 165 frames are considered
              -  Time Series has erased the last 15 points
    '''
    clean_last_points = kwargs.get('clean_last_points', True)
    clean_last_points = False
    data_conc_folder = os.listdir(data_folder)
    Columns = ['Label', 'QUALITY', 'POSITION_X', 'POSITION_Y', 'POSITION_Z', 'POSITION_T', 'MANUAL_COLOR', 'RADIUS',
               'VISIBILITY', 'CONTRAST', 'SNR', 'MEDIAN_INTENSITY', 'MIN_INTENSITY', 'MAX_INTENSITY',
               'TOTAL_INTENSITY', 'ESTIMATED_DIAMETER']


    conc = {}

    doc_cols = ['Cell_Number', 'First_spot_ID', 'Track_ID', 'Starting_frame', 'Ending_frame', 'Condition']
    doc = pd.DataFrame(columns=doc_cols)
    for j in data_conc_folder:
        folders = os.listdir(data_folder + j)

        df2 = pd.DataFrame()
        aux,aux2 = [],[]
        i = 0

        for k in folders:
            file = data_folder + j + '/' + k + '/' + os.listdir(data_folder + j + '/' + k)[0]
            df = pd.read_csv(file,sep='\t')
            df = df.set_index(['TRACK_ID', 'FRAME'], drop=False)
            for algo, data in df.groupby(level=0):
                if clean_last_points:
                    if len(data.iloc[0:-15]) > 150:
                        description = [i, data.ID.iloc[0], data.TRACK_ID.iloc[0]
                            , data.FRAME.iloc[0], data.FRAME.iloc[-1], k + '/' + os.listdir(data_folder + j + '/' + k)[0]]

                        # doc = doc.append(pd.DataFrame([description], columns = doc_cols),ignore_index=True,sort = False)
                        doc = doc.append(pd.DataFrame([description], columns=doc_cols), ignore_index=True)

                        df2 = pd.concat([df2, data.iloc[0:-15]],sort=True)
                        aux = aux + [i] * len(data.iloc[0:-15])
                        aux2 = aux2 + (data.iloc[0:-15]).index.get_level_values(1).tolist()
                        i = i + 1
                else:
                    if len(data) > 150:
                        description = [i, data.ID.iloc[0], data.TRACK_ID.iloc[0]
                            , data.FRAME.iloc[0], data.FRAME.iloc[-1],
                                       k + '/' + os.listdir(data_folder + j + '/' + k)[0]]

                        # doc = doc.append(pd.DataFrame([description], columns = doc_cols),ignore_index=True,sort = False)
                        doc = doc.append(pd.DataFrame([description], columns=doc_cols), ignore_index=True)

                        df2 = pd.concat([df2, data], sort=True)
                        aux = aux + [i] * len(data)
                        aux2 = aux2 + (data).index.get_level_values(1).tolist()
                        i = i + 1

        index = [aux, aux2]
        index_tuples = list(zip(*index))
        index = pd.MultiIndex.from_tuples(index_tuples, names=['cell', 'frames'])

        df2 = df2.set_index(index, drop=False)
        df2.drop(Columns, axis=1, inplace=True)

        conc[j] = df2


    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        writer = pd.ExcelWriter(str(save_folder) + str(save_name)+ '.xlsx', engine='xlsxwriter')
        doc.to_excel(writer, sheet_name='ID_Data')
        writer.save()

    return (conc)


def drop_cells(conc,drop_cells):
    "drop cells is a dic. Its keys are the concentration name and the values are a list of the cells you want to drop"
    len_cells = []
    for c in drop_cells:
        df = conc[c]
        df.drop(drop_cells[c], axis=0,level='cell', inplace=True)
        conc[c] = df
    return (conc)


def clean_cells(conc,clean_cells):
    "clean cells is a dic. Its keys are the concentration name and the values are a list of the cells you want to clean (last 10 points)"
    len_cells = []
    for c in clean_cells:
        df = conc[c]
        df2 = pd.DataFrame()
        for cells, data in df.groupby(level='cell'):
            if cells in clean_cells[c]:
                print(len(data))
                data.drop(data.tail(20).index,inplace=True)
                print(len(data))
                df2 = pd.concat([df2, data], sort=True)
            else:
                df2 = pd.concat([df2, data], sort=True)
        conc[c] = df2
    return (conc)