from scipy.signal import butter, filtfilt, freqs


def butter_filter(conc,low,high):
    ''' Applies a Butter Filter over the Time Series, using the zero phase filtering strategy.

    RETURN - conc with aditional columns with the Filtered Time Series [named, MEAN_INTENSITY_HP, LP and BP respectevely]
    '''

    for c in conc:
        df = conc[c]
        auxHP, auxLP, auxBP = [], [], []
        for cells, data in df.groupby(level='cell'):

            if len(data['MEAN_INTENSITY']) > 20:

                # Highpass filter
                b_HP, a_HP = butter(4, low, btype='highpass')
                auxHP.extend(filtfilt(b_HP, a_HP, data['MEAN_INTENSITY'], padtype='odd'))

                # Lowpass filter
                b_LP, a_LP = butter(4, high, btype='lowpass')
                auxLP.extend(filtfilt(b_LP, a_LP, data['MEAN_INTENSITY'], padtype='odd'))

                # Bandpass filter
                auxBP.extend(filtfilt(b_LP, a_LP, filtfilt(b_HP, a_HP, data['MEAN_INTENSITY'], padtype='odd')
                                      , padtype='odd'))

            else:
                auxHP.extend([None] * len(data['MEAN_INTENSITY']))
                auxLP.extend([None] * len(data['MEAN_INTENSITY']))
                auxBP.extend([None] * len(data['MEAN_INTENSITY']))
        df['MEAN_INTENSITY_HP'] = auxHP
        df['MEAN_INTENSITY_LP'] = auxLP
        df['MEAN_INTENSITY_BP'] = auxBP
    return conc

