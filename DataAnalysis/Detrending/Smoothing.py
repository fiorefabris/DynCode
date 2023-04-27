import pandas as pd
from scipy.signal import butter, filtfilt, freqs,freqz
import matplotlib.pyplot as plt
import numpy as np



def smooth_signal(conc,traze_value,len_tw):
    for c in conc:
        df = conc[c]; sm_data = pd.DataFrame([])
        for cells, data in df.groupby(level='cell'):
            sm_data = pd.concat([sm_data,(data[traze_value].rolling(window=len_tw,center=True)).mean()])

        sm_data = pd.DataFrame(sm_data.values, index=pd.MultiIndex.from_tuples(sm_data.index, names=['cell', 'frames']),
                                   columns=[traze_value])
        mask= sm_data[traze_value].isna()
        sm_data.loc[mask,traze_value] = df[traze_value][mask].values
        df['sm_'+traze_value] = sm_data

#%%
def frq_smooth(conc,traze_value,highcut, fs,plot_filt=False, order=5):
    ''' fs: sampling rate
        highcut: desired frequency to filter
                    '''
    for c in conc:
        df = conc[c]
        auxLP = []
        nyq = 0.5 * fs
        high = highcut / nyq
        b_LP, a_LP = butter(order, high, btype='lowpass')

        for cells, data in df.groupby(level='cell'):
                auxLP.extend(filtfilt(b_LP, a_LP, data['MEAN_INTENSITY'], padtype='odd'))
        df[traze_value+'_LP'] = auxLP

        if plot_filt:
            plt.figure(1)
            plt.clf()

            w, h = freqz(b_LP, a_LP, worN=2000)
            plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
            plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
                     '--', label='sqrt(0.5)')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Gain')
            plt.grid(True)
            plt.legend(loc='best')