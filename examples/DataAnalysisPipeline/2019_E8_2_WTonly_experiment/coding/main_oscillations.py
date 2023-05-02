''' This is a file for generating the data analysis post-review commons revision'''

#%%
''' Importings'''
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 


from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)


import pickle
import sys


sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
colors = sns.color_palette('muted')
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
#%matplotlib inline
matplotlib.rcParams['savefig.transparent'] = True


def load_file(filename):
    infile = open(filename,'rb')
    output = pickle.load(infile)
    infile.close()
    return(output)

def save_file(variable,filename):
    outfile= open(filename,'wb')
    pickle.dump(variable,outfile)
    outfile.close()

def get_linear_detrend(conc):
    '''
    Calculates the linear detrend of each df.MEAN_INTENSITY 
    in conc.
    
    -------------------------------------------------
    INPUT:
        - conc (dict): dictionary of dataframes
    
    -------------------------------------------------
    OUTPUT:
        - conc with two aditional columns: 'trend' and 'detrended'
    '''
    for key in conc.keys():
        df = conc[key]
        trend_aux = []
        detrended_aux= []
        
        for cells, data in df.groupby(level='cell'):
            
        # fit linear model
            X = data.index.get_level_values(1)
            X = np.reshape(X, (len(X), 1))
            y = data.MEAN_INTENSITY.values
            model = LinearRegression()
            model.fit(X, y)
        # calculate trend
            trend = model.predict(X)
            detrended = [y[i]-trend[i] for i in range(0, len(data.MEAN_INTENSITY))]
            trend_aux.extend(list(trend))
            detrended_aux.extend(list(detrended))
        conc[key]['trend'] = trend_aux
        conc[key]['detrended'] = detrended_aux

def ax_fontsize(fontsize,ax,rotation=None):
    plt.setp(ax.get_xticklabels(), rotation=rotation, fontsize=fontsize)


#%%
################################################
#### Fast Fourier Transform computation
################################################  


def norm_fft_yf(data, Dt, max_freq = None):
    '''
    norm_fft_yf(data, Dt, max_freq = None)
    Computes a one-dimentional discrete Fourier Transform. 
    
    Computes the one-dimensional discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform (FFT)
    algorithm. Returns the FFT on a orthonomal basis so Parseval equality is satisfied. 
    
    Parameters
    -----------
    data : list of floats
        amplitude time serie.
    Dt : float
        sampling rate.
    max_freq: float, optional --- actually is not working
        maximum frequency of the FFT. 
        
    Returns
    --------
        yf : complex ndarray 
        The Fast Fourier Transfrom values, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_xf: for obtenining the frequencies involved on the FFT transform.

    '''
    
    N = data.shape[0]
    Nf = N // 2 #if max_freq is None else int(max_freq * T)
    yf =  np.fft.fft(data,norm='ortho')
    return (yf[:Nf])[1:][0::]

def norm_fft_xf(data, Dt, max_freq = None):
    '''
    norm_fft_xf(data, Dt, max_freq = None)
    Gives the frequencies involved on the one-dimentional discrete Fourier Transform. 

    
    Parameters
    -----------
    data : list of floats
        amplitude time serie.
    Dt : float
        sampling rate.
    max_freq: float, optional --- actually is not working
        maximum frequency of the FFT. 
        
    Returns
    --------
        xf : complex ndarray 
        The angular frequencies involved on the one-dimentional discrete Fourier Transfor, in which only the positive frequencies are considered and the main frequency is not considered.
        The maximum frequency resolved is the nyquist frequency.
        The frequency resolution is 2*Dt / T, where T is the time series lenght.  
        
    See Also
    ---------
    norm_fft_yf: for computing he one-dimentional discrete Fourier Transform. 

    '''
    N = data.shape[0]
    Nf = N // 2 #if max_freq is None else int(max_freq * T)
    xf = np.linspace(0.0, 0.5 / Dt, N // 2)
    return (xf[:Nf]*(2*np.pi))[1:][0::]

def autocorr(x):
    ''' 
    autocorr(x)
    reurns the autocorrelation function of the time traze x
    
    Parameters
    ----------
    x: 1-d vector
        time series
    
    Returns
    -------
    out : 1-d vector
        autocorrelation vector output with non-negative timelags
    '''
    result = np.correlate(x, x, mode='full')
    return result[result.size//2:]

def autocorr_normed(x):
    ''' 
    autocorr_normed(x)
    reurns the normalized autocorrelation function of the time traze x
    
    Parameters
    ----------
    x: 1-d vector
        time series
    
    Returns
    -------
    out : 1-d vector
        autocorrelation vector output with non-negative timelags
    '''
    y = x - np.mean(x)
    norm = np.sum(y ** 2)
    
    result = np.correlate(y, y, mode='full')/norm
    return result[result.size//2:]

#%%
