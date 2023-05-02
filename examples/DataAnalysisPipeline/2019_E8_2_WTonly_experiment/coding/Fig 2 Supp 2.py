# Esto esta inspirado en una figura que es de /home/fabris/Documents/Dyncode/Dyncode_Figures/05-DATA_ANALYSIS-2019_E8_2_WTonly.py

import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)

sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
%matplotlib inline

def my_shuffle(array):
    random.shuffle(array)
    return array

#%%
    
from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import pickle
import sys



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

#%%    
filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/2019_E8_2_WTonly.pkl'
conc = load_file(filename)

main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
conc_data = {}
key = ['an_WT_ESL_PD03','an_WT_ESL']
keywords = ['NC','WT']
TH=('187.17241379310343', '2068.9655172413795')
slope_th = 187
amp_th = 2068
for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)

#%%
ipi_conc = ['an_WT_ESL']
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.95, hspace=0.0)
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[0], wspace=0.3, hspace=0.3)

colors =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))
colors = [colors[i] for i in [15,25]]


col1 = 'amp_peaks'
col2= 'dt_peaks'
col3 = 'IPI'

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    col4= (df['amp_peaks'].dropna() + df['amp_peaks'].dropna().groupby(level='cell').shift(1))/2


    ax1 = plt.subplot(gs0[k,0]); 
    ax1.plot(my_shuffle(df['dt_peaks'].dropna().values),my_shuffle(df['amp_peaks'].dropna().values),'o',markersize=2,alpha=1,color = 'grey',linewidth = 0)
    ax1.plot(df['dt_peaks'].dropna().values,df['amp_peaks'].dropna().values,'o',markersize=3,alpha=1,color = colors[k+1],linewidth = 0)
    ax1.set_xlim([0, 52]);
    ax1.set_ylim( [0, 13000])
    ax1.axhspan(0,amp_th, color='darkgoldenrod',alpha = 0.1,linewidth=0.0)
    X = np.arange(0,60)
    Y = slope_th/2 * X
    ax1.axvspan(0,5,color='darkgoldenrod',alpha = 0.1,linewidth=0.0)
    ax1.fill_between(X,X*0,Y,color='darkgoldenrod',alpha = 0.1,linewidth=0.0)
    
    set_scale(ax1,np.arange(0,100,15),np.arange(0,16000,2000))
    ax1.set_xticklabels([int(x) for x in np.arange(0,100,15)*20/60],fontsize=6)
    ax1.set_yticklabels([int(x) for x in np.arange(0,16000,2000)/100],fontsize=6)
    ax1.xaxis.set_label_coords(0.5, -0.10);ax1.yaxis.set_label_coords(-0.20,0.5)

    ax1.set_ylabel('amplitude (a.u.)',fontsize=8,rotation = 90); ax1.set_xlabel('duration (min)',fontsize=8);


plt.savefig(save_folder+'Fig. 2 Supp. 2.pdf', format='pdf')