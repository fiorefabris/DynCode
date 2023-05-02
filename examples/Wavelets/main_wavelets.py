#/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2019/PycharmProjects/DynCode

from wavelets import WaveletAnalysis,wavelets
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import sys

import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)

def download_data(filename):
    infile =  open(filename,'rb')
    results = pickle.load(infile)
    infile.close()
    return(results)
    
def save_data(data, filename):
    outfile= open(filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()
    
#%%

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')
from plotting_main import load_file,save_file,get_linear_detrend,set_scale_bars,joint_duration,is_consecutive


filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/2019_E8_2_WTonly.pkl'
conc = load_file(filename)

main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
conc_data = {}
key = ['an_WT_ESL_PD03','an_WT_ESL']
keywords = ['NC','WT']
TH=('187.17241379310343', '2068.9655172413795')
amp_th = 2068
slope_th = 187
for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)


#for key in list(conc_data.keys()):
#    dm =find_dm(conc_data[key])

def ax_fontsize(fontsize,ax,rotation=None):
    plt.setp(ax.get_xticklabels(), rotation=rotation, fontsize=fontsize)

#%%
w

# given a signal x(t) and a sample spacing dt
#conc = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/local_quantifiers/conc.pkl')
conc = download_data('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/new_cdc/conc_data.pkl')
#x = conc['an_ESC_FAX'].MEAN_INTENSITY.values
x = conc['an_WT_ESL'].MEAN_INTENSITY.values
dt = 1/3
#%%
x = np.sin(2*np.pi/20*np.arange(0,200,1/3))
wa = WaveletAnalysis(x, dt=dt,dj=dj,wavelet=wavelets.Morlet(w0=3))
power = wa.wavelet_power # wavelet power spectrum
scales = wa.scales # scales 
t = wa.time # associated time vector
rx = wa.reconstruction() # reconstruction of the original data

#%%
from DataAnalysis.Detrending.MinPolDT import Search_min,Fit_pol_min


filter_name = 'sm'
traze_value = 'sm_MEAN_INTENSITY'
dist_bet_peaks_order = 1
order_mins=1

Search_min(conc_data, [1],traze_value = 'sm_MEAN_INTENSITY' )
Fit_pol_min(conc_data, 5, [1])


#Search_min(conc, [1],traze_value = 'sm_MEAN_INTENSITY' ) busca los minimos de orden [1] 
#Fit_pol_min(conc, 5, [1])
#%%

fig, ax = plt.subplots()
T, S = np.meshgrid(t, scales)
ax.contourf(T, S, power, 100)
a,b=wa.coi
ax.fill_between(a,b,1000,color = 'gray',alpha = 0.5)
ax.set_yscale('log')
#fig.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/local_quantifiers/figures_wavelets/wavelets.png')
#%%
plt.close('all')
#Eq 5. Es la resoloción temporal del paso de la frecuencua en el espacio de funciones wavelet. 
dt = 1/3

dj = 1/32
#Esto tiene que ver con la resolución de 

#for cell,x in conc['an_WT_ESL'].groupby(level='cell'):
for cell,x in conc_data['an_WT_ESL'].groupby(level='cell'):

   # time_series = x.sm_MEAN_INTENSITY.values - x.FITM_O1.values
    time_series = x.sm_MEAN_INTENSITY.values - np.mean(x.sm_MEAN_INTENSITY.values)
    wa = WaveletAnalysis(time_series, dt=dt,dj=dj,wavelet=wavelets.DOG(m=4))#wavelets.Morlet(w0=2)
    power = wa.wavelet_power # wavelet power spectrum
    scales = wa.scales # scales 
    t = wa.time # associated time vector
    rx = wa.reconstruction() # reconstruction of the original data
   
    fig, ax = plt.subplots(2,2,gridspec_kw={'height_ratios': [1.5, 0.5],'width_ratios': [0.95,0.05]})
    plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.90)
   
    ax00 = ax[0,0]
    ax01 = ax[0,1]
    ax1 = ax[1,0]
    fig.delaxes(ax[1,1])

   
    T, S = np.meshgrid(t, scales)
    levels = np.linspace(100*(100**2), 3000*(100**2),100)
    img=ax00.contourf(T, S, power, 100,vmin = 100*(100**2),vmax = 3000*(100**2),levels = levels,extend='both') #max = np.quantile(power,0.80)
    ax00.set_ylim([0,30])
    #ax00.set_yscale('log')
    ax00.set_ylabel('period (min)')
    #ax00.set_xlabel('time (min)')
    ax00.set_xlim([0,100])
    
    cbar_ticks = np.linspace(100*(100**2), 3000*(100**2),5)
    cbar = plt.colorbar(img, cax=ax01,ticks=cbar_ticks)
    cbar.ax.set_yticklabels(np.round(i/(100**2),2) for i in cbar_ticks)
    cbar.set_label('power (a.u.)', rotation=270)
    
    #ax1.plot(x.sm_MEAN_INTENSITY.values)
    ax1.plot(time_series)
    #ax1.plot(rx,color='red',linestyle = '-',alpha = 0.1)
    ax1.set_xlabel('time (min)')
    ax1.set_ylabel('KTR intensity (a.u.)')
    
    ax1.set_xlim([ax00.get_xticks()[0]*3,ax00.get_xticks()[-1]*3])
    ax1.set_xticks([i*3 for i in ax00.get_xticks()])
    ax1.set_xticklabels([np.round(i/3,1) for i in ax1.get_xticks()])
    
    #ax1.set_ylim([58000,65000])
    ax1.set_ylim([-5000,10000])
    #ax1.set_yticks([58000,65000])
    #ax1.set_ylim([43000,63000])
    #ax1.set_yticks([43000,63000])
    ax1.set_yticklabels([int(np.round(i/100,0)) for i in ax1.get_yticks()])

    name = 'wavelets_cell_'+str(cell)
 #%%   fig.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/local_quantifiers/figures_wavelets/'+name+'.png')

#%%
    plt.plot(x.sm_MEAN_INTENSITY.values)
    plt.plot(rx)
    
#%%
for cell,x in conc['an_WT_ESL_PD03'].groupby(level='cell'):
    
    time_series = x.MEAN_INTENSITY.values #- np.mean(x.MEAN_INTENSITY.values)#x.sm_MEAN_INTENSITY.values - x.FITM_O1.values
    wa = WaveletAnalysis(time_series, dt=dt,dj=dj)
    power = wa.wavelet_power # wavelet power spectrum
    scales = wa.scales # scales 
    t = wa.time # associated time vector
    rx = wa.reconstruction() # reconstruction of the original data
   
    fig, ax = plt.subplots(2,2,gridspec_kw={'height_ratios': [1.5, 0.5],'width_ratios': [0.95,0.05]})
    plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.90)
   
    ax00 = ax[0,0]
    ax01 = ax[0,1]
    ax1 = ax[1,0]
    fig.delaxes(ax[1,1])

   
    T, S = np.meshgrid(t, scales)
    levels = np.linspace(100*(100**2),3000*(100**2),100)
    img=ax00.contourf(T, S, power, 100,vmin = 100*(100**2),vmax = 3000*(100**2),levels = levels,extend='both') #max = np.quantile(power,0.80)
    ax00.set_ylim([0,30])
    #ax00.set_yscale('log')
    ax00.set_ylabel('period (min)')
    #ax00.set_xlabel('time (min)')
    ax00.set_xlim([0,100])
    
    cbar_ticks = np.linspace(100*(100**2),3000*(100**2),5)
    cbar = plt.colorbar(img, cax=ax01,ticks=cbar_ticks)
    cbar.ax.set_yticklabels(np.round(i/(100**2),2) for i in cbar_ticks)
    cbar.set_label('power (a.u.)', rotation=270)
    
    #ax1.plot(x.sm_MEAN_INTENSITY.values)
    ax1.plot(time_series)
    #ax1.plot(rx,color='red',linestyle = '-',alpha = 0.1)
    ax1.set_xlabel('time (min)')
    ax1.set_ylabel('KTR intensity (a.u.)')
    
    ax1.set_xlim([ax00.get_xticks()[0]*3,ax00.get_xticks()[-1]*3])
    ax1.set_xticks([i*3 for i in ax00.get_xticks()])
    ax1.set_xticklabels([np.round(i/3,1) for i in ax1.get_xticks()])
    
    #ax1.set_ylim([58000,65000])
    ax1.set_ylim([-5000,10000])
    #ax1.set_yticks([58000,65000])
    #ax1.set_ylim([43000,63000])
    #ax1.set_yticks([43000,63000])
    ax1.set_yticklabels([int(np.round(i/100,0)) for i in ax1.get_yticks()])

    name = 'wavelets_cell_'+str(cell)