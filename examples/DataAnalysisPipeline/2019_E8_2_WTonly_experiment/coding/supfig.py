import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 
#from supfig_peaks.py import peaks####
import pickle5 as pickle


from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.1
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)


import pickle
import sys
#sys.path.insert(0, '/home/fiore/Documents/DO2019/analysis_2019_E8_2_WTonly/coding/')
#from th_analysis_plotting import peaks,main_WT


sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
colors = sns.color_palette('muted')
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
%matplotlib inline
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.pyplot.tick_params(direction='out', length=2)


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

#%% Function for intermidiate step

def peaks(df_dic,params):
    
    """ Identifies peaks of a time serie & returns the corresponding 
    characterization: amp, dt, max, min, raise, fall, IPI.
    
    -------------------------------------------------------
    INPUTS:
            - params: amp_th and pendiente_th of a corresponding level.
            - df_dic (with identified peaks) 
    -------------------------------------------------------
    OUTPUT:
            - df: a DataFrame with the old and new information.
            
    """
           
    c,amp_th,pendiente_th  = params[0],params[1][1],params[1][0]
    #c,amp_th,pendiente_th  = params[0],params[1],params[2]
    df = df_dic[c]

    
    main_aux_amp,main_aux_dt = [],[]
    main_aux_max,main_aux_min = [],[]
    main_aux_raise,main_aux_fall = [],[]
    
    for cells, data in df.groupby(level='cell'):
        #Auxiliary variables for storing the calculations of the loop
        L = data.shape[0]
 
        aux_amp,aux_dt = [None]*L,[None]*L
        aux_max,aux_min = [None]*L,[None]*L       
        aux_raise,aux_fall = [None]*L,[None]*L
 
    
        col_mins = 'MIN_O1'
        mask_mins = data[col_mins].isna()
        mins = data[~ mask_mins][col_mins]
        
        col_max = 'sm_MEAN_INTENSITY_PEAKS_O1'
        mask_max = data[col_max].isna()
        maxs = data[~ mask_max][col_max]
        
        maxs_index = maxs.index.get_level_values(1)
        mins_index = mins.index.get_level_values(1)
        
        # Ensures that there are minima at the borders of the TS
        test1= True
        while test1:
            if maxs_index[0] <= mins_index[0]:
                maxs.drop([maxs_index[0]],level=1,inplace=True)
                maxs_index = maxs.index.get_level_values(1)
            else:
                test1=False
                
                
        test2=True
        while test2:
            if mins_index[-1] < maxs_index[-1]:
                maxs.drop([maxs_index[-1]],level=1,inplace=True)
                maxs_index = maxs.index.get_level_values(1)
            else:
                test2 = False
        
        # Searchs for the peaks
        amp_ant,idx_ant,idx_mins_der_ant,idx_mins_izq_ant = 0,0,0,0
        flag_share_ix = False
        
        for (i,j) in enumerate(maxs_index):
            flag_max = True;flag_len=True
            maximo = maxs[cells].loc[j]
            min_izq_index = -1
            min_der_index = 0
            
            # cheks left amplitude for the candidate maxima
            while ((maximo - mins[mins_index<j][cells].iloc[min_izq_index]) < amp_th) or ((maxs[cells].loc[j] - mins[mins_index<j][cells].iloc[min_izq_index])/(j- mins[mins_index<j].index.get_level_values(1)[min_izq_index])< pendiente_th) :
                min_izq_index = min_izq_index -1
                if abs(min_izq_index) > len(mins[mins_index<j][cells]):
                    min_izq_index = min_izq_index + 1
                    flag_len = False
                    break

            # cheks right amplitude for the candidate maxima   
            while ((maximo -mins[mins_index>j][cells].iloc[min_der_index]) < amp_th) or ((maxs[cells].loc[j] -  mins[mins_index>j][cells].iloc[min_der_index])/(mins[mins_index>j].index.get_level_values(1)[min_der_index]-j )< pendiente_th):
                min_der_index = min_der_index + 1
                if (abs(min_der_index) >= len(mins[mins_index>j][cells])):
                    min_der_index = min_der_index - 1
                    flag_len=False
                    break
                
            if flag_len:
                # Cheks than between two minima there is a maxima
                idx_mins_izq = mins_index[mins_index<j][min_izq_index] -data.index.get_level_values(1)[0]
                
                amp_der = (maximo -  mins[mins_index>j][cells].iloc[min_der_index])
                amp_izq = (maximo - mins[mins_index<j][cells].iloc[min_izq_index])
                amp = (amp_izq+amp_der)/2
                

#                if  idx_mins_der_ant > idx_mins_izq:
#                    if amp > amp_ant:
#                        # delete previous maxima and go on
#                        aux_amp[idx_ant] = None
#                        aux_dt[idx_ant] = None
#                        aux_max[idx_ant] = None
#                        if flag_share_ix==False: 
#                            aux_min[idx_mins_izq_ant] = None# puede ser que se comparta
#                        aux_min[idx_mins_der_ant] = None     
#                        aux_raise[idx_ant] = None
#                        aux_fall[idx_ant] = None
#         
#               
#                    if amp <= amp_ant:
#                        # go to next maxima
#                        flag_max = False
                        
         
                dt_izq = j- mins[mins_index<j].index.get_level_values(1)[min_izq_index]
                dt_der = mins[mins_index>j].index.get_level_values(1)[min_der_index]-j
                    
                
                if flag_max: 
                        amp = (amp_izq + amp_der)/2
                        dt =  dt_der + dt_izq
                        raise_ = amp_izq/dt_izq
                        fall_ = amp_der/dt_der
                        
                        idx = j-data.index.get_level_values(1)[0]
                        aux_amp[idx] = amp
                        aux_dt[idx] = dt
                        aux_max[idx] = maximo
                        idx_mins_izq = mins_index[mins_index<j][min_izq_index] -data.index.get_level_values(1)[0]
                        aux_min[idx_mins_izq] = mins[mins_index<j][cells].iloc[min_izq_index]
                        
                        idx_mins_der = mins_index[mins_index>j][min_der_index] -data.index.get_level_values(1)[0]
                        aux_min[idx_mins_der] = mins[mins_index>j][cells].iloc[min_der_index]            
             
                        aux_raise[idx] = raise_
                        aux_fall[idx] = fall_
                        
                        if idx_mins_der_ant == idx_mins_izq: 
                            flag_share_ix = True
                        else:
                            flag_share_ix = False
                        
                        amp_ant = amp
                        idx_ant = idx
                        idx_mins_der_ant = idx_mins_der
                        idx_mins_izq_ant = idx_mins_izq
                    
        main_aux_amp = main_aux_amp + (aux_amp)    
        main_aux_dt = main_aux_dt+ (aux_dt)
        main_aux_max = main_aux_max + (aux_max)
        main_aux_min = main_aux_min + (aux_min)
        main_aux_raise = main_aux_raise + (aux_raise)
        main_aux_fall = main_aux_fall + (aux_fall)
        
    df['amp_peaks']= main_aux_amp

    df['dt_peaks']= main_aux_dt
    df['max_']= main_aux_max
    df['min_']= main_aux_min
    df['raise_']= main_aux_raise
    df['fall_']= main_aux_fall

    mask = df['max_'].notna()
    df['IPI'] = df[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-

    return(df)
    

#%%

def set_scale_bars(ax, x_bar=None, y_bar=None, xunits='', yunits='', x_scale=1, y_scale=1, round_x=False,round_y=False):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if x_bar is not None:
        x0, xf = x_bar[0], x_bar[1]
        ax.spines['bottom'].set_bounds(x0, xf)
        ax.set_xticks([(xf + x0) / 2])
        scale_x = (xf - x0) // x_scale
        #if round_x:       ax.set_xticklabels([str(int(scale_x)) + xunits])#, rotation='vertical', va='center')
        #else:         ax.set_xticklabels([str(scale_x) + xunits])#, rotation='vertical', va='center')
        ax.set_xticklabels([])
        ax.xaxis.set_tick_params(length=0)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
    if y_bar is not None:
        y0, yf = y_bar[0], y_bar[1]
        ax.spines['left'].set_bounds(y0, yf)
        ax.set_yticks([(yf + y0) / 2])
        scale = (yf - y0) // y_scale
        if round_y:
            scale = int(scale)
        #ax.set_yticklabels([str(scale) + yunits], rotation='vertical', va='center')
        ax.set_yticklabels([])
        ax.yaxis.set_tick_params(length=0)
    else:
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
        
        
        
        
        
# =============================================================================
#  Download data
# =============================================================================
        
        #%%    DOWNLOAD WILD TYPE DATA
filename= '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/2019_E8_2_WTonly.pkl'
conc = load_file(filename)

main_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
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


def ax_fontsize(fontsize,ax,rotation=None):
    plt.setp(ax.get_xticklabels(), rotation=rotation, fontsize=fontsize)
    
#%% charge intermidiate steps
th = [187.17241379310343, 2068.9655172413795]
df_dic_ = {}
df_dic_['an_WT_ESL'] = peaks(conc_data.copy(),['an_WT_ESL',th])
df_dic_['an_WT_ESL_PD03'] = peaks(conc_data.copy(),['an_WT_ESL_PD03',th])

save_file(df_dic_,main_folder + 'df_dic_.pkl')

#%%
# =============================================================================
#  Figure
# =============================================================================
        
#%%

# EL GRIDSPECT
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'

#Figura principal
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
gs_main = gridspec.GridSpec(nrows=1, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.90, hspace=0.0)

#3 filas en toda la hoja
gs_main_cols = gridspec.GridSpecFromSubplotSpec(nrows=7, ncols=1, subplot_spec=gs_main[0], wspace=0.4, hspace=0.0, height_ratios=[0.1,0.1,0.1,0.1,0.1,1,0.5])

#%Para cada fila, 3 columnas
gs_row0  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=5, subplot_spec=gs_main_cols[0], wspace=0.2, hspace=0.2,width_ratios = [1,1,1,1,0.5])
gs_row1  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=5, subplot_spec=gs_main_cols[1], wspace=0.2, hspace=0.2,width_ratios = [1,1,1,1,0.5])
gs_row2  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=5, subplot_spec=gs_main_cols[2], wspace=0.2, hspace=0.2,width_ratios = [1,1,1,1,0.5])
gs_row2b = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=5, subplot_spec=gs_main_cols[3], wspace=0.2, hspace=0.2,width_ratios = [1,1,1,1,0.5])
gs_row3  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=5, subplot_spec=gs_main_cols[4], wspace=0.2, hspace=0.2,width_ratios = [1,1,1,1,0.5])
gs_row4  = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main_cols[5], wspace=0.2, hspace=0.0,height_ratios = [0.3,0.7],width_ratios = [1,1,0.025])
gs_row5  = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_main_cols[6], wspace=0.0, hspace=0.0,height_ratios = [0.3,0.7])


#%%
# =============================================================================
# Selected traces
# =============================================================================
ylim=[43000,63000]

M = max([conc_data[key].max().FRAME for key in conc_data.keys()])
#l = ['*','**','***','****','*****','*','**','***','****','*****']
#traze_label = [0,1,2]
cells_conc = {'an_WT_ESL':[19,3,31,20]} 
cells_conc = {'an_WT_ESL_PD03':[37,66],'an_WT_ESL':[40,43]}

conc_labels = {'an_WT_ESL_PD03':'MEKi','an_WT_ESL': 'serum + LIF'}
filter_labels = ['raw data','smoothed data', 'max and min', 'after amp. threshold', 'detected pulses']

box=[];handles = [None]*2
colors =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))
colors = [colors[i] for i in [25,15]]


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
k = 0; c = 'an_WT_ESL'
df = conc_data[c];df.sort_index(inplace=True); 
color=colors[k]

ax = plt.subplot(gs_row0[-1])
ax.text(-0.4, 0.9,   filter_labels[0], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

numero = 0
for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row0[numero]); numero=numero+1
        ax.plot(df.loc[cells_conc[c][w]]['MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)

        if w == 0: 
            ax.text(0.05, 1.5,  conc_labels[c], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = colors[k])
       
        set_scale_bars(ax)




            
# =============================================================================
# # smooth
# =============================================================================

ax = plt.subplot(gs_row1[-1])
ax.text(-0.4, 0.9,   filter_labels[1], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

numero = 0
for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row1[numero]); numero=numero+1
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)

# =============================================================================
# # max and min
# =============================================================================

ax = plt.subplot(gs_row2[-1])
ax.text(-0.4, 0.9,   filter_labels[2], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)


numero = 0
for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2[numero]); numero=numero+1
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY_PEAKS_O1'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['MIN_O1'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)




# =============================================================================
# # Intermidate step
# =============================================================================

df_dic_ = load_file(main_folder + 'df_dic_.pkl')
ax = plt.subplot(gs_row2b[-1])
ax.text(-0.4, 0.9,   filter_labels[3], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)


numero = 0
for k,c in enumerate(cells_conc):
    df_ = df_dic_[c];df_.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2b[numero]); numero=numero+1
        ax.plot(df_.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        handles[0], = ax.plot(df_.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df_.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)

# =============================================================================
# # Maximo final
# =============================================================================

ax = plt.subplot(gs_row3[-1])
ax.text(-0.4, 0.9,   filter_labels[4], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)



numero = 0
for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row3[numero]); numero=numero+1
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
        if (w ==0) & (k ==0):
            ax.set_ylabel('KTR signal \n 200 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
            ax.xaxis.set_label_coords(0.2, -0.15);
        


##########################################
### CHARGING THE FILES ###################
##########################################

data_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/'
dataset_name='2019_E8_2_WTonly'

filename = data_folder + dataset_name + '.pkl'
conc_df = pd.read_pickle(filename)


filename = data_folder + 'df_results_WT_Norm_.pkl'
infile =  open(filename,'rb')
df_result_WT = pickle.load(infile)
infile.close()

filename = data_folder + 'df_results_CN_Norm_.pkl'
infile =  open(filename,'rb')
df_result_CN = pickle.load(infile)
infile.close()



#heatmap_values = df_result_WT#heatmap
heatmap_condition = ['WT','NC']
lcurves_values = df_result_CN.values #curvas de nivel
lcurves_condition = 'CN'
    

################
##plot levercurves
################


for w, heatmap_values_ in enumerate([df_result_CN,df_result_WT]):
    
    x = np.array(heatmap_values_.columns) #pendiente
    y = np.array(heatmap_values_.index) #amplitú

    xx,yy = np.meshgrid(x,y)
    Z=heatmap_values_.values
    
    ax = plt.subplot(gs_row4[1,w])

    plt.grid(0)
    levels= np.geomspace(1e-5,1e-1,num=50) #50
    levels=np.arange(1e-5,0.05+5e-5,5e-7)
    
    CSF = ax.contourf(xx, yy, Z,levels=levels,extend='both',cmap='Greys',vmin=1e-5,
                      vmax=0.05,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.05,linthresh=1e-5,linscale=0))
    c_levels= [1e-4]
    
    CS = ax.contour(xx, yy, lcurves_values,levels=c_levels,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.05,linthresh=1e-5), colors='green',linewidths=1)
    ax.plot([187], [2068],marker="v", markersize=4, alpha=1, color="red",linewidth = 0)

    ax.tick_params(axis='x',rotation=0,length=2)

    set_scale(ax,np.arange(0,500,100),[0,6000])
    ax.set_yticklabels([0/100,6000/100],fontsize=6)
    ax.set_xticklabels(np.arange(0,500,100)/100/20*60,fontsize=6)
    ax.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax.set_rasterization_zorder(-10)


    
    if w==0:
        ax.set_ylabel(r'$A_{th}$ (a.u.)',fontsize=8); 
        ax.set_xlabel(r'$v_{th}$ (a.u./min)',fontsize=8);
        ax.xaxis.set_label_coords(0.5,-0.17);
        ax.yaxis.set_label_coords(-0.05,0.5);

ax_cbar = plt.subplot(gs_row4[1,w+1])

cbar = plt.colorbar(CSF,  cax=ax_cbar, label= 'averaged pulse density ' + r'$(min^{-1})$',
                        ticks= [1e-5,1e-4,0.05],spacing='proportional',extendrect = True)
cbar.add_lines(CS)
cbar.ax.set_yticklabels([x for x in [3e-5,3e-4,0.15]])
cbar.ax.set_ylabel( 'averaged pulse density ' + r'(min$^{-1})$',rotation=270,fontsize=8)



cbar.ax.tick_params(labelsize=6,length=2)
cbar.ax.yaxis.label.set_size(8)

cbar.ax.yaxis.set_label_coords(10,0.5)
cbar.ax.tick_params(labelsize=6,direction='out', pad=1,length=2)
# =============================================================================
# boxplot
# =============================================================================

data_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/files/'
save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/'
dataset_name='2019_E8_2_WTonly'

target_conc = 'WT'

#sns.set(context='talk', style='white')
#plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
#level_description= ['Threshold chriteria 1e-5','Threshold chriteria 1e-4','Threshold chriteria 1e-3','Threshold chriteria 5e-3','Threshold chriteria 1e-2']
level_description= ['Threshold chriteria 1e-5','Threshold chriteria 1e-4','Threshold chriteria 1e-3','Threshold chriteria 5e-3','Threshold chriteria 1e-2']

y = [1e-5,1e-4,1e-3,5e-3,1e-2]
c_levels= [1e-4]

ax2 = plt.subplot(gs_row5[1]);

for k,l in enumerate(['level_1']):
    #ax2.text(0.85, 0.95, level_description[k], ha='center', va='center', transform=ax2.transAxes, fontsize=10)
    filename = data_folder+'df_results_curve_'+target_conc+'_'+l+'.pkl'
    infile = open(filename,'rb')
    results = pickle.load(infile)
    infile.close()

    filename = data_folder+lcurves_condition+'_'+l+'.pkl'
    infile = open(filename,'rb')
    level = pickle.load(infile)
    infile.close()

    label = [tuple([int(i[0]),int(i[1])]) for i in level]
    box= [results[i].amp_peaks.groupby(level='cell').count()/results[i].FRAME.groupby(level=0).count() for i in results]
    X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
        
    ax2.plot(26, 0.017, marker="v", markersize=4, alpha=1, color="red",linewidth = 0)
    p1 = ax2.bar(np.arange(len(box)),height = [np.mean(i) for i in box], yerr = [np.std(i)/np.sqrt(len(i)) for i in box],color=colors[1],alpha=0.5,linewidth=0.0,width=0.8)
    
    
    ax2.tick_params(axis='x', labelsize=6,direction='out', pad=1,length=2); ax2.tick_params(axis='y', labelsize=6,direction='out', pad=1,length=2)

 
    ax2.set_ylim([0,2e-2]); #ax2.set_yticks([0,0.01,0.03,0.05,0.1])
    ax2.set_xlim([-1,len(label)+1]) #[-1,len(label)+1]
    set_scale(ax2,np.arange(0,len(label),8),[0,2e-2])

    label2 = [(np.round(l[0]/100/20*60,2),np.round(l[1]/100,2)) for l in label]
    label3=[label2[i] for i in np.arange(0,len(label),8)]
    ax2.set_xticklabels(label3,rotation = 0)
    ax2.set_yticklabels([0,2e-2/20*60],rotation = 0)



    ax2.set_ylabel( 'averaged \n pulse density \n ' + r'(min$^{-1})$',fontsize=8,rotation = 90)
    ax2.set_xlabel(r'$(v_{th}$ (a.u./min) , $A_{th}$ (a.u.))',fontsize=8);
    ax2.xaxis.set_label_coords(0.5, -0.35);
    ax2.yaxis.set_label_coords(-0.03,0.5)
    ax2.invert_xaxis()
ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)

save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'
plt.savefig(save_folder+'Fig2Sup1.pdf', format='pdf')


    #%%

