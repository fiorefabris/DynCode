import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 
#from supfig_peaks.py import peaks####


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
matplotlib.rcParams['savefig.transparent'] = True

import pickle
import sys



sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
 


sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')

from plotting_main import load_file, save_file, get_linear_detrend, make_square_axes, set_scale_bars , peaks
        #%%    DOWNLOAD DATA
filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/2019_E8_6.pkl'
main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/'
conc = load_file(filename)

key = ['an_0ng','an_2-5ng','an_5ng','an_20ng']
keywords = ['NC','2-5ng','5ng','20ng']
TH =  ('172.41379310344828', '1847.3730489010406')
amp_th = 1847
slope_th=172

conc_data = {}
for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)

#%%

df_dic_ = {}
th = [172.41379310344828, 1847.3730489010406]
if False:
    for c in conc_data.keys():
        df_dic_[c] = peaks(conc_data.copy(),[c,th])
    save_file(df_dic_,main_folder + 'df_dic_.pkl')

#%%

# EL GRIDSPECT
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/figures/'

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


cells_conc = {'an_0ng':[1,32,53,55,59],'an_2-5ng':[0,8,19,33,15],'an_5ng':[1,4,49,53,55],'an_20ng':[2,6,17,26,64]} 
cells_conc = {'an_0ng':[59],'an_2-5ng':[19],'an_5ng':[1],'an_20ng':[2]} 

conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

filter_labels = ['raw data','smoothing', 'max and min', 'amplitude threshold', 'pulses']

box=[];handles = [None]*len(conc_labels)
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]


plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
#df = conc_data[c];df.sort_index(inplace=True); 


ax = plt.subplot(gs_row0[-1])
ax.text(-0.4, 0.9,  filter_labels[0], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row0[k])
        ax.plot(df.loc[cells_conc[c][w]]['MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)

        if w == 0: 
            ax.text(0.05, 1.5,  conc_labels[c], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = colors[k])


# =============================================================================
# # smooth
# =============================================================================

ax = plt.subplot(gs_row1[-1])
ax.text(-0.4, 0.9,   filter_labels[1], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row1[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)

        
# =============================================================================
# # max and min
# =============================================================================

ax = plt.subplot(gs_row2[-1])
ax.text(-0.4, 0.9,  filter_labels[2], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY_PEAKS_O1'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['MIN_O1'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')

        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        set_scale_bars(ax)



# =============================================================================
# # Intermidate step
# =============================================================================

ax = plt.subplot(gs_row2b[-1])
ax.text(-0.4, 0.9,   filter_labels[3], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)
df_dic_ = load_file(main_folder + 'df_dic_.pkl')

for k,c in enumerate(cells_conc):
    df_ = df_dic_[c];df_.sort_index(inplace=True); color=colors[k]

    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2b[k])
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
silent_ax(ax)
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row3[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=5, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=4, alpha=1,label = 'Peak minima')

        ax.set_xlim([0, M]);
        ax.set_ylim(ylim)
        
        
        set_scale_bars(ax, x_bar=(0,60), y_bar=ylim, xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
        if w==0 and k==0:
           ax.set_ylabel('KTR signal \n 200 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
           ax.xaxis.set_label_coords(0.3, -0.2);
           ax.yaxis.set_label_coords(-0.05,0.6);
           ax.tick_params(labelsize=6,direction='out', pad=1,length=2)


##########################################
### CHARGING THE FILES ###################
##########################################
           


data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/'
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/figures/'
dataset_name='2019_E8_6'


filename = data_folder + dataset_name + '.pkl'
conc_df = load_file(filename)

filename = data_folder + 'df_results_20ng_Norm_.pkl'
df_result_20ng = load_file(filename)

filename = data_folder + 'df_results_0ng_Norm_.pkl'
df_result_NC = load_file(filename)

results = [df_result_NC,df_result_20ng]
conditions = ['NC','20ng']

lcurves_condition = 'NC'
target_conc = '20ng'
lcurves_values = df_result_NC.values #curvas de nivel
    

# =============================================================================
# Plot levelcurves
# =============================================================================

for w, heatmap_values_ in enumerate([df_result_NC,df_result_20ng]):
    
    x = np.array(heatmap_values_.columns) #pendiente
    y = np.array(heatmap_values_.index) #amplitú

    xx,yy = np.meshgrid(x,y)
    Z=heatmap_values_.values
    
    ax = plt.subplot(gs_row4[1,w])

    plt.grid(0)
    levels= np.geomspace(1e-5,1e-1,num=50)
    levels=np.arange(1e-5,0.05+5e-5,5e-7)
    
    CSF = ax.contourf(xx, yy, Z,levels=levels,extend='both',cmap='Greys',vmin=1e-5,
                      vmax=0.05,norm=matplotlib.colors.SymLogNorm(vmin=1e-4,vmax=0.05,linthresh=1e-4,linscale=0))
    c_levels= [5e-3]
    
    CS = ax.contour(xx, yy, lcurves_values,levels=c_levels,norm=matplotlib.colors.SymLogNorm(vmin=5e-4,vmax=0.05,linthresh=5e-4), colors='lightblue',linewidths=1)
    ax.tick_params(axis='x',rotation=0,length=2)
    ax.plot([slope_th], [amp_th], marker="v", markersize=4, alpha=1, color="red",linewidth = 0)

    
    set_scale(ax,np.arange(0,500,100),[0,3000])
    ax.set_yticklabels([0/100,3000/100],fontsize=6)
    ax.set_xticklabels(np.arange(0,500,100)/100/20*60,fontsize=6)
    ax.tick_params(labelsize=6,direction='out', pad=1,length=2)


    
    if w==0:
        ax.set_ylabel(r'$A_{th}$ (a.u.)',fontsize=8); 
        ax.set_xlabel(r'$v_{th}$ (a.u./min)',fontsize=8);
        ax.xaxis.set_label_coords(0.5,-0.17);
        ax.yaxis.set_label_coords(-0.05,0.5);
    

ax_cbar = plt.subplot(gs_row4[1,w+1])
cbar = plt.colorbar(CSF, cax=ax_cbar,label= 'averaged pulse density ' + r'$(min^{-1})$',
                        ticks= [5e-4,5e-3,0.05],spacing='proportional',extendrect = True)
cbar.add_lines(CS)
cbar.ax.set_yticklabels([x for x in [15e-4, 0.015,0.15]])
cbar.ax.set_ylabel( 'averaged pulse density ' + r'(min$^{-1})$',rotation=270,fontsize=8)


cbar.ax.tick_params(labelsize=6,length=2)
cbar.ax.yaxis.label.set_size(8)

cbar.ax.yaxis.set_label_coords(10,0.5)
cbar.ax.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# boxplot
# =============================================================================


dataset_name='2019_E8_6'
target_conc = '20ng'

y = [1e-5,1e-4,1e-3,5e-3,1e-2]

c_levels= [5e-3]
ax2 = plt.subplot(gs_row5[1]);

for k,l in enumerate(['level_1']):
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
            
    ax2.plot(31, 0.026, marker="v", markersize=4, alpha=1, color="red",linewidth = 0)
    p1 = ax2.bar(np.arange(len(box)),height = [np.mean(i) for i in box], yerr = [np.std(i)/np.sqrt(len(i)) for i in box],color=colors[-1],alpha=0.5,linewidth=0.0,width=0.8)    

    ax2.tick_params(axis='x', labelsize=6,direction='out', pad=1); ax2.tick_params(axis='y', labelsize=6,direction='out', pad=1)

 
    ax2.set_ylim([0,3e-2]); 
    ax2.set_xlim([-1,len(label)+1])
    set_scale(ax2,np.arange(0,len(label),10),[0,3e-2])

    label2 = [(np.round(l[0]/100/20*60,2),np.round(l[1]/100,2)) for l in label]
    label3=[label2[i] for i in np.arange(0,len(label),10)]
    ax2.set_xticklabels(label3,rotation = 0)
    ax2.set_yticklabels([0,3e-2/20*60],rotation = 0)



    ax2.set_ylabel( 'averaged \n pulse density \n ' + r'(min$^{-1})$',fontsize=8,rotation = 90)
    ax2.set_xlabel(r'$(v_{th}$ (a.u./min) , $A_{th}$ (a.u.))',fontsize=8);
    ax2.xaxis.set_label_coords(0.5, -0.35);
    ax2.yaxis.set_label_coords(-0.03,0.5)
    ax2.invert_xaxis()

    
ax2.tick_params(labelsize=6,direction='out', pad=1)


save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/'
plt.savefig(save_folder+'F3Sup3.pdf', format='pdf')