import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import mode 
import sys


from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale
from sklearn.linear_model import LinearRegression
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import matplotlib

plt.rcdefaults()
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.5

matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.minor.size'] = 1

matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['xtick.minor.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5


def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z


sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
 

sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/')

from plotting_main import load_file, save_file, get_linear_detrend, make_square_axes, set_scale_bars , peaks,load_file5



#%% EPISC


epi_main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/EpiSC/data/'
epi_conc_data = {}
epi_color = sns.color_palette()[5]
epi_key =[ 'an_EpiSC_FAX_PD03', 'an_EpiSC_FAX']
keywords = ['CN','WT']
epi_TH = ('33.1315210671255', '305.08474576271186')
epi_slope_th = 33
epi_amp_th = 305

for j,k in enumerate(keywords):
    filename = epi_main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    epi_conc_data[epi_key[j]] = aux[epi_TH]
get_linear_detrend(epi_conc_data)

if False:
    epi_df_dic_ = {}
    epi_th = [33.1315210671255, 305.08474576271186]
    for c in epi_conc_data.keys():
        epi_df_dic_[c] = peaks(epi_conc_data.copy(),[c,epi_th])
    save_file(epi_df_dic_,epi_main_folder + 'epi_df_dic_pulses.pkl')



#%%ESC

esc_main_folder ='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/'
esc_conc_data = {}
esc_color = sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
esc_key =[ 'an_ESC_FAX_PD03', 'an_ESC_FAX']
keywords = ['CN','WT']
esc_TH =('57.25643730066835', '406.77966101694915')
esc_slope_th = 57
esc_amp_th = 406

for j,k in enumerate(keywords):
    filename = esc_main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    esc_conc_data[esc_key[j]] = aux[esc_TH]
get_linear_detrend(esc_conc_data)


if False:
    esc_df_dic_ = {}
    esc_th = [57.25643730066835, 406.77966101694915]
    for c in esc_conc_data.keys():
        esc_df_dic_[c] = peaks(esc_conc_data.copy(),[c,esc_th])
    save_file(esc_conc_data,esc_main_folder + 'esc_df_dic_pulses.pkl')


#%%

# EL GRIDSPECT
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/'

#Figura principal
fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
gs_main = gridspec.GridSpec(nrows=1, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.90, hspace=0.0)

#3 filas en toda la hoja
gs_main_cols = gridspec.GridSpecFromSubplotSpec(nrows=7, ncols=1, subplot_spec=gs_main[0], wspace=0.4, hspace=0.0, height_ratios=[0.1,0.1,0.1,0.1,0.1,1,0.5])

#%Para cada fila, 4 columnas
gs_row0  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=6, subplot_spec=gs_main_cols[0], wspace=0.2, hspace=0.2,width_ratios = [1,1,0.2,1,1,0.5])
gs_row1  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=6, subplot_spec=gs_main_cols[1], wspace=0.2, hspace=0.2,width_ratios = [1,1,0.2,1,1,0.5])
gs_row2  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=6, subplot_spec=gs_main_cols[2], wspace=0.2, hspace=0.2,width_ratios = [1,1,0.2,1,1,0.5])
gs_row2b = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=6, subplot_spec=gs_main_cols[3], wspace=0.2, hspace=0.2,width_ratios = [1,1,0.2,1,1,0.5])
gs_row3  = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=6, subplot_spec=gs_main_cols[4], wspace=0.2, hspace=0.2,width_ratios = [1,1,0.2,1,1,0.5])
gs_row4  = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=7, subplot_spec=gs_main_cols[5], wspace=0.4, hspace=0.2,height_ratios = [0.3,0.7],width_ratios = [1,1,0.1,0.1,1,1,0.1])
gs_row5  = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main_cols[6], wspace=0.2, hspace=0.2,height_ratios = [0.4,0.6])
for row in [gs_row0,gs_row1,gs_row2,gs_row2b,gs_row3]:
    ax = plt.subplot(row[2])
    ax.set_visible(False)

gs_row0  = [gs_row0[0],gs_row0[1],gs_row0[3],gs_row0[4],gs_row0[5]]
gs_row1  = [gs_row1[0],gs_row1[1],gs_row1[3],gs_row1[4],gs_row1[5]]
gs_row2  = [gs_row2[0],gs_row2[1],gs_row2[3],gs_row2[4],gs_row2[5]]
gs_row2b  = [gs_row2b[0],gs_row2b[1],gs_row2b[3],gs_row2b[4],gs_row2b[5]]
gs_row3  = [gs_row3[0],gs_row3[1],gs_row3[3],gs_row3[4],gs_row3[5]]
gs_row4  = [gs_row4[1,0],gs_row4[1,1],gs_row4[1,2],gs_row4[1,4],gs_row4[1,5],gs_row4[1,6]]
gs_row5 = [gs_row5[1,0],gs_row5[1,2]]

#%%


# =============================================================================
# Selected traces
# =============================================================================
cells_conc = {'an_EpiSC_FAX_PD03':[47],'an_EpiSC_FAX':[27],'an_ESC_FAX_PD03':[28],'an_ESC_FAX':[8]} 

conc_labels = {'an_EpiSC_FAX_PD03':'MEKi ','an_EpiSC_FAX': 'EpiSC ', 
               'an_ESC_FAX_PD03': 'MEKi','an_ESC_FAX': 'ESC'}
filter_labels = ['raw data','smoothing', 'max and min', 'amplitude threshold', 'pulses']

conc_data = merge_two_dicts(epi_conc_data,esc_conc_data)

colors = [epi_color,esc_color]
Y_lim = [[62000,65000],[60000,64000]]
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
ax = plt.subplot(gs_row0[-1])
ax.text(-0.2, 0.9,  filter_labels[0], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

box=[];handles = [None]*len(conc_labels)


#for conc_data in [epi_conc_data,esc_conc_data]:
    
M = max([conc_data[key].max().FRAME for key in conc_data.keys()])

for k,c in enumerate(cells_conc):
    print(c)
    df = conc_data[c];df.sort_index(inplace=True); #color=colors[k]
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row0[k])
        ax.plot(df.loc[cells_conc[c][w]]['MEAN_INTENSITY'].values, color=colors[k//2], linewidth=0.75);
        ax.set_xlim([0, M]);
        ax.set_ylim(Y_lim[k//2])
        set_scale_bars(ax)

        if w == 0: 
            ax.text(0.05, 1.5,  conc_labels[c], ha='left', va='center', transform=ax.transAxes, fontsize=8,color = colors[k//2])
# =============================================================================
# # smooth
# =============================================================================

ax = plt.subplot(gs_row1[-1])
ax.text(-0.2, 0.9,   filter_labels[1], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); 
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row1[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k//2], linewidth=0.75);
        ax.set_xlim([0, M]);
        ax.set_ylim(Y_lim[k//2])
        set_scale_bars(ax)

        
# =============================================================================
# # max and min
# =============================================================================

ax = plt.subplot(gs_row2[-1])
ax.text(-0.2, 0.9,  filter_labels[2], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); 
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k//2], linewidth=0.75);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY_PEAKS_O1'].values,linestyle="None", marker = "." , color='blue', markersize=1, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['MIN_O1'].values,linestyle="None", marker = ".", color='black', markersize=1, alpha=1,label = 'Peak minima')

        ax.set_xlim([0, M]);
        ax.set_ylim(Y_lim[k//2])
        set_scale_bars(ax)



# =============================================================================
 # Intermidate step
# =============================================================================

ax = plt.subplot(gs_row2b[-1])
ax.text(-0.2, 0.9,   filter_labels[3], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
set_scale_bars(ax)

epi_conc_data_pulses = load_file(epi_main_folder + 'epi_df_dic_pulses.pkl')
esc_conc_data_pulses =  load_file(esc_main_folder + 'esc_df_dic_pulses.pkl')
conc_data_pulses = merge_two_dicts(epi_conc_data,esc_conc_data)


for k,c in enumerate(conc_data_pulses):
    df_ = conc_data_pulses[c];df_.sort_index(inplace=True); color=colors[k//2]

    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row2b[k])
        ax.plot(df_.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k//2], linewidth=0.75);
        handles[0], = ax.plot(df_.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=1, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df_.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=1, alpha=1,label = 'Peak minima')

        ax.set_xlim([0, M]);
        ax.set_ylim(Y_lim[k//2])
        set_scale_bars(ax)




# =============================================================================
# # Maximo final
# =============================================================================
ax = plt.subplot(gs_row3[-1])
ax.text(-0.2, 0.9,   filter_labels[4], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = 'k')
silent_ax(ax)
set_scale_bars(ax)

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); 
    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(gs_row3[k])
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k//2], linewidth=0.75);
        handles[0], = ax.plot(df.loc[cells_conc[c][w]]['max_'].values,linestyle="None", marker = "." , color='blue', markersize=1, alpha=1,label = 'Peak maxima',markeredgecolor='blue')
        handles[1], = ax.plot(df.loc[cells_conc[c][w]]['min_'].values,linestyle="None", marker = ".", color='black', markersize=1, alpha=1,label = 'Peak minima')

        ax.set_xlim([0, M]);
        ax.set_ylim(Y_lim[k//2])
        
        
        set_scale_bars(ax, x_bar=(0,60), y_bar=Y_lim[k//2], xunits='', yunits=r'$\times 10^2$', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
        if w==0 and k==0:
           ax.set_ylabel('KTR signal \n 30 a.u.',fontsize=10); ax.set_xlabel('20 minutes',fontsize=10)
           ax.xaxis.set_label_coords(0.3, -0.2);
           ax.yaxis.set_label_coords(-0.05,0.6);
           ax.tick_params(labelsize=8,direction='out', pad=1,length=2)

        if k==2:
           ax.set_ylabel('KTR signal \n 40 a.u.',fontsize=10); ax.set_xlabel('20 minutes',fontsize=10)
           ax.xaxis.set_label_coords(0.3, -0.2);
           ax.yaxis.set_label_coords(-0.05,0.6);
           ax.tick_params(labelsize=8,direction='out', pad=1,length=2)




           
# =============================================================================
# Plot levelcurves EPI
# =============================================================================
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.minor.size'] = 1

matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['xtick.minor.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
           
           
##########################################
### CHARGING THE FILES ###################
##########################################

epi_data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/EpiSC/data/'

filename = epi_main_folder + 'df_results_WT_Norm_.pkl'
epi_df_result_WT = load_file(filename)


filename = epi_main_folder + 'df_results_CN_Norm_.pkl'
epi_df_result_NC = load_file(filename)


lcurves_condition = 'NC'
target_conc = 'WT'
epi_lcurves_values = epi_df_result_NC.values #curvas de nivel




# =============================================================================
# Plot levelcurves
# =============================================================================

for w, heatmap_values_ in enumerate([epi_df_result_NC,epi_df_result_WT]):
    
    x = np.array(heatmap_values_.columns) #pendiente
    y = np.array(heatmap_values_.index) #amplitú

    xx,yy = np.meshgrid(x,y)
    Z=heatmap_values_.values
    
    ax = plt.subplot(gs_row4[w])

    plt.grid(0)
    #levels= np.geomspace(1e-5,1e-1,num=50)
    levels=np.arange(5e-6,9e-3,5e-8)
    
    CSF = ax.contourf(xx, yy, Z,levels=levels,extend='both',cmap='Greys',vmin=1e-5,
                      vmax=0.005,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.005,linthresh=1e-5,linscale=0))
    c_levels= [1e-3]
    
    CS = ax.contour(xx, yy, epi_lcurves_values,levels=c_levels,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.005,linthresh=1e-5), colors='brown',linewidths=0.75)
    ax.tick_params(axis='x',rotation=0,length=2)
    ax.plot([epi_slope_th], [epi_amp_th], marker="v", markersize=4, alpha=1, color="red",linewidth = 0)

    ax.set_ylim([0,1500]);ax.set_xlim([0,150])
    set_scale(ax,np.arange(0,500,50),[0,1500])
    ax.set_yticklabels([0/100,1500/100],fontsize=8)
    ax.set_xticklabels(np.arange(0,500,50)/100/20*60,fontsize=8)
    ax.tick_params(labelsize=8,direction='out', pad=1,length=2)


    
   # if w==0:
    if False:
        ax.set_ylabel(r'$A_{th}$ (a.u.)',fontsize=10); 
        ax.set_xlabel(r'$v_{th}$ (a.u./min)',fontsize=10);
        ax.xaxis.set_label_coords(0.5,-0.17);
        ax.yaxis.set_label_coords(-0.05,0.5);
    

ax_cbar = plt.subplot(gs_row4[w+1])
cbar = plt.colorbar(CSF, cax=ax_cbar,label= '',
                        ticks= [5e-5,5e-4,0.005],spacing='proportional',extendrect = True)
cbar.add_lines(CS)
cbar.ax.set_yticklabels([x for x in [15e-5, 15e-4,0.015]])
#cbar.ax.set_ylabel( 'averaged pulse density ' + r'(min$^{-1})$',rotation=270,fontsize=10)


cbar.ax.tick_params(labelsize=6,length=2)
cbar.ax.yaxis.label.set_size(8)

cbar.ax.yaxis.set_label_coords(10,0.5)
cbar.ax.tick_params(labelsize=8,direction='out', pad=1,length=2)
silent_ax(ax_cbar)

# =============================================================================
# boxplot
# =============================================================================
lcurves_condition = 'CN'


y = [1e-5,1e-4,1e-3,5e-3,1e-2]
ax2 = plt.subplot(gs_row5[0]);

for k,l in enumerate(['level_1']):
    filename = epi_data_folder+'df_results_curve_'+target_conc+'_'+l+'.pkl'
    epi_results = load_file(filename)


    filename = epi_data_folder+lcurves_condition+'_'+l+'.pkl'
    epi_level = load_file(filename)


    label = [tuple([int(i[0]),int(i[1])]) for i in epi_level]
    box= [epi_results[i].amp_peaks.groupby(level='cell').count()/epi_results[i].FRAME.groupby(level=0).count() for i in epi_results]
    X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
            
    ax2.plot(10, 0.026, marker="v", markersize=4, alpha=1, color="red",linewidth = 0)
    p1 = ax2.bar(np.arange(len(box)),height = [np.mean(i) for i in box], yerr = [np.std(i)/np.sqrt(len(i)) for i in box],color=esc_color,alpha=0.5,linewidth=0.0,width=0.8)    

    ax2.tick_params(axis='x', labelsize=8,direction='out', pad=1); ax2.tick_params(axis='y', labelsize=8,direction='out', pad=1)

 
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

    
ax2.tick_params(labelsize=8,direction='out', pad=1)


# =============================================================================
# Plot levelcurves ESC
# =============================================================================
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.minor.size'] = 1

matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['xtick.minor.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
           
##########################################
### CHARGING THE FILES ###################
##########################################

esc_main_folder ='/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/'
esc_conc_data = {}
esc_color = sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]


filename = esc_main_folder + 'df_results_WT_Norm_.pkl'
esc_df_result_WT = load_file5(filename)


filename = esc_main_folder + 'df_results_CN_Norm_.pkl'
esc_df_result_NC = load_file5(filename)


lcurves_condition = 'NC'
target_conc = 'WT'
esc_lcurves_values = esc_df_result_NC.values #curvas de nivel




# =============================================================================
# Plot levelcurves
# =============================================================================

for w, heatmap_values_ in enumerate([esc_df_result_NC,esc_df_result_WT]):
    
    x = np.array(heatmap_values_.columns) #pendiente
    y = np.array(heatmap_values_.index) #amplitú

    xx,yy = np.meshgrid(x,y)
    Z=heatmap_values_.values
    
    ax = plt.subplot(gs_row4[w+3])

    plt.grid(0)
    #levels= np.geomspace(1e-5,1e-1,num=50)
    levels=np.arange(5e-6,9e-2,5e-7)
    
    CSF = ax.contourf(xx, yy, Z,levels=levels,extend='both',cmap='Greys',vmin=1e-5,
                      vmax=0.05,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.05,linthresh=1e-5,linscale=0))
    c_levels= [1e-3]
    
    CS = ax.contour(xx, yy, esc_lcurves_values,levels=c_levels,norm=matplotlib.colors.SymLogNorm(vmin=1e-5,vmax=0.005,linthresh=1e-5), colors='green',linewidths=0.75)
    ax.tick_params(axis='x',rotation=0,length=2)
    ax.plot([esc_slope_th], [esc_amp_th], marker="v", markersize=4, alpha=1, color="red",linewidth = 0)

    ax.set_ylim([0,1500]);ax.set_xlim([0,150])
    set_scale(ax,np.arange(0,500,50),[0,1500])
    ax.set_yticklabels([0/100,1500/100],fontsize=8)
    ax.set_xticklabels(np.arange(0,500,50)/100/20*60,fontsize=8)
    ax.tick_params(labelsize=8,direction='out', pad=1,length=2)

    
    #if w==0:
    if False:
        ax.set_ylabel(r'$A_{th}$ (a.u.)',fontsize=10); 
        ax.set_xlabel(r'$v_{th}$ (a.u./min)',fontsize=10);
        ax.xaxis.set_label_coords(0.5,-0.17);
        ax.yaxis.set_label_coords(-0.05,0.5);
    

ax_cbar = plt.subplot(gs_row4[w+4])
cbar = plt.colorbar(CSF, cax=ax_cbar,label= '',
                        ticks= [5e-5,5e-4,0.005],spacing='proportional',extendrect = True) #averaged pulse density ' + r'$(min^{-1})$
cbar.add_lines(CS)
cbar.ax.set_yticklabels([x for x in [15e-5, 15e-4,0.015]])
#cbar.ax.set_ylabel( 'averaged pulse density ' + r'(min$^{-1})$',rotation=270,fontsize=10)


cbar.ax.tick_params(labelsize=6,length=2)
cbar.ax.yaxis.label.set_size(8)

cbar.ax.yaxis.set_label_coords(10,0.5)
cbar.ax.tick_params(labelsize=8,direction='out', pad=1,length=2)





# =============================================================================
# boxplot
# =============================================================================
lcurves_condition = 'CN'


y = [1e-5,1e-4,1e-3,5e-3,1e-2]
ax2 = plt.subplot(gs_row5[1]);

for k,l in enumerate(['level_1']):
    filename = esc_main_folder+'df_results_curve_'+target_conc+'_'+l+'.pkl'
    esc_results = load_file(filename)


    filename = esc_main_folder+lcurves_condition+'_'+l+'.pkl'
    esc_level = load_file(filename)


    label = [tuple([int(i[0]),int(i[1])]) for i in esc_level]
    box= [esc_results[i].amp_peaks.groupby(level='cell').count()/esc_results[i].FRAME.groupby(level=0).count() for i in esc_results]
    X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
            
    ax2.plot(10, 0.026, marker="v", markersize=4, alpha=1, color="red",linewidth = 0)
    p1 = ax2.bar(np.arange(len(box)),height = [np.mean(i) for i in box], yerr = [np.std(i)/np.sqrt(len(i)) for i in box],color=esc_color,alpha=0.5,linewidth=0.0,width=0.8)    

    ax2.tick_params(axis='x', labelsize=8,direction='out', pad=1); ax2.tick_params(axis='y', labelsize=8,direction='out', pad=1)

 
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

    
ax2.tick_params(labelsize=8,direction='out', pad=1)

plt.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/trial.pdf', format='pdf')
