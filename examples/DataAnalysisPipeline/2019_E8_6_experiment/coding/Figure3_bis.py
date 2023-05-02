import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker


from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

from sklearn.linear_model import LinearRegression
from scipy.stats import mode 
from DataAnalysis.Detrending.PlotDetrend import plot_detrend
from DataAnalysis.Preprocesing.PlotData import grid_plot

import pickle
import sys

import matplotlib
plt.rcdefaults()

matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.75

matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.major.width'] = 0.75
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['xtick.minor.width'] = 0.75


sys.path.append('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from consecutive_main import consecutive_non_cumulative,makeColours,consecutive_cumulative,mask_low_variance
from matplotlib import cm

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
%matplotlib inline

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
        
        
# te da 1 si es consecutivo y 0 si no, con el criterio de medio dt
def is_consecutive(df):
    ''' Function that determines if an IPI is between two consecutive pulses (return 1) or not (return 0).
    An IPI is between two consecutive pulses if the silent region between pulses is smaller than 0.5* mixed dt
    -------------------------
    IPUNT:
        df: dataframe
        
    -------------------------
    OUTPUT:
        consecutive: a list with the IPI len'''
    
    y = df['IPI'].dropna().values-mixed_dt(df) #silencios entre picos
    x = mixed_dt(df) 
    consecutive = np.zeros(len(x))
    for i in range(len(x)):
        if x[i]*0.5 >= y[i]:
            consecutive[i] = 1
    return(consecutive)
    
#Time correlation 3: con el dt miti miti

def mixed_dt(df):
    values_dt_mixed = []
    for cells, data in df.groupby(level='cell'):
        cell_t_M = data[data['max_'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
            #print(dt_raise,dt_fall)
    return(values_dt_mixed)

#%% 

#%%    
def is_consecutive(df):
    # discrimina los intervalos de consecutividad vs no consecutividad
    y = df['IPI'].dropna().values-mixed_dt(df)
    x = mixed_dt(df)
    consecutive = np.zeros(len(x))
    for i in range(len(x)):
        if x[i]*0.5 >= y[i]:
            consecutive[i] = 1
    return(consecutive)


def is_consecutive_cell(df):
    #solo pasar datafranme de una celula (creo)
    #le agrega el primero a cada tren de pulsos
    #PROBLEMA: tiene en cuenta las celulas que tienen solo un pulso
    y = df['IPI'].dropna().values-mixed_dt(df)
    x = mixed_dt(df)
    consecutive = []
    for i in range(len(x)):
        if x[i]*0.5 >= y[i]:
            if (len(consecutive) ==0 or consecutive[-1] ==0): consecutive.append(1)
            consecutive.append(1)
        else:
            consecutive.append(0) 
    return(consecutive)

    
def is_isolated_cell(df):
    #solo pasar datafranme de una celula (creo)
    #calcula el numero de picos no consecutivis x concentracion
    # tiene en cuenta las celulas de un solo pulso
    total_peaks = df[df['amp_peaks'].notna()].amp_peaks.count()
    consecutive_ = sum(is_consecutive_cell(df))
    isolated_peaks = total_peaks - consecutive_
    return(isolated_peaks)
    

    
#%%    
filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/2019_E8_6.pkl'
conc = load_file(filename)
filename_aux= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/df_results_curve_NC_level_1.pkl'
aux = load_file(filename_aux)


#%%
main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/files/'
conc_data = {}
key = ['an_0ng','an_2-5ng','an_5ng','an_20ng']
keywords = ['NC','2-5ng','5ng','20ng']
TH =  ('172.41379310344828', '1847.3730489010406')
amp_th = 1847
slope_th=172
for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_level_1.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)

#for key in list(conc_data.keys()):
#    dm =find_dm(conc_data[key])
    #%%
    
    
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/analysis_2019_E8_6/figures/' 



fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
#Main 6 filas

gs_main0 = gridspec.GridSpec(nrows=2, ncols=1, figure=fig)#,width_ratios= [1,4],height_ratios=[1,4]); 
gs_main0.update(left=0.1, right=0.95, bottom=0.1, top=1, hspace=0.0,wspace = 0.0)

gs_main = gridspec.GridSpecFromSubplotSpec(nrows=6, ncols=1, subplot_spec=gs_main0[0], wspace=0.0, hspace=0.5)

#gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); 
#gs_main.update(left=0.1, right=0.95, bottom=0.1, top=1, hspace=0.5,wspace = 0.0)
gs_main_ = [gs for gs in gs_main]


# =============================================================================
# 
# =============================================================================
#Lugar 0: esquema 
#lugar 1 traza2
#gs_0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main_[1], wspace=0.5, hspace=0.5)





ylim = [43500, 62800]
M = max([conc_data[key].max().FRAME for key in conc_data.keys()])

cells_conc = {'an_0ng':[1,32,53,55,59],'an_2-5ng':[0,8,19,33,15],'an_5ng':[1,4,49,53,55],'an_20ng':[2,6,17,26,64]} 
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

box=[];handles = [None]*len(conc_labels)
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]



# =============================================================================
# trazas
# =============================================================================
inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=len(cells_conc), subplot_spec=gs_main_[1], wspace=0.2, hspace=0.4)

plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 

for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]

    inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc[c]), ncols=1, subplot_spec=inner_gs0[k], wspace=0.0, hspace=0.0)

    for w,z in enumerate(cells_conc[c]):
        ax = plt.subplot(inner_inner_gs0[w]);       
        ax.patch.set_alpha(0.0)
        ax.patch.set_linewidth(0.0)
    
        ax.plot(df.loc[cells_conc[c][w]]['sm_MEAN_INTENSITY'].values, color=colors[k], linewidth=0.7);
        ax.set_xlim([0, M]);
        ylim_=[43000,63000]
        ax.set_ylim(ylim_)

        #if (k==len(cells_conc)-1) and (w == len(cells_conc[c])-1) :
        if (k==0) and (w == len(cells_conc[c])-1):   
           set_scale_bars(ax, x_bar=(0,60), y_bar=ylim_, xunits='', yunits='', x_scale=60/20, y_scale=100, round_x = True,round_y=True)
           ax.set_ylabel('     KTR signal \n 200 a.u.',fontsize=8); ax.set_xlabel('20 minutes',fontsize=8)
           ax.xaxis.set_label_coords(0.3, -0.5);
           ax.yaxis.set_label_coords(-0.08,2.2);
           
        elif (k!=0) and (w == len(cells_conc[c])-1): 
            set_scale_bars(ax, x_bar=(0,60), y_bar=ylim_, xunits='', yunits='', x_scale=60/20, y_scale=100, round_x = True,round_y=True)

        else:
            set_scale_bars(ax)

        if w == 0: 
            ax.text(0.05, 1.5,  conc_labels[c], ha='left', va='center', transform=ax.transAxes, fontsize=6,color = colors[k])
        


# =============================================================================
# statitics
# =============================================================================

ipi_conc = ['an_2-5ng','an_5ng','an_20ng']
gs3_inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[3:6], wspace=0.4, hspace=0.0)
gs0_inner = gridspec.GridSpecFromSubplotSpec(nrows=len(ipi_conc), ncols=3, subplot_spec=gs3_inner[0:3], wspace=0.4, hspace=0.4)
col3 = 'IPI'#ax1
col2='dt_peaks'#ax3
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 



for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax1 = plt.subplot(gs0_inner[k,1]); ax3 = plt.subplot(gs0_inner[k,0]); #ax4 = plt.subplot(gs0[k,2]);
    bins1= ax1.hist(df[col3],bins=np.arange(2,80,4),range=(2,None), density=1, facecolor=colors[k+1], alpha=0.8,linewidth=0); 
    bins3= ax3.hist(df[col2],bins=np.arange(2,46,2),range=(2,None), density=1, facecolor=colors[k+1], alpha=0.8,linewidth=0); 


    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax3.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.set_xlim([0,80]);ax1.set_ylim([0,0.07])
    ax3.set_xlim([0,46]);ax3.set_ylim([0,0.12]);


    #mode_IPI = mode(df[col3])[0][0]*20/60; mode_IPI = 'mode: '+str(np.round(mode_IPI,2))+' \n'
    mode_IPI =  (bins1[1][np.argmax(bins1[0])] + bins1[1][np.argmax(bins1[0])+1])/2 ; mode_IPI = 'mode: '+str(np.round(mode_IPI*20/60,2))+' \n'
    ax1.text(1, 0.80, mode_IPI, ha='right', va='center', transform=ax1.transAxes, fontsize=6)
    ax1.text(1, 1.05,r'$Q$: '+str(np.round(df[col3].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col3].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col3].quantile(0.75)*20/60,2)) , ha='right', va='center', transform=ax1.transAxes, fontsize=6)

    #mode_dt = mode(df[col2])[0][0]*20/60 ; mode_dt = 'mode: '+str(np.round(mode_dt,2))+' \n'
    mode_dt =  (bins3[1][np.argmax(bins3[0])] + bins3[1][np.argmax(bins3[0])+1])/2; mode_dt = 'mode: '+str(np.round(mode_dt*20/60,2))+' \n'
    ax3.text(1, 0.80, mode_dt, ha='right', va='center', transform=ax3.transAxes, fontsize=6)
    ax3.text(1, 1.05,r'$Q$: '+str(np.round(df[col2].quantile(0.25)*20/60,2))+' ; '+str(np.round(df[col2].quantile(0.50)*20/60,2))+' ; '
                                  +str(np.round(df[col2].quantile(0.75)*20/60,2)) , ha='right', va='center', transform=ax3.transAxes, fontsize=6)

    if k != len(ipi_conc)-1:
        set_scale(ax1,np.arange(0,80,30),[0,0.07]) ;ax1.set_xticklabels([None,None]);ax1.set_yticklabels([None,None]);
        set_scale(ax3,np.arange(0,46,15),[0,0.12]);ax3.set_xticklabels([None,None]);ax3.set_yticklabels([None,None]);
        
    
    ax3.text(0.05, 0.9,  conc_labels[c], ha='left', va='center', transform=ax3.transAxes, fontsize=6,color = colors[k+1])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
set_scale(ax1,np.arange(0,80,30),[0,0.07]); ax1.set_xlabel('IPI (min)',fontsize=8); ax1.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.5, -0.3);ax1.yaxis.set_label_coords(-0.01,0.5)
ax1.set_xticklabels([int(x) for x in np.arange(0,80,30)*20/60],fontsize=6)
ax1.set_yticklabels([np.round(x*3,2) for x in [0,0.07]],fontsize=8)


set_scale(ax3,np.arange(0,46,15),[0,0.12]); ax3.set_xlabel('duration (min)',fontsize=8); ax3.set_ylabel('probability\ndensity (1/min)',fontsize=8,rotation = 90)
ax3.set_xticklabels([int(x) for x in np.arange(0,46,15)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.3);ax3.yaxis.set_label_coords(-0.01,0.5)
ax3.set_yticklabels([np.round(x*3,2) for x in [0,0.12]],fontsize=8)

for ax_ in [ax1,ax3]: ax_.tick_params(labelsize=6,direction='out', pad=1,length=2)


    

# =============================================================================
# silence vs mixed dt
# =============================================================================


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    
    consecutive_ = is_consecutive(df)
    total_consecutive     = consecutive_.sum()
    total_non_consecutive = len(consecutive_) - total_consecutive

    ax3 = plt.subplot(gs0_inner[k,2]); 
    ax3.plot(df['IPI'].dropna().values-mixed_dt(df),mixed_dt(df),'o',markersize=1.5,alpha=0.5,color = colors[k+1])
    ax3.set_xlim([-5, 60]);
    ax3.set_ylim( [-5, 60])
    outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-mixed_dt(df),mixed_dt(df))]
    outliers_sum = np.sum([1 for i,j in outliers if i>60 or j >60])
    
    #print(c,sum(1* (df['IPI'].dropna().values-mixed_dt(df) < 0))) 
    
    print('outliers : ' + str(outliers_sum))
    print('total_data: ' + str(len(outliers)))
    
    X = np.arange(0,200)
    ax3.plot(X,X*2,0,color = 'black',alpha=0.8,linewidth=0.8,linestyle='--') 
   
    
    x_label = 'joint duration (min)'
    if k ==2: ax3.set_xlabel('silence (min)',fontsize=8); 
    if k ==2: ax3.set_ylabel(x_label,fontsize=8); 
    if k != len(ipi_conc)-1:
            set_scale(ax3,np.arange(0,240,30),np.arange(0,240,30));ax3.set_xticklabels([None,None]);ax3.set_yticklabels([None,None]);


    text_ = 'consecutive : '+ str(np.round(100*total_consecutive/len(consecutive_),1))+ ' %'
    ax3.text(1.15, 1.05, text_, ha='right', va='center', transform=ax3.transAxes, fontsize=6,color = 'black')
    ax3.set_aspect(1)
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
set_scale(ax3,np.arange(0,240,30),np.arange(0,240,30))
ax3.set_xticklabels([int(x) for x in np.arange(0,240,30)*20/60],fontsize=6)
ax3.set_yticklabels([int(x) for x in np.arange(0,240,30)*20/60],fontsize=6)
ax3.xaxis.set_label_coords(0.5, -0.3); ax3.yaxis.set_label_coords(-0.30,0.5)

ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)

# =============================================================================
# consecutive
# =============================================================================
#####

conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}
hypothesis_list_df = [conc_data[c] for c in conc_labels]
labels_hypothesis_list_df = [conc_labels[c] for c in conc_labels]

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
gs3_inner_inner = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs3_inner[3], height_ratios=[1.2,0.9,0.9])#, wspace=0.4, hspace=0.4)
ax2 = plt.subplot(gs3_inner_inner[1:3]);

for i,df_consecutive in enumerate(hypothesis_list_df):
    label = labels_hypothesis_list_df[i]
    
    consecutive_cumulative_obj = consecutive_cumulative(df_consecutive)
    box_plot_consecutive_cumulative = consecutive_cumulative_obj.get_consecutive_trains_of_pulses()
    Norm = len(df_consecutive.index.get_level_values(0).unique())
    
    ax2.plot(np.arange(1,len(box_plot_consecutive_cumulative)+1),[i/Norm for i in box_plot_consecutive_cumulative], color=colors[i], linewidth=0.75, marker = "." , markersize=4, alpha=1,label=labels_hypothesis_list_df[i])
    ax2.set_yscale('log')
    ax2.set_xlim([0, 10.5]); ax2.set_ylim([0.01,10])
    ax2.set_xticks([1,4,7,10])
    ax2.set_ylabel('counts x cell',fontsize=8); ax2.set_xlabel('length of sequence \n of consecutive pulses ',fontsize=8)
    ax2.xaxis.set_label_coords(0.5, -0.09);
    ax2.yaxis.set_label_coords(-0.3,0.4);
ax2.tick_params(axis='x', labelsize=6,length=2); ax2.tick_params(axis='y', labelsize=6,length=2)
ax2.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)


# =============================================================================
# pulse density
# =============================================================================
#####

#gs3_inner_inner[0].remove()
gs1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=4, subplot_spec=gs_main[2:4], wspace=0.6, hspace=0.5,width_ratios=[1.5,0.2,0.6,0.9])


ylim = ylim
M = max([conc_data[key].max().FRAME for key in conc_data.keys()])
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
box=[];handles = [None]*2


for k,c in enumerate(cells_conc):
    df = conc_data[c];df.sort_index(inplace=True); color=colors[k]
    box.append((df['amp_peaks'].groupby(level='cell').count() / df['amp_peaks'].groupby(level='cell').size()).to_list()) 

#label = [conc_labels[c] for c in cells_conc]
label = [0,2.5,5,20]
X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
gs2_ = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs1[0:2,3], wspace=0.0, hspace=0.0,height_ratios=[1.9,0.1])

ax2 = plt.subplot(gs2_[0]);

for i in range(len(X)):
    xA = np.random.normal(0, 0.1, len(box[i])), 
    ax2.scatter(xA+X[i],box[i],  alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
bp = ax2.boxplot(box,whis=[5, 95],vert=True,patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp['boxes']):
     box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp['whiskers']):
    whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp['caps']):
    cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp['medians']):
    median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp['fliers']):
    flyer.set(markeredgecolor=colors[i])## change color and linewidth of the medians


ax2.tick_params(axis='x', labelsize=6,length=2); ax2.tick_params(axis='y', labelsize=6,length=2)
ax2.set_xticklabels(label,rotation = 0)
label_peaks = r'$\rm{\left(min^{-1}\right)}$'
ax2.set_ylabel('pulse density ' + label_peaks,fontsize=8)
ax2.set_xlabel('[FGF4] (ng/ml)',fontsize=8)
ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.01,0.55)


yticks = ax2.get_yticks(); yticks[-1] = 0.06 
minutes_ylim = yticks[-1]/(20/60)
ax2.set_yticklabels([0,minutes_ylim],fontsize=6)
ax2.set_ylim([-0.005,yticks[-1]]); ax2.set_yticks([0,yticks[-1]])


ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
#ax2.invert_yaxis()

# =============================================================================
# activity
# =============================================================================
plt.rc('axes.spines', top=False, bottom=True, left=False, right=False); 

gs0_inner_ = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc), ncols=1, subplot_spec=gs1[0,0:2], wspace=0.0, hspace=0.0)

conc_labels= {'an_0ng': '0ng/ml of [FGF]', 'an_2-5ng': '2.5ng/ml of [FGF]',
 'an_5ng': '5ng/ml of [FGF]','an_20ng': '20ng/ml of [FGF]'}

for k,c in enumerate(cells_conc):
    ax1 = plt.subplot(gs0_inner_[k]);
    df = conc_data[c]
    activity = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count() *  100
    activity = np.sort(activity)[::-1]
    silent = np.ones(len(activity)) * 100 - activity
    p1 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),silent,width=0.8,color='darkgray',alpha=0.5,linewidth=0.0)
    p2 = ax1.bar(np.arange(0,len(df.index.get_level_values(0).unique())),activity,bottom=silent,width=0.8,color=colors[k],alpha=0.8,linewidth=0.0)
    #if k == 0: plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 2.6), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)
    ax1.set_xlim([-1,df.index.get_level_values(0)[-1]+1]);ax1.set_ylim([0,100])
    ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
    
    if k != len(ipi_conc):
        silent_ax(ax1);
    else:
        ax1.set_yticks([0,100]); ax1.set_xticks([]);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2);


    #ax1.set_xlabel( conc_labels[c] + ' concentration' ,fontsize=10); 
    #ax1.set_xticks([0,df.index.get_level_values(0)[-1]])
ax1.set_ylabel('% of cell track' ,fontsize=8); 
ax1.yaxis.set_label_coords(-0.05,1.5)
ax1.set_xlabel('individual cells per condition' ,fontsize=8); 
ax1.xaxis.set_label_coords(0.5,-0.15)


# =============================================================================
# mean activity
# =============================================================================

ax1 = plt.subplot(gs1[0,2]);

activity = []; silent = []
activity_err = []; silent_err = []

for c in cells_conc:
    df = conc_data[c]
    activity_ = df['dt_peaks'].groupby(level='cell').sum() / df['FRAME'].groupby(level='cell').count()
    silent_ = np.ones(len(activity_)) - activity_
    activity.append(activity_.mean()*100)
    silent.append(silent_.mean()*100)
    silent_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))
    activity_err.append( silent_.std()*100/ np.sqrt(len(df.index.get_level_values(0).unique())))



p1 = ax1.barh(y = np.arange(len(cells_conc)),width = silent,xerr=silent_err,left = 0, color='darkgray',alpha=0.5,linewidth=0.0 )
p2 = ax1.barh(y = np.arange(len(cells_conc)),width = activity,xerr=activity_err,left = silent, color=colors,alpha=0.8,linewidth=0.0 )
plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(-0.1, 1.43), loc='upper left',frameon=False,fontsize=6,markerscale=0.2,handletextpad=0.1,labelspacing=0.1)



#p1 = ax1.bar(np.arange(len(cells_conc)),silent,yerr=silent_err,width=0.8,color='darkgray',alpha=0.5,linewidth=0.4)
#p2 = ax1.bar(np.arange(len(cells_conc)),activity,bottom=silent,yerr = activity_err,width=0.8,color=colors,alpha=0.8,linewidth=0.4)
#plt.legend((p1[0], p2[0]), ('silence', 'pulsing'),fancybox=False,bbox_to_anchor=(0.55, 1.4), loc='upper left',frameon=False,fontsize=6,markerscale=0.2)

labels = [0,2.5,5,20]
plt.xlabel('% of cell track' ,fontsize=8); 
plt.ylabel('[FGF4] (ng/ml)',fontsize=8); 
plt.xticks([0,50,100]);plt.xlim([0,100]);plt.ylim([-0.55,3.5])
plt.yticks(np.arange(len(cells_conc)), labels,fontsize=6);ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
ax1.xaxis.set_label_coords(0.5, -0.2);ax1.yaxis.set_label_coords(-0.24,0.5)
ax1.invert_yaxis()


plt.savefig(save_folder+'Figure3.pdf', format='pdf',transparent=True)

#%%

for c in data_conc:
    df = data_conc[c]
    

#%%
fig = plt.figure(constrained_layout=False, figsize=(8.27/4, 11.69/4))
gs_main = gridspec.GridSpec(nrows=3, ncols=1, figure=fig); 

ipi_conc = ['an_2-5ng','an_5ng','an_20ng']
col = 'dt_peaks'#ax1
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 


for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ix = len(conc_data[c].index.get_level_values(0).unique())
    colors =  sns.color_palette(sns.color_palette("Blues_d",ix))
    ax1 = plt.subplot(gs_main[k]); 
    for cell, data in df.groupby(level='cell'): 
        ax1.plot(data[col].dropna().values[0:-1],data[col].dropna().values[1:],'o',color = colors[cell],label=cell,markersize=2,alpha=0.5); 

    ax1.set_xlim([0,60]);ax1.set_ylim([0,60])
    ax1.axvspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)
    ax1.axhspan(0,2, color='darkgoldenrod',alpha = 0.4,linewidth=0.0)

    
    
    if k != len(ipi_conc)-1:
        set_scale(ax1,[0,60],[0,60]) ;ax1.set_xticklabels([None,None]);ax1.set_yticklabels([None,None]);
        
    ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
set_scale(ax1,[0,60],[0,60]); ax1.set_xlabel('i-th duration (min)',fontsize=8); ax1.set_ylabel('i+1-th duration (min)',fontsize=8,rotation = 90)
ax1.xaxis.set_label_coords(0.5, -0.05);ax1.yaxis.set_label_coords(-0.1,0.5)
ax1.set_xticklabels([int(x*20/60) for x in [0,60]],fontsize=6)
ax1.set_yticklabels([int(x*20/60) for x in [0,60]],fontsize=6)

ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)



plt.savefig(save_folder+'Figure_timecorr.pdf', format='pdf',transparent=True)

#%%
# =============================================================================
# boxplot consecutiveness
# =============================================================================

#si sumas box1 y box2, te puede dar uno (si hay pico) o cero (si no hay picos)
    # Filtramos la celulas que tienen cero picos en orden de poder hacer estadistica solo de picos

gs3_inner_inner = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs3_inner[3], height_ratios=[1.2,0.9,0.9])#, wspace=0.4, hspace=0.4)

plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 
handles = [None]*2

box1 = []
box2 = []

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    c_aux=[]
    i_aux=[]
    for cells, data in df.groupby(level='cell'):
        consecutive_ = is_consecutive_cell(data);c_aux.append(np.nan_to_num(sum(consecutive_)/data.amp_peaks.count(),nan=-1))
        isolated_ = is_isolated_cell(data);i_aux.append(np.nan_to_num(isolated_/data.amp_peaks.count(),nan=-1))
    box1.append(c_aux)
    box2.append(i_aux)

def delete_non_pulsatile_cells(box):
    new_box=[]
    for condition in box:
        new_box.append([i for i in condition if i>-0.5])
    return(new_box)


box = delete_non_pulsatile_cells(box1)
label = [2.5,5,20]
X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
ax2 = plt.subplot(gs3_inner_inner[1:3]);

bp = ax2.boxplot(box,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

for i,box_ in enumerate(bp['boxes']):
     box_.set( color=colors[i+1], linewidth=0.0,facecolor=colors[i+1],alpha = 0.1)# change outline color
for i,whisker in enumerate(bp['whiskers']):
    whisker.set(color=colors[i//2+1],linestyle = '-', linewidth=1,alpha=0.3)
for i,cap in enumerate(bp['caps']):
    cap.set(color=colors[i//2+1],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
for i,median in enumerate(bp['medians']):
    median.set(color=colors[i+1],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
for i,flyer in enumerate(bp['fliers']):
    flyer.set(markeredgecolor='black')## change color and linewidth of the medians


for i in range(len(X)):
    xA = np.random.normal(0, 0.1, len(box[i])), 
    ax2.scatter(xA+X[i],box[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

ax2.tick_params(axis='x', labelsize=6,length=2); ax2.tick_params(axis='y', labelsize=6,length=2)
ax2.set_xticklabels(label,rotation = 0)
#label_peaks = r'$\rm{\left(\frac{cons. peaks \, \times \, cell}{Total peaks \, \times \,cell}\right)}$'

ax2.set_ylabel('% consecutive pulses',fontsize=8)
    
ax2.set_xlabel('[FGF4] (ng/ml)',fontsize=8)
ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.01,0.5)
xticks = ax2.get_yticks(); ax2.set_ylim([-0.05,1.05]); ax2.set_yticks([0,1])
ax2.set_yticklabels([0,100],fontsize=6) #porque es porcentaje


ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)


    

