
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

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
# Figure 01
#%%

sns.despine();sns.set(context='paper', style='ticks');plt.grid(0);

colors_ = sns.cubehelix_palette(500, start=.30, rot=-.40)[250]

#plt.rc('axes.spines', top=False, bottom=False, left=False, right=False); 
cells_conc = {'an_WT_ESL': [7,8,10,29,44],'an_WT_N2chHep':[23,24,8,1,0]}
conc_labels = {'an_0ng': 'KO-0ng of FGF4','an_WT_ESL':
                  'WT-ESL','an_WT_N2chHep':'WT-N2B27 medium' }

    
conc_aux = TimeSerieDC.TimeSerieDC(data_folder,dataset_name,save_folder,clean_last_points = True)
conc_order = ['an_WT_ESL','an_WT_N2chHep','an_0ng']
conc_aux.def_conc_order(conc_order)
conc_aux.get_values();    
    
ylim = [25000,64000]
M = 1200/1.75

fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69/3))
gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.95, wspace=0.0,hspace=0.5)

# =============================================================================
# 
# =============================================================================


gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[0], wspace=0.3, hspace=0.0)


inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_ESL']), ncols=1, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_ESL']), ncols=1, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)
inner_gs2 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_ESL']), ncols=1, subplot_spec=gs0[2], wspace=0.0, hspace=0.0)

conc = conc_aux.data
traze_value = conc_aux.traze_value
#th = 2500; order = 2

df = conc['an_WT_ESL'];df.sort_index(inplace=True);
for k,c in enumerate(cells_conc['an_WT_ESL']):
    ax = plt.subplot(inner_gs0[k]); color=colors_
    time = len(df.loc[c][traze_value].values)
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.set_xlim([0, 450]);
    ax.set_ylim(ylim)
    #ax.text(0.8, 0.1, conc_labels['an_WT_ESL'], ha='center', va='center', transform=ax.transAxes, fontsize=5)


    if k==len(cells_conc['an_WT_ESL'])-1:
        set_scale_bars(ax, x_bar=(0,34.3), y_bar=[25000,50000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
        ax.set_ylabel('KTR signal \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
        ax.xaxis.set_label_coords(0.1, -0.15);
        ax.yaxis.set_label_coords(-0.05,0.5);
        #set_scale(ax,[0, 450],ylim); ax.set_ylabel('Intensity (a.u.)',fontsize=8); ax.set_xlabel('Time (min)',fontsize=8)
        #ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
        #ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
        #ax.set_yticklabels([i/1000 for i in ylim], fontsize=6)
    else:
        set_scale_bars(ax)



conc = filtmin_conc.data
peaks_conc = filtmin_conc.peaks
traze_value =filtmin_conc.traze_value
th = 2500; order = 2
#ylim = filtmin_conc.ylim;

M = filtmin_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_ESL'];df.sort_index(inplace=True);
for k,c in enumerate(cells_conc['an_WT_ESL']):
    ax = plt.subplot(inner_gs1[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)

    ax.set_xlim([0, 450]);
    ax.set_ylim( [-5000, 25000])

    if k==len(cells_conc['an_WT_ESL'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-5000,20000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('KTR signal \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
#        set_scale(ax,[0, 450],[-5000, 25000]); 
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-5000, 25000]], fontsize=6)

    else:
        set_scale_bars(ax)



conc = buterfilt_conc.data
peaks_conc = buterfilt_conc.peaks
traze_value =buterfilt_conc.traze_value
th = 2400; order = 2
#ylim = buterfilt_conc.ylim
M = buterfilt_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_ESL'];df.sort_index(inplace=True)
for k,c in enumerate(cells_conc['an_WT_ESL']):
    ax = plt.subplot(inner_gs2[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)
    ax.set_xlim([0, 450]);
    ax.set_ylim( [-15000, 15000])

    if k==len(cells_conc['an_WT_ESL'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-15000,10000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('KTR signal \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
        
#        set_scale(ax,[0, 450],[-15000, 15000]); #ax.set_ylabel('Intensity  (a.u.)',fontsize=10); ax.set_xlabel('Frames',fontsize=10)
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-15000, 15000]], fontsize=6)
    else:
        set_scale_bars(ax)
################################
        ################################
        ##################################
        ##########################################
    
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[1], wspace=0.3, hspace=0.0)

inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)
inner_gs2 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[2], wspace=0.0, hspace=0.0)


conc = conc_aux.data
traze_value = conc_aux.traze_value
#th = 2500; order = 2


conc = filtmin_conc.data
peaks_conc = filtmin_conc.peaks
traze_value =filtmin_conc.traze_value
th = 2500; order = 2
#ylim = filtmin_conc.ylim;

M = filtmin_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_ESL'];df.sort_index(inplace=True);
for k,c in enumerate(cells_conc['an_WT_ESL']):
    ax = plt.subplot(inner_gs1[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    #ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    #ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)

    ax.set_xlim([0, 450]);
    ax.set_ylim( [-5000, 25000])

    if k==len(cells_conc['an_WT_ESL'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-5000,20000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('KTR signal \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
#        set_scale(ax,[0, 450],[-5000, 25000]); 
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-5000, 25000]], fontsize=6)

    else:
        set_scale_bars(ax)



conc = buterfilt_conc.data
peaks_conc = buterfilt_conc.peaks
traze_value =buterfilt_conc.traze_value
th = 2400; order = 2
#ylim = buterfilt_conc.ylim
M = buterfilt_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_ESL'];df.sort_index(inplace=True)
for k,c in enumerate(cells_conc['an_WT_ESL']):
    ax = plt.subplot(inner_gs2[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    #ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    #ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)
    ax.set_xlim([0, 450]);
    ax.set_ylim( [-15000, 15000])

    if k==len(cells_conc['an_WT_ESL'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-15000,10000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('KTR signal \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
        
#        set_scale(ax,[0, 450],[-15000, 15000]); #ax.set_ylabel('Intensity  (a.u.)',fontsize=10); ax.set_xlabel('Frames',fontsize=10)
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-15000, 15000]], fontsize=6)
    else:
        set_scale_bars(ax)
        
fig.savefig('/home/fabris/Documents/Dyncode/low_res_datasets/2018_E1_28/Sup_Figs/Figure1.pdf', format='pdf')


#%%
# =============================================================================
# el otro wt
# =============================================================================
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[1], wspace=0., hspace=0.0)

inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[0], wspace=0.0, hspace=0.0)
inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[1], wspace=0.0, hspace=0.0)
inner_gs2 = gridspec.GridSpecFromSubplotSpec(nrows=len(cells_conc['an_WT_N2chHep']), ncols=1, subplot_spec=gs0[2], wspace=0.0, hspace=0.0)

conc = conc_aux.data
traze_value = conc_aux.traze_value
#th = 2500; order = 2
ylim = [25000,64000]
M = 1200/1.75



df = conc['an_WT_N2chHep'];df.sort_index(inplace=True);
for k,c in enumerate(cells_conc['an_WT_N2chHep']):
    ax = plt.subplot(inner_gs0[k]); colors_='mediumslateblue'
    time = len(df.loc[c][traze_value].values)
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.set_xlim([0, 450]);
    ax.set_ylim(ylim)
    #ax.text(0.8, 0.1, conc_labels['an_WT_ESL'], ha='center', va='center', transform=ax.transAxes, fontsize=5)


    if k==len(cells_conc['an_WT_ESL'])-1:
        set_scale_bars(ax, x_bar=(0,34.3), y_bar=[25000,50000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
        ax.set_ylabel('intensity \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
        ax.xaxis.set_label_coords(0.1, -0.15);
        ax.yaxis.set_label_coords(-0.05,0.5);
        #set_scale(ax,[0, 450],ylim); ax.set_ylabel('Intensity (a.u.)',fontsize=8); ax.set_xlabel('Time (min)',fontsize=8)
        #ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
        #ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
        #ax.set_yticklabels([i/1000 for i in ylim], fontsize=6)
    else:
        set_scale_bars(ax)



conc = filtmin_conc.data
peaks_conc = filtmin_conc.peaks
traze_value =filtmin_conc.traze_value
th = 2500; order = 2
ylim = filtmin_conc.ylim
M = filtmin_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_N2chHep'];df.sort_index(inplace=True);
for k,c in enumerate(cells_conc['an_WT_N2chHep']):
    ax = plt.subplot(inner_gs1[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)

    ax.set_xlim([0, 450]);
    ax.set_ylim( [-5000, 25000])

    if k==len(cells_conc['an_WT_N2chHep'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-5000,20000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('intensity \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
#        set_scale(ax,[0, 450],[-5000, 25000]); 
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-5000, 25000]], fontsize=6)

    else:
        set_scale_bars(ax)



conc = buterfilt_conc.data
peaks_conc = buterfilt_conc.peaks
traze_value =buterfilt_conc.traze_value
th = 2400; order = 2
ylim = buterfilt_conc.ylim
M = buterfilt_conc.traze_max_values
col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

df = conc['an_WT_N2chHep'];df.sort_index(inplace=True)
for k,c in enumerate(cells_conc['an_WT_N2chHep']):
    ax = plt.subplot(inner_gs2[k]); color=colors_
    
    time = len(df.loc[c][traze_value].values)
    ax.axhline(y = th, linewidth=0.5, color='black',alpha = 1,linestyle=':')
    ax.plot(np.arange(time),df.loc[c][traze_value].values, color=colors_, linewidth=1);
    ax.plot(np.arange(time),df.loc[c][col], 'o', color='k', markersize=1, alpha=1)
    ax.set_xlim([0, 450]);
    ax.set_ylim( [-15000, 15000])

    if k==len(cells_conc['an_WT_N2chHep'])-1:
         set_scale_bars(ax, x_bar=(0,34.3), y_bar=[-15000,10000], xunits='', yunits=r'$\times 10^2$', x_scale=1.75, y_scale=1000, round_x = True,round_y=True)
         ax.set_ylabel('intensity \n 250 a.u.',fontsize=8); ax.set_xlabel('60 min',fontsize=8)
         ax.xaxis.set_label_coords(0.1, -0.15);
         ax.yaxis.set_label_coords(-0.05,0.5);
        
#        set_scale(ax,[0, 450],[-15000, 15000]); #ax.set_ylabel('Intensity  (a.u.)',fontsize=10); ax.set_xlabel('Frames',fontsize=10)
#        ax.xaxis.set_label_coords(0.5, -0.1);ax.yaxis.set_label_coords(-0.15,0.5)
#        ax.set_xticklabels([0, int(450 * 1.75)], fontsize=6)
#        ax.set_yticklabels([i/1000 for i in  [-15000, 15000]], fontsize=6)
    else:
        set_scale_bars(ax)



fig.savefig('/home/fabris/Documents/Dyncode/low_res_datasets/2018_E1_28/Sup_Figs/Figure1.pdf', format='pdf')

