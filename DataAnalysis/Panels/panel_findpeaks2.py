import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import pandas as pd

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale

def plot_find_peaks_panel_2(conc,peaks_conc_dt,th_dt,vmax_dt,x_range_dt,dist_bet_peaks_conc_dt,num_bins,total_peaks_conc_dt,
                          order,cells_conc,conc_labels,filter_labels,traze_values, M, ylim,
                          save=False, **kwargs):

    plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); colors = sns.color_palette('colorblind');sns.set(context='talk', style='white')


    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
    gs_main = gridspec.GridSpec(nrows=1, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.9, hspace=0.0)
    gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs_main[0], wspace=0.1, hspace=0.1)
    # uno para cada filtro o detrending

    for n, traze_value in enumerate(traze_values):
        col = str(traze_value) + '_' + str(th_dt[n]) + 'PEAKS_O' + str(order);

        inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(conc), ncols=2, subplot_spec=gs0[0,n], wspace=0.0, hspace=0.0)


        for k,c in enumerate(conc):
            df = conc[c];df.sort_index(inplace=True);
            inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=5, ncols=1, subplot_spec=inner_gs0[k,0], wspace=0.3, hspace=0.0)
            ax = plt.subplot(inner_inner_gs0[0]);ax2= plt.subplot(inner_inner_gs0[1:4])

            ax.plot(df.loc[cells_conc[c]][traze_value], color=colors[k], linewidth=0.5);
            #ax.set_xlim([0, M]);ax.set_ylim(ylim[n + 1])
            ax.plot(df.loc[cells_conc[c]][col], 'o', color='k', markersize=0.8, alpha=1)
            silent_ax(ax)



            aux, aux2 = [], []
            for cell, data in df.groupby(level='cell'):
                aux.extend(np.arange(0, len(data)));aux2.extend([cell] * len(data))
            index = [np.array(aux2), np.array(aux)];index_tuples = list(zip(*index))
            index = pd.MultiIndex.from_tuples(index_tuples, names=['cell', 'frames'])
            df = df.set_index(index);df.sort_index(inplace=True)
            mask = df[col].fillna(0).unstack().isnull()
            cmap = plt.cm.get_cmap('cool');cmap.set_under(color='linen')
            im = sns.heatmap((df[col].fillna(0).unstack()), mask=mask, cmap=cmap, vmin=th_dt[n], vmax=vmax_dt[n]
                             , cbar_kws={'label': 'Peaks Amplitude'}, cbar=False, ax=ax2)
            ax2.set_ylabel('');ax2.set_xlabel('');ax2.set_xticks([]);ax2.set_yticks([])


        #mappable = im.get_children()[0]
        #cb_ax = ax.add_axes([0.47, 0.14, 0.2, 0.005])
        #cbar = ax.colorbar(mappable, cax=cb_ax, orientation='horizontal')
        #cbar.set_ticks([th_dt[n], vmax_dt[n]])
        #cbar.ax.tick_params(labelsize=10)

        # histogramas
        peaks_conc = peaks_conc_dt[n]; dist_bet_peaks_conc = dist_bet_peaks_conc_dt[n]
        for k,c in enumerate(conc):
            df = conc[c]; df.sort_index(inplace=True);peaks_df = peaks_conc[c]; peaks_df.sort_index(inplace=True)
            ax = plt.subplot(inner_gs0[k,1]); ax.set_xlim(list(x_range_dt[n]))
            ax.hist(dist_bet_peaks_conc[c].loc[:,'time_dist_bet_peaks_'+col], num_bins, density=True,range=x_range_dt[n],
                                       stacked= False,facecolor=colors[k],  alpha=1,label = 'Number of peaks = '+str(total_peaks_conc_dt[n][c]))
            ax.tick_params(axis='x', labelsize=15)
            ax.tick_params(axis='y', labelsize=15)
            ax.legend(loc="upper right", fontsize=8)
            silent_ax(ax)


    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(str(save_folder) + str(save_name) + '.pdf', format='pdf')

    plt.show()

