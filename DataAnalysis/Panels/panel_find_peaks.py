import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.gridspec as gridspec
from DataAnalysis.FindPeaks.Threshold_analysis import find_peaks_dens_th

from DataAnalysis.Preprocesing.PlotData import silent_ax
from DataAnalysis.Panels.panel_detrending import set_scale


def plot_find_peaks_panel(conc,peaks_conc_dt,test_conditions,traze_values_dt,th_values_dt,Thresholds,
                          order,cells_conc,conc_labels,filter_labels,traze_values, M, ylim,
                          save=False, **kwargs):

    M = np.max(M);plt.rcdefaults();matplotlib.rc('axes',edgecolor='silver')

    fig = plt.figure(constrained_layout=False,figsize=(8.27, 11.69))

    gs_main= gridspec.GridSpec(nrows = 1, ncols =1,figure = fig);gs_main.update( left=0.1, right=0.9,bottom = 0.5,top = 0.9, hspace=0.0)
    gs0 = gridspec.GridSpecFromSubplotSpec(nrows = 1, ncols =2,subplot_spec=gs_main[0], wspace=0.3,hspace=0.0) #uno para cada filtro o detrending

    for n,detrend in enumerate(traze_values):
        th_values = th_values_dt[n]
        inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs0[n], wspace=0.3, hspace=0.0)
        peaks_conc = peaks_conc_dt[n]

        peaks_dens_th, colors = find_peaks_dens_th(peaks_conc, traze_values_dt[n], Thresholds, order,test_conditions)
        peaks_dens_th = [j / np.max(j) for j in peaks_dens_th]  # normalizacion

        for j, peaks_dens_th_condition in enumerate(peaks_dens_th):
            ax = plt.subplot(inner_gs0[1])
            set_scale(ax,th_values,[0,0.1,1]); ax.set_xlabel('Threshold', fontsize=10);ax.set_ylabel('Density of peaks', fontsize=10)
            ax.yaxis.set_label_coords(-0.01, 0.5)
            ax.plot(Thresholds, peaks_dens_th_condition, linestyle=":", marker='o', markersize=6, color=colors[j],
                    alpha=0.5, label=conc_labels[j])
            ax.grid(0);ax.legend(bbox_to_anchor=(0.5,1), loc="upper left",fontsize=8)#ax.legend(loc="upper right", fontsize=8)


        ### Plot cell panel
        inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(th_values) + 1 , ncols=1, subplot_spec=inner_gs0[0],
                                                           wspace=0.0, hspace=0.1)#Uno para cada threshold que quiero mostrar y el histograma

        for k,th in enumerate(th_values):
            inner_inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(test_conditions), ncols=1, subplot_spec=inner_inner_gs0[k], wspace=0.0, hspace=0.0)
            col = str(detrend) + '_' + str(th) + 'PEAKS_O' + str(order)

            for z,c in enumerate(test_conditions):
                color = colors[z]; df = conc[c]; ax = plt.subplot(inner_inner_inner_gs0[z])


                ax.plot(df.loc[cells_conc[c]][traze_values_dt[n]], color=color, linewidth=0.5);
                ax.set_xlim([0, M]);ax.set_ylim(ylim[n + 1])

                #ax.axhline(th, color='black', linestyle=':', linewidth=1)

                ax.plot(df.loc[cells_conc[c]][col], 'o', color='k', markersize=0.8, alpha=1)

                if k ==0 and z==0:     ax.set_title(filter_labels[n] ,horizontalalignment='left',verticalalignment='bottom')
                if (k == len(th_values) - 1 and z == (len(test_conditions) - 1)):
                    set_scale(ax, [0, M], ylim[n + 1])
                    ax.set_ylabel('Intensity \n (a.u.)', fontsize=10);ax.set_xlabel('Frames', fontsize=10)
                    ax.xaxis.set_label_coords(0.5, -0.01);ax.yaxis.set_label_coords(-0.01, 0.5)
                else:
                    silent_ax(ax)

        ### histograma
        inner_inner_inner_gs0 = gridspec.GridSpecFromSubplotSpec(nrows=len(test_conditions), ncols=1,
                                                                 subplot_spec=inner_inner_gs0[len(th_values)], wspace=0.0,
                                                                 hspace=0.0)
        for z, c in enumerate(test_conditions):
            color = colors[z];ax = plt.subplot(inner_inner_inner_gs0[z]);
            # ax.axvline(0, color='black', linestyle=':', linewidth=1)

            ax.hist(conc[c][traze_values[n]], 50, range=tuple(ylim[n + 1]), density=1,facecolor=color, alpha=0.8)
            ax.set_xlim(tuple(ylim[n+1]));ax.set_ylim([0,0.004])
            if z ==0: silent_ax(ax)
            else: set_scale(ax,th_values_dt[n],[0,0.0040]); ax.set_xlabel('Intensity \n (a.u.)',fontsize=10)





    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    #plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.tight_layout()

    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None);
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()
