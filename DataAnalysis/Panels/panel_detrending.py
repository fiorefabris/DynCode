import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

import matplotlib.gridspec as gridspec

from DataAnalysis.Preprocesing.PlotData import silent_ax



def set_scale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([]))
    #ax.yaxis.set_tick_params(grid_alpha = 0.5,labelleft=True)

    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.tick_params(axis="y",direction="out", width=1,labelsize=8)
    ax.tick_params(axis="x",direction="out", width=1,labelsize=8)


def plot_detrend_panel(conc,cells_conc,conc_labels,plot_value_filter,traze_values, M, ylim,plot_value, save=False, **kwargs):
    colors = sns.color_palette('colorblind'); M = np.max(M); plt.rcdefaults()


    fig = plt.figure(constrained_layout=False,figsize=(8.27, 11.69))

    #One gridspec for each filter + non fliter
    gs0 = gridspec.GridSpec(nrows = len(conc), ncols =1,figure = fig);gs0.update( left=0.1, right=0.35, hspace=0.05)
    gs1_main = gridspec.GridSpec(nrows = len(conc), ncols = 1,figure = fig);gs1_main.update( left=0.40, right=0.65, hspace=0.05)
    gs2_main = gridspec.GridSpec(nrows = len(conc), ncols = 1,figure = fig);gs2_main.update( left=0.70, right=0.95, hspace=0.05)



    for k,C in enumerate(conc):
        color = colors[k]; df = conc[C]

        gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs1_main[k], wspace=0.0, hspace=0.0)
        gs2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs2_main[k], wspace=0.0, hspace=0.0)


        inner_gs1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs1[0], wspace=0.0, hspace=0.0)
        inner_gs2 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs2[0], wspace=0.0, hspace=0.0)


        ax0 = plt.subplot(gs0[k]);ax1 = plt.subplot(inner_gs1[0]); ax1b = plt.subplot(inner_gs1[1])
        ax1c = plt.subplot(gs1[1]); ax2c = plt.subplot(gs2[1])
        ax2 = plt.subplot(inner_gs2[0]); ax2b = plt.subplot(inner_gs2[1])



        ax0.text(0.99, 0.01, conc_labels[C], ha='right', va='bottom', transform=ax0.transAxes, fontsize=10,color = color)
        ax0.plot(df.loc[cells_conc[C]][plot_value], color=color, linewidth=0.5); ax0.set_ylim(ylim[0]); ax0.set_xlim([0, M])

        ax1.plot(df.loc[cells_conc[C]][plot_value], color=color, linewidth=0.5);ax1.set_xlim([0, M]);ax1.set_ylim(ylim[0])
        ax1.plot(df.loc[cells_conc[C]][plot_value_filter[0]], alpha=0.8, color='black', linewidth=0.5);ax1b.set_xlim([0, M]);ax1b.set_ylim(ylim[1])
        ax1b.plot(df.loc[cells_conc[C]][traze_values[0]], color=color, linewidth=0.5)

        ax2.plot(df.loc[cells_conc[C]][plot_value], color=color, linewidth=0.5);ax2.set_xlim([0, M]);ax2.set_ylim(ylim[0])
        ax2.plot(df.loc[cells_conc[C]][plot_value_filter[1]], alpha=0.8, color='black', linewidth=0.5)
        ax2b.plot(df.loc[cells_conc[C]][traze_values[1]], color=color, linewidth=0.5);ax2b.set_xlim([0, M]);ax2b.set_ylim(ylim[2])

        ax1c.hist(df[traze_values[0]], 50, density=1, facecolor=color, alpha=0.8); ax1c.set_xlim(ylim[1]);ax1c.set_ylim([0,0.0025])
        ax2c.hist(df[traze_values[1]], 50, density=1, facecolor=color, alpha=0.8); ax2c.set_xlim(ylim[2]);ax2c.set_ylim([0,0.0025])


        if k ==0:
            ax0.set_title('Raw data',fontsize=10) #horizontalalignment
            ax1.set_title('Minima polynomial detrending',fontsize=10)
            ax2.set_title('Butter filter detrending',fontsize=10)

        if k==len(conc)-1:
            set_scale(ax0,[0, M],ylim[0]); ax0.set_ylabel('Intensity \n (a.u.)',fontsize=10); ax0.set_xlabel('Frames',fontsize=10)
            ax0.xaxis.set_label_coords(0.5, -0.01);ax0.yaxis.set_label_coords(-0.01,0.5)
            silent_ax(ax1)
            set_scale(ax1b,[0, M],ylim[1]);ax1b.set_xlabel('Frames',fontsize=10);ax1b.xaxis.set_label_coords(0.5, -0.01)
            silent_ax(ax2)
            set_scale(ax2b,[0, M],ylim[2]);ax2b.set_xlabel('Frames',fontsize=10);ax2b.xaxis.set_label_coords(0.5, -0.01)
            ax1c.yaxis.tick_right();set_scale(ax1c,[0,ylim[1][-1]],[0,0.0025]);ax1c.set_xlabel('Intensity \n (a.u.)',fontsize=10)
            ax2c.yaxis.tick_right();set_scale(ax2c,[0,ylim[2][-1]],[0,0.0025]);ax2c.set_xlabel('Intensity \n (a.u.)',fontsize=10)
            ax2c.set_ylabel('Ocurrences',fontsize=10)


        else:
            AX = [ax0,ax1,ax1b,ax2,ax2b]
            for ax in AX:
                silent_ax(ax)
                ax.set_xlim([0, M])
                silent_ax(ax1c);silent_ax(ax2c)


    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    #plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()


    if save:
        save_folder = kwargs.get('save_folder', None);
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()