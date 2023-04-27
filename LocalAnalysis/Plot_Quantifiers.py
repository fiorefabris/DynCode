import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

from DataAnalysis.Preprocesing.PlotData import silent_ax
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
plt.rcdefaults()


def silent_ax(ax):
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

def set_scale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([]))
    #ax.yaxis.set_tick_params(grid_alpha = 0.5,labelleft=True)

    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.tick_params(axis="y",direction="out", width=1,labelsize=8)
    ax.tick_params(axis="x",direction="out", width=1,labelsize=8)

def plot_local_value(conc,traze_value,plot_value, len_cells, M, ylim ,plot_value_ylim, save=False,Cols=5, **kwargs):

    if 'PE' in plot_value:            title = 'Shannon Permutation Entropy '
    elif 'F' in plot_value:           title = 'Fisher Permutation Entropy'
    #elif 'C' in plot_value:           title = 'CJS '
    elif 'LOCAL_VAR' in plot_value:   title = 'Local Variance'
    elif '_LOCAL_MEAN' in plot_value: title = 'Local Mean'
    elif '_CV2' in plot_value:        title = 'Local Cv2'


    plot_value = traze_value + '_'+plot_value
    colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
    colors = [colors[i] for i in [15, 5]]

    plt.rcdefaults(); plt.rc('axes.spines', top=True, bottom=True, left=True, right=True);
    fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69))
    gs_row = gridspec.GridSpec(nrows=len(conc)+2, ncols=1, figure=fig, wspace=0.0, hspace=0.3)

    fig.suptitle(title,va='baseline', fontsize=12)
    for k,C in enumerate(conc):

        color = colors[k]; df = conc[C]; df.sort_index(inplace=True)
        cells = df.index.get_level_values(0).unique()

        Rows = len_cells[k] // Cols;
        if len_cells[k] % Cols != 0:
            Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

        for z,cell in enumerate(cells):
            #ax = axs[z]; ax.grid(False)
            ax = plt.subplot(inner_gs[z])
            ax1 = ax.twinx(); #ax1.grid(False)
            ax.set_ylim(ylim);
            ax.set_xlim([0, M])
            ax1.set_ylim(plot_value_ylim)

            if cells[-1] >= cell:
                ax.plot(df.loc[cells[cell]][traze_value], color=color, linewidth=0.5)
                ax.text(0.9, 0.1, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
                X = df.loc[cells[cell]][plot_value].index.get_level_values(0)
                ax1.fill_between(X, 0, df.loc[cells[cell]][plot_value].fillna(0), linewidth=0.0, alpha=0.3,
                                 color='gray')
                #print(ax1.get_ylim())

            if z==((Rows-1)*Cols):
                set_scale(ax, [0, M], ylim);
                ax.set_ylabel('KTR signal (a.u.)', fontsize=8);
                ax.set_xlabel('time (min) ', fontsize=8)
                ax.xaxis.set_label_coords(0.5, -0.25);
                ax.yaxis.set_label_coords(-0.2, 1)
                ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
                ax.set_yticklabels([int(i/100) for i in ylim], fontsize=6)
            else:
                silent_ax(ax)

            if z == 0:
                #ax1.yaxis.set_label_position("right")
                #ax1.yaxis.tick_right()
                #ax.yaxis.set_ticks_position('both')
                set_scale(ax1, [0, M], plot_value_ylim);
                ax1.set_ylabel('local quant', fontsize=8);
                #ax1.set_xlabel('time (min) ', fontsize=8)
                #ax1.xaxis.set_label_coords(0.5, -0.25);
                ax1.yaxis.set_label_coords(-0.2, 1)
                #ax1.set_xticklabels([0, int(M * 20/60)], fontsize=6)
                ax1.set_yticklabels(plot_value_ylim, fontsize=6,rotation=90);
            else:
                silent_ax(ax1)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)

    if save:
        save_folder_name = kwargs.get('save_folder_name', None)
        fig.savefig(save_folder_name+ '.pdf', format='pdf')
        plt.close()


def plot_local_value_all(conc,traze_value,plot_value_list, len_cells, M, ylim ,plot_value_ylim_list, save=False,Cols=5, **kwargs):

    colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
    colors = [colors[i] for i in [15, 5]]
    colors2 = sns.color_palette("hls", len(plot_value_list))

    plt.rcdefaults(); plt.rc('axes.spines', top=True, bottom=True, left=True, right=True);
    fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69))
    gs_row = gridspec.GridSpec(nrows=len(conc), ncols=1, figure=fig, wspace=0.0, hspace=0.3)

    for k,C in enumerate(conc):

        color = colors[k]; df = conc[C]; df.sort_index(inplace=True)
        cells = df.index.get_level_values(0).unique()

        Rows = len_cells[k] // Cols;
        if len_cells[k] % Cols != 0:
            Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

        for z,cell in enumerate(cells):
            ax = plt.subplot(inner_gs[z])
            ax.set_ylim(ylim);
            ax.set_xlim([0, M])
            if cells[-1] >= cell:
                ax.plot(df.loc[cells[cell]][traze_value], color=color, linewidth=0.5)
                ax.text(0.9, 0.1, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)

            for j,plot_value in enumerate(plot_value_list):
                plot_value = traze_value + plot_value
                ax1 = ax.twinx();
                ax1.set_ylim(plot_value_ylim_list[j])
                if cells[-1] >= cell:
                    X = df.loc[cells[cell]][plot_value].index.get_level_values(0)
                    ax1.fill_between(X, 0, df.loc[cells[cell]][plot_value].fillna(0), linewidth=0.1, alpha=0.3,
                                 color=colors2[j],label = plot_value)
                    silent_ax(ax1)
                    #ax1.legend()
            if z==((Rows-1)*Cols):
                set_scale(ax, [0, M], ylim);
                ax.set_ylabel('KTR signal (a.u.)', fontsize=8);
                ax.set_xlabel('time (min) ', fontsize=8)
                ax.xaxis.set_label_coords(0.5, -0.25);
                ax.yaxis.set_label_coords(-0.2, 1)
                ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
                ax.set_yticklabels([int(i/100) for i in ylim], fontsize=6)
            else:
                silent_ax(ax)

           # if z == 0:
                #ax1.yaxis.set_label_position("right")
                #ax1.yaxis.tick_right()
                #ax.yaxis.set_ticks_position('both')
               # set_scale(ax1, [0, M], plot_value_ylim);
               # ax1.set_ylabel('local quant', fontsize=8);
                #ax1.set_xlabel('time (min) ', fontsize=8)
                #ax1.xaxis.set_label_coords(0.5, -0.25);
                #ax1.yaxis.set_label_coords(-0.2, 1)
                #ax1.set_xticklabels([0, int(M * 20/60)], fontsize=6)
               # ax1.set_yticklabels(plot_value_ylim, fontsize=6,rotation=90);
           # else:


    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)

    if save:
        save_folder_name = kwargs.get('save_folder_name', None)
        fig.savefig(save_folder_name+ '.pdf', format='pdf')
        plt.close()
