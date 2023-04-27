import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


from DataAnalysis.Preprocesing.PlotData import silent_ax,set_scale

def quality_tracking_test(conc, len_cells, M,ylim,var_th,var_ylim, description, save=False,long=1,cols=5, **kwargs):
    colors = sns.color_palette('muted');    plt.rcdefaults()

    k = 0;
    z = 0;

    Tot = np.sum(len_cells);
    Cols = cols;
    Rows = Tot // Cols;
    if Tot % Cols != 0:
        Rows = Rows + 1

    M = np.max(M)

    fig, axs = plt.subplots(Rows, Cols, sharex=False, sharey=False, figsize=(8.27, 11.69*long))
    plt.rcdefaults()
    axs = axs.ravel()
    handles = [None] * len(conc)

    for C in conc:
        color = colors[k];
        k = k + 1;
        df = conc[C];
        df.sort_index(inplace=True);
        cells = df.index.get_level_values(0).unique()

        for cell in cells:
            ax = axs[z];
            z = z + 1;
            ax.grid(False)
            ax1 = ax.twinx();
            ax1.grid(False)

            ax1.axhline(y = var_th,linestyle = ':',alpha = 0.7, color = 'gray',linewidth=0.5)
            ax.set_ylim(ylim);
            ax.set_xlim([0, M]);
            ax1.set_ylim(var_ylim)

            if cells[-1] >= cell:
                handles[k-1], = ax.plot(df.loc[cells[cell]]['MEAN_INTENSITY'], color=color, linewidth=0.5,label=C)
                ax.text(0.9, 0.1, 'Cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
                X = df.loc[cells[cell]]['STANDARD_DEVIATION'].index.get_level_values(0)
                ax1.fill_between(X, 0, df.loc[cells[cell]]['STANDARD_DEVIATION'], linewidth=0.0, alpha=0.3,
                                 color='gray')

            if z == (len(axs) - Cols+1):
                if (ylim[0] < 0):
                    set_scale(ax, [0, M], [ylim[0], 0, ylim[1]]);
                else:
                    set_scale(ax, [0, M], ylim);
                ax.set_ylabel('Intensity \n (a.u.)', fontsize=10);
                ax.set_xlabel('Frames', fontsize=10)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.15, 0.5)
                ax1.yaxis.set_major_locator(ticker.NullLocator())
                ax1.yaxis.set_minor_locator(ticker.NullLocator())
            else:
                silent_ax(ax);silent_ax(ax1)
                #if z == len(axs) - 1:
                #    set_scale(ax1, [None, None], var_ylim);
                #    ax1.set_ylabel('Tracking \n Variance', fontsize=10);
                #else:
                #

    for ax in axs[z:]: silent_ax(ax)
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    axs.flatten()[-1].legend(handles=handles, loc='upper right', bbox_to_anchor=(2.4, 1.5), ncol=1, markerscale=20,
                             title='Datasets',fontsize = 10,title_fontsize=10,frameon=False)
    description_label = kwargs.get('description_label', False);
    if description_label:
        plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None);
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()