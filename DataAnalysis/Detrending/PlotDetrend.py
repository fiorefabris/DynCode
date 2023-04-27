import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from DataAnalysis.Preprocesing.PlotData import silent_ax,set_scale


def plot_detrend(conc, len_cells, M, plot_value, plot_value_filter, description, ylim, save=False,cols = 5,long=1, **kwargs):
    colors = sns.color_palette('muted')
    # No entiendo que graficar arriba de que


    Tot = np.sum(len_cells);
    Cols = cols;
    Rows = Tot // Cols;
    if Tot % Cols != 0:
        Rows = Rows + 1

    M = np.max(M)

    fig, axs = plt.subplots(Rows, Cols, sharex=False, sharey=False, figsize=(8.27, 11.69*long))
    axs = axs.ravel()
    zo=0;handles=[None]*len(conc)
    for k,C in enumerate(conc):
        color = colors[k];
        df = conc[C];
        df.sort_index(inplace=True);
        cells = df.index.get_level_values(0).unique()
        plt.rcdefaults()

        for z,cell in enumerate(cells):
            ax = axs[zo];zo=zo+1
            # ax.tick_params('y', colors=color)
            ax.set_ylim(ylim);
            ax.set_xlim([0, M])

            ax.grid(False);

            if cells[-1] >= cell:
                handles[k], = ax.plot(df.loc[cell][plot_value], color=color, linewidth=0.5,label=C)
                ax.plot(df.loc[cell][plot_value_filter], alpha=0.8, color='black', linewidth=0.5)
                ax.text(0.9, 0.1, 'Cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)

    for i,ax in enumerate(axs):
        if i == len(axs)-Cols:
            if (ylim[0] < 0):
                set_scale(ax, [0, M], [ylim[0], 0, ylim[1]]);
            else:
                set_scale(ax, [0, M], ylim);
            ax.set_ylabel('Intensity \n (a.u.)', fontsize=10);
            ax.set_xlabel('Frames', fontsize=10)
            ax.xaxis.set_label_coords(0.5, -0.1);
            ax.yaxis.set_label_coords(-0.15, 0.5)
        else:
            silent_ax(ax)

    description_label = kwargs.get('description_label', False);
    if description_label:
        plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()




    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    axs.flatten()[-1].legend(handles=handles, loc='upper right', bbox_to_anchor=(2.4, 1.5), ncol=1, markerscale=20,
                             title='Datasets', fontsize=10, title_fontsize=10, frameon=False)

    plt.figtext(0.1, 0.05, 'Figure: ' + description)

    plt.show();

    if save:
        save_folder = kwargs.get('save_folder', None);
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()

def heatmap_plot(conc, plot_value, ylim, description, save=False, **kwargs):
            ## ponerle los bordes. Poner las celulas del mismo ancho. ver tema escalas
            plt.rc('axes.spines', top=True, bottom=True, left=True, right=True)

            Tot = len(conc);
            Cols = 1;
            Rows = Tot // Cols;
            Rows += Tot % Cols
            cmap = sns.diverging_palette(220, 20, sep=3, as_cmap=True, s=99)
            fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27, 11.69))
            plt.grid(0)

            # caption = 'Time Series for '+str(C)
            # ax.text(0.0,-0.1,caption,transform=ax.transAxes,fontsize = 40)

            C = [c for c in conc]
            k = 0
            for i in range(Rows):
                if k < len(conc):
                    df = conc[C[k]]
                    im = sns.heatmap((df[plot_value]).unstack(), vmin=ylim[0], vmax=ylim[1], cmap=cmap,
                                     cbar_kws={'label': 'Intensity'},
                                     cbar=False, ax=axs[i])
                    im.set_ylabel('')
                    im.set_xlabel('')
                    axs[i].set_title('Concentration {}'.format(C[k]), fontsize=8)
                    im.set_xticks([])
                    im.set_yticks([])
                    k = k + 1
                    # im.tick_params(colors='steelblue')
                    # im.splines.set_color('k')
            fig.subplots_adjust(bottom=0.1, top=0.9, left=0.15, right=0.7, wspace=0.07, hspace=0.07)
            mappable = im.get_children()[0]
            cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
            cbar = fig.colorbar(mappable, cax=cb_ax)

            # set the colorbar ticks and tick labels
            # cbar.set_ticks(np.arange(-1, 5, 2))
            # cbar.set_ticklabels(['low', 'medium', 'high'])

            plt.figtext(0.1, 0.05, 'Figure: ' + description)

            plt.show();

            if save:
                save_folder = kwargs.get('save_folder', None);
                save_name = kwargs.get('save_name', None)
                fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
                plt.close()