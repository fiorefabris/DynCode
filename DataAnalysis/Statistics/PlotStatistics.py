import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


from DataAnalysis.Preprocesing.PlotData import silent_ax

def amplitud_histograms(conc, description,cumulative=False, save=False, *argv, **kwargs):

    """Plot the amplitud histograms of the *argv columns of each df in conc"""

    plot_values = [plot_value for plot_value in argv]
    range_values = kwargs.get('range_values', None)

    color = sns.color_palette('colorblind')
    Cols = len(plot_values)
    Rows = len(conc)


    fig, axs = plt.subplots(Rows, Cols, sharex=False, sharey=True, figsize=(8.27, 11.69))
    j = 0
    for plot_value in plot_values:
        i = 0
        axs[i, j].set_title(plot_values[j])
        for C in conc:
            ax = axs[i, j]
            num_bins = 50

            if range_values[j][0] < 0 and range_values[j][1] > 0:
                ax.axvline(0, color='black', linestyle=':', linewidth=1)

            n, bins, patches = ax.hist(conc[C][plot_value], num_bins, range=range_values[j], density=1,
                                       facecolor=color[i], alpha=0.8, cumulative=cumulative)
            ax.set_xlim(range_values[j])
            ax.grid(False)
            silent_ax(ax)
            i = i + 1
        j = j + 1

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()
# Falta que return las scales

def plot_statistics(doc,percentiles_dict,description,save=False,**kwargs):

    plot_values = doc.value.unique()
    statistical = np.array(doc.columns[3:])

    Cols = len(plot_values)
    Rows = len(doc.columns) - 3 + 1

    fig, axs = plt.subplots(Rows, Cols, sharex=False, sharey='row', figsize=(8.27, 11.69))

    j = 0
    for plot_value in plot_values:
        df = doc.loc[doc.value == plot_value]
        df.set_index(['Concentration'], inplace=True)
        i = 0
        axs[i, j].set_title(plot_value)
        for st in statistical:
            ax = axs[i, j]
            silent_ax(ax)
            if j == 0:
                ax.set_ylabel(st, rotation=90, size='large')
            aux_st = [aux_df[st].values for value, aux_df in df.groupby(level='Concentration')]
            sns.boxplot(data=aux_st, orient='v', ax=ax, palette='colorblind', linewidth=0.7, width=0.5, fliersize=1.5)
            # sns.swarmplot(data=aux_st,orient='v' , ax=ax,color='black')
            i = i + 1
        j = j + 1
    j = 0
    color = sns.color_palette('colorblind')
    cmap = LinearSegmentedColormap.from_list("my_colormap", color)
    axs[i, j].set_ylabel('Values of percentiles', rotation=90, size='large')
    for P in percentiles_dict:
        ax = axs[i, j]
        ax.grid(False)
        percentiles_dict[P].T.plot(grid=False, cmap=cmap, ax=ax, legend=None, linewidth=0.5, marker='o',
                                 markersize=2)
        j = j + 1
        silent_ax(ax)
        if P == 'MEAN_INTENSITY_BP': ax.legend(loc='upper right',bbox_to_anchor=(0, -0.5, 1.02, 0.5),ncol=3,fontsize =5)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()

