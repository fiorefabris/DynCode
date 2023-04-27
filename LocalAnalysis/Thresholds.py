def local_quantifiers_histograms(conc, description, cumulative=False, save=False, *argv, **kwargs):
    """Plot the amplitud histograms of the *argv columns of each df in conc"""

    for c in conc:
        df = conc[c]
        Data = df[traze_value + '_PE']
        plot_values = [plot_value for plot_value in argv]
        range_values = kwargs.get('range_values', None)




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





#%%

def set_threshold(Data, Thresholds):
    ''' Searchs for values over an Thresholds list of Data.
    Returns
    - Data_Peaks : a Pandas Series with the values and positions of the local maxima
    '''

    for c in conc:
        df = conc[c]
        for th in Thresholds:
            Data_Over_Th = Data.where(Data > Abs_Thr, 0)


    '''
    -Thresholds  : un array de los valores que quiera del Threshold
    -Traze value : La columna del Df que quiera analizar
    -orders = lista de ordenes que quiero ver

    RETURNS:
        - conc with peaks information
        - peaks_df: dictionary of df. Each df has the information of Number_Peaks, Frames, Dens_Peaks for each concentration.
    '''

    peaks_conc = {}

    for c in conc:
        df = conc[c]
        peaks_df = pd.DataFrame()

        for order in orders:
            for th in Thresholds:
                auxNP, auxF, auxDP, auxCell = [], [], [], []
                aux_peaks = []

                for cells, data in df.groupby(level='cell'):
                    Peaks, Number_Peaks, Frames, Dens_Peaks = Find_Peaks(data[traze_value], th, order)

                    aux_peaks.extend(Peaks)
                    auxNP.append(Number_Peaks)
                    auxF.append(Frames)
                    auxDP.append(Dens_Peaks)
                    auxCell.append(cells)

                df[str(traze_value) + '_' + str(th) + 'PEAKS_O' + str(order)] = np.array(aux_peaks)
                peaks_df[str(traze_value) + '_' + str(th) + 'NP_O' + str(order)] = auxNP
                peaks_df[str(traze_value) + '_' + str(th) + 'F_O' + str(order)] = auxF
                peaks_df[str(traze_value) + '_' + str(th) + 'DP_O' + str(order)] = auxDP

        peaks_df['Cell'] = auxCell
        peaks_df = peaks_df.set_index(['Cell'], drop=True)
        peaks_conc[c] = peaks_df