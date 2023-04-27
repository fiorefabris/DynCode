


def cell_cycle_peaks(conc,peaks_conc,dist_bet_peaks_conc,total_peaks_conc,
                               th, traze_value, order, vmax, x_range, num_bins, save=False,**kwargs):

    col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)
    sns.set(context='talk', style='white')
    plt.rc('axes.spines',top=True,bottom=True,left=True,right=True)
    colors = sns.color_palette('colorblind')

    fig, axs = plt.subplots(len(conc),2, sharex='col', sharey='col',figsize = (8.27, 11.69))
    z = 0

    for i,c in enumerate(conc):
        df       = conc[c][col]; df.sort_index(inplace=True)
        aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
        for cell,data in df.groupby(level = 'cell'):
            cell_count += 1
            if len(data) > len(aux_peaks):
                aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
                cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(0,len(data)-len(cell_count_array)),'constant', constant_values=(0,0))

            else:
                aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
                cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(0,len(cell_count_array)-len(data)),'constant', constant_values=(0,0))

        cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'cell_number':cell_count_array,'beg_aligned':np.arange(0,len(aux_peaks))})
        cell_cycle.set_index('beg_aligned')











            aux.extend(np.arange(0,len(data))); aux2.extend([cell] * len(data) )
    #Esto claramente se pyede mejorar

        index = [np.array(aux2),np.array(aux)];  index_tuples = list(zip(*index))
        index = pd.MultiIndex.from_tuples(index_tuples, names=['cell', 'frames'])

        df = df.set_index(index); df.sort_index(inplace=True)

        mask = df[col].fillna(0).unstack().isnull()
        cmap=plt.cm.get_cmap('cool'); cmap.set_under(color='linen')

        im = sns.heatmap((df[col].fillna(0).unstack()),mask =mask,cmap = cmap,vmin = th, vmax = vmax
                         ,cbar_kws={'label': 'Peaks Amplitude'},cbar = False,ax = ax)
        ax.set_ylabel(''); ax.set_xlabel('')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title('Concentration {}'.format(c),fontsize = 8)

        z = z+1

    ylim = axs[z-1,0].get_ylim()
    description = 'Find peaks statistics for ' + str(traze_value) + ' th: '+ str(th) + '& time window : '+ str(order*2+1) + ' frames\n Histograms of time between peaks in minutes \n bins = ' + str(bins) +'& ylim = '+str(ylim)


    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
    mappable = im.get_children()[0]
    #cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    cb_ax = fig.add_axes([0.47, 0.14, 0.2, 0.005])
    cbar = fig.colorbar(mappable ,cax=cb_ax,orientation = 'horizontal')

    cbar.set_ticks([th, vmax])
    cbar.ax.tick_params(labelsize = 10)
    plt.figtext(0.01,0.05,'Figure: ' + description,fontsize=11)

    if save:
        save_folder= kwargs.get('save_folder', None) ; save_name= kwargs.get('save_name', None)
        fig.savefig(str(save_folder)+str(save_name)+'.pdf', format='pdf')
        print(str(save_folder)+str(save_name)+'.pdf')