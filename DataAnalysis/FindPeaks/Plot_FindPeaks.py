# th = 11.0
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from DataAnalysis.Detrending.MinPolDT import Find_Min
from DataAnalysis.Preprocesing.PlotData import silent_ax
import matplotlib.ticker as ticker


def peaks_grid_plot(conc, peaks_conc, len_cells, M, ylim, th, order, traze_value, save=False,cols = 5,long=1, **kwargs):

    col = str(traze_value) + '_' + str(th) + 'PEAKS_O' + str(order)

    colors = sns.color_palette('muted')
    k = 0; z = 0
    Tot = np.sum(len_cells); Cols = cols;  Rows = Tot // Cols;

    if Tot % Cols != 0:
        Rows = Rows + 1

    M = np.max(M)

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27, 11.69*long))
    axs = axs.ravel()
    plot_mins = kwargs.get('plot_mins', False)
    order_mins = kwargs.get('order_mins', 2)
    if plot_mins:
        for C in conc:
            df = conc[C]
            aux = []
            for cells, data in df.groupby(level='cell'):
                Mins = Find_Min(data[traze_value], order_mins)
                aux.extend(Mins)
            df['MIN_O' + str(order_mins)] = np.array(aux)


    for c in conc:
        df = conc[c]
        df.sort_index(inplace=True)
        peaks_df = peaks_conc[c]
        peaks_df.sort_index(inplace=True)
        color = colors[k];  k = k + 1
        cells = df.index.get_level_values(0).unique()
        plt.rcdefaults()

        for cell in cells:
            ax = axs[z];z = z + 1
            ax.set_ylim(ylim)
            ax.set_xlim([0, M])

            ax.grid(False)
            silent_ax(ax)


            if cells[-1] >= cell:
                ax.plot(df.loc[cell][traze_value], color=color, linewidth=0.5)
                ax.text(0.9, 0.1, 'Cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
                ax.plot(df.loc[cell ][col], 'o', color='r', markersize=0.3, alpha=1)
                if plot_mins:
                    col_mins = 'MIN_O' + str(order_mins)
                    ax.plot(df.loc[cell][col_mins], 'o', color='k', markersize=0.3, alpha=1)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    description = kwargs.get('description', None)
    plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')


def plot_time_between_peaks_st(conc,peaks_conc,dist_bet_peaks_conc,total_peaks_conc,
                               th, traze_value, order, vmax, x_range, num_bins, save=False,**kwargs):

    col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)
    sns.set(context='talk', style='white')
    plt.rc('axes.spines',top=True,bottom=True,left=True,right=True)
    colors = sns.color_palette('muted')

    fig, axs = plt.subplots(len(conc),2, sharex='col', sharey='col',figsize = (8.27, 11.69))
    z = 0

    for c in conc:
        df       = conc[c]; df.sort_index(inplace=True)
        peaks_df = peaks_conc[c]; peaks_df.sort_index(inplace=True)

        #histogramas
        ax =  axs[z,0]; silent_ax(axs[z,0])
        ax.set_xlim(list(x_range))
        n, bins, patches = ax.hist(dist_bet_peaks_conc[c].loc[:,'time_dist_bet_peaks_'+col], num_bins, density=True,range=x_range,
                                   stacked= False,facecolor=colors[z],  alpha=1)
        #,label = 'Number of peaks = '+str(total_peaks_conc[c])
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        #ax.set_title("\n".join(wrap( str(c)+'- Data = '+str(np.sum(n)) ,55)),fontsize=15)
        #ax.grid(False)
        #ax.legend(loc="upper right", fontsize=8)


        ax =  axs[z,1]
        aux,aux2 = [],[]
        for cell,data in df.groupby(level = 'cell'):
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

def plot_time_between_peaks_st2(conc,peaks_conc,dist_bet_peaks_conc,th,traze_value,order,x_range,num_bins,save=False,**kwargs):

    col = str(traze_value)+'_'+str(th)+'DP_O'+str(order)
    box = [peaks_conc[c][col].tolist() for c in peaks_conc]
    label = [c for c in conc]
    colors = sns.color_palette('colorblind')


    X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]


    fig, axs = plt.subplots(2,1, sharex=False, sharey= False,figsize = (8.27, 11.69))
    sns.despine()
    sns.set_style("ticks")
    plt.grid(0)

    for i in range(len(X)):
        axs[0].scatter(X[i], box[i], alpha=1,s = 2)

    bp = axs[0].boxplot(box,vert=True,patch_artist=False,showmeans=False,meanline=True )
    for box in bp['boxes']:
         box.set( color='k', linewidth=1)# change outline color

    for whisker in bp['whiskers']:
        whisker.set(color='k',linestyle = '-.', linewidth=1)

    for cap in bp['caps']:
        cap.set(color='k',linestyle = '-.', linewidth=1)## change color and linewidth of the caps

    for median in bp['medians']:
        median.set(color='k',linestyle = ':', linewidth=1)## change color and linewidth of the medians

    axs[0].tick_params(axis='x', labelsize=10); axs[0].tick_params(axis='y', labelsize=10)
    axs[0].set_xticklabels(label,rotation = 45)
    axs[0].text(0.3, 0.9,'Density of peaks - TH: '+str(th)+' - Time window: '+str(order*2+1)
            , ha='center', va='center', transform=axs[0].transAxes,fontsize = 10)

    yticks = axs[0].get_yticks(); axs[0].set_ylim([0,yticks[-1]]); axs[0].set_yticks([0,yticks[-1]])


    axs[1].set_xlim(list(x_range))
    col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)

    for i,c in enumerate(conc):
        aux = []
        for cell, data in conc[c].groupby(level='cell'):
            aux.append(len(data))
        norm = np.sum(aux)

        hist, bin_edges = np.histogram(dist_bet_peaks_conc[c].loc[:,'time_dist_bet_peaks_'+col], num_bins,range=x_range, density=False)
        bin_edges_center = 0.5*(bin_edges[1:]+bin_edges[:-1])
        axs[1].plot(bin_edges_center,np.divide(hist,norm),color=colors[i],linewidth=0.4,marker='.',markersize = 10,label=c)


    yticks = axs[1].get_yticks(); axs[1].set_ylim([0,yticks[-1]]); axs[1].set_yticks([0,yticks[-1]])
    axs[1].tick_params(axis='x', labelsize=10); axs[1].tick_params(axis='y', labelsize=10)
    axs[1].legend(loc="upper right", fontsize=8)

    if save:
        save_folder= kwargs.get('save_folder', None) ; save_name= kwargs.get('save_name', None)
        fig.savefig(str(save_folder)+str(save_name)+'.pdf', format='pdf')


def plot_dist_bet_peaks_hist_all(conc,dist_bet_peaks_conc,total_peaks_conc,th,traze_value,order,x_range,num_bins,cumulative = False,save=False,**kwargs):

    sns.set(context='talk', style='white')
    plt.rc('axes.spines',top=True,bottom=True,left=True,right=True)
    colors = sns.color_palette('muted')

    fig, axs = plt.subplots(1, sharex=False, sharey= False,figsize = (8.27, 11.69/2))
    sns.despine()
    sns.set_style("ticks")
    plt.grid(0)

    axs.set_xlim(list(x_range))
    col = str(traze_value) + '_' + str(th) + 'PEAKS_O' + str(order)

    for i, c in enumerate(conc):
        hist, bin_edges = np.histogram(dist_bet_peaks_conc[c].loc[:, 'time_dist_bet_peaks_' + col], num_bins,
                                       range=x_range, density=False)
        bin_edges_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        if cumulative == False:


            aux = []
            for cell, data in conc[c].groupby(level='cell'):
                aux.append(len(data))
            norm = np.sum(aux)
            axs.plot(bin_edges_center, np.divide(hist,norm), color=colors[i], linewidth=0.4, marker='.', markersize=10,label=c)
            axs.legend(loc="upper right", fontsize=8)

        else:
            aux = []
            for cell, data in conc[c].groupby(level='cell'):
                aux.append(len(data))
            norm = np.sum(aux)
            axs.plot(bin_edges_center, np.cumsum(np.divide(hist,norm)), color=colors[i], linewidth=0.4, marker='.', markersize=10,label=c)
            axs.legend(loc="upper left", fontsize=8)

    yticks = axs.get_yticks()
    axs.set_ylim([0, yticks[-1]])
    axs.set_yticks([0, yticks[-1]])
    axs.tick_params(axis='x', labelsize=10)
    axs.tick_params(axis='y', labelsize=10)


    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7, wspace=0.1, hspace=0.5)
    description = 'Histograms of time between peaks in minutes \n The "ocurrencies" are devided by the number of frames'
    #plt.figtext(0.1, 0.05, 'Figure: ' + description, fontsize=11)
    plt.show()
    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(str(save_folder) + str(save_name) + '.pdf', format='pdf')

def plot_time_between_peaks_hist(conc,peaks_conc,dist_bet_peaks_conc,total_peaks_conc, th, traze_value, order
                               ,x_range,num_bins,cumulative = False,save=False,**kwargs):

    col = str(traze_value)+'_'+str(th)+'PEAKS_O'+str(order)
    sns.set(context='talk', style='white')
    plt.rc('axes.spines',top=True,bottom=True,left=True,right=True)
    colors = sns.color_palette('muted')

    #fig, axs = plt.subplots(len(conc)-1, sharex=True, sharey=True,figsize = (8.27, 11.69))
    fig, axs = plt.subplots(len(conc), sharex=False, sharey=True, figsize=(8.27, 11.69))
    z = 0

    conc_ = [c for c in conc]
    #for c in conc_[1:]:
    for c in conc_:
        df       = conc[c]; df.sort_index(inplace=True)
        peaks_df = peaks_conc[c]; peaks_df.sort_index(inplace=True)

        #histogramas
        if len(conc) > 1:
            ax =  axs[z]; z = z+1
            if z != (len(conc)):
                ax.xaxis.set_major_locator(ticker.NullLocator())
                ax.xaxis.set_minor_locator(ticker.NullLocator())
            else:
                ax.set_xlabel('Time(minutes)')
        else:
            ax = axs
        ax.set_xlim(list(x_range))
        if cumulative == False:
            #ax.set_ylim([0,0.2])
            n, bins, patches = ax.hist(dist_bet_peaks_conc[c].loc[:,'time_dist_bet_peaks_'+col], num_bins, density=False,range=x_range,
                                       stacked= False,facecolor=colors[z],  alpha=1)
            #label = 'Number of peaks = ' + str(total_peaks_conc[c]
        else:
            ax.set_ylim([0,1])
            n, bins, patches = ax.hist(dist_bet_peaks_conc[c].loc[:,'time_dist_bet_peaks_'+col], num_bins, density=True,cumulative = True,range=x_range,
                                       stacked= False,facecolor=colors[z],  alpha=1,label = 'Number of peaks = '+str(total_peaks_conc[c]))

        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_title('Concentration {}'.format(c),fontsize = 8)
        #ax.legend(loc="upper right", fontsize=8)
    print(bins)
    plt.show()
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7, wspace=0.1, hspace=0.5)
    description = 'Histograms of time between peaks in minutes \n bins = ' + str(bins)
    plt.figtext(0.1, 0.05, 'Figure: ' + description,fontsize = 11)

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(str(save_folder) + str(save_name) + '.pdf', format='pdf')


