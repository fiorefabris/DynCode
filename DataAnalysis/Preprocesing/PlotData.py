import matplotlib.ticker as ticker
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib

matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['axes.linewidth'] = 0.75

matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['xtick.major.width'] = 0.75
matplotlib.rcParams['xtick.minor.size'] = 1
matplotlib.rcParams['xtick.minor.width'] = 0.75
matplotlib.rcParams['ytick.major.width'] = 0.75
matplotlib.rcParams['ytick.major.width'] = 0.75
plt.rcdefaults()


#from DataAnalysis.Panels.panel_detrending import set_scale


def silent_ax(ax):
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

def set_yscale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.NullLocator())

def set_scale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([]))
    #ax.yaxis.set_tick_params(grid_alpha = 0.5,labelleft=True)

    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.tick_params(axis="y",direction="out", width=1,labelsize=8)
    ax.tick_params(axis="x",direction="out", width=1,labelsize=8)


def grid_plot(conc, len_cells, M, ylim,plot_value, description, save=False,long=1,cols=5, **kwargs):
    # En main version 5 está algunas alternativas de como intentar que el indice se muestre ensolo un grafico.

    colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
    colors = [colors[i] for i in [15, 5]]
    colors = sns.color_palette(sns.dark_palette("red", 30, reverse=True))
    colors = [colors[i] for i in [15, 5]]

    colors = [sns.color_palette()[5],sns.color_palette()[5]]
    esc_color = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=False))[15]
    colors = [esc_color,esc_color]

    #colors = sns.color_palette(sns.dark_palette("#3498db", 100, reverse=True))
   # colors = [colors[i] for i in [1, 30, 50, 80]]
    #colors = ['mediumslateblue',sns.cubehelix_palette(500, start=.30, rot=-.40)[250],sns.cubehelix_palette(500, start=.30, rot=-.40)[400]]


    #conc_labels = {'an_WT_N2chHep':'WT-N2B27 medium','an_WT_ESL':  'WT-ESL',"an_0ng":"KO-0ng FGF4",}
   # conc_labels = {'an_0ng': '0ng/ml of FGF4', 'an_2-5ng': '2.5ng/ml of FGF4',
   #                'an_5ng': '5ng/ml of FGF4', 'an_20ng': '20ng/ml of FGF4'}
    #conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}
    conc_labels = {'an_ESC_FAX': 'serum + LIF', 'an_ESC_FAX_PD03':'Meki'}
    #conc_labels = {'an_EpiSC_FAX': 'EPISCs in FAX', 'an_EpiSC_FAX_PD03':'EPISCs in FAX and MEKi'}



    #Tot = np.sum(len_cells);



    Cols = cols;
    M = np.max(M);  M = 330 ;     #
    M = 1200/1.75
    M = 230*3

    fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69))
    #para el only_Wt: #gs_row = gridspec.GridSpec(nrows=len(conc_labels)*2, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
    plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);
    gs_row = gridspec.GridSpec(nrows=len(conc_labels)+2, ncols=1, figure=fig, wspace=0.0, hspace=0.3)


    #axs = axs.ravel();z=0
    # last = len(axs) - 1
    #handles = [None] * len(conc)

    for k,C in enumerate(conc_labels):
        Rows = len_cells[k] // Cols;
        if len_cells[k] % Cols != 0:
            Rows = Rows + 1

        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

        label = C;
        color = colors[k];
        df = conc[C];
        df.sort_index(inplace=True);
        cells = df.index.get_level_values(0).unique()
        plt.rcdefaults()

        for zo,cell in enumerate(cells):
            #ax = axs[z];z = z+1
            ax = plt.subplot(inner_gs[zo])
            ax.set_ylim(ylim);
            ax.set_xlim([0, M])

            pointplot = kwargs.get('pointplot', False);
            if zo == 0:  # (w == len(cells_conc[c])-1):
                ax.text(0, 1.9, conc_labels[C], ha='left', va='top', transform=ax.transAxes, fontsize=8,
                        color=color)

            if pointplot:
                if cells[-1] >= cell:
                    # ax.plot(df.loc[cells[cell]][plot_value], color=color, linewidth=0.5)
                    ax.plot(df.loc[cell][plot_value],'o', color=color, markersize=0.5,label=label)
                    ax.text(0.9, 0.1, 'cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
            else:
                if cells[-1] >= cell:
                    #ax.plot(df.loc[cells[cell]][plot_value], color=color, linewidth=0.5)
                    ax.axhline(0,linestyle = ':',color='gray',linewidth=0.25)
                    ax.plot(df.loc[cell][plot_value], color=color, linewidth=0.4,label=label)
                    ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
            ax.grid(False);
            ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

            if zo==((Rows-1)*Cols):
                set_scale(ax, [0, M], ylim);
                ax.set_ylabel('KTR signal (a.u.)', fontsize=8);
                ax.set_xlabel('time (min) ', fontsize=8)
                ax.xaxis.set_label_coords(0.5, -0.25);
                ax.yaxis.set_label_coords(-0.2, 1)
                #ax.set_xticklabels([0, int(M * 1.75)], fontsize=6)
                ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
                ax.set_yticklabels([int(i/100) for i in ylim], fontsize=6)

            else:
                silent_ax(ax)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
    #axs.flatten()[-2].legend(handles = handles,loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3,markerscale = 20)
    #axs.flatten()[-1].legend(handles=handles, loc='upper right', bbox_to_anchor=(1.8, 1.5), ncol=1, markerscale=20,
    #                         title='Datasets',fontsize = 8,title_fontsize=10,frameon=False)
    description_label = kwargs.get('description_label', False);
    if False:
    #if description_label:
        plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None);
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
        plt.close()