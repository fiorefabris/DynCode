import sys

sys.path.append('/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_2_WTonly/coding/data_analysis_post_review_commons')
from main_oscillations import autocorr_normed
#%% Esto es para la autocorrelacion normalizada
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#

ylim = [-1,1]
xlim = [-1,330] #110 min
long=1;cols=5
save_folder =  '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/after_RC/figures_05_2022/'
save_name ='new_AC_single_cells_normed'

colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]
conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))


plt.rcdefaults()
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=len(conc_labels)*2, ncols=1, figure=fig, wspace=0.0, hspace=0.3)
plt.rcdefaults(); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);


handles = [None] * len(conc_data)

for k,C in enumerate(conc_labels):
    Rows = len_cells[k] // Cols;
    if len_cells[k] % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[k], wspace=0, hspace=0.0)

    label = C;
    color = colors[k];
    df = conc_data[C];
    df.sort_index(inplace=True);
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        data = df.query('cell ==' + str(cell))
        AC = autocorr_normed(data.detrended.values)
    
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
    
        if cells[-1] >= cell:
            handles[k],= ax.plot(AC, color=color, linewidth=0.4,label=label)
            ax.text(0.95, 0.2, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)

    
        if zo==((Rows-1)*Cols):
            #set_scale(ax, [], ylim);
            ax.set_ylabel('autocorrelation', fontsize=8);
            ax.set_xlabel('time lag (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -1);
            ax.yaxis.set_label_coords(-0.2, 1)
            ax.set_xticks([0,xlim[-1]])
            ax.set_yticks(ylim)
            ax.set_xticklabels([i/3 for i in [0,xlim[-1]]], fontsize=6)
            #ax.set_xticklabels([0, int(M * 20/60)], fontsize=6)
            ax.set_yticklabels([-1,1], fontsize=6)

        else:
            silent_ax(ax)

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')
plt.close()