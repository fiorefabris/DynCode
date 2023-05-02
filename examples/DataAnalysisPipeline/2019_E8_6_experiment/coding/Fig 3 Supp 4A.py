#%%

save_folder = '/home/fabris/Documents/Dyncode/DO2020/analysis_2019_E8_6/figures/'

import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)

conc_labels = {'an_0ng':'0ng/ml ','an_2-5ng': '2.5ng/ml ',
               'an_5ng': '5ng/ml','an_20ng': '20ng/ml'}
ipi_conc = ['an_2-5ng','an_5ng','an_20ng']

fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); gs_main.update(left=0.2, right=0.9, bottom=0.2, top=0.95, hspace=0.0)
gs0 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=len(ipi_conc), subplot_spec=gs_main[0], wspace=0.4, hspace=0.6)

col1 = 'amp_peaks'
colors =  sns.color_palette(sns.dark_palette("#3498db",100,reverse=True))
colors = [colors[i] for i in [1,30,50,80]]

for k,c in enumerate(ipi_conc):
    df = conc_data[c]
    ax2 = plt.subplot(gs0[k]);
    bins2 = ax2.hist(df[col1],bins=15,range =(amp_th,12000), density=1, facecolor=colors[k+1],  alpha=0.8,linewidth=0); 
    ax2.axvspan(0,amp_th, color='darkgoldenrod',alpha = 0.1,linewidth=0.0)

    ax2.set_xlim([0,12000]);ax2.set_ylim([0,0.0004]);
    set_scale(ax2,np.arange(0,12000,2000),[0,0.0004]) ;#ax2.set_xticklabels([None,None,None]);ax2.set_yticklabels([None,None]);
    ax2.set_xticklabels([int(x) for x in np.arange(0,12000,2000)/100],fontsize=6)

    n, bins = np.histogram(df[col1],bins=15,range =(amp_th,12000), density=1)

    mode_amp =  (bins2[1][np.argmax(bins2[0])] + bins2[1][np.argmax(bins2[0])+1])/2; mode_amp = 'mode: '+str(np.round(mode_amp/100,1))+' \n'
    ax2.text(1, 0.90, mode_amp, ha='right', va='center', transform=ax2.transAxes, fontsize=6)
    ax2.text(1, 1.05,r'$Q$: '+str(np.round(df[col1].quantile(0.25)/100,2))+' ; '+str(np.round(df[col1].quantile(0.50)/100,2))+' ; '
                                  +str(np.round(df[col1].quantile(0.75)/100,2)) , ha='right', va='center', transform=ax2.transAxes, fontsize=6)
    ax2.xaxis.set_label_coords(0.5, -0.1);ax2.yaxis.set_label_coords(-0.01,0.5)


    ax2.text(0.05, 1.05,  conc_labels[c], ha='left', va='center', transform=ax2.transAxes, fontsize=6,color = colors[k+1])
    if k == 0: ax2.set_xlabel('amplitude (a.u.)',fontsize=8); ax2.set_ylabel('frequency',fontsize=8,rotation = 90)
    ax2.xaxis.set_label_coords(0.5, -0.15);ax2.yaxis.set_label_coords(-0.05,0.5)
#ax2.set_xticklabels([0,amp_th/100,7000/100])
    ax2.set_yticklabels([0,0.04])
    ax2.tick_params(labelsize=6,direction='out', pad=1)


plt.savefig(save_folder + 'Fig. 3 Supp. 4.pdf', format='pdf') 

