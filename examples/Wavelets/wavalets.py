import matplotlib
plt.rcdefaults()

matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0.0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)


def silent_ax(ax):
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
#%%
plt.close('all')
dt = 1/3
dj = 1/32

long=1;Cols=3

conc_labels = {'an_WT_ESL': 'serum + LIF', 'an_WT_ESL_PD03':'Meki'}

len_cells = []
for c in conc_data:
    len_cells.append(len(conc_data[c].index.get_level_values(0).unique()))

handles = [None] * len(conc_data)

#%%
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5) #11.69
gs_row = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, wspace=0.0, hspace=0.1,height_ratios=[1.98, 0.02])
plt.rc('axes.spines', top=False, bottom=True, left=True, right=False);

c = 'an_WT_ESL_PD03'
len_cells_ = len_cells[0]

for k,C in enumerate([c]):
    Rows = len_cells_ // Cols;
    if len_cells_ % Cols != 0:
        Rows = Rows + 1
        inner_gs = gridspec.GridSpecFromSubplotSpec(nrows=Rows, ncols=Cols, subplot_spec=gs_row[0], wspace=0.05, hspace=0.1)

    label = C;
    df = conc_data[C];
    df.sort_index(inplace=True);
    cells = df.index.get_level_values(0).unique()


    for zo,cell in enumerate(cells):
        ax = plt.subplot(inner_gs[zo])
        data = df.query('cell ==' + str(cell))
        
        #time_series = data.MEAN_INTENSITY.values - data.FITM_O1.values
        time_series = data.MEAN_INTENSITY.values - np.mean(data.MEAN_INTENSITY.values)
        #time_series = data.sm_MEAN_INTENSITY.values
        
        wa = WaveletAnalysis(time_series, dt=dt,dj=dj,wavelet=wavelets.Morlet(w0=6))
        power = wa.wavelet_power # wavelet power spectrum
        scales = wa.scales # scales 
        t = wa.time # associated time vector
        rx = wa.reconstruction() # reconstruction of the original data

        T, S = np.meshgrid(t, scales)
        levels = np.linspace(100*(100**2), 3000*(100**2),100)
        img=ax.contourf(T, S, power, 100,vmin = 100*(100**2),vmax = 3000*(100**2),levels = levels,extend='both')
        a,b=wa.coi
        ax.fill_between(a,b,30,color = 'gray',alpha = 0.5,linewidth = 0)

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.5)
        ax.tick_params(width=0.5)
    
        if cells[-1] >= cell:
            #handles[k],= ax.plot(FFT_omega,FFT_plot, color=color, linewidth=0.7,label=label)
            ax.text(0.95, 0.9, str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.grid(False);
        
        ax.tick_params(labelsize=6, direction='out', pad=1, length=2)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.5)
        ax.tick_params(width=0.5)
        ax.set_xlim([0,110])
        ax.set_ylim([0,20])

    
        if zo==((Rows-1)*Cols):
            ax.set_ylabel('period (min)', fontsize=8);
            ax.set_xlabel('time (min)', fontsize=8)
            ax.xaxis.set_label_coords(0.5, -0.5);
            ax.yaxis.set_label_coords(-0.1, 1)
            ax.set_xticks([0,110])
            ax.set_yticks([0,20])
#            ax.set_xticklabels([0,1,2,3], fontsize=6)
#            ax.set_yticklabels([0,ylim[-1]/1e8], fontsize=6)
#
        else:
            silent_ax(ax)
    
    cbar_ax = plt.subplot(gs_row[1])
    cbar_ticks = np.linspace(100*(100**2), 3000*(100**2),5)
    cbar = plt.colorbar(img, cax=cbar_ax,ticks=cbar_ticks, orientation="horizontal")
    cbar.ax.set_xticklabels(np.round(i/(100**2),2) for i in cbar_ticks)
    cbar.set_label('power (a.u.)', rotation=0)


#fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig('/home/fabris/Documents/Trabajo/Doctorado/Dyncode/Dyncode_Figures/wavelets/NEW/PDO3.pdf', format='pdf')
plt.close()
#%%

