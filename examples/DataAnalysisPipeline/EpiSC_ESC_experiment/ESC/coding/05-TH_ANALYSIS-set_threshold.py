# =============================================================================
# boxplot for finding the threshold value
# =============================================================================

data_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/data_explore_th/'
save_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/explore_th_figures/'
dataset_name='ESC'

TH = [ 89.37424044, 508.47457627] # (89, 508) #el primer try 

# c_levels = [1e-3,5e-3,1e-2]
# level 1 (th_analysis) : (57, 406) + ('57.25643730066835', '406.77966101694915')
# level 1B (th_analysis) : (67, 438) +  ('67.79661016949152', '438.5496707041748')
# level 2 (th_analysis) : (47, 305) + ('47.984553993628516', '305.08474576271186')
# level 3 (th_analysis) : (30, 305) + ('30.52109447851462', '305.08474576271186')

target_conc = 'WT'
lcurves_condition = 'CN'

level_description= ['Threshold chriteria 1e-3','Threshold chriteria 5e-3','Threshold chriteria 1e-2']

c_levels = [1e-3,5e-3,1e-2]
c_levels= [1e-3]

fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.69))
gs_main = gridspec.GridSpec(nrows=3, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.2, top=0.90, hspace=0.0)


for k,l in enumerate(['level_1B']):
    ax2 = plt.subplot(gs_main[k]);
    filename = data_folder+'df_results_curve_'+target_conc+'_'+l+'.pkl'
    infile = open(filename,'rb')
    results = pickle.load(infile)
    infile.close()

    filename = data_folder+lcurves_condition+'_'+l+'.pkl'
    infile = open(filename,'rb')
    level = pickle.load(infile)
    infile.close()

    label = [tuple([int(i[0]),int(i[1])]) for i in level]
    box= [results[i].amp_peaks.groupby(level='cell').count()/results[i].FRAME.groupby(level=0).count() for i in results]
    X = [np.ones(len(box[i]))*(i+1) for i in range(0,len(box))]
        
    ax2.plot(11, 0.025, marker="v", markersize=4, alpha=1, color="red",linewidth = 0)
    p1 = ax2.bar(np.arange(len(box)),height = [np.mean(i) for i in box], yerr = [np.std(i)/np.sqrt(len(i)) for i in box],color=colors[1],alpha=0.5,linewidth=0.0,width=0.8)
    
    
    ax2.tick_params(axis='x', labelsize=6,direction='out', pad=1,length=2); ax2.tick_params(axis='y', labelsize=6,direction='out', pad=1,length=2)

 
    ax2.set_ylim([0,5e-2]); #ax2.set_yticks([0,0.01,0.03,0.05,0.1])
    ax2.set_xlim([-1,len(label)+1]) #[-1,len(label)+1]
    set_scale(ax2,np.arange(0,len(label),8),[0,5e-2])

    label2 = [(np.round(l[0]/100/20*60,2),np.round(l[1]/100,2)) for l in label]
    label3=[label2[i] for i in np.arange(0,len(label),8)]
    ax2.set_xticklabels(label3,rotation = 0)
    ax2.set_yticklabels([0,5e-2/20*60],rotation = 0)



    ax2.set_ylabel( 'averaged \n pulse density \n ' + r'(min$^{-1})$',fontsize=8,rotation = 90)
    ax2.set_xlabel(r'$(v_{th}$ (a.u./min) , $A_{th}$ (a.u.))',fontsize=8);
    ax2.xaxis.set_label_coords(0.5, -0.35);
    ax2.yaxis.set_label_coords(-0.03,0.5)
    ax2.invert_xaxis()
ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
fig.savefig(save_folder + l + '.pdf', format='pdf')

#%%

def get_linear_detrend(conc):
    '''
    Calculates the linear detrend of each df.MEAN_INTENSITY 
    in conc.
    
    -------------------------------------------------
    INPUT:
        - conc (dict): dictionary of dataframes
    
    -------------------------------------------------
    OUTPUT:
        - conc with two aditional columns: 'trend' and 'detrended'
    '''
    for key in conc.keys():
        df = conc[key]
        trend_aux = []
        detrended_aux= []
        
        for cells, data in df.groupby(level='cell'):
            
        # fit linear model
            X = data.index.get_level_values(1)
            X = np.reshape(X, (len(X), 1))
            y = data.MEAN_INTENSITY.values
            model = LinearRegression()
            model.fit(X, y)
        # calculate trend
            trend = model.predict(X)
            detrended = [y[i]-trend[i] for i in range(0, len(data.MEAN_INTENSITY))]
            trend_aux.extend(list(trend))
            detrended_aux.extend(list(detrended))
        conc[key]['trend'] = trend_aux
        conc[key]['detrended'] = detrended_aux

#%% 

filename= '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/ESC.pkl'
conc = load_file(filename)

main_folder = '/home/fabris/Documents/Trabajo/Doctorado/Dyncode/DO2020/EpiSC_ESC_Experiment/ESC/data/data_explore_th/'
conc_data = {}
key =[ 'an_ESC_FAX_PD03', 'an_ESC_FAX']
keywords = ['CN','WT']

TH = ('67.79661016949152', '438.5496707041748')
slope_th = 67
amp_th = 438
level = 'level_1B'

# before: TH = ('89.37424043658498', '508.4745762711865') # (89, 508)
# c_levels = [1e-3,5e-3,1e-2] level1B = 5e-4
# level 1B (th_analysis) : (67, 438) +  ('67.79661016949152', '438.5496707041748')
# level 1 (th_analysis) : (57, 406) + ('57.25643730066835', '406.77966101694915')
# level 2 (th_analysis) : (47, 305) + ('47.984553993628516', '305.08474576271186')
# level 3 (th_analysis) : (30, 305) + ('30.52109447851462', '305.08474576271186')



for j,k in enumerate(keywords):
    filename = main_folder + 'df_results_curve_'+k+'_'+level+'.pkl'
    aux = load_file(filename)
    conc_data[key[j]] = aux[TH]
get_linear_detrend(conc_data)





#%%
figure_name = 'pulses_'+level
    
len_cells=[54,52]
ylim = [58000,65000]

colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
colors = [colors[i] for i in [15, 5]]

z = 0
Tot = np.sum(len_cells); Cols = 4;  Rows = Tot // Cols;
if Tot % Cols != 0: Rows = Rows + 1
M = max([conc_data[key].max().FRAME for key in conc_data.keys()])


fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*2, 11.69*2))
axs = axs.ravel()

for k,c in enumerate([ 'an_ESC_FAX_PD03', 'an_ESC_FAX']):
    df = conc_data[c]
    df.sort_index(inplace=True)
    #peaks_df = peaks_conc[c]
    #    peaks_df.sort_index(inplace=True)
    color = colors[k];  
    cells = df.index.get_level_values(0).unique()
    plt.rcdefaults()

        
    for cell in cells:
        ax = axs[z];z = z + 1
        ax.set_ylim(ylim)
        ax.set_xlim([0-5, M+5])

        ax.grid(False)
        silent_ax(ax)


        if cells[-1] >= cell:
            ax.plot(df.loc[cell]['sm_MEAN_INTENSITY'].values, color=color, linewidth=0.7,alpha=0.7)
            ax.text(0.9, 0.1, 'Cell ' + str(cell), ha='center', va='center', transform=ax.transAxes, fontsize=5)
            ax.plot(df.loc[cell]['max_'].values, 'o', color='r', markersize=0.3, alpha=1)
            ax.plot(df.loc[cell]['min_'].values, 'o', color='k', markersize=0.3, alpha=1)
#           ax.plot(df.loc[cell]['sm_MEAN_INTENSITY_PEAKS_O1'].values, 'o', color='r', markersize=0.1, alpha=1)
#           #col_mins = 'MIN_O1' + str(order_mins)
#           ax.plot(df.loc[cell]['MIN_O1'].values, 'o', color='k', markersize=0.1, alpha=1)
                
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)

fig.savefig(save_folder + figure_name + '.pdf', format='pdf')
#%%