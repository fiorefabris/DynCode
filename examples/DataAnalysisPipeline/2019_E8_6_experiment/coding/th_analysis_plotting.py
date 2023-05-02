import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import SymLogNorm

import pickle
import seaborn as sns
import numpy as np

def set_scale(ax,xlim,ylim):
    ax.yaxis.set_major_locator(ticker.FixedLocator(ylim))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([]))
    #ax.yaxis.set_tick_params(grid_alpha = 0.5,labelleft=True)

    ax.xaxis.set_major_locator(ticker.FixedLocator(xlim))
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.tick_params(axis="y",direction="out", width=1,labelsize=8)
    ax.tick_params(axis="x",direction="out", width=1,labelsize=8)


#%% 

def get_heatmap(df_result,condition,save=False,save_folder=None):
    ''' Squared seaborn heatmap of the threshold analysis. 
    
    NOTE: The function is specific for the amp_th= np.linspace(0,4000,100)
        & pendiente_th =  np.linspace(0,1000,100) exploration.
    -----------------------------------------------
    INPUT:
        - df_result       : DataFrame with corresponding x,y and z axis
        - condition (str) : The condition of the expermient

    ------------------------------------------------
    OUTPUT DESCRIPTION:
        file name : 01_threshold_sq_heatmap_**condition**
        x axis    : slope threshold 
        y axis    : Amplitude Threshold
        z axis    : Number of peaks x frame (averaged x cell) [log scale]
    '''

    fig, axs = plt.subplots(1,1, sharex='row', sharey='row' ,figsize = (8.27, 11.69/2.5))
    ax=axs
    hm = sns.heatmap(df_result.fillna(-1),annot=False,annot_kws={'fontsize':5},
                cbar=True,ax = ax,vmin=0.001,vmax=0.06,cbar_kws={'label': r'Number of peaks x frame '+
                                            ' \n (averaged x cells)' ,'extend': 'both','ticks':[0.001, 0.01,0.03,0.05, 0.06]},cmap="RdYlBu_r",
                norm=SymLogNorm(vmin=0.001,vmax=0.06,linthresh=0.001,linscale=0),square=True)
    ax.set_ylabel('Amplitude threshold (a.u.)',fontsize=15); ax.set_xlabel('Slope threshold (a.u./frames)',fontsize=15);
    
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)
    
    ax.invert_yaxis()

    
    set_scale(ax,np.linspace(0,100,11), np.linspace(0,100,9))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter([str(x) for x in np.arange(0,1100,100)]))
    axs.yaxis.set_major_formatter(ticker.FixedFormatter([str(x) for x in np.arange(0,4500,500)]))


    axs.xaxis.set_label_coords(0.5, -0.20);axs.yaxis.set_label_coords(-0.20,0.5)
    fig.subplots_adjust(bottom=0.20, top=0.9, left=0.15, right=0.9,wspace=0.3, hspace=0.7)
    axs.tick_params(labelsize=12)
    axs.tick_params(axis='x',rotation=45)


    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.yaxis.label.set_size(15)
    cbar.ax.set_yticklabels([x for x in cbar.get_ticks()])
    
    max_value= round(max(df_result.fillna(-1).max()),4)
    ax.text(1.1, 1.05,'Max value: '+str(max_value), ha='center', va='center', transform=ax.transAxes, fontsize=10)
    fig.subplots_adjust(bottom=0.20, top=0.9, left=0.15, right=0.9,wspace=0.0, hspace=0.0)

    if save:
        plt.savefig(save_folder + '01_threshold_sq_heatmap_'+condition+'.pdf', format='pdf')
    else:
        plt.show()
#%%
        
def get_level_curves(heatmap_values,heatmap_condition,
                     lcurves_values,lcurves_condition,
                     save=False,save_folder=None,
                     plot_lcurves=False,save_lcurves = False):
    ''' Matplotlib heatmap of the threshold analysis. The specified 
        levels curves are plotted. 
    
    NOTE: The function is specific for the amp_th= np.linspace(0,4000,100)
        & pendiente_th =  np.linspace(0,1000,100) exploration.
    -----------------------------------------------
    INPUT:
        - heatmap_values          : DataFrame with corresponding x,y and z axis
        - heatmap_condition (str) : The condition of the expermient
        - lcurves_values          : DataFrame with corresponding x,y and z axis, 
        corresponding to the level curve condition
        - lcurves_condition (str) : The condition of the level curve
        
        ** [plot_lcurves = True] plots the individual level curves 
        ** [save_lcurves = True] saves a .pkl file with the specified level curves
    ------------------------------------------------
    OUTPUT DESCRIPTION:
        file name : 02_threshold_levelcurves_**condition**_**level_curves_condition**
        x axis    : slope threshold 
        y axis    : Amplitude Threshold
        z axis    : Number of peaks x frame (averaged x cell) [log scale]
        lines     : level curves
    ''' 
    x = np.array(heatmap_values.columns) #pendiente
    y = np.array(heatmap_values.index) #amplitú
    
    xx,yy = np.meshgrid(x,y)
    Z=heatmap_values.values
    
    fig, axs = plt.subplots(1,1, sharex='row', sharey='row' ,figsize = (8.27, 11.69/3))
    ax=axs
    plt.grid(0)
    levels= np.geomspace(1e-5,1e-1,num=50)
    levels=np.arange(1e-5,0.05+5e-5,5e-5)
    CSF = ax.contourf(xx, yy, Z,levels=levels,extend='both',cmap='RdYlBu_r',vmin=1e-5,
                      vmax=0.05,norm=SymLogNorm(vmin=1e-5,vmax=0.05,linthresh=1e-5,linscale=0))
    c_levels= [5e-3]
    #c_levels= [1e-5,1e-4,1e-3,5e-3,1e-2]
    c_levels= [1e-4]

    CS = ax.contour(xx, yy, lcurves_values,levels=c_levels,norm=SymLogNorm(vmin=1e-5,vmax=0.05,linthresh=1e-5), colors='green')
    
    ax.clabel(CS, inline=1, fontsize=10,fmt='%1.5f' )#
    
    ax.set_ylabel('Amplitude threshold (a.u.)',fontsize=15); 
    ax.set_xlabel('Slope threshold (a.u./frames)',fontsize=15);
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)
    
    axs.xaxis.set_label_coords(0.5, -0.19);
    axs.yaxis.set_label_coords(-0.19,0.5)
    fig.subplots_adjust(bottom=0.20, top=0.9, left=0.15, right=0.9,wspace=0.3, hspace=0.7)
    axs.tick_params(labelsize=12)
    axs.tick_params(axis='x',rotation=45)
    
    
    cbar_ticks=[1e-5,1e-4,1e-3,5e-3,1e-2,0.05]
    cbar = fig.colorbar(CSF, ax=ax,label= r'Number of peaks x frame'+' \n (averaged x cells)', 
                        ticks= cbar_ticks,spacing='proportional')
    cbar.add_lines(CS)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.yaxis.label.set_size(15)
    cbar.ax.set_yticklabels([x for x in cbar_ticks])
    
    ax.set_title('Colormap: ' + heatmap_condition + ' - Lines: ' + lcurves_condition,fontsize=15)
    fig.subplots_adjust(bottom=0.20, top=0.9, left=0.15, right=0.9,wspace=0.0, hspace=0.0)
        
    if save:
        plt.savefig(save_folder+ '02_threshold_levelcurves_'+heatmap_condition+'_'+lcurves_condition+'.pdf', format='pdf')
    plt.show()
    
    if plot_lcurves:
        plt.figure()
        for i,level in enumerate(range(len(CS.allsegs))):
            print(level)
            if CS.allsegs[level] == []:
                pass
            else:
                dat0= CS.allsegs[level][0]
                plt.plot(dat0[:,0],dat0[:,1],label=str(i+1))
                plt.grid(0)
                legend = plt.legend(title='Level curve '+lcurves_condition,fontsize=10)
                plt.setp(legend.get_title(),fontsize=12)
        plt.show()
        
    if save_lcurves:
        for i,level in enumerate(range(len(CS.allsegs))):
            if CS.allsegs[level] == []:
                pass
            else:
                dat0= CS.allsegs[level][0]
                filename = save_folder+lcurves_condition+'_level_'+str(level+1)+'.pkl'
                outfile = open(filename,'wb')
                pickle.dump(dat0,outfile)
                outfile.close()    