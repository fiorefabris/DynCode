#### PARA EL 28
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


#%%%
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(0,len(data)-len(cell_count_array)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(0,len(cell_count_array)-len(data)),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_beg':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_beg')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_beg", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks');  axs[i,0].set_xlabel('Time (aligned at the beg)') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(0,len(data)-len(cell_count_array)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(0,len(cell_count_array)-len(data)),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_beg':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_beg')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_beg", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1])
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_BEG.pdf', format='pdf')

plt.show()

    


#%%

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(0,len(data)-len(cell_count_array)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(0,len(cell_count_array)-len(data)),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_beg':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_beg')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks/ cells');  axs[i,0].set_xlabel('Time (aligned at the beg)') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(0,len(data)-len(cell_count_array)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(0,len(cell_count_array)-len(data)),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_beg':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_beg')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1],alpha=1)
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_23_cell_cycle_density_BEG.pdf', format='pdf')

plt.show()

 
#%% con alineacion al final

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(len(data)-len(aux_peaks),0),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(len(data)-len(cell_count_array),0),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(len(aux_peaks)-len(data),0),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(len(cell_count_array)-len(data),0),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_end')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks');  axs[i,0].set_xlabel('Time (aligned at the end)') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(len(data)-len(aux_peaks),0),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(len(data)-len(cell_count_array),0),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(len(aux_peaks)-len(data),0),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(len(cell_count_array)-len(data),0),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_end')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1])
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_END.pdf', format='pdf')

plt.show()

#%%   
col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(len(data)-len(aux_peaks),0),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(len(data)-len(cell_count_array),0),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(len(aux_peaks)-len(data),0),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(len(cell_count_array)-len(data),0),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_end')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks/cells');  axs[i,0].set_xlabel('Time (aligned at the end)') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(len(data)-len(aux_peaks),0),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(len(data)-len(cell_count_array),0),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(len(aux_peaks)-len(data),0),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(len(cell_count_array)-len(data),0),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_end')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1])
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_density_END.pdf', format='pdf')

plt.show()

#%% con alineacion intermedia

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(cell_count_array)
            aux_peaks= np.array(data.notna()*1) + np.pad(aux_peaks ,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(cell_count_array)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_end')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks');  axs[i,0].set_xlabel('Time ') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(cell_count_array)
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(cell_count_array)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_end')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="total_peaks",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1])
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_MID.pdf', format='pdf')

plt.show()

#%%   
col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(cell_count_array)
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(cell_count_array)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    cell_cycle = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle.set_index('Time_aligned_end')
    cell_cycle['peaks_per_cell']=cell_cycle['total_peaks']/cell_cycle['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle,legend=False,ax = axs[i,0])
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')

axs[i,0].set_ylabel('Total number of peaks/cells');  axs[i,0].set_xlabel('Time ') 
axs[0,0].set_title('Filtmin detrend')
   
col='MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks = np.array([]);cell_count = 0;cell_count_array = np.array([])
    for cell,data in df.groupby(level = 'cell'):
        cell_count+=1
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(cell_count_array)
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = np.ones(len(data)) + np.pad(cell_count_array,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(cell_count_array)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            cell_count_array = cell_count_array + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    cell_cycle_bf = pd.DataFrame({'total_peaks':aux_peaks,'number_of_cells':cell_count_array,'Time_aligned_end':np.arange(0,len(aux_peaks))})
    cell_cycle_bf.set_index('Time_aligned_end')
    cell_cycle_bf['peaks_per_cell']=cell_cycle_bf['total_peaks']/cell_cycle_bf['number_of_cells']
    
    sns.scatterplot(x="Time_aligned_end", y="peaks_per_cell",hue = 'number_of_cells', 
                    sizes=(1, 100),palette="cool",data=cell_cycle_bf,legend=False,ax = axs[i,1])
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')
axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_density_MID.pdf', format='pdf')

plt.show()



#%% ahora hist


save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
for i,c in enumerate(filtmin_conc.conc_order):

    df       = filtmin_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []
    ax =  axs[i,0]
    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
        
    cell_cycle_df = pd.DataFrame({'normed_time':normed_time_aux,'peaks':aux_peaks})
    cell_cycle_df.set_index('normed_time',inplace=True)
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    n, bins, patches = ax.hist(plot_value, bins = np.arange(0,1,0.05), cumulative = False,density=True,stacked= False,facecolor=colors[i],  alpha=1)

axs[i,0].set_ylabel('Ocurrences');  axs[i,0].set_xlabel('Location of the peak \n (relative to the cell cycle)')
axs[0,1].set_title('Buterfilt detrend')
axs[0,0].set_title('Filtermin detrend')

col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):  
    df       = buterfilt_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []
    ax =  axs[i,1]
    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
    
    cell_cycle_df = pd.DataFrame({'normed_time':normed_time_aux,'peaks':aux_peaks})
    cell_cycle_df.set_index('normed_time',inplace=True)
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    n, bins, patches = ax.hist(plot_value, bins=np.arange(0,1,0.05), cumulative=False,density=False,stacked= False,facecolor=colors[i],  alpha=1)
    
    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_hist_normed.pdf', format='pdf')

plt.show()

#%% cum normed

save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order)-1,2, sharex=True, sharey=True,figsize = (8.27, 11.69))

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
for i,c in enumerate(filtmin_conc.conc_order):

    df       = filtmin_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []
    ax =  axs[i-1,0]
    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
        
    cell_cycle_df = pd.DataFrame({'normed_time':normed_time_aux,'peaks':aux_peaks})
    cell_cycle_df.set_index('normed_time',inplace=True)
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    if i !=0: n, bins, patches = ax.hist(plot_value, bins = np.arange(0,1,0.05),cumulative = True, density=True,stacked= False,facecolor=colors[i],  alpha=1)
    if i!=0: ax.plot(bins,bins,color = "black", linewidth=0.4, linestyle='-.',label = "random")

axs[i-1,0].set_ylabel('Density of ocurrences');  axs[i-1,0].set_xlabel('Location of the peak \n (relative to the cell cycle)')
axs[0,1].set_title('Buterfilt detrend')
axs[0,0].set_title('Filtermin detrend')

col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):  
    df       = buterfilt_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []
    ax =  axs[i-1,1]
    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
    
    cell_cycle_df = pd.DataFrame({'normed_time':normed_time_aux,'peaks':aux_peaks})
    cell_cycle_df.set_index('normed_time',inplace=True)
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    if i !=0: n, bins, patches = ax.hist(plot_value, bins=np.arange(0,1,0.05),cumulative = True, density=True,stacked= False,facecolor=colors[i],  alpha=1)
    if i!=0: ax.plot(bins,bins,color = "black", linewidth=0.4, linestyle='-.',label = "random")

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_hist_cum_normed.pdf', format='pdf')

plt.show()
   
   
#%% plotea los histogramas all 
#%%

colors = sns.color_palette('colorblind');
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'
sns.set(context='paper', style='white')
col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'


fig, axs = plt.subplots(2,2, sharex=True, sharey='row',figsize = (8.27, 11.69))

for i,c in enumerate(filtmin_conc.conc_order):
    ax =  axs[0,0]

    df       = filtmin_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []
    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
        
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    n, bin_edges=np.histogram(plot_value, bins = np.arange(0,1,0.05), density=True)
    bin_edges_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    if i!=0: ax.plot(bin_edges_center, n, color=colors[i], linewidth=0.4,  marker='.', markersize=10,label=c)
    
    ax=axs[1,0]
    
    if i!=0: ax.plot(bin_edges_center, np.cumsum(n), color=colors[i], linewidth=0.4, marker='.', markersize=10,label=c)
    
ax.plot(bin_edges_center,bin_edges_center*(20),color = "black", linewidth=0.4, linestyle='-.',label = "random")
    
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'

for i,c in enumerate(buterfilt_conc.conc_order):  
    df       = buterfilt_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_peaks = []
    normed_time_aux = []

    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        normed_time_aux.extend(aux)
        aux_peaks.extend(data[col].notna()*1)
    
    cell_cycle_df = pd.DataFrame({'normed_time':normed_time_aux,'peaks':aux_peaks})
    cell_cycle_df.set_index('normed_time',inplace=True)
    plot_value = np.multiply(normed_time_aux,aux_peaks); plot_value = plot_value[plot_value>0]
    n, bin_edges=np.histogram(plot_value, bins=np.arange(0,1,0.05),density=True)
    bin_edges_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    ax =  axs[0,1]

    if i!=0: ax.plot(bin_edges_center, n, color=colors[i], linewidth=0.4, marker='.', markersize=10,label=c)
    
    ax = axs[1,1]
    if i!=0: ax.plot(bin_edges_center, np.cumsum(n), color=colors[i], linewidth=0.4, marker='.', markersize=10,label=c)
    
ax.plot(bin_edges_center,bin_edges_center*20,color = "black", linewidth=0.4, linestyle='-.',label = "random")
axs[1,0].set_ylabel('Density of ocurrences');  axs[1,0].set_xlabel('Location of the peak \n (relative to the cell cycle)')
axs[0,1].set_title('Buterfilt detrend')
axs[0,0].set_title('Filtermin detrend')
axs[1,1].legend(loc="lower left", fontsize=8,bbox_to_anchor=(0.5, -0.00),
      ncol=2)




fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_hist_all_density.pdf', format='pdf')

plt.show()

#%%


save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

colors = sns.color_palette('colorblind')
sns.set(context='paper', style='dark')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)
cmap = sns.color_palette("YlGnBu", 5)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'

for i,c in enumerate(filtmin_conc.conc_order):
    ax=axs[i,0]
    df       = filtmin_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_int_peaks,time_interval,cell_aux = [],[],[]

    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        prev_limit = 0
        for limit in np.arange(0.05,1.05,0.05):
            #no toma el ultimo punyo
            mask = (aux < limit) & (aux>prev_limit)
            aux_int_peaks.append((data[col].notna()*1)[mask].sum())
            prev_limit = limit
           
            
        cell_aux.extend([cell] * len(np.arange(0.05,1.05,0.05)) )
        time_interval.extend(np.arange(0.05,1.05,0.05))
    
    index = [np.array(cell_aux),np.array(time_interval)];  index_tuples = list(zip(*index))
    index = pd.MultiIndex.from_tuples(index_tuples, names=['cell', 'time_interval'])
    
    cell_cycle_df = pd.DataFrame({'peaks':aux_int_peaks})
    cell_cycle_df = cell_cycle_df.set_index(index); df.sort_index(inplace=True)

    im = sns.heatmap(cell_cycle_df.unstack(),cbar_kws={'label': 'Peaks'},cbar = False,cmap = cmap,ax=ax,xticklabels = np.arange(0.05,1.05,0.05).round(2))
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel(''); #ax.set_xticks([])
    ax.set_yticks([])
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    if i == (len(filtmin_conc.conc_order)-1): ax.set_xticks(np.arange(0.05,1.05,0.05).round(2)); 

axs[i,0].set_ylabel('Cells');  axs[i,0].set_xlabel('Time \n (in frames and relative to the cell cycle)')
axs[0,0].set_title('Filtermin detrend \n Concentration an_0ng')
axs[0,1].set_title('Buterfilt detrend \n Concentration an_0ng')


col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):  
    ax=axs[i,1]
    df       = buterfilt_conc.data[c][[col,'FRAME']]; df.sort_index(inplace=True)  
    aux_int_peaks,time_interval,cell_aux = [],[],[]

    for cell,data in df.groupby(level = 'cell'):
        aux = data.FRAME-data.FRAME.iloc[0]
        aux = aux / aux.iloc[-1]
        prev_limit = 0
        for limit in np.arange(0.05,1.05,0.05):
            #no toma el ultimo punyo
            mask = (aux < limit) & (aux>prev_limit)
            aux_int_peaks.append((data[col].notna()*1)[mask].sum())
            prev_limit = limit
           
            
        cell_aux.extend([cell] * len(np.arange(0.05,1.05,0.05)) )
        time_interval.extend(np.arange(0.05,1.05,0.05))
    
    index = [np.array(cell_aux),np.array(time_interval)];  index_tuples = list(zip(*index))
    index = pd.MultiIndex.from_tuples(index_tuples, names=['cell', 'time_interval'])
    
    cell_cycle_df = pd.DataFrame({'peaks':aux_int_peaks})
    cell_cycle_df = cell_cycle_df.set_index(index); df.sort_index(inplace=True)

    im = sns.heatmap(cell_cycle_df.unstack(),cbar_kws={'label': 'Peaks'},cbar = False,cmap = cmap,ax=ax,xticklabels = np.arange(0.05,1.05,0.05).round(2))
    ax.set_ylabel(''); ax.set_xlabel(''); ax.set_xticks([])
    ax.set_yticks([])
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.2)
mappable = im.get_children()[0]
#cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cb_ax = fig.add_axes([0.47, 0.14, 0.2, 0.005])
cbar = fig.colorbar(mappable ,cax=cb_ax,orientation = 'horizontal')

cbar.set_ticks([1,4])
cbar.set_label(label='Number of peaks')
cbar.ax.tick_params(labelsize = 10)

fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_TS.pdf', format='pdf')

plt.show()
#%%
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks,frames_count = np.array([]),np.array([])

    for cell,data in df.groupby(level = 'cell'):
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            frames_count = np.ones(len(data)) + np.pad(frames_count,(0,len(data)-len(frames_count)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            frames_count = frames_count + np.pad(np.ones(len(data)),(0,len(frames_count)-len(data)),'constant', constant_values=(0,0))
    
    peaks_sum =  np.add.reduceat(aux_peaks, np.arange(0, len(aux_peaks), 10))
    frames_sum = np.add.reduceat(frames_count, np.arange(0, len(frames_count), 10))
    
    peaks_dens = peaks_sum / frames_sum
    
    cell_cycle = pd.DataFrame({'peaks_dens':peaks_dens[:-1],'peaks_sum':peaks_sum[:-1], 'frames_sum':frames_sum[:-1],'error':( np.sqrt(((1-peaks_dens)*peaks_dens)/frames_sum))[:-1],'Time_aligned_beg':np.arange(0,len(aux_peaks),10)[:-1]})
    cell_cycle.set_index('Time_aligned_beg',inplace=True,drop=False)
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_dens",data=cell_cycle,legend=False,ax = axs[i,0])    
    axs[i,0].errorbar( x= list(cell_cycle["Time_aligned_beg"]), y= list(cell_cycle["peaks_dens"]), 
       yerr=list(cell_cycle['error']), fmt=' ',alpha=0.5)
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')
    #if i!= len(filtmin_conc.conc_order)-1:    axs[i,0].set_yticks([]) and  axs[i,0].set_xticks([])



axs[i,0].set_ylabel('Density of peaks \n (peaks/frames)');  axs[i,0].set_xlabel('Time (aligned at the beg)') 
axs[0,0].set_title('Filtmin detrend')
   
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks,frames_count = np.array([]),np.array([])
    
    for cell,data in df.groupby(level = 'cell'):
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            frames_count = np.ones(len(data)) + np.pad(frames_count,(0,len(data)-len(frames_count)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            frames_count = frames_count + np.pad(np.ones(len(data)),(0,len(frames_count)-len(data)),'constant', constant_values=(0,0))
    
    peaks_sum =  np.add.reduceat(aux_peaks, np.arange(0, len(aux_peaks), 10))
    frames_sum = np.add.reduceat(frames_count, np.arange(0, len(frames_count), 10))
    
    peaks_dens = peaks_sum / frames_sum
    
    
    cell_cycle = pd.DataFrame({'peaks_dens':peaks_dens[:-1],'peaks_sum':peaks_sum[:-1], 'frames_sum':frames_sum[:-1],'error':( np.sqrt(((1-peaks_dens)*peaks_dens)/frames_sum))[:-1],'Time_aligned_beg':np.arange(0,len(aux_peaks),10)[:-1]})
    cell_cycle.set_index('Time_aligned_beg',inplace=True,drop=False)
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_dens",data=cell_cycle,legend=False,ax = axs[i,1])    
    axs[i,1].errorbar( x= list(cell_cycle["Time_aligned_beg"]), y= list(cell_cycle["peaks_dens"]), 
       yerr=list(cell_cycle['error']), fmt=' ',alpha=0.5)
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')


axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_norm_BEG.pdf', format='pdf')

plt.show()
#%%poisson o binomial
    
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))
fig2, axs2 = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))
fig3, axs3 = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))
cmap = sns.color_palette("YlGnBu", 5)


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    M = filtmin_conc.traze_max_values[i]
    M=651
    cell_index,frame_index,peaks_values = [],[],[]
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data.notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(np.pad(peaks_sum,(0,M//10-len(peaks_sum)),'constant', constant_values=(-100,-100)))
        frame_index.extend(np.arange(0,M,10)[:-1])
        cell_index.extend(np.ones(M//10) * [cell])
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'peaks_values':peaks_values})
    cell_cycle.set_index(['cell','frames_int'],inplace = True,drop=True)
    cell_cycle.peaks_values[cell_cycle.peaks_values < 0] = None
    
    
    
    ax = axs[i,0]                  #np.arange(0,M,10)[:-1].round(2)
#    im = sns.heatmap(cell_cycle.unstack(),cbar_kws={'label': 'Peaks'},cbar = False,cmap = cmap,ax=ax,xticklabels = 10)    
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel(''); ax.set_xticks([]); 
    #if i == (len(filtmin_conc.conc_order)-1) : ax.set_xticks(list(np.linspace(0,1,10)*ax.get_xlim()[1]))
    ax.set_yticks([]); 
    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Cells'); ax.set_xlabel('Frames')
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)




    MEAN = cell_cycle.mean(level='frames_int').values[:,0];    mean = np.linspace(0,max(MEAN),10)
    VAR = cell_cycle.std(level='frames_int').values[:,0]
    
    ax=axs2[i,0]
    ax.plot(mean,mean,color='red', linestyle=':',label='poisson')
    ax.plot(mean,mean*(1-mean/10),color='green', linestyle=':',label='binomial')
    ax.plot(MEAN,VAR,'o')
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel(''); ax.set_xticks([]); 
    ax.set_yticks([0,2])
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Variance'); ax.set_xlabel('Mean');ax.legend(); ax.set_xticks([0,2]); ax.set_yticks([0,1])



    ax= axs3[i,0]
    cell_cycle.mean(level='frames_int').plot(marker='o',ax= ax,label=False,legend=False)
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel(''); ax.set_xticks([]); 
    ax.set_yticks([0,2]);
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Mean'); ax.set_xlabel('Time');ax.set_xticks([0,M]);ax.set_yticks([0,2]);

    
    #axs3[i,0].errorbar( x= list(cell_cycle.mean(level='frames_int').index), y= list(cell_cycle.mean(level='frames_int').values[:,0]), yerr=list(cell_cycle.std(level='frames_int').values[:,0]), fmt='',alpha=0.5)
    
    ######
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    M = buterfilt_conc.traze_max_values[i]
    M = 651
    cell_index,frame_index,peaks_values = [],[],[]
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data.notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(np.pad(peaks_sum,(0,M//10-len(peaks_sum)),'constant', constant_values=(-100,-100)))
        frame_index.extend(np.arange(0,M,10)[:-1])
        cell_index.extend(np.ones(M//10) * [cell])
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'peaks_values':peaks_values})
    cell_cycle.set_index(['cell','frames_int'],inplace = True,drop=True)
    cell_cycle.peaks_values[cell_cycle.peaks_values < 0] = None
    
    ax = axs[i,1]
    im = sns.heatmap(cell_cycle.unstack(),cbar_kws={'label': 'Peaks'},cbar = False,cmap = cmap,ax=ax,
                xticklabels = 10)    
    if i != (len(buterfilt_conc.conc_order)-1) : ax.set_xlabel('');ax.set_xticks([]); 
    ax.set_ylabel(''); ax.set_yticks([]); 
    if i == (len(buterfilt_conc.conc_order)-1) : ax.set_xlabel('Frames')
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    

    MEAN = cell_cycle.mean(level='frames_int').values[:,0];    mean = np.linspace(0,max(MEAN),10)
    VAR = cell_cycle.std(level='frames_int').values[:,0]
    
    ax=axs2[i,1]
    ax.plot(mean,mean,color='red', linestyle=':',label='poisson')
    ax.plot(mean,mean*(1-mean/10),color='green', linestyle=':',label='binomial')
    ax.plot(MEAN,VAR,'o')
    if i != (len(buterfilt_conc.conc_order)-1) :ax.set_xlabel(''); ax.set_xticks([]); 
    ax.set_yticks([]); ax.set_ylabel('');
    if i == (len(buterfilt_conc.conc_order)-1) :  ax.set_xlabel('Mean');ax.set_xticks([0,2]); 
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)


    ax= axs3[i,1]
    cell_cycle.mean(level='frames_int').plot(marker='o',ax= ax,label=False,legend=False)
    if i != (len(buterfilt_conc.conc_order)-1) :ax.set_xlabel(''); ax.set_xticks([]); 
    ax.set_yticks([]);ax.set_ylabel('')
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    if i == (len(buterfilt_conc.conc_order)-1) :  ax.set_xlabel('Time');ax.set_xticks([0,M]);
    
    #axs3[i,1].errorbar( x= list(cell_cycle.mean(level='frames_int').index), y= list(cell_cycle.mean(level='frames_int').values[:,0]), yerr=list(cell_cycle.std(level='frames_int').values[:,0]), fmt='',alpha=0.5)

  
axs[0,1].set_title('Buterfilt detrend  \n an_0ng')
axs[0,0].set_title('Filtermin detrend \n an_0ng')
axs2[0,1].set_title('Buterfilt detrend  \n an_0ng')
axs2[0,0].set_title('Filtermin detrend \n an_0ng')
axs3[0,1].set_title('Buterfilt detrend  \n an_0ng')
axs3[0,0].set_title('Filtermin detrend \n an_0ng')
    
    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.2)
fig.savefig(str(save_folder)+'2018_E1_28_TS_bin.pdf', format='pdf')

fig2.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.2)
fig2.savefig(str(save_folder)+'2018_E1_28_Dist_plot.pdf', format='pdf')

fig3.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.2)
fig3.savefig(str(save_folder)+'2018_E1_28_mean_value.pdf', format='pdf')

plt.show()
#%%
    #plt.scatter(x=cell_cycle.mean(level='frames_int'),y=cell_cycle.mean(level='frames_int'))
    #plt.scatter(x=cell_cycle.mean(level='frames_int'),y=cell_cycle.mean(level='frames_int')*10*(1-cell_cycle.mean(level='frames_int')));plt.show()
    
  
    
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')
    #if i!= len(filtmin_conc.conc_order)-1:    axs[i,0].set_yticks([]) and  axs[i,0].set_xticks([])



axs[i,0].set_ylabel('Density of peaks \n (peaks/frames)');  axs[i,0].set_xlabel('Time (aligned at the beg)') 
axs[0,0].set_title('Filtmin detrend')
   
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks,frames_count = np.array([]),np.array([])
    
    for cell,data in df.groupby(level = 'cell'):
        if len(data) > len(aux_peaks):
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(0,len(data)-len(aux_peaks)),'constant', constant_values=(0,0))
            frames_count = np.ones(len(data)) + np.pad(frames_count,(0,len(data)-len(frames_count)),'constant', constant_values=(0,0))

        else:
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(0,len(aux_peaks)-len(data)),'constant', constant_values=(0,0))
            frames_count = frames_count + np.pad(np.ones(len(data)),(0,len(frames_count)-len(data)),'constant', constant_values=(0,0))
    
    peaks_sum =  np.add.reduceat(aux_peaks, np.arange(0, len(aux_peaks), 10))
    frames_sum = np.add.reduceat(frames_count, np.arange(0, len(frames_count), 10))
    
    peaks_dens = peaks_sum / frames_sum
    
    
    cell_cycle = pd.DataFrame({'peaks_dens':peaks_dens[:-1],'peaks_sum':peaks_sum[:-1], 'frames_sum':frames_sum[:-1],'error':( np.sqrt(((1-peaks_dens)*peaks_dens)/frames_sum))[:-1],'Time_aligned_beg':np.arange(0,len(aux_peaks),10)[:-1]})
    cell_cycle.set_index('Time_aligned_beg',inplace=True,drop=False)
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_dens",data=cell_cycle,legend=False,ax = axs[i,1])    
    axs[i,1].errorbar( x= list(cell_cycle["Time_aligned_beg"]), y= list(cell_cycle["peaks_dens"]), 
       yerr=list(cell_cycle['error']), fmt=' ',alpha=0.5)
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')


axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
#fig.savefig(str(save_folder)+'2018_E1_23_cell_cycle_norm_BEG.pdf', format='pdf')

plt.show()


#%%



save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks,frames_count = np.array([]),np.array([])

    for cell,data in df.groupby(level = 'cell'):
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(frames_count)
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            frames_count = np.ones(len(data)) + np.pad(frames_count,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(frames_count)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            frames_count = frames_count + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    
    peaks_sum =  np.add.reduceat(aux_peaks, np.arange(0, len(aux_peaks), 10))
    frames_sum = np.add.reduceat(frames_count, np.arange(0, len(frames_count), 10))
    
    peaks_dens = peaks_sum / frames_sum
    
    cell_cycle = pd.DataFrame({'peaks_dens':peaks_dens[:-1],'peaks_sum':peaks_sum[:-1], 'frames_sum':frames_sum[:-1],'error':( np.sqrt(((1-peaks_dens)*peaks_dens)/frames_sum))[:-1],'Time_aligned_beg':np.arange(0,len(aux_peaks),10)[:-1]})
    cell_cycle.set_index('Time_aligned_beg',inplace=True,drop=False)
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_dens",data=cell_cycle,legend=False,ax = axs[i,0])    
    axs[i,0].errorbar( x= list(cell_cycle["Time_aligned_beg"]), y= list(cell_cycle["peaks_dens"]), 
       yerr=list(cell_cycle['error']), fmt=' ',alpha=0.5)
    axs[i,0].set_ylabel('');  axs[i,0].set_xlabel('')
    #if i!= len(filtmin_conc.conc_order)-1:    axs[i,0].set_yticks([]) and  axs[i,0].set_xticks([])



axs[i,0].set_ylabel('Density of peaks \n (peaks/frames)');  axs[i,0].set_xlabel('Time (aligned at the mid)') 
axs[0,0].set_title('Filtmin detrend')
   
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c][col]; df.sort_index(inplace=True)  
    aux_peaks,frames_count = np.array([]),np.array([])
       
    
    for cell,data in df.groupby(level = 'cell'):
        if len(data) > len(aux_peaks):
            N = len(data)-len(aux_peaks); n = len(data)-len(frames_count)
            aux_peaks = np.array(data.notna()*1) + np.pad(aux_peaks,(N//2,N//2+N%2),'constant', constant_values=(0,0))
            frames_count = np.ones(len(data)) + np.pad(frames_count,(n//2,n//2+n%2),'constant', constant_values=(0,0))

        else:
            N = len(aux_peaks)-len(data);n = len(frames_count)-len(data)
            aux_peaks = aux_peaks + np.pad(np.array(data.notna()*1),(N//2,N//2+N%2),'constant', constant_values=(0,0))
            frames_count = frames_count + np.pad(np.ones(len(data)),(n//2,n//2+n%2),'constant', constant_values=(0,0))
    
    peaks_sum =  np.add.reduceat(aux_peaks, np.arange(0, len(aux_peaks), 10))
    frames_sum = np.add.reduceat(frames_count, np.arange(0, len(frames_count), 10))
    
    peaks_dens = peaks_sum / frames_sum
    
    
    cell_cycle = pd.DataFrame({'peaks_dens':peaks_dens[:-1],'peaks_sum':peaks_sum[:-1], 'frames_sum':frames_sum[:-1],'error':( np.sqrt(((1-peaks_dens)*peaks_dens)/frames_sum))[:-1],'Time_aligned_beg':np.arange(0,len(aux_peaks),10)[:-1]})
    cell_cycle.set_index('Time_aligned_beg',inplace=True,drop=False)
    
    sns.scatterplot(x="Time_aligned_beg", y="peaks_dens",data=cell_cycle,legend=False,ax = axs[i,1])    
    axs[i,1].errorbar( x= list(cell_cycle["Time_aligned_beg"]), y= list(cell_cycle["peaks_dens"]), 
       yerr=list(cell_cycle['error']), fmt=' ',alpha=0.5)
    axs[i,1].set_ylabel('');  axs[i,1].set_xlabel('')


axs[0,1].set_title('Buterfilt detrend')

    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.07)
fig.savefig(str(save_folder)+'2018_E1_28_cell_cycle_norm_MID.pdf', format='pdf')

plt.show()

#%%

    
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='ticks')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex=True, sharey=True,figsize = (8.27, 11.69))
cmap = sns.color_palette("YlGnBu", 5)


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c]; df.sort_index(inplace=True)  
    M = filtmin_conc.traze_max_values[i]
    cell_index,frame_index,peaks_values = [],[],[]
    frame_exp_index = []
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data[col].notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(peaks_sum)
        frame_index.extend(np.arange(0,len(data),10)[:-1])
        cell_index.extend(np.ones(len(peaks_sum)) * [cell])
        frame_exp_index.extend(np.arange(0,len(data),10)[:-1] + data.FRAME.iloc[0] )
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'frames_exp_int':frame_exp_index,'peaks_values':peaks_values})
    #cell_cycle.set_index(['cell','frames_int'],inplace = True,drop=True)
    #cell_cycle.peaks_values[cell_cycle.peaks_values < 0] = None
    
    
    
    ax = axs[i,0]                  #np.arange(0,M,10)[:-1].round(2)
    legend=False
    sns.scatterplot(x="frames_int", y="frames_exp_int",hue = 'peaks_values', 
                    s=5,palette="YlGnBu",data=cell_cycle,legend=legend,ax = ax,alpha=0.7,marker = '+')
    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Exp time (10 frames)'); ax.set_xlabel('Time rel to cell cycle \n (10 frames)')
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel(''); 
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    ax.set_xlim(0,M)
    

col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'

for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c]; df.sort_index(inplace=True)  
    M = filtmin_conc.traze_max_values[i]
    cell_index,frame_index,peaks_values = [],[],[]
    frame_exp_index = []
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data[col].notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(peaks_sum)
        frame_index.extend(np.arange(0,len(data),10)[:-1])
        cell_index.extend(np.ones(len(peaks_sum)) * [cell])
        frame_exp_index.extend(np.arange(0,len(data),10)[:-1] + data.FRAME.iloc[0] )
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'frames_exp_int':frame_exp_index,'peaks_values':peaks_values})
    
    ax = axs[i,1]      
    legend=False
    if i == len(buterfilt_conc.conc_order)-1: legend='brief' ;    
    sns.scatterplot(x="frames_int", y="frames_exp_int",hue = 'peaks_values', 
                    s=5,palette="YlGnBu",data=cell_cycle,legend=legend,ax = ax,alpha=0.7,marker = '+')
    ax.set_ylabel(''); ax.set_xlabel(''); 
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    ax.set_xlim(0,M);ax.set_ylim(0,M)
    if i == len(buterfilt_conc.conc_order)-1:ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. ,ncol=1,frameon=False)

fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.05, hspace=0.3)
axs[0,1].set_title('Buterfilt detrend  \n an_0ng')
axs[0,0].set_title('Filtermin detrend \n an_0ng')
fig.savefig(str(save_folder)+'2018_E1_28_TvsT.pdf', format='pdf')


plt.show()

#%% nuevi

save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

col = 'DT_MEAN_INTENSITY_2500PEAKS_O2'
colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),2, sharex='row', sharey='row' ,figsize = (8.27, 11.69))
cmap = sns.color_palette("YlGnBu", 5)


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c]; df.sort_index(inplace=True)  
    cell_index,frame_index,peaks_values = [],[],[]
    frame_exp_index = []
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data[col].notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(peaks_sum)
        frame_index.extend(np.arange(0,len(data),10)[:-1])
        cell_index.extend(np.ones(len(peaks_sum)) * [cell])
        frame_exp_index.extend(np.arange(0,len(data),10)[:-1] + data.FRAME.iloc[0] )
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'frames_exp_int':frame_exp_index,'peaks_values':peaks_values})

    aux_mean,aux_var=[],[]
    aux_frames,aux_frames_exp=[],[]
    max_x, max_y= 0,0
    if max_x<cell_cycle.frames_int.max()+50: max_x = cell_cycle.frames_int.max()+50
    if max_y<cell_cycle.frames_exp_int.max()+50: max_y = cell_cycle.frames_exp_int.max()+50

    for frames_int_interval in np.arange(0,cell_cycle.frames_int.max()+50,50):
        mask = (cell_cycle.frames_int < frames_int_interval+50) & (cell_cycle.frames_int > frames_int_interval)
        df = cell_cycle[mask]
        for frames_exp_int_interval in np.arange(cell_cycle.frames_exp_int.min(),cell_cycle.frames_exp_int.max()+50,50):
            mask2 = (df.frames_exp_int<frames_exp_int_interval+50) & (df.frames_exp_int>frames_exp_int_interval)
            aux_mean.append(df[mask2].peaks_values.mean())
            aux_var.append(df[mask2].peaks_values.var())
            aux_frames.append(frames_int_interval)
            aux_frames_exp.append(frames_exp_int_interval)
        
    sq_cell_cycle_mean = pd.DataFrame({'frames_int': aux_frames,'frames_int_exp':aux_frames_exp,' ':aux_mean})
    sq_cell_cycle_mean.set_index(['frames_int','frames_int_exp'],inplace = True,drop=True)
    sq_cell_cycle_var = pd.DataFrame({'frames_int': aux_frames,'frames_int_exp':aux_frames_exp,' ':aux_var})
    sq_cell_cycle_var.set_index(['frames_int','frames_int_exp'],inplace = True,drop=True)
 
    #lo mudo es peak_mean y peak_var
    
    ax = axs[i,0]                  
    sns.heatmap(sq_cell_cycle_mean.unstack().T,annot=True,annot_kws={'fontsize':5},
                cbar=False,ax = ax,cbar_kws={'label': 'Peaks'},cmap="YlGnBu",vmin=0,vmax=2)

    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel('');
    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Experimental Time'); ax.set_xlabel('Time \n (in frames and relative to the cell cycle)');
    
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)

    #ax.set_xlim(0,600)
    #ax.set_ylim(0,max_y)
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    ax.invert_yaxis()
    ax = axs[i,1]                  
    sns.heatmap(sq_cell_cycle_var.unstack().T,annot=True,annot_kws={'fontsize':5},
                cbar=False,ax = ax,cbar_kws={'label': 'Peaks'},cmap="YlGnBu",vmin=0,vmax=2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)

    ax.set_ylabel(''); ax.set_xlabel('');
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)

    #ax.set_xlim(0,600)
    #ax.set_ylim(0,max_y)
    ax.invert_yaxis()
axs[0,0].set_title('Mean Peaks - Filtermin detrend \n Concentration an_0ng')
axs[0,1].set_title('Variance Peaks - Filtermin detrend \n Concentration an_0ng')
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.9,wspace=0.3, hspace=0.7)
fig.savefig(str(save_folder)+'2018_E1_28_TvsT_heatmap_filtmin.pdf', format='pdf')


#%%

#Falta poner que sea menor tam,bien
col = 'MEAN_INTENSITY_BP_2400PEAKS_O2'
 
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(buterfilt_conc.conc_order),2, sharex='row', sharey='row' ,figsize = (8.27, 11.69))
cmap = sns.color_palette("YlGnBu", 5)


for i,c in enumerate(buterfilt_conc.conc_order):
    df       = buterfilt_conc.data[c]; df.sort_index(inplace=True)  
    cell_index,frame_index,peaks_values = [],[],[]
    frame_exp_index = []
    
    for cell,data in df.groupby(level = 'cell'):
        peaks_sum   =   np.add.reduceat(np.array(data[col].notna()*1), np.arange(0, len(data), 10))[:-1]
        peaks_values.extend(peaks_sum)
        frame_index.extend(np.arange(0,len(data),10)[:-1])
        cell_index.extend(np.ones(len(peaks_sum)) * [cell])
        frame_exp_index.extend(np.arange(0,len(data),10)[:-1] + data.FRAME.iloc[0] )
   
    
    cell_cycle = pd.DataFrame({'cell': cell_index,'frames_int':frame_index,'frames_exp_int':frame_exp_index,'peaks_values':peaks_values})

    aux_mean,aux_var=[],[]
    aux_frames,aux_frames_exp=[],[]
    max_x, max_y= 0,0
    if max_x<cell_cycle.frames_int.max()+50: max_x = cell_cycle.frames_int.max()+50
    if max_y<cell_cycle.frames_exp_int.max()+50: max_y = cell_cycle.frames_exp_int.max()+50

    for frames_int_interval in np.arange(0,cell_cycle.frames_int.max()+50,50):
        mask = (cell_cycle.frames_int < frames_int_interval+50) & (cell_cycle.frames_int > frames_int_interval)
        df = cell_cycle[mask]
        for frames_exp_int_interval in np.arange(cell_cycle.frames_exp_int.min(),cell_cycle.frames_exp_int.max()+50,50):
            mask2 = (df.frames_exp_int<frames_exp_int_interval+50) & (df.frames_exp_int>frames_exp_int_interval)
            aux_mean.append(df[mask2].peaks_values.mean())
            aux_var.append(df[mask2].peaks_values.var())
            aux_frames.append(frames_int_interval)
            aux_frames_exp.append(frames_exp_int_interval)
        
    sq_cell_cycle_mean = pd.DataFrame({'frames_int': aux_frames,'frames_int_exp':aux_frames_exp,' ':aux_mean})
    sq_cell_cycle_mean.set_index(['frames_int','frames_int_exp'],inplace = True,drop=True)
    sq_cell_cycle_var = pd.DataFrame({'frames_int': aux_frames,'frames_int_exp':aux_frames_exp,' ':aux_var})
    sq_cell_cycle_var.set_index(['frames_int','frames_int_exp'],inplace = True,drop=True)
 
    #lo mudo es peak_mean y peak_var
    
    ax = axs[i,0]                  
    sns.heatmap(sq_cell_cycle_mean.unstack().T,annot=True,annot_kws={'fontsize':5},
                cbar=False,ax = ax,cbar_kws={'label': 'Peaks'},cmap="YlGnBu",vmin=0,vmax=1)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)
    
    if i != (len(buterfilt_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel('');
    if i == (len(buterfilt_conc.conc_order)-1) : ax.set_ylabel('Experimental Time'); ax.set_xlabel('Time \n (in frames and relative to the cell cycle)');

    #ax.set_xlim(0,max_x)
    #ax.set_ylim(0,max_y)
    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)
    ax.invert_yaxis()
    ax = axs[i,1]                  
    sns.heatmap(sq_cell_cycle_var.unstack().T,annot=True,annot_kws={'fontsize':5},
                cbar=False,ax = ax,cbar_kws={'label': 'Peaks'},cmap="YlGnBu",vmin=0,vmax=1)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)
    ax.set_ylabel(''); ax.set_xlabel('');

    if i!=0: ax.set_title('Concentration {}'.format(c),fontsize = 8)

    #ax.set_xlim(0,max_x)
    #ax.set_ylim(0,max_y)
    ax.invert_yaxis()
axs[0,0].set_title('Mean Peaks - Buterfilt detrend \n Concentration an_0ng')
axs[0,1].set_title('Variance Peaks - Buterfilt detrend \n Concentration an_0ng')
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.9,wspace=0.3, hspace=0.7)
fig.savefig(str(save_folder)+'2018_E1_28_TvsT_heatmap_buterfilt.pdf', format='pdf')


plt.show()


#%%
#%%number of cells
    
save_folder = '/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Cell_cycle_analysis/'

colors = sns.color_palette('colorblind')
sns.set(context='paper', style='white')
plt.rc('axes.spines',top=False,bottom=True,left=True,right=False)

fig, axs = plt.subplots(len(filtmin_conc.conc_order),1, sharex=True, sharey=True,figsize = (8.27, 11.69))
cmap = sns.color_palette("YlGnBu", 5)


for i,c in enumerate(filtmin_conc.conc_order):
    df       = filtmin_conc.data[c]; df.sort_index(inplace=True)  
    M = np.max(filtmin_conc.traze_max_values)
    cell_sum = np.zeros(M+1)
    
    for cell,data in df.groupby(level = 'cell'):
        cell_sum_aux =  np.array(list(np.zeros(data.FRAME.iloc[0]))+list(np.ones(len(data))) + list(np.zeros(M-data.FRAME[-1])))
        cell_sum = cell_sum + cell_sum_aux

    
    
    ax = axs[i]                  #np.arange(0,M,10)[:-1].round(2)
    ax.plot(cell_sum,'+',markersize = 1)
    if i != (len(filtmin_conc.conc_order)-1) : ax.set_ylabel(''); ax.set_xlabel('');

    if i == (len(filtmin_conc.conc_order)-1) : ax.set_ylabel('Number of cells'); ax.set_xlabel('Frames')
    ax.set_title('Concentration {}'.format(c),fontsize = 8)
    ax.set_xlim(0,M)


   
    
fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.7,wspace=0.07, hspace=0.2)
fig.savefig(str(save_folder)+'2018_E1_28_NumberOfCells.pdf', format='pdf')


plt.show()
    