
#sanity check
for c in cells_conc:
    df = conc_data[c];df.sort_index(inplace=True);
    for cells, data in df.groupby(level='cell'):
        res = data['IPI'].dropna().values-mixed_dt(data)
        for i,j in enumerate(res):
            if j< 0: print(c,cells,i,j)
        

     #%% Cell analysis
#Total cells
for c in conc_data:
    df = conc_data[c]
    c_c = 0
    for cells,data in df.groupby(level='cell'):
        c_c=c_c+1
    print(c,c_c)
     
#Total cells with pulses    
for c in conc_data:
    df = conc_data[c]
    c_c = 0
    for cells,data in df.groupby(level='cell'):
        if data.amp_peaks.dropna().count() !=0:
            c_c = c_c+1
    print(c,c_c)
    
#Cells with isolated pulses
for c in conc_data:
    df = conc_data[c]
    c_c = 0
    for cells,data in df.groupby(level='cell'):
        if is_isolated_cell(data)!=0:
            c_c=c_c+1
    print(c,c_c)

#Cells with consecutive pulses
for c in conc_data:
    df = conc_data[c]
    c_c = 0
    for cells,data in df.groupby(level='cell'):
        if sum(is_consecutive_cell(data))!=0:
            c_c=c_c+1
    print(c,c_c)

#%%
    #Number of pulses
for c in conc_data:
    df = conc_data[c]
    print(c,df.amp_peaks.dropna().count())

    #Number of pulses
for c in conc_data:
    df = conc_data[c]
    print(c,df.IPI.dropna().count())
    
    #outliers
for c in conc_data:
    df = conc_data[c]
    outliers = [[i,j] for i,j in zip(df['IPI'].dropna().values-mixed_dt(df),mixed_dt(df))]
    print(c,np.sum([1 for i,j in outliers if i>60 or j >60]))
    
#%%
    
         
def get_consecutive_data(conc_data):
    cons_box = {}; isol_box = {}
    for k,c in enumerate(conc_data):
        df = conc_data[c] 
        c_aux=[]; i_aux=[]

        for cells, data in df.groupby(level='cell'):
            consecutive_ = is_consecutive_cell(data);c_aux.append(np.nan_to_num(sum(consecutive_)/data.amp_peaks.count(),nan=-1))
            isolated_ = is_isolated_cell(data);i_aux.append(np.nan_to_num(isolated_/data.amp_peaks.count(),nan=-1))
        cons_box[c] = c_aux
        isol_box[c] = i_aux
    cons_box = delete_non_pulsatile_cells(cons_box)
    isol_box = delete_non_pulsatile_cells(isol_box)
    return(cons_box,isol_box)

def delete_non_pulsatile_cells(box):
    new_box={}
    for condition in box:
        new_box[condition]= [i for i in box[condition] if i>-0.5]
    return(new_box)

    #%%

from scipy import stats

with pd.ExcelWriter(save_folder+'tests.xlsx') as writer:
     results_list = {}  
     sheet_names = ['amplitude','duration','IPI']
     labels = ['amp_peaks','dt_peaks','IPI']
     for l,shn in zip(labels,sheet_names):         
         Results = pd.DataFrame(index=conc_data.keys(),columns=conc_data.keys()) 
         for m,key1 in enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']) :
                 df = conc_data[key1]
                 data1 = df[l].dropna().values
                 for n,key2 in  enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']):
                     df = conc_data[key2]
                     data2 = df[l].dropna().values
                     Results[key1][key2] = stats.ks_2samp(data1,data2)[1]
         results_list[shn] = Results
         Results.to_excel(writer, sheet_name=shn)
     
     sheet_names = 'pulse rate'
     Results = pd.DataFrame(index=conc_data.keys(),columns=conc_data.keys()) 
     for m,key1 in enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']) :
             df = conc_data[key1]
             data1 = (df['amp_peaks'].groupby(level='cell').count() / df['amp_peaks'].groupby(level='cell').size()).replace([None],[0])
             for n,key2 in  enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']):
                 df = conc_data[key2]
                 data2 =  (df['amp_peaks'].groupby(level='cell').count() / df['amp_peaks'].groupby(level='cell').size()).replace([None],[0])
                 Results[key1][key2] = stats.ks_2samp(data1,data2)[1]
     results_list[shn] = Results
     Results.to_excel(writer, sheet_name=sheet_names)
  
    
     sheet_names = 'consecutive pulses'
     cons, isol = get_consecutive_data(conc_data)
     Results = pd.DataFrame(index=conc_data.keys(),columns=conc_data.keys()) 
     for m,key1 in enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']) :
             data1 = cons[key1]
             for n,key2 in  enumerate(['an_0ng','an_2-5ng', 'an_5ng', 'an_20ng']):
                 df = conc_data[key2]
                 data2 =  cons[key2]
                 Results[key1][key2] = stats.ks_2samp(data1,data2)[1]
     results_list[shn] = Results
     Results.to_excel(writer, sheet_name=sheet_names)            

