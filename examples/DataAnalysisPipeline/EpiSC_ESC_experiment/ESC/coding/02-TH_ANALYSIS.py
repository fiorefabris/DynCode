import pickle
from time import time
import itertools
import multiprocessing as mp
from functools import partial
import pandas as pd
import numpy as np
#%%
def peak_finding(df_dic,params):
    """ Identifies peaks of a time serie:
            - amp_th and pendiente_th
            - df_dic (with identified peaks) """
    
    
    c,amp_th,pendiente_th  = params[0],params[1],params[2]
    df = df_dic[c]
    
    main_aux_amp,main_aux_dt = [],[]
    main_aux_max,main_aux_min = [],[]
    main_aux_raise,main_aux_fall = [],[]
    
    for cells, data in df.groupby(level='cell'):
        #Auxiliary variables for storing the calculations of the loop
        L = data.shape[0]
 
        aux_amp,aux_dt = [None]*L,[None]*L
        aux_max,aux_min = [None]*L,[None]*L       
        aux_raise,aux_fall = [None]*L,[None]*L
 
    
        col_mins = 'MIN_O1'
        mask_mins = data[col_mins].isna()
        mins = data[~ mask_mins][col_mins]
        
        col_max = 'sm_MEAN_INTENSITY_PEAKS_O1'
        mask_max = data[col_max].isna()
        maxs = data[~ mask_max][col_max]
        
        maxs_index = maxs.index.get_level_values(1)
        mins_index = mins.index.get_level_values(1)
        
        # Ensures that there are minima at the borders of the TS
        test1= True
        while test1:
            if maxs_index[0] <= mins_index[0]:
                maxs.drop([maxs_index[0]],level=1,inplace=True)
                maxs_index = maxs.index.get_level_values(1)
            else:
                test1=False
                
                
        test2=True
        while test2:
            if mins_index[-1] < maxs_index[-1]:
                maxs.drop([maxs_index[-1]],level=1,inplace=True)
                maxs_index = maxs.index.get_level_values(1)
            else:
                test2 = False
        
        # Searchs for the peaks
        amp_ant,idx_ant,idx_mins_der_ant,idx_mins_izq_ant = 0,0,0,0
        flag_share_ix = False
        
        for (i,j) in enumerate(maxs_index):
            flag_max = True;flag_len=True
            maximo = maxs[cells].loc[j]
            min_izq_index = -1
            min_der_index = 0
            
            # cheks left amplitude for the candidate maxima
            while ((maximo - mins[mins_index<j][cells].iloc[min_izq_index]) < amp_th) or ((maxs[cells].loc[j] - mins[mins_index<j][cells].iloc[min_izq_index])/(j- mins[mins_index<j].index.get_level_values(1)[min_izq_index])< pendiente_th) :
                min_izq_index = min_izq_index -1
                if abs(min_izq_index) > len(mins[mins_index<j][cells]):
                    min_izq_index = min_izq_index + 1
                    flag_len = False
                    break

            # cheks right amplitude for the candidate maxima   
            while ((maximo -mins[mins_index>j][cells].iloc[min_der_index]) < amp_th) or ((maxs[cells].loc[j] -  mins[mins_index>j][cells].iloc[min_der_index])/(mins[mins_index>j].index.get_level_values(1)[min_der_index]-j )< pendiente_th):
                min_der_index = min_der_index + 1
                if (abs(min_der_index) >= len(mins[mins_index>j][cells])):
                    min_der_index = min_der_index - 1
                    flag_len=False
                    break
                
            if flag_len:
                # Cheks than between two minima there is a maxima
                idx_mins_izq = mins_index[mins_index<j][min_izq_index] -data.index.get_level_values(1)[0]
                
                amp_der = (maximo -  mins[mins_index>j][cells].iloc[min_der_index])
                amp_izq = (maximo - mins[mins_index<j][cells].iloc[min_izq_index])
                amp = (amp_izq+amp_der)/2
                

                if  idx_mins_der_ant > idx_mins_izq:
                    if amp > amp_ant:
                        # delete previous maxima and go on
                        aux_amp[idx_ant] = None
                        aux_dt[idx_ant] = None
                        aux_max[idx_ant] = None
                        if flag_share_ix==False: 
                            aux_min[idx_mins_izq_ant] = None# puede ser que se comparta
                        aux_min[idx_mins_der_ant] = None     
                        aux_raise[idx_ant] = None
                        aux_fall[idx_ant] = None
         
               
                    if amp <= amp_ant:
                        # go to next maxima
                        flag_max = False
                        
         
                dt_izq = j- mins[mins_index<j].index.get_level_values(1)[min_izq_index]
                dt_der = mins[mins_index>j].index.get_level_values(1)[min_der_index]-j
                    
                
                if flag_max: 
                        amp = (amp_izq + amp_der)/2
                        dt =  dt_der + dt_izq
                        raise_ = amp_izq/dt_izq
                        fall_ = amp_der/dt_der
                        
                        idx = j-data.index.get_level_values(1)[0]
                        aux_amp[idx] = amp
                        aux_dt[idx] = dt
                        aux_max[idx] = maximo
                        idx_mins_izq = mins_index[mins_index<j][min_izq_index] -data.index.get_level_values(1)[0]
                        aux_min[idx_mins_izq] = mins[mins_index<j][cells].iloc[min_izq_index]
                        
                        idx_mins_der = mins_index[mins_index>j][min_der_index] -data.index.get_level_values(1)[0]
                        aux_min[idx_mins_der] = mins[mins_index>j][cells].iloc[min_der_index]            
             
                        aux_raise[idx] = raise_
                        aux_fall[idx] = fall_
                        
                        if idx_mins_der_ant == idx_mins_izq: 
                            flag_share_ix = True
                        else:
                            flag_share_ix = False
                        
                        amp_ant = amp
                        idx_ant = idx
                        idx_mins_der_ant = idx_mins_der
                        idx_mins_izq_ant = idx_mins_izq
                    
        main_aux_amp = main_aux_amp + (aux_amp)    
        main_aux_dt = main_aux_dt+ (aux_dt)
        main_aux_max = main_aux_max + (aux_max)
        main_aux_min = main_aux_min + (aux_min)
        main_aux_raise = main_aux_raise + (aux_raise)
        main_aux_fall = main_aux_fall + (aux_fall)
        
        
    df['amp_peaks']= main_aux_amp
    return((df.amp_peaks.groupby(level=0).count()/df.FRAME.groupby(level=0).count()).mean())

#%%

def main(amp_th,pendiente_th,nprocs,df_dic): 
    df_name = [c for c in df_dic ]
    start = time()
    pool = mp.Pool(processes=nprocs) 
    params = itertools.product(df_name,amp_th,pendiente_th)
    
    func = partial(peak_finding,df_dic)
    results= pool.map(func,params) 

    pool.close()
    pool.join() 
    print(np.round(time()-start, 3))  
    
    values = list(itertools.product(df_name,amp_th,pendiente_th))
    df_results_WT = pd.DataFrame(index=amp_th,columns = pendiente_th)
    df_results_CN = pd.DataFrame(index=amp_th,columns = pendiente_th)

    for (v,res) in zip(values,results):
        if v[0] == 'an_ESC_FAX':
            df_results_WT.loc[v[1],v[2]] = res
        if v[0] == 'an_ESC_FAX_PD03':
            df_results_CN.loc[v[1],v[2]] = res
            
    filename_WT = 'df_results_WT_Norm_.pkl'
    filename_CN = 'df_results_CN_Norm_.pkl'


    df_results_WT.to_pickle(filename_WT)
#    outfile_WT= open(filename_WT,'wb')
#    pickle.dump(df_results_WT,outfile_WT)
#    outfile_WT.close()
    
    df_results_CN.to_pickle(filename_CN)

#    outfile_CN= open(filename_CN,'wb')
#    pickle.dump(df_results_CN,outfile_CN)
#    outfile_CN.close()
    
    return(0)

#%%

    
filename= 'ESC.pkl'
amp_th= np.linspace(0,6000,60)
pendiente_th =  np.linspace(0,500,60)
nprocs = 80

infile = open(filename,'rb')
df_dic = pd.read_pickle(infile)
infile.close()

main(amp_th,pendiente_th,nprocs,df_dic)

