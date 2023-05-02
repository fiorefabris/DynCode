'''Time series forecasting module
   Tries to reproduce a stochastic pulsing dynamics with experimental parameters
'''


from scipy.optimize import curve_fit
import pickle5 as pickle #pickle5
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def download_data(filename):
    infile =  open(filename,'rb')
    results = pickle.load(infile)
    infile.close()
    return(results)

def save_data(data, filename):
    outfile= open(filename,'wb')
    pickle.dump(data,outfile,protocol=-1)
    outfile.close()
#%% 
# =============================================================================
# define silence intervals
# =============================================================================
        
def joint_duration(data):
    values_dt_mixed = []
    cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
    cell_t_m = data[data['min_'].notna()].FRAME.values
    for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
        dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
        dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
        mixed_dt = dt_raise + dt_fall
        values_dt_mixed.append(mixed_dt)
    return(values_dt_mixed)

def silence(data):
    return(data['IPI'].dropna().values-joint_duration(data))

def find_dm(df, flag_contour  = False, flag_isolated  = False):
    """ find_dm(df, flag_contour  = False, flag_isolated  = False)
    Detects the silent regions (i.e. consecutive minima) and compute their duration.
    
   ---------------------------------------------------------------------------------------
    INPUT:
        - df: DataFrame 
            dataframe containing time series
        - flag_contour : bool.
            if True, considers infinite contour conditions. Remember that if True, 
            flag_isolated will be True also. 
        - flag_isolated : bool
            if True, also considers silent regions of time series with isolated pulses or 
            without pulses.

    --------------------------------------------------------------------------------------
    OUTPUT:
        - modified df
            generates a new column on df called dm_ with the silent regions duration
    """
    if flag_contour:
        print('flag isolated is True')
        flag_isolated  = True
   
    dm_   =  pd.DataFrame(np.array([np.nan]* len(df)),index = df.index)
    
    for cells, data in df.groupby(level='cell'):        
        values_dt_mixed = joint_duration(data)
 
        if len(data.amp_peaks.dropna().values) < 2:
            #si tiene uno o cero pulsos
            if flag_isolated:
                
                a = data.FRAME.values[0]; b = data.FRAME.values[-1]
                cell_t_m = data[data['min_'].notna()].FRAME.values

                if len(cell_t_m) == 2:
                    dm_.at[cells,a] = cell_t_m[0] - a
                    dm_.at[cells,b] = b - cell_t_m[1]
                    
                elif len(cell_t_m) == 0:
                    if flag_contour:
                        dm_.at[cells,a] = b - a
                    else: 
                        pass
                else:
                    print('ERROR')
        
        if flag_contour:

            a = data.FRAME.values[0]; b = data.FRAME.values[-1]
            cell_t_m = data[data['min_'].notna()].FRAME.values
            
            if len(cell_t_m) > 2:
                #si tiene mas de dos minimos
        
                dm_.at[cells,a] = cell_t_m[0] - a
                dm_.at[cells,b] = b - cell_t_m[-1]

        for i,(ix,IPI) in enumerate(zip(data['IPI'].dropna().index,data['IPI'].dropna().values)):
            dm_.at[ix] =  IPI - values_dt_mixed[i]
    df['dm_']= dm_
    return(0)

#%%
    
def aux_find_dm_moving_boundaries(data,it,inc):
    """ aux_find_dm_moving_boundaries(data,it,inc)
    Detects the silent regions (i.e. consecutive minima) and compute their duration.
    At the boundaries, we assume there is a pulse. The boundaries are mooving boundaries.
    The moving boundaries is inc % frames per interaction
   ---------------------------------------------------------------------------------------
    INPUT:
        - df: DataFrame 
            dataframe containing time series
        - it : int
            number of iteration 
        - inc : int
            the porcentage of length of the time serues that we extend to each size 
            in each iteration
 
    --------------------------------------------------------------------------------------
    OUTPUT:
        - dm_ = list
            list of the the silent regions duration
    """
    D = np.round(len(data.FRAME)*inc/100,0); assert (D > 0) #inc porciento de incremento
    
    a = data.FRAME.values[0] - D * it 
    b = data.FRAME.values[-1] + D * it 
    dm_ = []
    

    values_dt_mixed = joint_duration(data)    
    assert len(data['IPI'].dropna().values) == len(values_dt_mixed)
    cell_t_m = data[data['min_'].notna()].FRAME.values
    
    if len(data.amp_peaks.dropna().values) == 0:
        #si no tengo pulsos
        dm_.append(b-a)
    elif len(data.amp_peaks.dropna().values) > 0:
        dm_.append(cell_t_m[0] - a)
        dm_.append(b - cell_t_m[-1])
    else:
        assert False,'ERROR: pulsos negativos'
    
    for i,(ix,IPI) in enumerate(zip(data['IPI'].dropna().index,data['IPI'].dropna().values)):
        dm_.append(IPI - values_dt_mixed[i])
    
    assert len(dm_) == (len(data.amp_peaks.dropna().values) + 1), dm_
    return(dm_)


#%%
# =============================================================================
# population fitting functions 
# =============================================================================

def linear(x,lambda_):
    '''
    linear(x,lambda_)
    linearized formula from the exponential distribution
    
    ---------------------------------------------------------------------------------------
    INPUT:
        - x: vector
            x vector
        - lambda_: float
            parameter of the exponential distribution

    ---------------------------------------------------------------------------------------
    OUT:            
        - out: vector
            the y values of the linearized formula from the exponential distribution with 
            lambda_ coefficient computed to the x values
            
    ---------------------------------------------------------------------------------------
    SEE ALSO:
        - get_linear_fit
    '''
    return np.log(lambda_) - lambda_ * x


def get_linear_fit(hist,slice_max = False):
    ''' get_linear_fit (hist,slice_ = False)
        gets an exponential fitting by linearizing the input data and fitting it with a linear
        fit. 
    
        Gets a linear fitting from a normalized histogram. Y axis is in log scale. The linear 
        fitting corresponds to the linearized formula from the exponential distribution. Uses
        non-linear least squares to fit the linear function to data.
        
        ---------------------------------------------------------------------------------------
        INPUT:
            - hist: output from matplotlib hist
                x,y fitting values. The histogram has to be normed. For fitting, the function takes
                the  log scale of hist amplitude values. The fitting ignores the -inf values of this
                log scale. 
                
            - slice (optional): int (in minutes units!)
                consideres the x values until the slice value. 
            
        
        ---------------------------------------------------------------------------------------
        OUTPUT:
            - 
        ---------------------------------------------------------------------------------------
        SEE ALSO:
            - linear
    '''
    if slice_max == False:
        (counts,bins,_) = hist
        X = (bins[:-1]+ bins[1:])/2 ;  
        y = np.log(counts) 
        mask = ~ np.isinf(y)
        return curve_fit(linear,X[mask], y[mask])
    
    else:
        (counts,bins,_) = hist
        X = (bins[:-1]+ bins[1:])/2 ;  
        slice_ = X*20/60 < slice_max; X = X[slice_]
        y = np.log(counts[slice_]) ; 
        mask = ~ np.isinf(y)
        return curve_fit(linear,X[mask], y[mask])
#%%
# función para generar las sieries temporales population
def get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name):
    ''' 
    get_population_time_series(conc_data,ipi_conc,exponential,save_folder_save_name)
    generates time series responding to stochastic pulsing dynamics.
    
    For generating this time series, picks random values of exponential. As pulse duration,
    picks random values from the experimental pulse duration (dt_peaks) distribution. All
    choices are with replacement. Time series can not start with a allready started pulse,
    but they can end without finishing a pulse (so the pulse is not ending in this time 
    serie). The total numer of time series generated is the number of experimental time series
    of the data set. The duration of each of these generated time series is the duration of each 
    of these experimental time series. 
    
    ---------------------------------------------------------------------------------------
    INPUT:
        - conc_data: dic
            experimental data
        - ipi_conc : list
            list of labels of the experimental data you want to consider
        - exponential: list
            exponential distribution values. They are going to be the time between pulses of the
            generated data. The sampling of this distribution is with replacement, so no minimal
            length is needed. 
        -save_folder_save_name: str
            path to save the generated dataframe as a pkl file
        
    ---------------------------------------------------------------------------------------
    OUTPUT:
        saved file: DataFrame as a pickle file.
        Silences has 0 amplitude. Beggining and ending points of the pulses (minimum values of the pulse)
        has 0.2 amplitude, and the rest of the pulse has amplitude 1. The columns of this dataframe are
            * 'exponential_IPI' : time series amplitud values 
            * 'dt_peaks' : pulses duration. Its position is at the middle of the pulse
            * 'amp_peaks' : amplitude value of the pulse. It has a 1 at the middle of the pulse and in practice
            is only for knowing where is the peak of the pulse.
            * 'min_': minimum values of the pulse. It has a 0.2 when the pulse is starting and ending.
            * 'FRAME' frame values of the experiment that these time series are inspired on.
            * 'IPI' : interpulse interval. Its position is at the middle of the second of the pair of consecutive pulses. 

    '''
        
    for c in ipi_conc:
        df = conc_data[c]; 
        
        df_aux = np.array([])
        main_aux_amp,main_aux_dt,main_aux_min = np.array([]),np.array([]),np.array([])
        
        
        for cell,data in df.groupby(level= 'cell'):
            time = data.index.get_level_values(1).values; L = len(time)
            
            t = 0; 
            time_series = np.zeros(L)
            aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
            
            while (t < L):
                dm = int(np.round(np.random.choice(exponential),0))# Esto es en frames y antes tenia un multiplicado x 3
                pulse = int(np.random.choice(df.dt_peaks.dropna().values))
                
                min_izq_index = t + dm
                max_index = t+ dm + pulse//2
                min_der_index = t + dm + pulse
                flag_start = True 
                
                if not pulse == (min_der_index-min_izq_index): print('ERROR')
                
                if L >= min_der_index:
                    
                    time_series[min_izq_index:min_der_index] = np.ones(pulse)
                    if flag_start : time_series[min_izq_index] = 0.2
                    time_series[min_der_index-1] = 0.2
                    
                    
                    aux_amp[max_index] = 1
                    aux_dt[max_index] = pulse
                    aux_min[min_izq_index] = 0.2
                    aux_min[min_der_index-1] = 0.2
    
                else:
                    
                    if (min_izq_index + 3) < L :
                        print('finished a pulse after sampling : cell ',cell); 
                        time_series[min_izq_index:] = np.ones(len(time_series[min_izq_index:] ))
                        time_series[min_izq_index] = 1
                        
                        if max_index < L: aux_amp[max_index] = None; aux_dt[max_index] = None
                        else: aux_amp[-1] = None; aux_dt[-1] = None
                        aux_min[min_izq_index] = None
                        aux_min[-1] = None
                        
    #                    if max_index < L: aux_amp[max_index] = 1; aux_dt[max_index] = L - min_izq_index
    #                    else: aux_amp[-1] = 1; aux_dt[-1] = L - min_izq_index
    #                    aux_min[min_izq_index] = 0.2
    #                    aux_min[-1] = 0.2
                        
                t = min_der_index
             
            df_aux = np.concatenate((df_aux,time_series),axis=0)
            
            main_aux_amp = np.concatenate((main_aux_amp,aux_amp),axis=0)
            main_aux_dt = np.concatenate((main_aux_dt,aux_dt),axis=0)
            main_aux_min = np.concatenate((main_aux_min,aux_min),axis=0)
     
    
    df_forecast = pd.DataFrame({'exponential_IPI': df_aux,'dt_peaks':main_aux_dt,'amp_peaks':main_aux_amp,'min_':main_aux_min,'FRAME':conc_data[c].index.get_level_values(1)},index=conc_data[c].index)
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    return(0)
    
# función para generar las sieries temporales single cell

def get_single_cell_time_series(conc_data,ipi_conc,save_folder_save_name,flag_contour = False):
    ''' 
    get_single_cell_time_series(conc_data,ipi_conc,save_folder_save_name)
    generates time series responding to stochastic pulsing dynamics.
    
    For generating this time series, picks random values of exponential silence distribution. 
    The random choice is with reposition. For each cell, it generates an exponential dirstribution 
    of silent intervals by using as the parameter of this exponential distribution the single 
    cells mean value of silent intervals. For eacj cell, every pulse has the same pulse duration,
    estimates as the mean value of the pulse duration of the experimental single cell we are 
    using for taking the parameters in this time series. Time series can not start with a allready started pulse,
    but they can end without finishing a pulse (so the pulse is not ending in this time 
    serie). The total numer of time series generated is the number of experimental time series
    of the data set. The duration of each of these generated time series is the duration of each 
    of these experimental time series. 
    
    ---------------------------------------------------------------------------------------
    INPUT:
        - conc_data: dic
            experimental data
        - ipi_conc : list
            list of labels of the experimental data you want to consider
        -save_folder_save_name: str
            path to save the generated dataframe as a pkl file
        - flag_contour : bool 
            if True, infinite contour conditions are considered
        
    ---------------------------------------------------------------------------------------
    OUTPUT:
        saved file: DataFrame as a pickle file.
        Silences has 0 amplitude. Beggining and ending points of the pulses (minimum values of the pulse)
        has 0.2 amplitude, and the rest of the pulse has amplitude 1. The columns of this dataframe are
            * 'exponential_IPI' : time series amplitud values 
            * 'dt_peaks' : pulses duration. Its position is at the middle of the pulse
            * 'amp_peaks' : amplitude value of the pulse. It has a 1 at the middle of the pulse and in practice
            is only for knowing where is the peak of the pulse.
            * 'min_': minimum values of the pulse. It has a 0.2 when the pulse is starting and ending.
            * 'FRAME' frame values of the experiment that these time series are inspired on.
            * 'IPI' : interpulse interval. Its position is at the middle of the second of the pair of consecutive pulses. 

    '''


    
    lambda_box = []
    
    for c in ipi_conc:
        df = conc_data[c]; find_dm(df,flag_contour,True)
        df_aux = np.array([])
        main_aux_amp,main_aux_dt,main_aux_min = np.array([]),np.array([]),np.array([])
        
        for cell,data in df.groupby(level= 'cell'):
            time = data.index.get_level_values(1).values; L = len(time)
            
            t = 0; time_series = np.zeros(L)
            aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
            
            lambda_ = np.mean(data.dm_.dropna().values)
            lambda_box.append(lambda_)
            if lambda_ > 0:
                #print(1/lambda_,'------------')
                exponential_cell = np.random.exponential(lambda_,1000)
            else: 
                 exponential_cell = np.array([])
                 
            if exponential_cell.size > 0:
                while (t < L):
                    dm = int(np.round(np.random.choice(exponential_cell),0))
                    
                    pulse = np.round(data.dt_peaks.mean(),0)
                    if np.isnan(pulse):
                        pulse = int(np.round(df.dt_peaks.mean(),0))
                    else: 
                        pulse = int(pulse)
                    
                    min_izq_index = t + dm
                    max_index = t+ dm + pulse//2
                    min_der_index = t + dm + pulse
                    flag_start = True 
                    
                    if not pulse == (min_der_index-min_izq_index): print(00000)
                    
                    if L >= min_der_index:
                        
                        time_series[min_izq_index:min_der_index] = np.ones(pulse)
                        if flag_start : time_series[min_izq_index] = 0.2
                        time_series[min_der_index-1] = 0.2
                        
                        
                        aux_amp[max_index] = 1
                        aux_dt[max_index] = pulse
                        aux_min[min_izq_index] = 0.2
                        aux_min[min_der_index-1] = 0.2
        
                    else:
                        
                        if (min_izq_index + 3) < L :
                            print('finished a pulse after sampling : cell ',cell); 
                            time_series[min_izq_index:] = np.ones(len(time_series[min_izq_index:] ))
                            time_series[min_izq_index] = 1
                            
                            if max_index < L: aux_amp[max_index] = None; aux_dt[max_index] = None
                            else: aux_amp[-1] = None; aux_dt[-1] = None
                            aux_min[min_izq_index] = None
                            aux_min[-1] = None
                            
    
                    t = min_der_index
            else:
                pass
    
            df_aux = np.concatenate((df_aux,time_series),axis=0)
            main_aux_amp = np.concatenate((main_aux_amp,aux_amp),axis=0)
            main_aux_dt = np.concatenate((main_aux_dt,aux_dt),axis=0)
            main_aux_min = np.concatenate((main_aux_min,aux_min),axis=0)
            
    df_forecast = pd.DataFrame({'exponential_IPI': df_aux,'dt_peaks':main_aux_dt,'amp_peaks':main_aux_amp,'min_':main_aux_min,'FRAME':conc_data[c].index.get_level_values(1)},index=conc_data[c].index)
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    
#%%
    
def generate_exponential(lambda_):
    if lambda_ > 0:
        exponential_cell = np.random.exponential(lambda_,1000)
    else: 
        exponential_cell = np.array([])
    return(exponential_cell)
    
def get_pulse(data,df):
    pulse = np.round(data.dt_peaks.mean(),0)
    if np.isnan(pulse):
        pulse = int(np.round(df.dt_peaks.mean(),0))
    else: 
        pulse = int(pulse)
    return(pulse)
        
    
def get_single_cell_time_series_moving_boundaries(conc_data,ipi_conc,save_folder_save_name):
 
    lambda_box = []
    inc = 1
    
    for c in ipi_conc:
        df = conc_data[c]; 
        df_aux = np.array([])
        main_aux_amp,main_aux_dt,main_aux_min = np.array([]),np.array([]),np.array([])
        
        dm_cells = []
        for cell,data in df.groupby(level= 'cell'):
            
            n_pulse = np.inf ; it = 0; 
            time = data.index.get_level_values(1).values
            L = len(time)

            while (n_pulse - data.amp_peaks.count()) > 0  :
                #n_pulse: cuenta los pulsos de la serie simulada. Es una función decreciente de los moving bd.
                dm_ = aux_find_dm_moving_boundaries(data,it,inc)
                
                t = 0; time_series = np.zeros(L)
                aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
                lambda_ = np.mean(dm_)
                exponential_cell = generate_exponential(lambda_)
                                
                assert exponential_cell.size > 0, ('ERROR: lambda = ',lambda_)
                
                
                while (t < L):
                    dm = int(np.round(np.random.choice(exponential_cell),0))                        
                    pulse = get_pulse(data,df) #dt
                    
                    min_izq_index = t + dm
                    max_index = t+ dm + pulse//2
                    min_der_index = t + dm + pulse
                    
                    assert pulse == (min_der_index-min_izq_index)
                    
                    if L >= min_der_index:
                        
                        time_series[min_izq_index:min_der_index] = np.ones(pulse)
                        time_series[min_izq_index] = 0.2
                        time_series[min_der_index-1] = 0.2
                        
                        
                        aux_amp[max_index] = 1
                        aux_dt[max_index] = pulse
                        aux_min[min_izq_index] = 0.2
                        aux_min[min_der_index-1] = 0.2
        
                    else:
                        
                        if min_izq_index  < L :
                            #finished a pulse after sampling 
                            #no contamos los pulsos que no terminan!
                            
                            time_series[min_izq_index:] = np.ones(len(time_series[min_izq_index:] ))
                            time_series[min_izq_index] = 0.2
                            
    
                    t = min_der_index
                
                n_pulse_ant = n_pulse
                n_pulse = np.sum([1 for i in aux_amp if i == 1]) ; 
                it = it + 1
                
                
            assert ( (((n_pulse -  data.amp_peaks.count()) < 0) and  ((it-1)==0)) or ((n_pulse -  data.amp_peaks.count()) <= 0)),(it,(n_pulse -  data.amp_peaks.count()),(n_pulse_ant -  data.amp_peaks.count()))
            #la cantidad de pulsos estimados es menor que la que hay en los datos si y solo si estamos en la primera interacción.
            # De lo contrario, los pulsos estimados son la misma cantidad que los datos
            
            
            dm_cells.append(dm_);lambda_box.append(lambda_)
            df_aux = np.concatenate((df_aux,time_series),axis=0)
            main_aux_amp = np.concatenate((main_aux_amp,aux_amp),axis=0)
            main_aux_dt = np.concatenate((main_aux_dt,aux_dt),axis=0)
            main_aux_min = np.concatenate((main_aux_min,aux_min),axis=0)
            
    df_forecast = pd.DataFrame({'exponential_IPI': df_aux,'dt_peaks':main_aux_dt,'amp_peaks':main_aux_amp,'min_':main_aux_min,'FRAME':conc_data[c].index.get_level_values(1)},index=conc_data[c].index)
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    
    return(dm_cells,lambda_box)
    
   


#%%

def get_single_cells_pulses_shuffle(conc_data,ipi_conc,save_folder_save_name,seed):
    np.random.seed(seed)
    for c in ipi_conc:
        df = conc_data[c]; 
        df_aux = np.array([])
        main_aux_amp,main_aux_dt,main_aux_min = np.array([]),np.array([]),np.array([])
        main_aux_dm = np.array([])
        
        for cell,data in df.groupby(level= 'cell'):
            n_pulse = data.amp_peaks.count()
            time = data.index.get_level_values(1).values            
            dt = data.dt_peaks.mean()
            
            if np.isnan(dt):
                dt = 0
                assert n_pulse == 0
            else:
                dt = int(np.round(data.dt_peaks.mean(),0))
                assert n_pulse > 0

            L = len(time)
            time_series = np.zeros(L)
            aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
            min_der_ix = []; min_izq_ix = []
            
            n = 0; k =0
            index = np.arange((time[-1]-dt+1)-time[0])
            time_pulses = set()
            
            while n < n_pulse:
                k = k + 1
                if k > 1e7:
                    print(cell)
                    break
                else:
                    #la idea es elegir de manera random el principio de un pulso. 
                    #Los pulsos no pueden estar cortados ni al principio ni al final
                    t0 = np.random.choice(index) 
                    time_pulses_cand = time[t0:t0+dt] # asi dejamos que el pulso pueda estar al final 
                    
                    if not time_pulses.intersection(time_pulses_cand):
                        time_series[t0:t0+dt] = np.ones(dt)    
                        n = n + 1 
                        time_pulses.update(time_pulses_cand)
                        
                        time_series[t0] = 0.2
                        time_series[t0+dt-1] = 0.2
                        
                        max_index = t0+dt//2
                        aux_amp[max_index] = 1
                        aux_dt[max_index] = dt
                        aux_min[t0] = 0.2
                        aux_min[t0+dt-1] = 0.2
                        min_izq_ix.append(t0)
                        min_der_ix.append(t0+dt-1)
            
            min_der_ix.sort() ; min_izq_ix.sort()
            assert all([i > j for i,j in zip(min_der_ix,min_izq_ix)])
            
            dm_cells= [j-i for i,j in zip(min_izq_ix[0:-1],min_der_ix[1:])]
            main_aux_dm = np.concatenate((main_aux_dm,dm_cells),axis=0)
            df_aux = np.concatenate((df_aux,time_series),axis=0)
            main_aux_amp = np.concatenate((main_aux_amp,aux_amp),axis=0)
            main_aux_dt = np.concatenate((main_aux_dt,aux_dt),axis=0)
            main_aux_min = np.concatenate((main_aux_min,aux_min),axis=0)
            
    df_forecast = pd.DataFrame({'exponential_IPI': df_aux,'dt_peaks':main_aux_dt,'amp_peaks':main_aux_amp,'min_':main_aux_min,'FRAME':conc_data[c].index.get_level_values(1)},index=conc_data[c].index)
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    
    return(0)

    

#%%
def generar_poisson(df,L,save_folder_save_name):
    '''
    Generates a time serie with a poissonean time interval distribution
    toma el dm de el valor medio de df
    toma lambda_ de la poissoneana como la que ajusta a la distribucion de dm de df
    df es conc_data['an_WT_ESL']
    
    '''
    
    #calcular un lambda para usar de parámetro
    find_dm(df)
    bins1 = plt.hist(df['dm_'],bins=np.arange(2,200,6), density=1);
    popt,pcov = get_linear_fit(bins1)
    lambda_ = popt[0]   
    exponential = np.random.exponential(1/lambda_,len(df['dt_peaks']))
    dt = int(np.round(df.dt_peaks.dropna().values.mean(),0))
    #generar las time series condicion inicial    
    
    time = np.arange(L)
    t = 0; 
    time_series = np.zeros(L)
    aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
    
    
    #contruir la serie temporal
    while (t < L):
        dm = int(np.round(np.random.choice(exponential),0))
        pulse = dt
                
        min_izq_index = t + dm
        max_index = t+ dm + pulse//2
        min_der_index = t + dm + pulse
        flag_start = True 
                
        assert pulse == (min_der_index-min_izq_index)
                
        if L >= min_der_index:
                    
            time_series[min_izq_index:min_der_index] = np.ones(pulse)
            if flag_start : time_series[min_izq_index] = 0.2
            time_series[min_der_index-1] = 0.2
                    
                    
            aux_amp[max_index] = 1
            aux_dt[max_index] = pulse
            aux_min[min_izq_index] = 0.2
            aux_min[min_der_index-1] = 0.2
    
        else:
                    
            if (min_izq_index + 3) < L :
                time_series[min_izq_index:] = np.ones(len(time_series[min_izq_index:] ))
                time_series[min_izq_index] = 1
                        
                if max_index < L: aux_amp[max_index] = None; aux_dt[max_index] = None
                else: aux_amp[-1] = None; aux_dt[-1] = None
                aux_min[min_izq_index] = None
                aux_min[-1] = None
                        
        
        t = min_der_index
             
    
    
    df_forecast = pd.DataFrame({'exponential_IPI': time_series,'dt_peaks':aux_dt,'amp_peaks':aux_amp,'min_':aux_min,'FRAME':time})
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    return(lambda_)
#%%
def shuffle_poison(original_data_folder_name,save_folder_save_name):
    
    data = download_data(original_data_folder_name+'.pkl')
    n_pulse = data.amp_peaks.count()
    time = data.index.values            
    dt = int(np.round(data.dt_peaks.mean(),0))

    L = len(time)
    time_series = np.zeros(L)
    aux_amp,aux_dt,aux_min = [None]*L,[None]*L,[None]*L
    min_der_ix = []; min_izq_ix = []
            
    n = 0; k =0
    index = np.arange((time[-1]-dt+1)-time[0])
    time_pulses = set()
            
    while n < n_pulse:
        k = k + 1
        if k > 1e7:
            break
        else:
            #la idea es elegir de manera random el principio de un pulso. 
            #Los pulsos no pueden estar cortados ni al principio ni al final
            t0 = np.random.choice(index) 
            time_pulses_cand = time[t0:t0+dt] # asi dejamos que el pulso pueda estar al final 
                    
            if not time_pulses.intersection(time_pulses_cand):
                time_series[t0:t0+dt] = np.ones(dt)    
                n = n + 1 
                time_pulses.update(time_pulses_cand)
                        
                time_series[t0] = 0.2
                time_series[t0+dt-1] = 0.2
                        
                max_index = t0+dt//2
                aux_amp[max_index] = 1
                aux_dt[max_index] = dt
                aux_min[t0] = 0.2
                aux_min[t0+dt-1] = 0.2
                min_izq_ix.append(t0)
                min_der_ix.append(t0+dt-1)
            
    min_der_ix.sort() ; min_izq_ix.sort()
    assert all([i > j for i,j in zip(min_der_ix,min_izq_ix)])
            
    dm_cells= [j-i for i,j in zip(min_izq_ix[0:-1],min_der_ix[1:])]
            
    df_forecast = pd.DataFrame({'exponential_IPI': time_series,'dt_peaks':aux_dt,'amp_peaks':aux_amp,'min_':aux_min,'FRAME':time})
    
    mask = df_forecast['amp_peaks'].notna()
    df_forecast['IPI'] = df_forecast[mask].FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-
    save_data(df_forecast,save_folder_save_name+'.pkl')
    
    return(dm_cells)

#%%
def get_silent_poisson(df_forecast):    
    dm_   =  pd.DataFrame(np.array([np.nan]* len(df_forecast)),index = df_forecast.index)
    values_dt_mixed = joint_duration(df_forecast) 
    assert len(df_forecast.amp_peaks.dropna().values) > 2
    for i,(ix,IPI) in enumerate(zip(df_forecast['IPI'].dropna().index,df_forecast['IPI'].dropna().values)):
            dm_.at[ix] =  IPI - values_dt_mixed[i]
    return(dm_)