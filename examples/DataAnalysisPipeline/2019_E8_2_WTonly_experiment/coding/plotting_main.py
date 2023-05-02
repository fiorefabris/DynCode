# =============================================================================
# FIGURE 2 & 4
# =============================================================================
from sklearn.linear_model import LinearRegression
import pickle
import pickle5
import numpy as np
import matplotlib.ticker as ticker


def load_file5(filename):
    infile = open(filename,'rb')
    output = pickle5.load(infile)
    infile.close()
    return(output)

def load_file(filename):
    infile = open(filename,'rb')
    output = pickle.load(infile)
    infile.close()
    return(output)
    
    
def save_file(variable,filename):
    outfile= open(filename,'wb')
    pickle.dump(variable,outfile)
    outfile.close()


def silent_ax(ax):
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

def make_square_axes(ax):
    """Make an axes square in screen units.

    Should be called after plotting.
    """
    ax.set_aspect(1 / ax.get_data_ratio())

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
        trend_aux = []; sm_trend_aux = []
        detrended_aux= []; sm_detrended_aux = []
        
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
            
        # fit linear model sm data
            X = data.index.get_level_values(1)
            X = np.reshape(X, (len(X), 1))
            y = data.sm_MEAN_INTENSITY.values
            model = LinearRegression()
            model.fit(X, y)
        # calculate trend
            trend = model.predict(X)
            detrended = [y[i]-trend[i] for i in range(0, len(data.sm_MEAN_INTENSITY))]
            sm_trend_aux.extend(list(trend))
            sm_detrended_aux.extend(list(detrended))
        
        conc[key]['sm_trend'] = sm_trend_aux
        conc[key]['sm_detrended'] = sm_detrended_aux

        conc[key]['trend'] = trend_aux
        conc[key]['detrended'] = detrended_aux
        
        
def set_scale_bars(ax, x_bar=None, y_bar=None, xunits='', yunits='', x_scale=1, y_scale=1, round_x=False,round_y=False):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if x_bar is not None:
        x0, xf = x_bar[0], x_bar[1]
        ax.spines['bottom'].set_bounds(x0, xf)
        ax.set_xticks([(xf + x0) / 2])
        scale_x = (xf - x0) // x_scale
        #if round_x:       ax.set_xticklabels([str(int(scale_x)) + xunits])#, rotation='vertical', va='center')
        #else:         ax.set_xticklabels([str(scale_x) + xunits])#, rotation='vertical', va='center')
        ax.set_xticklabels([])
        ax.xaxis.set_tick_params(length=0)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
    if y_bar is not None:
        y0, yf = y_bar[0], y_bar[1]
        ax.spines['left'].set_bounds(y0, yf)
        ax.set_yticks([(yf + y0) / 2])
        scale = (yf - y0) // y_scale
        if round_y:
            scale = int(scale)
        #ax.set_yticklabels([str(scale) + yunits], rotation='vertical', va='center')
        ax.set_yticklabels([])
        ax.yaxis.set_tick_params(length=0)
    else:
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
        
        
def set_scale_bars_twinx(ax, ax2, x_bar=None, y_bar=None, y_bar_twinx=None, xunits='', yunits='',yunits_twinx = '', x_scale=1, y_scale=1, y_scale_twinx = 1, round_x=False,round_y=False, round_y_twinx = False):
    ax.spines['top'].set_visible(False)

    if x_bar is not None:
        x0, xf = x_bar[0], x_bar[1]
        ax.spines['bottom'].set_bounds(x0, xf)
        ax.set_xticks([(xf + x0) / 2])
        scale_x = (xf - x0) // x_scale
        ax.set_xticklabels([])
        ax.xaxis.set_tick_params(length=0)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
        
    if y_bar is not None:
        y0, yf = y_bar[0], y_bar[1]
        ax.spines['left'].set_bounds(y0, yf)
        ax.set_yticks([(yf + y0) / 2])
        scale = (yf - y0) // y_scale
        if round_y:
            scale = int(scale)
        #ax.set_yticklabels([str(scale) + yunits], rotation='vertical', va='center')
        ax.set_yticklabels([])
        ax.yaxis.set_tick_params(length=0)
    else:
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)
      
    if y_bar_twinx is not None:
        y0, yf = y_bar_twinx[0], y_bar_twinx[1]
        ax2.spines['right'].set_bounds(y0, yf)
        ax2.set_yticks([(yf + y0) / 2])
        scale = (yf - y0) // y_scale_twinx
        if round_y_twinx:
            scale = int(scale)

        ax2.set_yticklabels([])
        ax2.yaxis.set_tick_params(length=0)
    else:
        ax2.spines['right'].set_visible(False)
        ax2.yaxis.set_visible(False)

        
def joint_duration(df):
    values_dt_mixed = []
    for cells, data in df.groupby(level='cell'):
        cell_t_M = data[data['max_'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
    return(values_dt_mixed)



def is_consecutive(df):
    # discrimina los intervalos de consecutividad vs no consecutividad
    y = df['IPI'].dropna().values-joint_duration(df)
    x = joint_duration(df)
    consecutive = np.zeros(len(x))
    for i in range(len(x)):
        if x[i]*0.5 >= y[i]:
            consecutive[i] = 1
    return(consecutive)
 
def sanity_check():
    for c in conc_data.keys():
        df = conc_data[c];df.sort_index(inplace=True);
        for cells, data in df.groupby(level='cell'):
            res = data['IPI'].dropna().values-joint_duration(data)
            for i,j in enumerate(res):
                if j< 0: print(c,cells,i,j); print(data['IPI'].dropna().values,joint_duration(data)); 
                
                
#%% Function for intermidiate step


def peaks(df_dic,params):
    
    """ 
    Intermidiate function for pulse detection histogram noseque. Originalmente Fig. 3 Supp 3 
    Identifies peaks of a time serie & returns the corresponding 
    characterization: amp, dt, max, min, raise, fall, IPI.
    
    -------------------------------------------------------
    INPUTS:
            - params: amp_th and pendiente_th of a corresponding level.
            - df_dic (with identified peaks) 
    -------------------------------------------------------
    OUTPUT:
            - df: a DataFrame with the old and new information.
            
    """
           
    c,amp_th,pendiente_th  = params[0],params[1][1],params[1][0]
    #c,amp_th,pendiente_th  = params[0],params[1],params[2]
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
                

#                if  idx_mins_der_ant > idx_mins_izq:
#                    if amp > amp_ant:
#                        # delete previous maxima and go on
#                        aux_amp[idx_ant] = None
#                        aux_dt[idx_ant] = None
#                        aux_max[idx_ant] = None
#                        if flag_share_ix==False: 
#                            aux_min[idx_mins_izq_ant] = None# puede ser que se comparta
#                        aux_min[idx_mins_der_ant] = None     
#                        aux_raise[idx_ant] = None
#                        aux_fall[idx_ant] = None
#         
#               
#                    if amp <= amp_ant:
#                        # go to next maxima
#                        flag_max = False
                        
         
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

    df['dt_peaks']= main_aux_dt
    df['max_']= main_aux_max
    df['min_']= main_aux_min
    df['raise_']= main_aux_raise
    df['fall_']= main_aux_fall

    mask = df['max_'].notna()
    df['IPI'] = df[mask].groupby(level='cell').FRAME.diff() ### Le resta a un valor su anterior. El primero es NaN-

    return(df)