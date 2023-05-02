from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
from itertools import groupby


def makeColours( vals ):
    colours = np.zeros( (len(vals),3) )
    norm = Normalize( vmin=min(vals), vmax=max(vals) )

    #Can put any colormap you like here.
    colours = [cm.ScalarMappable( norm=norm, cmap='viridis').to_rgba( val ) for val in vals]

    return colours

def count_iterable(iterator):
    return sum(1 for i in iterator)
#%%
    
class consecutive_cumulative:
    def __init__(self,df):
        self.df = df
        
    def joint_duration(self,data):
        values_dt_mixed = []
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
        return(values_dt_mixed)

    def is_consecutive_cell(self,data):
        #te devuelve par de pulsos consecutivos, nivel celula
        y = data['IPI'].dropna().values-self.joint_duration(data) #silence 
        x = self.joint_duration(data) #joint duration
        consecutive = []
        for i in range(len(x)):
            if x[i]*0.5 >= y[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        
    def is_consecutive_box(self):
        df = self.df
        #Te da estadística de pulsos en toda la poblacion PARES DE PULSOS, nivel poblacion
        box = []
        for cells, data in df.groupby(level='cell'):
            box.append(self.is_consecutive_cell(data))
        return(box)
    
    def raise_order_consecutiveness(self,box):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad
        new_box = []
        for consecutive_ in box:
            new_consecutive_ = []
            for i,j in zip(consecutive_[:-1],consecutive_[1:]):
                new_consecutive_.append(i*j)
            new_box.append(new_consecutive_)
        return(new_box)
                
    def get_consecutive_trains_of_pulses(self):
        #te devuelve un array en donde cada elemento tiene la cantidad de trenes de pulsos sobre el total de la población
        df  = self.df
        box = self.is_consecutive_box() #pares de pulsos
        box_plot_consecutive = [df.amp_peaks.count(),np.sum([np.sum(l) for l in box])]
        
        while box_plot_consecutive[-1] > 0:
            box= self.raise_order_consecutiveness(box)
            box_plot_consecutive.append(np.sum([np.sum(l) for l in box])) #cada elemento de la lista new_box(box) es una célula
        return(box_plot_consecutive[:-1]) # el último elemento es el que se va a cero, por eso lo sacamos
    
    def get_consecutive_trains_of_pulses_cells(self):
        #te devuelve un array en donde cada elemento tiene una lista donde en cada 
        # lugarcito tiene la cantidad de trenes de la célula
        df = self.df
        #cada lugarcito es la suma en una celula
        box = self.is_consecutive_box() #pares de pulsos
        box_plot_consecutive = [df.amp_peaks.groupby(level='cell').count().values,[np.sum(l) for l in box]]
        while np.sum(box_plot_consecutive[-1]) > 0:
            box = self.raise_order_consecutiveness(box)
            box_plot_consecutive.append([np.sum(l) for l in box])
            box = box
        return(box_plot_consecutive)
#%%
# =============================================================================
# ESTE MODULO ES PARA COMPUTAR CONSECUTIVIDAD CON LA IDEA DE LUIS
# =============================================================================
''' module for computing consecutive trains of pulses 
    como concepto general, la box tiene en cada lugar una celula
'''
class consecutive_non_cumulative:

    def __init__(self,df):
        self.df = df
        
    def joint_duration(self,data):
        values_dt_mixed = []
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
        return(values_dt_mixed)
        

    def silence(self,data):
        return(data['IPI'].dropna().values-self.joint_duration(data))


    def is_isolated_cell(self,MAX,joint_duration,dm):
        # solo para una celula
        # calcula el numero de picos no consecutivis (isolated)
        # tiene en cuenta las celulas de un solo pulso
        # te devuelve la suma del total de picos menos los que son parte de intervalos consecutivos de pulsos
        consecutive = self.is_consecutive_cell(joint_duration,dm)
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(MAX.count()  - count)
                    
    def is_isolated_box(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos aislados
        box = []
        for cell,data in df.groupby(level='cell'):
            box.append(self.is_isolated_cell(data.amp_peaks,self.joint_duration(data),self.silence(data)))
        return(box)
    

    def is_consecutive_cell(self,joint_duration,dm):
        # para cada celula (data), te devuelve un array con 1 representando pares de pulsos consecutivos y
        # 0 pares de pulsos no consecutivos
        consecutive = []
        for i in range(len(joint_duration)):
            if joint_duration[i]*0.5 >= dm[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        
    
    def is_consecutive_box(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos consecutivos en una lista
        box = []
        for cell,data in df.groupby(level='cell'):
                box.append(self.is_consecutive_cell(self.joint_duration(data),self.silence(data)))
        return(box)
        
    
    def raise_order_consecutiveness(self,box):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad
        new_box = []
        for consecutive_ in box:
            new_consecutive_ = []
            for i,j in zip(consecutive_[:-1],consecutive_[1:]):
                new_consecutive_.append(i*j)
            new_box.append(new_consecutive_)
        return(new_box)
    
    
    def count_consecutive_pulses(self,box):
        #te devuelve una lista.- En cada lugar de la lista está el número de unos aislados en cada célula .
        count_box = []
        for consecutive in box:
            count = 0
            for i, group in groupby(consecutive):
                if i == 1:
                    if count_iterable(group) == 1:
                        count = count + 1
                    else:
                        pass
            count_box.append(count) 
        return(count_box)
    
    
    def get_consecutive_trains_of_pulses(self):
        #te da la suma sobre todo el dataset
        
        isolated = self.is_isolated_box()
        box = self.is_consecutive_box()
        
        box_plot_consecutive = [np.sum(isolated)]
        
        while np.sum([np.sum(l) for l in box]) > 0:
            box_plot_consecutive.append(sum(self.count_consecutive_pulses(box)))
            box = self.raise_order_consecutiveness(box)
            
        return(box_plot_consecutive)
    
    def get_number_of_pulses(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos aislados
        pulses = []
        for cell,data in df.groupby(level='cell'):
            pulses.append(data.amp_peaks.count())
        return(pulses)
 
        
    
    def get_consecutive_trains_of_pulses_cells(self):
        # te da la estadística por célula
        isolated = self.is_isolated_box()
        box = self.is_consecutive_box()
        
        box_plot_consecutive = [isolated]
        
        while np.sum([np.sum(l) for l in box]) > 0:
            box_plot_consecutive.append(self.count_consecutive_pulses(box,False))
            box = self.raise_order_consecutiveness(box)
        return(box_plot_consecutive)
    
    
    
    def count_consecutive_pulses_cell_number(self,MAX,joint_duration,dm):
        # solo para una celula
        # calcula el numero de picos no consecutivis (isolated)
        # tiene en cuenta las celulas de un solo pulso
        # te devuelve la suma del total de picos menos los que son parte de intervalos consecutivos de pulsos
        consecutive = self.is_consecutive_cell(joint_duration,dm)
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(count)
                    
    def count_consecutive_pulses_number(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos consecutivos
        box = []
        for cell,data in df.groupby(level='cell'):
            box.append(self.count_consecutive_pulses_cell_number(data.amp_peaks,self.joint_duration(data),self.silence(data)))
        return(box)
#%%
        
    def count_consecutive_pulses_old(self,box, population = False):
        #cuenta la cantidad de unos aislados que hay en un vector dentro de box
        #si population es True , te devuelve la suma directamente sobretoda la poblacion de box
        count_box = []
        for consecutive_ in box:
            count = 0
            for j,n in enumerate(consecutive_):
                if n == 1: #si tenemos un conjunto de pulsos consecutivos
                    
                    # para el primer conjunto de pulsos
                    #si hay mas conjuntos de pulsos y el siguiente no es consecutivo, contalo
                    if j == 0: 
                        if len(consecutive_) > 1: 
                            if consecutive_[j+1] == 0: 
                                count = count + 1 
                    # para el primer conjunto de pulsos
                    #si hay un solo conjunto de pulsos, contalo
                        else:
                            count = count + 1 #este else no estoy segura! que pasa cuando hay solo un uno en el array?
                   
                    # si el conjunto de pulsos es el ultimo y el anteultimo no es consecutivo, contalo
                    elif j == len(consecutive_)-1 :
                        if consecutive_[j-1] == 0:
                            count = count + 1
                    else:
                    #si no es el primer ni el ultimo conjunto de pulsos, y los de alrededor no tienen pusos, contalo
                        if (consecutive_[j-1] + consecutive_[j+1]) == 0:
                            count = count + 1
                else: 
                    pass
            count_box.append(count)
        if population:
            return(np.sum(count_box)) #te tira el total sobre todo el dataset
        else:
            return(count_box) # te lo tira por celula
    

    
#%%
#Esto es para podes plotear el promedio de muchas consecutividades (tenes muchas datasets, y queres hacerun promediopoblCIONAL 
# de consecutividad  sobe todos esos datasets )
def mean_consecutive_value(trials):
    arr_aux = []
    for j in range(np.max([len(i) for i in trials])):
        aux = []
       
        for i in trials:
            if len(i) > j : 
                aux.append(i[j])
            else:
                aux.append(0)
            
        arr_aux.append(aux)
    
    return(np.array([np.mean(k) for k in arr_aux]),np.array([np.std(k) for k in arr_aux]))

# Esto es para calcular los pulsos aislados y totales en todos los datasets deö  mundo

def first_element_box_old(trials):
    aux = []
    for trial in trials:
        aux.append(trial[0])
    return(aux)

def second_element_box_old(trials):
    aux = []
    for trial in trials:
        aux.append(trial[1])
    return(aux)
    
#Esto es para podes plotear el promedio de muchas consecutividades 
# como lo de arriba, pero solamente tomando el 80% de la poblacion de cada dataset
#seleccionandolos por la actividad (las del medio se quedan las otras no) 
    
def filter_activity_cells(df_consecutive):

    activity = df_consecutive['dt_peaks'].groupby(level='cell').sum() / df_consecutive['FRAME'].groupby(level='cell').count() *  100   
    activity_index = np.argsort(activity.values)[::-1]
    activity = [activity[j] for j in activity_index]

    y = len(activity)*2/10 // 2 #lo que le quiero sacar de cada lado
    activity_filter_cells = []
    
    for n,(ix,cell) in enumerate(zip(activity_index,activity)):
        if (n > y) and (n<= len(activity) - y):
            activity_filter_cells.append(ix)
    
    return(activity_filter_cells)
    
    

def mask_low_variance(box_plot_consecutive_cumulative_trials_mean,box_plot_consecutive_cumulative_trials_std):
    #filtra los valores negativos de la varianza y los reemplaza con un 1e-10 para quese pueda plottear
    arr = (box_plot_consecutive_cumulative_trials_mean-box_plot_consecutive_cumulative_trials_std) 
    new_arr = []
    for i in arr:
        if i > 0:
            new_arr.append(i)
        else:
            new_arr.append(1e-10)
    return(np.array(new_arr))