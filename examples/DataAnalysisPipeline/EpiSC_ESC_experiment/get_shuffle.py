import multiprocessing as mp
from functools import partial 
from exponential_dm_main import get_single_cells_pulses_shuffle
import pandas as pd

def download_data(filename):
    return(pd.read_pickle(filename))


def get_single_cells_pulses_shuffle_aux(conc_data,ipi_conc,save_folder,i):
    save_name = 'shuffle_'
    save_folder_save_name = save_folder + save_name + str(i)
    get_single_cells_pulses_shuffle(conc_data,ipi_conc,save_folder_save_name,i)
    return(0)


if __name__ == '__main__':
    pool = mp.Pool(processes= mp.cpu_count())
    tuple_ = range(200)
    
    epi_conc_data = download_data('/home/fiore/DO2020/consecutive_shuffle/epi_conc_data.pkl')
    save_folder = '/home/fiore/DO2020/consecutive_shuffle/epi_data/'
    ipi_conc = ['an_EpiSC_FAX']
    get_single_cells_pulses_shuffle_ = partial(get_single_cells_pulses_shuffle_aux,epi_conc_data,ipi_conc,save_folder)
    pool.map(get_single_cells_pulses_shuffle_,tuple_)
    print('EPI shuffle data end')
    
    epi_conc_data = download_data('/home/fiore/DO2020/consecutive_shuffle/esc_conc_data.pkl')
    save_folder = '/home/fiore/DO2020/consecutive_shuffle/esc_data/'
    ipi_conc = ['an_ESC_FAX']
    get_single_cells_pulses_shuffle_ = partial(get_single_cells_pulses_shuffle_aux,epi_conc_data,ipi_conc,save_folder)
    pool.map(get_single_cells_pulses_shuffle_,tuple_)
    
    pool.close()
    pool.join()
    print('end shuffle!')