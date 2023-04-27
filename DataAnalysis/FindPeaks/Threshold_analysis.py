import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def find_peaks_dens_th(peaks_conc,traze_value,Thresholds,order,test_conditions):

    #colors = sns.color_palette('colorblind')
    #colors = sns.color_palette(sns.dark_palette("#2ecc71", 30, reverse=True))
    colors = sns.cubehelix_palette(500, start=.30, rot=-.40)
    colors = [colors[i] for i in [490,250]]
    peaks_dens_th=[]
    keys = peaks_conc.keys()

    aux_color = []
    for i, key in enumerate(keys):
        for condition in test_conditions:
            if key == condition:
                aux_color.append(colors[i])
    colors = aux_color

    for condition in test_conditions:
        peaks_conc_df=peaks_conc[condition]
        aux_peaks = []
        for th in Thresholds:
            col = str(traze_value) + '_' + str(th) + 'DP_O' + str(order) #hay otra forma de normalizar!
            aux_peaks.append(peaks_conc_df[col].sum())

        peaks_dens_th.append(aux_peaks)

    return(peaks_dens_th,colors)

def plot_threshold_analysis(peaks_conc_dt,traze_value_dt,Thresholds,th_parameter,order,test_conditions,save=False,**kwargs):

    plt.rcdefaults()
    fig, axs = plt.subplots(1,len(peaks_conc_dt), sharex=True, sharey='row', figsize=(8.27, 11.69/4))

    conc_labels = {'an_WT_ESL': 'Serum + LIF', 'an_0ng': 'KO - 0ng of FGF4'}

    for i,peaks_conc in enumerate(peaks_conc_dt):

        peaks_dens_th, colors = find_peaks_dens_th(peaks_conc, traze_value_dt[i], Thresholds, order, test_conditions)
        peaks_dens_th = [j / np.max(j) for j in peaks_dens_th] #normalizacion

        obja=[]
        for j,peaks_dens_th_condition in enumerate(peaks_dens_th):
            #axs[i].set_title(traze_value_dt[i])
            obja.append(axs[i].plot(Thresholds, peaks_dens_th_condition, linestyle="-", marker='o', markersize=4, color=colors[j],
                       alpha=1, label=conc_labels[test_conditions[j]]))
            axs[i].grid(0); #axs[i].legend(loc="upper right", fontsize=6)

            if j == 0:
                th_index = (np.where(peaks_dens_th_condition < th_parameter))[0][0] #fijarse bien que estoy elijiendo
                obj0 = axs[i].plot(Thresholds[th_index],peaks_dens_th_condition[th_index], marker='^', markersize=12, color='red',
                        alpha=1, label=r'$\delta =$ '+ str(th_parameter))
                print('Threshold = '+str(Thresholds[th_index]))

        #ratio = [x / y for x, y in zip(peaks_dens_th[0], peaks_dens_th[1])]
        difference = [y - x for x, y in zip(peaks_dens_th[0], peaks_dens_th[1])]

        max_difference_index = difference.index(max(difference))
        #over_5_ratio = next(obj for obj in ratio if obj <= 0.05)
        #over_5_ratio_index= ratio.index(over_5_ratio)

        #obj = axs[1,i].plot(Thresholds, ratio, linestyle=":", marker='o', markersize=6, color='blue',
        #              alpha=0.5, label='ratio')
        #obj1 = axs.plot(Thresholds[over_5_ratio_index], over_5_ratio, marker='^', markersize=8, color='black',
        #                alpha=1, label='< 0.05 ')

        #ax1=axs[1,i].twinx()
        obj2 = axs[i].plot(Thresholds, difference, linestyle=":", markersize=2, color='gray',
                        alpha=1, label='difference')

        #obj3 = axs[i].plot(Thresholds[max_difference_index], max(difference), marker='v', markersize=4, color='gray',
        #                alpha=1, label='max difference')

        objs = obja[0]+obja[1] #+obj0
        labels = [obj.get_label() for obj in objs]
        axs[1].legend(objs, labels, loc="upper right", fontsize=8,frameon=False);
        axs[i].set_xlim([0, 3000]);
        axs[i].set_ylim([0,1])
        axs[i].set_yticks(np.linspace(0,1,6));
        if i ==0 : axs[0].set_yticklabels(np.round(np.linspace(0,1,6),1), fontsize=8)
        axs[i].set_xticks(np.arange(0,3500,500) ); axs[i].set_xticklabels([i/100 for i in np.arange(0,3500,500) ], fontsize=8)

        if i==0 : axs[i].set_xlabel('threshold (a.u.)',fontsize=8);axs[i].set_ylabel('normalized peak density',fontsize=8)
        axs[i].xaxis.set_label_coords(1, - 0.15);

    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9, wspace=0.15, hspace=0.0)
    #description = kwargs.get('description', " ")
    #plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()

    if save:
        save_folder = kwargs.get('save_folder', None)
        save_name = kwargs.get('save_name', None)
        fig.savefig(save_folder + str(save_name) + '.pdf', format='pdf')

def plot_threshold(peaks_conc_dt,traze_value_dt,Thresholds,order,test_conditions,**kwargs):

    plt.rcdefaults()
    fig, axs = plt.subplots(len(peaks_conc_dt), sharex=True, sharey=True, figsize=(8.27, 11.69))

    for i,peaks_conc in enumerate(peaks_conc_dt):
        ax=axs[i]; ax.grid(0)

        peaks_dens_th, colors = find_peaks_dens_th(peaks_conc, traze_value_dt[i], Thresholds, order, test_conditions)

        for j,peaks_dens_th_condition in enumerate(peaks_dens_th):
            ax.plot(Thresholds, peaks_dens_th_condition, linestyle=":", marker='o', markersize=6, color=colors[j],
                       alpha=0.5, label=test_conditions[j])

        ax.legend(loc="upper right", fontsize=8)

    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
    description = kwargs.get('description', " ")
    plt.figtext(0.1, 0.05, 'Figure: ' + description)
    plt.show()
 #falta save y ponerle cual es cual a cada titulo

def plot_ratio(peaks_conc_dt,traze_value_dt,Thresholds,order,test_conditions,**kwargs):
    if len(test_conditions) == 2:

        plt.rcdefaults()
        fig, axs = plt.subplots(len(peaks_conc_dt), sharex=True, sharey=True, figsize=(8.27, 11.69))

        for i,peaks_conc in enumerate(peaks_conc_dt):
            ax=axs[i]; ax.grid(0)
            ax1 = ax.twinx(); ax1.grid(0)

            peaks_dens_th, colors = find_peaks_dens_th(peaks_conc, traze_value_dt[i], Thresholds, order, test_conditions)
            ratio = [x/y for x, y in zip(peaks_dens_th[0], peaks_dens_th[1])]
            difference = [y-x for x, y in zip(peaks_dens_th[0], peaks_dens_th[1])]

            obj = ax.plot(Thresholds,ratio, linestyle=":", marker='o', markersize=6, color='blue',
                       alpha=0.5,label='ratio')
            obj1 = ax1.plot(Thresholds,difference, linestyle=":", marker='o', markersize=6, color='red',
                       alpha=0.5,label='difference')

            objs = obj + obj1
            labels = [obj.get_label() for obj in objs]
            ax.legend(objs, labels,loc="upper right", fontsize=8)

        fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.07, hspace=0.07)
        description = kwargs.get('description', " ")
        plt.figtext(0.1, 0.05, 'Figure: ' + description)
        plt.show()

    else:
        print('ERROR')
