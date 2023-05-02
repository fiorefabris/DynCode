# =============================================================================
# #Quantifiers (running finding peaks is not necesary)
# =============================================================================


values = [conc,filtmin_conc,buterfilt_conc]
conc.traze_value = 'MEAN_INTENSITY'
for value in values:
    value.set_local_quantifiers(30,1)
    value.set_local_quantifiers(15,1)


#%%
    
CN = 'an_0ng'
CP = 'an_WT_ESL'
quantifiers = ['DT_MEAN_INTENSITY','DT_MEAN_INTENSITY_N15_LOCAL_MEAN','DT_MEAN_INTENSITY_N15_LOCAL_VAR', 
                              'DT_MEAN_INTENSITY_N15_CV2','DT_MEAN_INTENSITY_N15_PE','DT_MEAN_INTENSITY_N15_F',
                              'DT_MEAN_INTENSITY_N30_LOCAL_MEAN','DT_MEAN_INTENSITY_N30_LOCAL_VAR', 
                              'DT_MEAN_INTENSITY_N30_CV2','DT_MEAN_INTENSITY_N30_PE','DT_MEAN_INTENSITY_N30_F']

#####

filtmin_CN = filtmin_conc.data[CN][quantifiers]
filtmin_CP = filtmin_conc.data[CP][quantifiers]

quantifiers = ['MEAN_INTENSITY_BP','MEAN_INTENSITY_BP_N15_LOCAL_MEAN','MEAN_INTENSITY_BP_N15_LOCAL_VAR', 
                              'MEAN_INTENSITY_BP_N15_CV2','MEAN_INTENSITY_BP_N15_PE','MEAN_INTENSITY_BP_N15_F',
                              'MEAN_INTENSITY_BP_N30_LOCAL_MEAN','MEAN_INTENSITY_BP_N30_LOCAL_VAR', 
                              'MEAN_INTENSITY_BP_N30_CV2','MEAN_INTENSITY_BP_N30_PE','MEAN_INTENSITY_BP_N30_F']


buterfilt_CN = buterfilt_conc.data[CN][quantifiers]
buterfilt_CP = buterfilt_conc.data[CP][quantifiers]

quantifiers = ['MEAN_INTENSITY','MEAN_INTENSITY_N15_LOCAL_MEAN','MEAN_INTENSITY_N15_LOCAL_VAR', 
                              'MEAN_INTENSITY_N15_CV2','MEAN_INTENSITY_N15_PE','MEAN_INTENSITY_N15_F',
                              'MEAN_INTENSITY_N30_LOCAL_MEAN','MEAN_INTENSITY_N30_LOCAL_VAR', 
                              'MEAN_INTENSITY_N30_CV2','MEAN_INTENSITY_N30_PE','MEAN_INTENSITY_N30_F']

nofilt_CN = conc.data[CN][quantifiers]
nofilt_CP = conc.data[CP][quantifiers]

#####



result_filtmin = pd.concat([filtmin_CN,filtmin_CP], keys=[CN, CP])
result_buterfilt = pd.concat([buterfilt_CN,buterfilt_CP], keys=[CN, CP])
result_nofilt = pd.concat([nofilt_CN,nofilt_CP], keys=[CN, CP])


main_result = pd.concat([result_nofilt,result_filtmin,result_buterfilt],axis=1)

labels = np.array(main_result.index.get_level_values(level = None))
main_result['label'] = labels

main_result.to_pickle("/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Quantifiers/data.pkl")
main_result.to_csv("/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/Data_Analysis/Quantifiers/data.csv")


####
pd.plotting.scatter_matrix(filtmin_CN.replace('NaN', 0),hist_kwds={'bins':20 ,'density': True})
pd.plotting.scatter_matrix(filtmin_CP.replace('NaN', 0),hist_kwds={'bins':20 ,'density': True})
plt.show()

#####

plot = sns.pairplot(main_result,hue = 'label' ,markers=".",dropna = True,hue_order= ['an_WT_ESL','an_0ng'])
plot.set_xticklabels(plot.get_xticklabels(), fontsize=7,rotation = 90)
#plt.savefig('/home/Fiore/Documents/Dyncode/Datasets/2018_E1_28/quetepasacalabaza12.png', format='png')
plt.show()


g = sns.PairGrid(main_result,hue = 'label')
g = g.map_diag(plt.hist)
g = g.map_offdiag(plt.scatter)
g = g.add_legend()
 
g = sns.FacetGrid(main_result.replace('NaN', 0), hue="label")
g.map(plt.hist,'MEAN_INTENSITY_BP',alpha=.7)
g.add_legend();
plt.show()

