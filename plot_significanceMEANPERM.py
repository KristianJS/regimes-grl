#### #### #### ... #### #### #### 


#all packages
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

metrics='sig_2to6clus_500nrsamp'
#metrics='patcor_4clus'

#MODEL='ECEARTH31'
#MODEL='HadGEM3-GA3'
MODEL='MRIAGCM32'


rootpath='/home/strommen/Python/NAO/Regimes/IrenePaper/GRL_3model_comparison/Outputs_from_WRtool'

ece_txt = '%s/OUTPUT_ECEARTH31/OUTPUT_ECE/OUTtxt/' % rootpath
had_txt = '%s/OUTPUT_HadGEM3GA3/OUTPUT_HadGEM/OUTtxt/' % rootpath
mri_txt = '%s/OUTPUT_MRIAGCM3/OUTPUT_MRI/OUTtxt/' % rootpath
mri_txt_hi = '%s/OUTPUT_MRIAGCM3/OUTPUT_MRI/OUTtxt/' % rootpath

ece_perm = '%s/OUTPUT_ECEARTH31/OUTPUT_ECEperm/OUTtxt/' % rootpath
had_perm = '%s/OUTPUT_HadGEM3GA3/OUTPUT_HadGEMperm/OUTtxt/' % rootpath
mri_perm = '%s/OUTPUT_MRIAGCM3/OUTPUT_MRIperm/OUTtxt/' % rootpath
mri_perm2 = '%s/OUTPUT_MRIAGCM3/sig_perm/' % rootpath

PATH1='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_ecehadgem/OUTtxt/'
PATH2='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_ecehadgem/OUTtxt/'
PATH3='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_MRI/sig_perm/'
PATH_perm1='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_T255T511base012_perm/OUTtxt/'
PATH_perm2='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_N216N512base012_perm/OUTtxt/'
PATH_perm3='/home/mavilia/WEATHER_REGIMEStool/results/JASMIN_ECEMETcomparison/OUTPUT_MRI/sig_perm/'


def plot_mod(MODEL, axnum):

    if MODEL=='ECEARTH31':
        ##FOR ECEARTH31
        ERAname='{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2008_4pcs'.format(ece_txt,metrics)
        NCEPname='{0}{1}_zg500_day_NCEPNCAR_obs_144x73_1ens_DJF_EAT_1979-2008_4pcs'.format(ece_txt,metrics)
        tit='Significance for EC-Earth3.1'
    elif MODEL=='HadGEM3-GA3':
        #FOR HADGEM3-GA3
        ERAname='{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1986-2011_4pcs'.format(had_txt,metrics)
        NCEPname='{0}{1}_zg500_day_NCEPNCAR_obs_144x73_1ens_DJF_EAT_1986-2011_4pcs'.format(had_txt,metrics)
        tit='Significance for HadGEM3-GA3'
    elif MODEL=='MRIAGCM32':
        #FOR MRI-AGCM3.2
        ERAname='{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2010_4pcs'.format(mri_txt,metrics)
        NCEPname='{0}{1}_zg500_day_NCEPNCAR_obs_144x73_1ens_DJF_EAT_1979-2010_4pcs'.format(mri_txt,metrics)
        tit='Significance for MRI-AGCM3.2'

    if metrics=='patcor_4clus':
        a=np.arange(0.5,1.1,0.1)
        name_metrics='Pattern correlation averaged over the 4 clusters'
        name_ylabel='pattern correlation'
        ylim=[0.5,1.1]
        NCEP=np.mean(np.loadtxt(NCEPname+'.txt'))
        ERA=np.mean(np.loadtxt(ERAname+'.txt'))
    elif metrics=='sig_2to6clus_500nrsamp':
        a=np.arange(25,105,5)
        name_metrics='Significance of cluster partition'
        name_ylabel='significance (%)'
        ylim=[25,105]
        NCEP=np.loadtxt(NCEPname+'.txt')[2]
        ERA=np.loadtxt(ERAname+'.txt')[2]

    data1_perm1=[]
    data2_perm1=[]
    data3_perm1=[]
    data4_perm1=[]
    data5_perm1=[]
    data6_perm1=[]
    for numperm in range(0,3):
        print('\nPermutation {0} of ensemble members by one (single ensemble)'.format(numperm))
        data1_name='{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs_{2}.txt'.format(ece_txt,metrics,numperm)
        data2_name='{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs_{2}.txt'.format(ece_txt,metrics,numperm)
        data3_name='{0}{1}_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs_{2}.txt'.format(had_txt,metrics,numperm)
        data4_name='{0}{1}_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs_{2}.txt'.format(had_txt,metrics,numperm)
        data5_name='{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs_{2}.txt'.format(mri_txt,metrics,numperm)
        data6_name='{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs_{2}.txt'.format(mri_txt,metrics,numperm)
        if metrics=='patcor_4clus':
            data1_perm1.append(np.mean(np.loadtxt(data1_name)))
            data2_perm1.append(np.mean(np.loadtxt(data2_name)))
            data3_perm1.append(np.mean(np.loadtxt(data3_name)))
            data4_perm1.append(np.mean(np.loadtxt(data4_name)))
        elif metrics=='sig_2to6clus_500nrsamp':
            data1_perm1.append(np.loadtxt(data1_name)[2])
            data2_perm1.append(np.loadtxt(data2_name)[2])
            data3_perm1.append(np.loadtxt(data3_name)[2])
            data4_perm1.append(np.loadtxt(data4_name)[2])
            data5_perm1.append(np.loadtxt(data5_name)[2])
            data6_perm1.append(np.loadtxt(data6_name)[2])
    data1_perm1mean=np.mean(data1_perm1)
    data1_perm1std=np.std(data1_perm1)
    data2_perm1mean=np.mean(data2_perm1)
    data2_perm1std=np.std(data2_perm1)
    data3_perm1mean=np.mean(data3_perm1)
    data3_perm1std=np.std(data3_perm1)
    data4_perm1mean=np.mean(data4_perm1)
    data4_perm1std=np.std(data4_perm1)
    data5_perm1mean=np.mean(data5_perm1)
    data5_perm1std=np.std(data5_perm1)
    data6_perm1mean=np.mean(data6_perm1)
    data6_perm1std=np.std(data6_perm1)

    data1_perm2=[]
    data2_perm2=[]
    data3_perm2=[]
    data4_perm2=[]
    data5_perm2=[]
    data6_perm2=[]
    for numperm in range(0,6):
        print('\nPermutation {0} of ensemble members by two'.format(numperm))
        data1_name='{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(ece_perm,metrics,numperm)
        data2_name='{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(ece_perm,metrics,numperm)
        data3_name='{0}{1}_zg500_day_METOFFICE_base_N216L85_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(had_perm,metrics,numperm)
        data4_name='{0}{1}_zg500_day_METOFFICE_base_N512L85_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(had_perm,metrics,numperm)
        data5_name='{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(mri_perm2,metrics,numperm)
        data6_name='{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_2perm{2}_4pcs.txt'.format(mri_perm2,metrics,numperm)
        if metrics=='patcor_4clus':
            data1_perm2.append(np.mean(np.loadtxt(data1_name)))
            data2_perm2.append(np.mean(np.loadtxt(data2_name)))
            data3_perm2.append(np.mean(np.loadtxt(data3_name)))
            data4_perm2.append(np.mean(np.loadtxt(data4_name)))
        elif metrics=='sig_2to6clus_500nrsamp':
            data1_perm2.append(np.loadtxt(data1_name)[2])
            data2_perm2.append(np.loadtxt(data2_name)[2])
            data3_perm2.append(np.loadtxt(data3_name)[2])
            data4_perm2.append(np.loadtxt(data4_name)[2])
            data5_perm2.append(np.loadtxt(data5_name)[2])
            data6_perm2.append(np.loadtxt(data6_name)[2])
    data1_perm2mean=np.mean(data1_perm2)
    data1_perm2std=np.std(data1_perm2)
    data2_perm2mean=np.mean(data2_perm2)
    data2_perm2std=np.std(data2_perm2)
    data3_perm2mean=np.mean(data3_perm2)
    data3_perm2std=np.std(data3_perm2)
    data4_perm2mean=np.mean(data4_perm2)
    data4_perm2std=np.std(data4_perm2)
    data5_perm2mean=np.mean(data5_perm2)
    data5_perm2std=np.std(data5_perm2)
    data6_perm2mean=np.mean(data6_perm2)
    data6_perm2std=np.std(data6_perm2)

    data1_perm3=[]
    data2_perm3=[]
    data3_perm3=[]
    data4_perm3=[]
    data5_perm3=[]
    data6_perm3=[]
    for numperm in range(0,6):
        print('\nPermutation {0} of ensemble members by three'.format(numperm))
        data1_name='{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_3perm{2}_4pcs.txt'.format(ece_perm,metrics,numperm)
        data2_name='{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_3perm{2}_4pcs.txt'.format(ece_perm,metrics,numperm)
        data3_name='{0}{1}_zg500_day_METOFFICE_base_N216L85_3ens_DJF_EAT_3perm{2}_4pcs.txt'.format(had_perm,metrics,numperm)
        data4_name='{0}{1}_zg500_day_METOFFICE_base_N512L85_3ens_DJF_EAT_3perm{2}_4pcs.txt'.format(had_perm,metrics,numperm)
        data5_name='{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_perm2,metrics)
        data6_name='{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_perm2,metrics)
        if metrics=='patcor_4clus':
            data1_perm3.append(np.mean(np.loadtxt(data1_name)))
            data2_perm3.append(np.mean(np.loadtxt(data2_name)))
            data3_perm3.append(np.mean(np.loadtxt(data3_name)))
            data4_perm3.append(np.mean(np.loadtxt(data4_name)))
        elif metrics=='sig_2to6clus_500nrsamp':
            data1_perm3.append(np.loadtxt(data1_name)[2])
            data2_perm3.append(np.loadtxt(data2_name)[2])
            data3_perm3.append(np.loadtxt(data3_name)[2])
            data4_perm3.append(np.loadtxt(data4_name)[2])
            data5_perm3.append(np.loadtxt(data5_name)[2])
            data6_perm3.append(np.loadtxt(data6_name)[2])
    data1_perm3mean=np.mean(data1_perm3)
    data1_perm3std=np.std(data1_perm3)
    data2_perm3mean=np.mean(data2_perm3)
    data2_perm3std=np.std(data2_perm3)
    data3_perm3mean=np.mean(data3_perm3)
    data3_perm3std=np.std(data3_perm3)
    data4_perm3mean=np.mean(data4_perm3)
    data4_perm3std=np.std(data4_perm3)
    data5_perm3mean=np.mean(data5_perm3)
    data5_perm3std=np.std(data5_perm3)
    data6_perm3mean=np.mean(data6_perm3)
    data6_perm3std=np.std(data6_perm3)

    print('METRICS={0}'.format(metrics))

    #===================================================ECEARTH31
    #exp_names=['NCEP','ERA','1ens','2ens','3ens']
    exp_names=['NCEP','1member','2members','3members']

    mtoplot_NCEP=[NCEP,float('NaN'),float('NaN'),float('NaN')]
    #mtoplot_ERA=[float('NaN'),ERA,float('NaN'),float('NaN'),float('NaN')]
    mtoplot_data1=[float('NaN'),data1_perm1mean,data1_perm2mean,data1_perm3mean]
    mtoplot_data2=[float('NaN'),data2_perm1mean,data2_perm2mean,data2_perm3mean]
    mtoplot_data3=[float('NaN'),data3_perm1mean,data3_perm2mean,data3_perm3mean]
    mtoplot_data4=[float('NaN'),data4_perm1mean,data4_perm2mean,data4_perm3mean]
    mtoplot_data5=[float('NaN'),data5_perm1mean,data5_perm2mean,data5_perm3mean]
    mtoplot_data6=[float('NaN'),data6_perm1mean,data6_perm2mean,data6_perm3mean]

    
    ax = fig.add_subplot('22%s' % axnum)
    x=list(range(4))
    plt.scatter(x,mtoplot_NCEP, c='k', s=500, marker='*')      #c
    #plt.scatter(x,mtoplot_ERA, c='grey', s=1000, marker='*')   #coral

    #plt.scatter(x,mtoplot_data1, c='b', s=100, marker='o',label='T255base')
    #plt.plot(x,mtoplot_data1, c='b',label='T255base')
    #plt.scatter(x,mtoplot_data2, c='r', s=100, marker='o',label='T255stoc')
    #plt.plot(x,mtoplot_data2, c='r',label='T255stoc')
    if MODEL=='ECEARTH31':
        plt.errorbar(x,mtoplot_data1, yerr=[0,data1_perm1std,data1_perm2std,data1_perm3std], c='b', fmt='o', markersize='15', ecolor='b',capsize=25, elinewidth=8,label='low res T255')
        plt.errorbar(x,mtoplot_data2, yerr=[0,data2_perm1std,data2_perm2std,data2_perm3std], c='lime', fmt='o', markersize='15', ecolor='lime',capsize=20, elinewidth=6,label='high res T511')
        plt.axhline(y=ERA, c='k',label='ERA-Interim')
    elif MODEL=='HadGEM3-GA3':
        plt.errorbar(x,mtoplot_data3, yerr=[0,data3_perm1std,data3_perm2std,data3_perm3std], c='b', fmt='o', markersize='15', ecolor='b',capsize=25, elinewidth=8,label='low res N216')
        plt.errorbar(x,mtoplot_data4, yerr=[0,data4_perm1std,data4_perm2std,data4_perm3std], c='lime', fmt='o', markersize='15', ecolor='lime',capsize=20, elinewidth=6,label='high res N512')
        plt.axhline(y=ERA, c='k',label='ERA-Interim')
    elif MODEL=='MRIAGCM32':
        plt.errorbar(x,mtoplot_data5, yerr=[0,data5_perm1std,data5_perm2std,data5_perm3std], c='b', fmt='o', markersize='15', ecolor='b',capsize=25, elinewidth=8,label='low res TL95')
        plt.errorbar(x,mtoplot_data6, yerr=[0,data6_perm1std,data6_perm2std,data6_perm3std], c='lime', fmt='o', markersize='15', ecolor='lime',capsize=20, elinewidth=6,label='high res TL319')
        plt.axhline(y=ERA, c='k',label='ERA-Interim')
    ax=gca()
    legend =ax.legend(loc='2', prop={'size':12})
    ax.set_ylim(ylim)
    for tickx in ax.xaxis.get_major_ticks():
        tickx.label.set_fontsize(16)
    for ticky in ax.yaxis.get_major_ticks():
        ticky.label.set_fontsize(16)
    ax.set_yticks(a)
    plt.ylabel('{0}'.format(name_ylabel), fontsize=16)
    plt.xticks(x, exp_names)#, rotation=65)
    #plt.grid()
    #tit='{0}'.format(name_metrics)
    plt.title(tit, fontsize=16, fontweight='bold')
    #top=0.932
    #bottom=0.077
    #left=0.092
    #right=0.95
    #hspace=0.2
    #wspace=0.2
    #plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    #plt.show(block=True)


mods = ['ECEARTH31', 'HadGEM3-GA3', 'MRIAGCM32']
fig = plt.figure(figsize=(12,8))
for num in [0,1,2]:
    plot_mod(mods[num], num+1)


top=0.932
bottom=0.077
left=0.092
right=0.95
hspace=0.2
wspace=0.2
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show(block=True)




