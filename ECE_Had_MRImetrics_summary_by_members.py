#### #### #### ... #### #### #### 


#all packages
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import matplotlib


#metrics='patcor_4clus'
#metrics='freq_ORDasREF_4clus'
metrics='sig_2to6clus_500nrsamp'
#metrics='varopt_2to6clus'
#metrics='varoptNUM_4clus'
#metrics='varoptDEN_4clus'
#metrics='PCscal1' #lag1 autocorrelation

#rootpath='/home/strommen/Python/NAO/Regimes/IrenePaper/GRL_3model_comparison/Outputs_from_WRtool'
rootpath='/network/home/aopp/strommen/Python/NAO/Regimes/IrenePaper/GRL_3model_comparison/Outputs_from_WRtool'

ece_txt = '%s/OUTPUT_ECEARTH31/OUTPUT_ECE/OUTtxt/' % rootpath
had_txt = '%s/OUTPUT_HadGEM3GA3/OUTPUT_HadGEM/OUTtxt/' % rootpath
mri_txt = '%s/OUTPUT_MRIAGCM3/OUTPUT_MRI/OUTtxt/' % rootpath
#mri_txt_hi = '%s/OUTPUT_MRIAGCM3/OUTPUT_MRI959/OUTtxt/' % rootpath
mri_txt_hi = mri_txt

ERA_ECE=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2008_4pcs.txt'.format(ece_txt,metrics))
NCEP_ECE=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2008_4pcs.txt'.format(ece_txt,metrics))
ERA_HadGEM=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1986-2011_4pcs.txt'.format(had_txt,metrics))
NCEP_HadGEM=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1986-2011_4pcs.txt'.format(had_txt,metrics))
ERA_MRI=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_txt,metrics))
NCEP_MRI=np.loadtxt('{0}{1}_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_txt,metrics))
ECEl=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs.txt'.format(ece_txt,metrics))
ECEh=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs.txt'.format(ece_txt,metrics))
HadGEMl=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs.txt'.format(had_txt,metrics))
HadGEMh=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs.txt'.format(had_txt,metrics))
MRIl=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_txt,metrics))
MRIh=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs.txt'.format(mri_txt,metrics))
ECEl0=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs_0.txt'.format(ece_txt,metrics))
ECEh0=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs_0.txt'.format(ece_txt,metrics))
HadGEMl0=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs_0.txt'.format(had_txt,metrics))
HadGEMh0=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs_0.txt'.format(had_txt,metrics))
MRIl0=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs_0.txt'.format(mri_txt,metrics))
MRIh0=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs_0.txt'.format(mri_txt,metrics))
ECEl1=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs_1.txt'.format(ece_txt,metrics))
ECEh1=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs_1.txt'.format(ece_txt,metrics))
HadGEMl1=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs_1.txt'.format(had_txt,metrics))
HadGEMh1=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs_1.txt'.format(had_txt,metrics))
MRIl1=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs_1.txt'.format(mri_txt,metrics))
MRIh1=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs_1.txt'.format(mri_txt,metrics))
ECEl2=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs_2.txt'.format(ece_txt,metrics))
ECEh2=np.loadtxt('{0}{1}_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs_2.txt'.format(ece_txt,metrics))
HadGEMl2=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs_2.txt'.format(had_txt,metrics))
HadGEMh2=np.loadtxt('{0}{1}_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs_2.txt'.format(had_txt,metrics))
MRIl2=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs_2.txt'.format(mri_txt,metrics))
MRIh2=np.loadtxt('{0}{1}_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs_2.txt'.format(mri_txt,metrics))


name_ref='ERAInterim'
numclus=4
name_clus=['NAO+','Blocking','Atlantic Ridge','NAO-']

if metrics=='patcor_4clus':
    a=np.arange(0,1.2,0.2)
    name_metrics='Pattern correlation'
    name_ylabel='pattern correlation'
    ylim=[0.1,1.05]
    tit='{0} relative to {1}'.format(name_metrics,name_ref)
elif metrics=='freq_ORDasREF_4clus':
    a=np.arange(16,37,4)
    name_metrics='Frequency of occurrence'
    name_ylabel='frequency (%)'
    ylim=[16,36]
    tit='{0}'.format(name_metrics)
elif metrics=='sig_2to6clus_500nrsamp':
    a=np.arange(25,105,5)
    a1 = np.arange(84,101,2)
    a2 = a
    name_metrics='Significance of cluster partition'
    name_ylabel='significance (%)'
    ylim=[25,105]
    ylim1 = [85,101]
    ylim2 = ylim
    tit='{0}'.format(name_metrics)
elif metrics=='varopt_2to6clus':
    a=np.arange(0.6,0.88,0.02)
    name_metrics='Optimal ratio'
    name_ylabel='optimal ratio'
    ylim=[0.6,0.86]
    tit='{0}'.format(name_metrics)
elif metrics=='varoptNUM_4clus':
    a=np.arange(1e+09,3.1e+09,0.5e+09)
    name_metrics='Inter-clusters variance'
    name_ylabel='NUMERATOR'
    ylim=[1e+09,3e+09]
    tit='{0}'.format(name_metrics)
elif metrics=='varoptDEN_4clus':
    a=np.arange(1e+09,3.1e+09,0.5e+09)
    name_metrics='Intra-cluster variance'
    name_ylabel='DENOMINATOR'
    ylim=[1e+09,3e+09]
    tit='{0}'.format(name_metrics)

exp_names=['NCEP(1979-2008)','NCEP(1986-2011)','NCEP(1979-2010)','ECEall low','ECEall high','HadGEMall low','HadGEMall high','MRIall low','MRIall high','ECE0','ECE1','ECE2','ECE0','ECE1','ECE2','HadGEM0','HadGEM1','HadGEM2','HadGEM0','HadGEM1','HadGEM2','MRI0','MRI1','MRI2','MRI0','MRI1','MRI2']
exp_names1=['NCEP(1979-2008)','NCEP(1986-2011)','NCEP(1979-2010)','ECEall low','ECEall high','HadGEMall low','HadGEMall high','MRIall low','MRIall high']
exp_names2=['NCEP(1979-2008)','NCEP(1986-2011)','NCEP(1979-2010)','ECE0','ECE1','ECE2','ECE0','ECE1','ECE2','HadGEM0','HadGEM1','HadGEM2','HadGEM0','HadGEM1','HadGEM2','MRI0','MRI1','MRI2','MRI0','MRI1','MRI2']
empty_names=['','','','','','','','','','','','','','','','','','','','','','','','','','','']


print('METRICS={0}'.format(metrics))

if metrics=='patcor_4clus' or metrics=='freq_ORDasREF_4clus' or metrics=='varoptNUM_4clus' or metrics=='varoptDEN_4clus':
    # FOR HadGEM model the order of the clusters is NAO+,Blocking,NAO-,Altantic Ridge so we need to switch the last two elements:
    HadGEMl2[2],HadGEMl2[3]=HadGEMl2[3],HadGEMl2[2]
    HadGEMh2[2],HadGEMh2[3]=HadGEMh2[3],HadGEMh2[2]

    for clus in range(numclus):

        mtoplot_NCEP_ECE    =0*[None]+[NCEP_ECE[clus]]+26*[None]
        mtoplot_NCEP_HadGEM =1*[None]+[NCEP_HadGEM[clus]]+25*[None]
        mtoplot_NCEP_MRI    =2*[None]+[NCEP_MRI[clus]]+24*[None]
        mtoplot_ECEl        =3*[None]+[ECEl[clus]]+23*[None]
        mtoplot_ECEh        =4*[None]+[ECEh[clus]]+22*[None]
        mtoplot_HadGEMl =5*[None]+[HadGEMl[clus]]+21*[None]
        mtoplot_HadGEMh =6*[None]+[HadGEMh[clus]]+20*[None]
        mtoplot_MRIl        =7*[None]+[MRIl[clus]]+19*[None]
        mtoplot_MRIh        =8*[None]+[MRIh[clus]]+18*[None]

        mtoplot_ECEl0        =9*[None]+[ECEl0[clus]]+17*[None]
        mtoplot_ECEl1        =10*[None]+[ECEl1[clus]]+16*[None]
        mtoplot_ECEl2        =11*[None]+[ECEl2[clus]]+15*[None]
        mtoplot_ECEh0        =12*[None]+[ECEh0[clus]]+14*[None]
        mtoplot_ECEh1        =13*[None]+[ECEh1[clus]]+13*[None]
        mtoplot_ECEh2        =14*[None]+[ECEh2[clus]]+12*[None]
        mtoplot_HadGEMl0 =15*[None]+[HadGEMl0[clus]]+11*[None]
        mtoplot_HadGEMl1 =16*[None]+[HadGEMl1[clus]]+10*[None]
        mtoplot_HadGEMl2 =17*[None]+[HadGEMl2[clus]]+9*[None]
        mtoplot_HadGEMh0 =18*[None]+[HadGEMh0[clus]]+8*[None]
        mtoplot_HadGEMh1 =19*[None]+[HadGEMh1[clus]]+7*[None]
        mtoplot_HadGEMh2 =20*[None]+[HadGEMh2[clus]]+6*[None]
        mtoplot_MRIl0        =21*[None]+[MRIl0[clus]]+5*[None]
        mtoplot_MRIl1        =22*[None]+[MRIl1[clus]]+4*[None]
        mtoplot_MRIl2        =23*[None]+[MRIl2[clus]]+3*[None]
        mtoplot_MRIh0        =24*[None]+[MRIh0[clus]]+2*[None]
        mtoplot_MRIh1        =25*[None]+[MRIh1[clus]]+1*[None]
        mtoplot_MRIh2        =26*[None]+[MRIh2[clus]]+0*[None]

        mean_ECEl=(ECEl0[clus]+ECEl1[clus]+ECEl2[clus])/3
        mean_ECEh=(ECEh0[clus]+ECEh1[clus]+ECEh2[clus])/3
        mean_HadGEMl=(HadGEMl0[clus]+HadGEMl1[clus]+HadGEMl2[clus])/3
        mean_HadGEMh=(HadGEMh0[clus]+HadGEMh1[clus]+HadGEMh2[clus])/3
        mean_MRIl=(MRIl0[clus]+MRIl1[clus]+MRIl2[clus])/3
        mean_MRIh=(MRIh0[clus]+MRIh1[clus]+MRIh2[clus])/3
        std_ECEl=np.std([ECEl0[clus],ECEl1[clus],ECEl2[clus]])
        std_ECEh=np.std([ECEh0[clus],ECEh1[clus],ECEh2[clus]])
        std_HadGEMl=np.std([HadGEMl0[clus],HadGEMl1[clus],HadGEMl2[clus]])
        std_HadGEMh=np.std([HadGEMh0[clus],HadGEMh1[clus],HadGEMh2[clus]])
        std_MRIl=np.std([MRIl0[clus],MRIl1[clus],MRIl2[clus]])
        std_MRIh=np.std([MRIh0[clus],MRIh1[clus],MRIh2[clus]])


        print "NUMCLUS=%s" % clus
        print "ECE_LOW_MEAN = %.3f" % mean_ECEl
        print "ECE_LOW_STD = %.3f" % (2*std_ECEl)
        print "ECE_LOW_ALL = %.3f" % (ECEl[clus])     
   
        print "HAD_LOW_MEAN = %.3f" % mean_HadGEMl
        print "HAD_LOW_STD = %.3f" % (2*std_HadGEMl)
        print "HAD_LOW_ALL = %.3f" % (HadGEMl[clus])

        print "MRI_LOW_MEAN = %.3f" % mean_MRIl
        print "MRI_LOW_STD = %.3f" % (2*std_MRIl)
        print "MRI_LOW_ALL = %.3f" % (MRIl[clus])
        print 40*'-'

        print "ECE_HI_MEAN = %.3f" % mean_ECEh
        print "ECE_HI_STD = %.3f" % (2*std_ECEh)
        print "ECE_HI_ALL = %.3f" % (ECEh[clus])

        print "HAD_HI_MEAN = %.3f" % mean_HadGEMh
        print "HAD_HI_STD = %.3f" % (2*std_HadGEMh)
        print "HAD_HI_ALL = %.3f" % (HadGEMh[clus])

        print "MRI_HI_MEAN = %.3f" % mean_MRIh
        print "MRI_HI_STD = %.3f" % (2*std_MRIh)
        print "MRI_HI_ALL = %.3f" % (MRIh[clus])  


        print 50*'-'


        #fig = plt.figure(figsize=(14,10))
        #ax = fig.add_subplot(2, 2, clus+1)
        ax = plt.subplot(2, 2, clus+1)
        x=list(range(len(exp_names)))
        plt.scatter(x,mtoplot_NCEP_ECE, c='rosybrown', s=200, marker='o')
        plt.scatter(x,mtoplot_NCEP_HadGEM, c='gray', s=200, marker='x')
        plt.scatter(x,mtoplot_NCEP_MRI, c='k', s=200, marker='^')
        plt.scatter(x,mtoplot_ECEl, c='b', s=200, marker='o')
        plt.scatter(x,mtoplot_ECEh, c='r', s=200, marker='o')
        plt.scatter(x,mtoplot_HadGEMl, c='b', s=200, marker='x')
        plt.scatter(x,mtoplot_HadGEMh, c='r', s=200, marker='x')
        plt.scatter(x,mtoplot_MRIl, c='b', s=200, marker='^')
        plt.scatter(x,mtoplot_MRIh, c='r', s=200, marker='^')

        plt.scatter(x,mtoplot_ECEl0, c='b', s=100, marker='o')
        plt.scatter(x,mtoplot_ECEh0, c='r', s=100, marker='o')
        plt.scatter(x,mtoplot_HadGEMl0, c='b', s=100, marker='x')
        plt.scatter(x,mtoplot_HadGEMh0, c='r', s=100, marker='x')
        plt.scatter(x,mtoplot_MRIl0, c='b', s=100, marker='^')
        plt.scatter(x,mtoplot_MRIh0, c='r', s=100, marker='^')

        plt.scatter(x,mtoplot_ECEl1, c='b', s=100, marker='o')
        plt.scatter(x,mtoplot_ECEh1, c='r', s=100, marker='o')
        plt.scatter(x,mtoplot_HadGEMl1, c='b', s=100, marker='x')
        plt.scatter(x,mtoplot_HadGEMh1, c='r', s=100, marker='x')
        plt.scatter(x,mtoplot_MRIl1, c='b', s=100, marker='^')
        plt.scatter(x,mtoplot_MRIh1, c='r', s=100, marker='^')

        plt.scatter(x,mtoplot_ECEl2, c='b', s=100, marker='o')
        plt.scatter(x,mtoplot_ECEh2, c='r', s=100, marker='o')
        plt.scatter(x,mtoplot_HadGEMl2, c='b', s=100, marker='x')
        plt.scatter(x,mtoplot_HadGEMh2, c='r', s=100, marker='x')
        plt.scatter(x,mtoplot_MRIl2, c='b', s=100, marker='^')
        plt.scatter(x,mtoplot_MRIh2, c='r', s=100, marker='^')

        plt.axhline(y=ERA_ECE[clus], c='rosybrown')#,label='ERA(1979-2008)')
        plt.axhline(y=ERA_HadGEM[clus], c='gray')#,label='ERA(1986-2011)')
        plt.axhline(y=ERA_MRI[clus], c='k',label='ERA(1979-2010)')

        plt.plot((9, 11), (mean_ECEl,mean_ECEl), c='b',linewidth=2,label='low')
        plt.plot((12, 14), (mean_ECEh,mean_ECEh), c='r',linewidth=2,label='high')
        plt.plot((15, 17), (mean_HadGEMl,mean_HadGEMl), c='b',linewidth=2)
        plt.plot((18, 20), (mean_HadGEMh,mean_HadGEMh), c='r',linewidth=2)
        plt.plot((21, 23), (mean_MRIl,mean_MRIl), c='b',linewidth=2)
        plt.plot((24, 26), (mean_MRIh,mean_MRIh), c='r',linewidth=2)
        if clus==0:
            legend =ax.legend(loc='best', prop={'size':20})
        suptit=name_clus[clus]
        ax.add_patch(Rectangle((9,mean_ECEl-std_ECEl),2,2*std_ECEl, facecolor='b', alpha=0.2))
        ax.add_patch(Rectangle((12,mean_ECEh-std_ECEh),2,2*std_ECEh, facecolor='r', alpha=0.2))
        ax.add_patch(Rectangle((15,mean_HadGEMl-std_HadGEMl),2,2*std_HadGEMl, facecolor='b', alpha=0.2))
        ax.add_patch(Rectangle((18,mean_HadGEMh-std_HadGEMh),2,2*std_HadGEMh, facecolor='r', alpha=0.2))
        ax.add_patch(Rectangle((21,mean_MRIl-std_MRIl),2,2*std_MRIl, facecolor='b', alpha=0.2))
        ax.add_patch(Rectangle((24,mean_MRIh-std_MRIh),2,2*std_MRIh, facecolor='r', alpha=0.2))
        ax.set_yticks(a)
        ax.set_ylim(ylim)
        ax.set_title(suptit, fontsize=40, fontweight='bold')
        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(40)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(40)
        if clus==0:# or clus==2:
            plt.ylabel('{0}'.format(name_ylabel), fontsize=30)
        if clus==2 or clus==3:
            plt.xticks(x, exp_names, rotation=85, fontsize=15)
        else:
            plt.xticks(x, empty_names)
        #plt.grid()
    plt.suptitle(tit, fontsize=45, fontweight='bold')
    plt.subplots_adjust(top=0.85)
    top    = 0.89  # the top of the subplots of the figure
    bottom = 0.12    # the bottom of the subplots of the figure
    left   = 0.07    # the left side of the subplots of the figure
    right  = 0.98  # the right side of the subplots of the figure
    hspace = 0.16   # the amount of height reserved for white space between subplots
    wspace = 0.14    # the amount of width reserved for blank space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    plt.show(block=True)
    #plt.savefig('%s/test.png' % rootpath)
    #plt.close()

elif metrics=='sig_2to6clus_500nrsamp' or metrics=='varopt_2to6clus' or metrics=='PCscal1':

    mtoplot1_NCEP_ECE    =0*[None]+[NCEP_ECE[2]]+8*[None]
    mtoplot1_NCEP_HadGEM =1*[None]+[NCEP_HadGEM[2]]+7*[None]
    mtoplot1_NCEP_MRI    =2*[None]+[NCEP_MRI[2]]+6*[None]
    mtoplot1_ECEl        =3*[None]+[ECEl[2]]+5*[None]
    mtoplot1_ECEh        =4*[None]+[ECEh[2]]+4*[None]
    mtoplot1_HadGEMl =5*[None]+[HadGEMl[2]]+3*[None]
    mtoplot1_HadGEMh =6*[None]+[HadGEMh[2]]+2*[None]
    mtoplot1_MRIl        =7*[None]+[MRIl[2]]+1*[None]
    mtoplot1_MRIh        =8*[None]+[MRIh[2]]+0*[None]

    mtoplot2_NCEP_ECE    =0*[None]+[NCEP_ECE[2]]+20*[None]
    mtoplot2_NCEP_HadGEM =1*[None]+[NCEP_HadGEM[2]]+19*[None]
    mtoplot2_NCEP_MRI    =2*[None]+[NCEP_MRI[2]]+18*[None]
    mtoplot2_ECEl0        =3*[None]+[ECEl0[2]]+17*[None]
    mtoplot2_ECEl1        =4*[None]+[ECEl1[2]]+16*[None]
    mtoplot2_ECEl2        =5*[None]+[ECEl2[2]]+15*[None]
    mtoplot2_ECEh0        =6*[None]+[ECEh0[2]]+14*[None]
    mtoplot2_ECEh1        =7*[None]+[ECEh1[2]]+13*[None]
    mtoplot2_ECEh2        =8*[None]+[ECEh2[2]]+12*[None]
    mtoplot2_HadGEMl0 =9*[None]+[HadGEMl0[2]]+11*[None]
    mtoplot2_HadGEMl1 =10*[None]+[HadGEMl1[2]]+10*[None]
    mtoplot2_HadGEMl2 =11*[None]+[HadGEMl2[2]]+9*[None]
    mtoplot2_HadGEMh0 =12*[None]+[HadGEMh0[2]]+8*[None]
    mtoplot2_HadGEMh1 =13*[None]+[HadGEMh1[2]]+7*[None]
    mtoplot2_HadGEMh2 =14*[None]+[HadGEMh2[2]]+6*[None]
    mtoplot2_MRIl0        =15*[None]+[MRIl0[2]]+5*[None]
    mtoplot2_MRIl1        =16*[None]+[MRIl1[2]]+4*[None]
    mtoplot2_MRIl2        =17*[None]+[MRIl2[2]]+3*[None]
    mtoplot2_MRIh0        =18*[None]+[MRIh0[2]]+2*[None]
    mtoplot2_MRIh1        =19*[None]+[MRIh1[2]]+1*[None]
    mtoplot2_MRIh2        =20*[None]+[MRIh2[2]]+0*[None]

    mean_ECEl=(ECEl0[2]+ECEl1[2]+ECEl2[2])/3
    mean_ECEh=(ECEh0[2]+ECEh1[2]+ECEh2[2])/3
    mean_HadGEMl=(HadGEMl0[2]+HadGEMl1[2]+HadGEMl2[2])/3
    mean_HadGEMh=(HadGEMh0[2]+HadGEMh1[2]+HadGEMh2[2])/3
    mean_MRIl=(MRIl0[2]+MRIl1[2]+MRIl2[2])/3
    mean_MRIh=(MRIh0[2]+MRIh1[2]+MRIh2[2])/3
    std_ECEl=np.std([ECEl0[2],ECEl1[2],ECEl2[2]])
    std_ECEh=np.std([ECEh0[2],ECEh1[2],ECEh2[2]])
    std_HadGEMl=np.std([HadGEMl0[2],HadGEMl1[2],HadGEMl2[2]])
    std_HadGEMh=np.std([HadGEMh0[2],HadGEMh1[2],HadGEMh2[2]])
    std_MRIl=np.std([MRIl0[2],MRIl1[2],MRIl2[2]])
    std_MRIh=np.std([MRIh0[2],MRIh1[2],MRIh2[2]])

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    
    #ax = gca()
    x1=list(range(len(exp_names1)))
    plt.scatter(x1,mtoplot1_NCEP_ECE, c='k', s=200, marker='o')
    plt.scatter(x1,mtoplot1_NCEP_HadGEM, c='k', s=200, marker='x', linewidth=3)
    plt.scatter(x1,mtoplot1_NCEP_MRI, c='k', s=200, marker='^')
    plt.scatter(x1,mtoplot1_ECEl, c='b', s=200, marker='o')
    plt.scatter(x1,mtoplot1_ECEh, c='r', s=200, marker='o')
    plt.scatter(x1,mtoplot1_HadGEMl, c='b', s=200, marker='x', linewidth=3)
    plt.scatter(x1,mtoplot1_HadGEMh, c='r', s=200, marker='x', linewidth=3)
    plt.scatter(x1,mtoplot1_MRIl, c='b', s=200, marker='^')
    plt.scatter(x1,mtoplot1_MRIh, c='r', s=200, marker='^')

    plt.axhline(y=ERA_MRI[2], c='k', linewidth=2.0, label='ERA(1979-2010)')
    legend =ax1.legend(loc='upper left', prop={'size':20})

    ax1.set_yticks(a1)
    ax1.set_ylim(ylim1)
    for tickx in ax1.xaxis.get_major_ticks():
        tickx.label.set_fontsize(40)
    for ticky in ax1.yaxis.get_major_ticks():
        ticky.label.set_fontsize(40)
    plt.ylabel('{0}'.format(name_ylabel), fontsize=30)
    plt.xticks(x1, exp_names1, rotation=85, fontsize=20)
    plt.title('(a)', fontsize=30, loc='left')
    
    ax2 = fig.add_subplot(122)
    x2=list(range(len(exp_names2)))

    plt.scatter(x2,mtoplot2_NCEP_ECE, c='k', s=200, marker='o')
    plt.scatter(x2,mtoplot2_NCEP_HadGEM, c='k', s=200, marker='x', linewidth=3)
    plt.scatter(x2,mtoplot2_NCEP_MRI, c='k', s=200, marker='^')
    plt.scatter(x2,mtoplot2_ECEl0, c='b', s=100, marker='o')
    plt.scatter(x2,mtoplot2_ECEh0, c='r', s=100, marker='o')
    plt.scatter(x2,mtoplot2_HadGEMl0, c='b', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_HadGEMh0, c='r', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_MRIl0, c='b', s=100, marker='^')
    plt.scatter(x2,mtoplot2_MRIh0, c='r', s=100, marker='^')

    print "LOWRES"
    print [ECEl0[2], ECEl1[2], ECEl2[2], HadGEMl0[2], HadGEMl1[2], HadGEMl2[2], MRIl0[2], MRIl1[2], MRIl2[2]]
    print "HIRES"
    print [ECEh0[2], ECEh1[2], ECEh2[2], HadGEMh0[2], HadGEMh1[2], HadGEMh2[2], MRIh0[2], MRIh1[2], MRIh2[2]]
   
    plt.scatter(x2,mtoplot2_ECEl1, c='b', s=100, marker='o')
    plt.scatter(x2,mtoplot2_ECEh1, c='r', s=100, marker='o')
    plt.scatter(x2,mtoplot2_HadGEMl1, c='b', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_HadGEMh1, c='r', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_MRIl1, c='b', s=100, marker='^')
    plt.scatter(x2,mtoplot2_MRIh1, c='r', s=100, marker='^')
    
    plt.scatter(x2,mtoplot2_ECEl2, c='b', s=100, marker='o')
    plt.scatter(x2,mtoplot2_ECEh2, c='r', s=100, marker='o')
    plt.scatter(x2,mtoplot2_HadGEMl2, c='b', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_HadGEMh2, c='r', s=100, marker='x', linewidth=2)
    plt.scatter(x2,mtoplot2_MRIl2, c='b', s=100, marker='^')
    plt.scatter(x2,mtoplot2_MRIh2, c='r', s=100, marker='^')

    #plt.axhline(y=ERA_ECE[2], c='rosybrown')#,label='ERA(1979-2008)')
    #plt.axhline(y=ERA_HadGEM[2], c='gray')#,label='ERA(1986-2011)')
    plt.axhline(y=ERA_MRI[2], c='k', linewidth=2.0, label='ERA(1979-2010)')
    
    plt.plot((3, 5), (mean_ECEl,mean_ECEl), c='b',linewidth=2,label='low')
    plt.plot((6, 8), (mean_ECEh,mean_ECEh), c='r',linewidth=2,label='high')
    plt.plot((9, 11), (mean_HadGEMl,mean_HadGEMl), c='b',linewidth=2)
    plt.plot((12, 14), (mean_HadGEMh,mean_HadGEMh), c='r',linewidth=2)
    plt.plot((15, 17), (mean_MRIl,mean_MRIl), c='b',linewidth=2)
    plt.plot((18, 20), (mean_MRIh,mean_MRIh), c='r',linewidth=2)

    legend =ax2.legend(loc='best', prop={'size':20})
    plt.title('(b)', fontsize=30, loc='left')
    
    ax2.add_patch(Rectangle((3,mean_ECEl-std_ECEl),2,2*std_ECEl, facecolor='b', alpha=0.2))
    ax2.add_patch(Rectangle((6,mean_ECEh-std_ECEh),2,2*std_ECEh, facecolor='r', alpha=0.2))
    ax2.add_patch(Rectangle((9,mean_HadGEMl-std_HadGEMl),2,2*std_HadGEMl, facecolor='b', alpha=0.2))
    ax2.add_patch(Rectangle((12,mean_HadGEMh-std_HadGEMh),2,2*std_HadGEMh, facecolor='r', alpha=0.2))
    ax2.add_patch(Rectangle((15,mean_MRIl-std_MRIl),2,2*std_MRIl, facecolor='b', alpha=0.2))
    ax2.add_patch(Rectangle((18,mean_MRIh-std_MRIh),2,2*std_MRIh, facecolor='r', alpha=0.2))
    ax2.set_yticks(a2)
    ax2.set_ylim(ylim)
    for tickx in ax2.xaxis.get_major_ticks():
        tickx.label.set_fontsize(40)
    for ticky in ax2.yaxis.get_major_ticks():
        ticky.label.set_fontsize(40)
    #plt.ylabel('{0}'.format(name_ylabel), fontsize=30)
    plt.xticks(x2, exp_names2, rotation=85, fontsize=20)
    #plt.grid()
    plt.suptitle(tit, fontsize=45, fontweight='bold')
    top    = 0.90  # the top of the subplots of the figure
    bottom = 0.20    # the bottom of the subplots of the figure
    left   = 0.07    # the left side of the subplots of the figure
    right  = 0.98  # the right side of the subplots of the figure
    hspace = 0.16   # the amount of height reserved for white space between subplots
    wspace = 0.14    # the amount of width reserved for blank space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    #plt.savefig('%s/test.png' % rootpath)
    plt.show(block=True)
    #plt.close()


elif metrics=='sig_2to6clus_500nrsampXXX' or metrics=='varopt_2to6clusXXX' or metrics=='PCscal1XXX':

    sys.exit()
    mtoplot_NCEP_ECE    =0*[None]+[NCEP_ECE[2]]+26*[None]
    mtoplot_NCEP_HadGEM =1*[None]+[NCEP_HadGEM[2]]+25*[None]
    mtoplot_NCEP_MRI    =2*[None]+[NCEP_MRI[2]]+24*[None]
    mtoplot_ECEl        =3*[None]+[ECEl[2]]+23*[None]
    mtoplot_ECEh        =4*[None]+[ECEh[2]]+22*[None]
    mtoplot_HadGEMl =5*[None]+[HadGEMl[2]]+21*[None]
    mtoplot_HadGEMh =6*[None]+[HadGEMh[2]]+20*[None]
    mtoplot_MRIl        =7*[None]+[MRIl[2]]+19*[None]
    mtoplot_MRIh        =8*[None]+[MRIh[2]]+18*[None]

    mtoplot_ECEl0        =9*[None]+[ECEl0[2]]+17*[None]
    mtoplot_ECEl1        =10*[None]+[ECEl1[2]]+16*[None]
    mtoplot_ECEl2        =11*[None]+[ECEl2[2]]+15*[None]
    mtoplot_ECEh0        =12*[None]+[ECEh0[2]]+14*[None]
    mtoplot_ECEh1        =13*[None]+[ECEh1[2]]+13*[None]
    mtoplot_ECEh2        =14*[None]+[ECEh2[2]]+12*[None]
    mtoplot_HadGEMl0 =15*[None]+[HadGEMl0[2]]+11*[None]
    mtoplot_HadGEMl1 =16*[None]+[HadGEMl1[2]]+10*[None]
    mtoplot_HadGEMl2 =17*[None]+[HadGEMl2[2]]+9*[None]
    mtoplot_HadGEMh0 =18*[None]+[HadGEMh0[2]]+8*[None]
    mtoplot_HadGEMh1 =19*[None]+[HadGEMh1[2]]+7*[None]
    mtoplot_HadGEMh2 =20*[None]+[HadGEMh2[2]]+6*[None]
    mtoplot_MRIl0        =21*[None]+[MRIl0[2]]+5*[None]
    mtoplot_MRIl1        =22*[None]+[MRIl1[2]]+4*[None]
    mtoplot_MRIl2        =23*[None]+[MRIl2[2]]+3*[None]
    mtoplot_MRIh0        =24*[None]+[MRIh0[2]]+2*[None]
    mtoplot_MRIh1        =25*[None]+[MRIh1[2]]+1*[None]
    mtoplot_MRIh2        =26*[None]+[MRIh2[2]]+0*[None]

    mean_ECEl=(ECEl0[2]+ECEl1[2]+ECEl2[2])/3
    mean_ECEh=(ECEh0[2]+ECEh1[2]+ECEh2[2])/3
    mean_HadGEMl=(HadGEMl0[2]+HadGEMl1[2]+HadGEMl2[2])/3
    mean_HadGEMh=(HadGEMh0[2]+HadGEMh1[2]+HadGEMh2[2])/3
    mean_MRIl=(MRIl0[2]+MRIl1[2]+MRIl2[2])/3
    mean_MRIh=(MRIh0[2]+MRIh1[2]+MRIh2[2])/3
    std_ECEl=np.std([ECEl0[2],ECEl1[2],ECEl2[2]])
    std_ECEh=np.std([ECEh0[2],ECEh1[2],ECEh2[2]])
    std_HadGEMl=np.std([HadGEMl0[2],HadGEMl1[2],HadGEMl2[2]])
    std_HadGEMh=np.std([HadGEMh0[2],HadGEMh1[2],HadGEMh2[2]])
    std_MRIl=np.std([MRIl0[2],MRIl1[2],MRIl2[2]])
    std_MRIh=np.std([MRIh0[2],MRIh1[2],MRIh2[2]])

    
    ax = gca()
    x=list(range(len(exp_names)))
    plt.scatter(x,mtoplot_NCEP_ECE, c='rosybrown', s=200, marker='o')
    plt.scatter(x,mtoplot_NCEP_HadGEM, c='gray', s=200, marker='x')
    plt.scatter(x,mtoplot_NCEP_MRI, c='k', s=200, marker='^')
    plt.scatter(x,mtoplot_ECEl, c='b', s=200, marker='o')
    plt.scatter(x,mtoplot_ECEh, c='r', s=200, marker='o')
    plt.scatter(x,mtoplot_HadGEMl, c='b', s=200, marker='x')
    plt.scatter(x,mtoplot_HadGEMh, c='r', s=200, marker='x')
    plt.scatter(x,mtoplot_MRIl, c='b', s=200, marker='^')
    plt.scatter(x,mtoplot_MRIh, c='r', s=200, marker='^')

    

    plt.scatter(x,mtoplot_ECEl0, c='b', s=100, marker='o')
    plt.scatter(x,mtoplot_ECEh0, c='r', s=100, marker='o')
    plt.scatter(x,mtoplot_HadGEMl0, c='b', s=100, marker='x')
    plt.scatter(x,mtoplot_HadGEMh0, c='r', s=100, marker='x')
    plt.scatter(x,mtoplot_MRIl0, c='b', s=100, marker='^')
    plt.scatter(x,mtoplot_MRIh0, c='r', s=100, marker='^')

    plt.scatter(x,mtoplot_ECEl1, c='b', s=100, marker='o')
    plt.scatter(x,mtoplot_ECEh1, c='r', s=100, marker='o')
    plt.scatter(x,mtoplot_HadGEMl1, c='b', s=100, marker='x')
    plt.scatter(x,mtoplot_HadGEMh1, c='r', s=100, marker='x')
    plt.scatter(x,mtoplot_MRIl1, c='b', s=100, marker='^')
    plt.scatter(x,mtoplot_MRIh1, c='r', s=100, marker='^')
    
    plt.scatter(x,mtoplot_ECEl2, c='b', s=100, marker='o')
    plt.scatter(x,mtoplot_ECEh2, c='r', s=100, marker='o')
    plt.scatter(x,mtoplot_HadGEMl2, c='b', s=100, marker='x')
    plt.scatter(x,mtoplot_HadGEMh2, c='r', s=100, marker='x')
    plt.scatter(x,mtoplot_MRIl2, c='b', s=100, marker='^')
    plt.scatter(x,mtoplot_MRIh2, c='r', s=100, marker='^')

    plt.axhline(y=ERA_ECE[2], c='rosybrown')#,label='ERA(1979-2008)')
    plt.axhline(y=ERA_HadGEM[2], c='gray')#,label='ERA(1986-2011)')
    plt.axhline(y=ERA_MRI[2], c='k',label='ERA(1979-2010)')
    
    plt.plot((9, 11), (mean_ECEl,mean_ECEl), c='b',linewidth=2,label='low')
    plt.plot((12, 14), (mean_ECEh,mean_ECEh), c='r',linewidth=2,label='high')
    plt.plot((15, 17), (mean_HadGEMl,mean_HadGEMl), c='b',linewidth=2)
    plt.plot((18, 20), (mean_HadGEMh,mean_HadGEMh), c='r',linewidth=2)
    plt.plot((21, 23), (mean_MRIl,mean_MRIl), c='b',linewidth=2)
    plt.plot((24, 26), (mean_MRIh,mean_MRIh), c='r',linewidth=2)

    legend =ax.legend(loc='best', prop={'size':20})

    ax.add_patch(Rectangle((9,mean_ECEl-std_ECEl),2,2*std_ECEl, facecolor='b', alpha=0.2))
    ax.add_patch(Rectangle((12,mean_ECEh-std_ECEh),2,2*std_ECEh, facecolor='r', alpha=0.2))
    ax.add_patch(Rectangle((15,mean_HadGEMl-std_HadGEMl),2,2*std_HadGEMl, facecolor='b', alpha=0.2))
    ax.add_patch(Rectangle((18,mean_HadGEMh-std_HadGEMh),2,2*std_HadGEMh, facecolor='r', alpha=0.2))
    ax.add_patch(Rectangle((21,mean_MRIl-std_MRIl),2,2*std_MRIl, facecolor='b', alpha=0.2))
    ax.add_patch(Rectangle((24,mean_MRIh-std_MRIh),2,2*std_MRIh, facecolor='r', alpha=0.2))
    ax.set_yticks(a)
    ax.set_ylim(ylim)
    for tickx in ax.xaxis.get_major_ticks():
        tickx.label.set_fontsize(40)
    for ticky in ax.yaxis.get_major_ticks():
        ticky.label.set_fontsize(40)
    plt.ylabel('{0}'.format(name_ylabel), fontsize=30)
    plt.xticks(x, exp_names, rotation=85, fontsize=25)
    #plt.grid()
    plt.title(tit, fontsize=45, fontweight='bold')
    top    = 0.89  # the top of the subplots of the figure
    bottom = 0.12    # the bottom of the subplots of the figure
    left   = 0.07    # the left side of the subplots of the figure
    right  = 0.98  # the right side of the subplots of the figure
    hspace = 0.16   # the amount of height reserved for white space between subplots
    wspace = 0.14    # the amount of width reserved for blank space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    #plt.savefig('%s/test.png' % rootpath)
    plt.show(block=True)
    #plt.close() 
 
