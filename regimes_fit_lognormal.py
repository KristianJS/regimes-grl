import sys
sys.path.append('/network/home/aopp/strommen/Python/UsefulFunctions')
sys.path.append('/network/home/aopp/strommen/Python/Cubes')
sys.path.append('/network/home/aopp/strommen/Python/NAO/Regimes')
import timing
import matplotlib
import pylab as pl
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, norm, lognorm
import cartopy.crs as ccrs
import matplotlib.path as mpath
import iris



datapath = './Data'
figpath = './Figures'





###############################################################################
#
#              COMPUTATION OF TRANSITION PROBABILITIES FOR
#              ALL THE FOUR DIFFERENT REGIMES
#
###############################################################################        


#Computing transition probs of some data.
#If logodds=True, returns log-odds of probabilities
#Choose a particular regime (naoplus, naominus, altr or blk)
#skip=True means omit any probabilities that are 0 or undefined (ie. no days in that regime)
#skip=False just sets those probabilities to 1%
def transition_probs(data, dataname, regime, logodds=False, skip=False, save=True):


    print "Computing transition probabilities for %s" % dataname
    print "Regime considered: %s" % regime
    
    trans_naoplus = []
    trans_naominus = []
    trans_atlr = []
    trans_blk = []
    
    N = len(data)
    num_years = int(N/90.)
    
    print "Using %s years worth of data" % num_years
    
    if dataname in ['HADGEM_Low', 'HADGEM_Hi', 'ERAInterim_1986-2011']:
        naoplus_num = 0
        naominus_num = 2
        atlr_num = 3
        blk_num = 1
    else:
        naoplus_num = 0
        naominus_num = 3
        atlr_num = 2
        blk_num = 1
            
    
    if regime == 'NAO+':
        regnum = naoplus_num
        regname = 'naoplus'
    elif regime == 'NAO-':
        regnum = naominus_num
        regname = 'naominus'
    elif regime == 'ATLR':
        regnum = atlr_num
        regname = 'atlr'
    elif regime == 'BLK':
        regnum = blk_num
        regname = 'blk'
    else:
        sys.exit("Error in transition_probs: invalid choice of regime")
        
    
    other_nums = list(set([1,2,3,4]).difference(set([regnum])))    
    
    
    if 'ERAInterim' in dataname:
        ensnum = 1
    elif 'NCEP' in dataname:
        ensnum = 1
    else:
        ensnum = 3
    
    print "%s effective ensemble members" % ensnum
    
    for k in range(num_years):
            
             
        #Subset data to a year
        year = data[k*90*ensnum:(k+1)*90*ensnum]
       
        
        #Count total number of regime days
        total_days = np.count_nonzero(year == regnum)
        
        
        #Estimate transition probs of regime by looking for any days with the associated
        #index and seeing what happens the day after
        trans_naoplus_days = 0
        trans_naominus_days = 0
        trans_atlr_days = 0
        trans_blk_days = 0
        
        
        for day in range(len(year)-1):
            if (year[day] == regnum) and (year[day+1] == atlr_num):
                trans_atlr_days += 1
            elif (year[day] == regnum) and (year[day+1] == naominus_num):
                trans_naominus_days += 1
            elif (year[day] == regnum) and (year[day+1] == naoplus_num):
                trans_naoplus_days += 1
            elif (year[day] == regnum) and (year[day+1] == blk_num):
                trans_blk_days += 1
            else:
                pass       
        
        
        #Compute probability estimate
        if total_days == 0:
            prob_naoplus = 0.01
            prob_naominus = 0.01
            prob_atlr = 0.01
            prob_blk = 0.01
            
            if skip:
                pass
            else:
                trans_naoplus.append(prob_naoplus)
                trans_naominus.append(prob_naominus)
                trans_atlr.append(prob_atlr)
                trans_blk.append(prob_blk)    
            
            
        else:
            if trans_naoplus_days == 0.:
                prob_naoplus = 0.01
                if skip:
                    pass
                else:
                    trans_naoplus.append(prob_naoplus)                  
            else:    
                prob_naoplus = float(trans_naoplus_days) / total_days
                trans_naoplus.append(prob_naoplus)
                
                
            if trans_naominus_days == 0.:
                prob_naominus = 0.01
                if skip:
                    pass
                else:
                    trans_naominus.append(prob_naominus)                  
            else:    
                prob_naominus = float(trans_naominus_days) / total_days
                trans_naominus.append(prob_naominus)  
        
    
            if trans_atlr_days == 0.:
                prob_atlr = 0.01
                if skip:
                    pass
                else:
                    trans_atlr.append(prob_atlr)                  
            else:    
                prob_atlr = float(trans_atlr_days) / total_days
                trans_atlr.append(prob_atlr)  
    
            if trans_blk_days == 0.:
                prob_blk = 0.01
                if skip:
                    pass
                else:
                    trans_blk.append(prob_blk)                  
            else:    
                prob_blk = float(trans_blk_days) / total_days
                trans_blk.append(prob_blk)    
        
        
        
    
    trans_naoplus = np.array(trans_naoplus)
    trans_naominus = np.array(trans_naominus)
    trans_atlr = np.array(trans_atlr)
    trans_blk = np.array(trans_blk)
    
    print "NAO+ shape =" + str(trans_naoplus.shape)
    print "NAO- shape =" + str(trans_naominus.shape)
    print "ATLR shape =" + str(trans_atlr.shape)
    print "BLK shape =" + str(trans_blk.shape)
    
       
    
    if logodds:
        #Ensure trans_naoplus is not 1
        trans_naoplus[trans_naoplus == 1.] = 0.95
        trans_naominus[trans_naominus == 1.] = 0.95
        trans_atlr[trans_atlr == 1.] = 0.95
        trans_blk[trans_blk == 1.] = 0.95
        
        
        naoplus_odds = np.log(trans_naoplus/(1-trans_naoplus))
        naominus_odds = np.log(trans_naominus/(1-trans_naominus))
        atlr_odds = np.log(trans_atlr/(1-trans_atlr))
        blk_odds = np.log(trans_blk/(1-trans_blk))
        
        if save == True:
            np.savetxt('%s/%s_%s_to_naoplus_logodds.txt' % (datapath, dataname, regname), naoplus_odds)
            np.savetxt('%s/%s_%s_to_naominus_logodds.txt' % (datapath, dataname, regname), naominus_odds)
            np.savetxt('%s/%s_%s_to_atlr_logodds.txt' % (datapath, dataname, regname), atlr_odds)
            np.savetxt('%s/%s_%s_to_blk_logodds.txt' % (datapath, dataname, regname), blk_odds)
        else:
            return naoplus_odds, naominus_odds, atlr_odds, blk_odds
    
    else:
        if save == True:
            print "Saving output as textfiles..."
            np.savetxt('%s/Probs/%s_%s_to_naoplus_probs.txt' % (datapath, dataname, regname), trans_naoplus) 
            np.savetxt('%s/Probs/%s_%s_to_naominus_probs.txt' % (datapath, dataname, regname), trans_naominus)
            np.savetxt('%s/Probs/%s_%s_to_atlr_probs.txt' % (datapath, dataname, regname), trans_atlr)
            np.savetxt('%s/Probs/%s_%s_to_blk_probs.txt' % (datapath, dataname, regname), trans_blk)
            #print '%s/%s_%s_to_blk_probs.txt' % (datapath, dataname, regname)
            
            
        else:
            return trans_naoplus, trans_naominus, trans_atlr, trans_blk
            
    print 50*'-'


#Function to load persistence probs for a dataset and a regime choice
def load_pers_probs(dataname, regime):

    if regime == 'NAO+':
        regname = 'naoplus'
    elif regime == 'NAO-':
        regname = 'naominus'
    elif regime == 'ATLR':
        regname = 'atlr'
    elif regime == 'BLK':
        regname = 'blk'
    else:
        sys.exit("Error in transition_probs: invalid choice of regime")

    probs = np.loadtxt('%s/Probs/%s_%s_to_%s_probs.txt' % (datapath, dataname, regname, regname))
    
    return probs
    


###############################################################################
#
#              FITTING OF LOG-NORMAL DISTRIBUTIONS TO
#              PERSISTENCE PROBS AND 
#
############################################################################### 
        
    

def fit_lognormal_pers(dataname, regime):

    pers_probs = load_pers_probs(dataname, regime)
    
    #print dataname
    #print pers_probs.shape
    #print 50*'-'
    
    params = lognorm.fit(1-pers_probs)
    x = np.linspace(0,1,500)
    
    fit = lognorm.pdf(x, params[0], loc=params[1], scale=params[2])
    area = np.trapz(fit, dx=0.01)
    area = 1.

    fit_flipped = np.flipud(fit)/area
    
    return fit_flipped
    
    


def plot_lognormal(dataname, regime):

    if regime == 'NAO+':
        regname = 'naoplus'
    elif regime == 'NAO-':
        regname = 'naominus'
    elif regime == 'ATLR':
        regname = 'atlr'
    elif regime == 'BLK':
        regname = 'blk'
    else:
        sys.exit("Error in transition_probs: invalid choice of regime")

    pers_probs_model_low = fit_lognormal_pers('%s_Low' % dataname, regime)
    pers_probs_model_hi = fit_lognormal_pers('%s_Hi' % dataname, regime)
    
    if 'HADGEM' in dataname:
        firstyear = 1986
        lastyear = 2011
    elif 'MRI' or 'NCEP' in dataname:
        firstyear = 1979
        lastyear = 2010
    else:
        firstyear = 1979
        lastyear = 2008        
    
    pers_probs_era = fit_lognormal_pers('ERAInterim_%s-%s' % (firstyear, lastyear), regime)
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    x = np.linspace(0,1,500)
    
    plt.plot(x, pers_probs_model_low, 'r--', linewidth=3.5, label='%s_Low' % dataname)
    plt.plot(x, pers_probs_model_hi, 'r-', linewidth=3.5, label='%s_Hi' % dataname)
    plt.plot(x, pers_probs_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.xlabel('Persistence probabilities', fontsize=18)
    plt.ylabel('Frequency', fontsize=18)
    plt.legend(loc='upper left', fontsize=20)
    ax.tick_params(labelsize=16)
    plt.title('Reverse log-normal fits of %s persistence' % regime, fontsize=20)
    #plt.show()
    plt.savefig('%s/%s_%s_persistence_pdf.jpg' % (figpath, dataname, regname), dpi=300, format='jpg')
    plt.close()
     


def plot_lognormal_regime(regime):

    if regime == 'NAO+':
        regname = 'naoplus'
    elif regime == 'NAO-':
        regname = 'naominus'
    elif regime == 'ATLR':
        regname = 'atlr'
    elif regime == 'BLK':
        regname = 'blk'
    else:
        sys.exit("Error in transition_probs: invalid choice of regime")

    pers_probs_ece_low = fit_lognormal_pers('ECE_Low', regime)
    pers_probs_ece_hi = fit_lognormal_pers('ECE_Hi', regime)
    pers_probs_ece_era = fit_lognormal_pers('ERAInterim_%s-%s' % (1979, 2008), regime) 
    
    pers_probs_hadgem_low = fit_lognormal_pers('HADGEM_Low', regime)
    pers_probs_hadgem_hi = fit_lognormal_pers('HADGEM_Hi', regime)
    pers_probs_hadgem_era = fit_lognormal_pers('ERAInterim_%s-%s' % (1986, 2011), regime)
    
    pers_probs_mri_low = fit_lognormal_pers('MRI_Low', regime)
    pers_probs_mri_hi = fit_lognormal_pers('MRI_Hi', regime)
    pers_probs_mri_era = fit_lognormal_pers('ERAInterim_%s-%s' % (1979, 2010), regime)
    
    pers_probs_ncep = fit_lognormal_pers('NCEP', regime)
    pers_probs_ncep_era = pers_probs_mri_era
    
    
    #print pers_probs_ece_era.shape
    #print pers_probs_hadgem_era.shape
    #print pers_probs_mri_era.shape
    #print pers_probs_ncep_era.shape
   
    
    fig = plt.figure(figsize=(12,8))
    x = np.linspace(0,1,500)
    
    ax1 = fig.add_subplot(221)
    plt.plot(x, pers_probs_ncep, 'b--', linewidth=3.0, label='NCEP')
    plt.plot(x, pers_probs_mri_era, 'b-', linewidth=3.0, label='ERA-Interim')
    #plt.xlabel('Persistence probabilities', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax1.tick_params(labelsize=12)
    plt.title('(a) ERA-Interim vs NCEP (1979-2010)', fontsize=17)
    
    ax2 = fig.add_subplot(222)
    plt.plot(x, pers_probs_ece_low, 'r--', linewidth=3.0, label='ECE_Low')
    plt.plot(x, pers_probs_ece_hi, 'r-', linewidth=3.0, label='ECE_Hi')
    plt.plot(x, pers_probs_ece_era, 'b-', linewidth=3.0, label='ERA-Interim')
    #plt.xlabel('Persistence probabilities', fontsize=14)
    #plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax2.tick_params(labelsize=12)
    plt.title('(b) ERA-Interim vs EC-Earth (1979-2008)', fontsize=17)
    
    ax3 = fig.add_subplot(223)
    plt.plot(x, pers_probs_hadgem_low, 'r--', linewidth=3.0, label='HADGEM_Low')
    plt.plot(x, pers_probs_hadgem_hi, 'r-', linewidth=3.0, label='HADGEM_Hi')
    plt.plot(x, pers_probs_hadgem_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.xlabel('Persistence probabilities', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax3.tick_params(labelsize=12)
    plt.title('(c) ERA-Interim vs HADGEM (1986-2011)', fontsize=17)
    
    ax4 = fig.add_subplot(224)
    plt.plot(x, pers_probs_mri_low, 'r--', linewidth=3.0, label='MRI_Low')
    plt.plot(x, pers_probs_mri_hi, 'r-', linewidth=3.0, label='MRI_Hi')
    plt.plot(x, pers_probs_mri_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.xlabel('Persistence probabilities', fontsize=14)
    #plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax3.tick_params(labelsize=12)
    plt.title('(d) ERA-Interim vs MRI (1979-2011)', fontsize=17)
    
    plt.tight_layout()
    #plt.suptitle('PDFs of %s persistence probabilities' % regime, fontsize=20)
    #plt.subplots_adjust(top=0.85)
    plt.savefig('%s/%s_all_persistence_pdf.jpg' % (figpath, regname), dpi=300, format='jpg')  
    plt.close()
    #plt.show()
    


def plot_histogram(regime):


    if regime == 'NAO+':
        regname = 'naoplus'
    elif regime == 'NAO-':
        regname = 'naominus'
    elif regime == 'ATLR':
        regname = 'atlr'
    elif regime == 'BLK':
        regname = 'blk'
    else:
        sys.exit("Error in transition_probs: invalid choice of regime")

    pers_probs_ece_low = load_pers_probs('ECE_Low', regime)
    pers_probs_ece_hi = load_pers_probs('ECE_Hi', regime)
    pers_probs_ece_era = load_pers_probs('ERAInterim_%s-%s' % (1979, 2008), regime)

    pers_probs_hadgem_low = load_pers_probs('HADGEM_Low', regime)
    pers_probs_hadgem_hi = load_pers_probs('HADGEM_Hi', regime)
    pers_probs_hadgem_era = load_pers_probs('ERAInterim_%s-%s' % (1986, 2011), regime)

    pers_probs_mri_low = load_pers_probs('MRI_Low', regime)
    pers_probs_mri_hi = load_pers_probs('MRI_Hi', regime)
    pers_probs_mri_era = load_pers_probs('ERAInterim_%s-%s' % (1979, 2010), regime)

    pers_probs_ncep = load_pers_probs('NCEP', regime)
    pers_probs_ncep_era = pers_probs_mri_era

    print pers_probs_ece_low.shape
    print pers_probs_ece_era.shape 


    fig = plt.figure(figsize=(12,8))
    x = np.linspace(0,1,500)

    ax1 = fig.add_subplot(221)
    #plt.hist(pers_probs_ncep, 'b--', label='NCEP')
    plt.hist(pers_probs_ece_low, facecolor='g', normed=True, label='ERA-Interim')
    #plt.xlabel('Persistence probabilities', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax1.tick_params(labelsize=12)
    plt.title('ERA-Interim vs NCEP', fontsize=18)

    ax2 = fig.add_subplot(222)
    #plt.plot(x, pers_probs_ece_low, 'r--', linewidth=3.0, label='ECE_Low')
    #plt.plot(x, pers_probs_ece_hi, 'r-', linewidth=3.0, label='ECE_Hi')
    #plt.plot(x, pers_probs_ece_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.hist(pers_probs_ece_era, facecolor='g', normed=True, label='ERA-Interim')
    #plt.xlabel('Persistence probabilities', fontsize=14)
    #plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax2.tick_params(labelsize=12)
    plt.title('ERA-Interim vs EC-Earth', fontsize=18)

    ax3 = fig.add_subplot(223)
    #plt.plot(x, pers_probs_hadgem_low, 'r--', linewidth=3.0, label='HADGEM_Low')
    #plt.plot(x, pers_probs_hadgem_hi, 'r-', linewidth=3.0, label='HADGEM_Hi')
    #plt.plot(x, pers_probs_hadgem_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.hist(pers_probs_hadgem_era, facecolor='g', normed=True, label='ERA-Interim')
    plt.xlabel('Persistence probabilities', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax3.tick_params(labelsize=12)
    plt.title('ERA-Interim vs HADGEM', fontsize=18)

    ax4 = fig.add_subplot(224)
    #plt.plot(x, pers_probs_mri_low, 'r--', linewidth=3.0, label='MRI_Low')
    #plt.plot(x, pers_probs_mri_hi, 'r-', linewidth=3.0, label='MRI_Hi')
    #plt.plot(x, pers_probs_mri_era, 'b-', linewidth=3.0, label='ERA-Interim')
    plt.hist(pers_probs_mri_era, facecolor='g', normed=True, label='ERA-Interim')
    plt.xlabel('Persistence probabilities', fontsize=14)
    #plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.ylim([0,10])
    ax3.tick_params(labelsize=12)
    plt.title('ERA-Interim vs MRI', fontsize=18)

    plt.tight_layout()
    #plt.suptitle('PDFs of %s persistence probabilities' % regime, fontsize=20)
    #plt.subplots_adjust(top=0.85)
    plt.savefig('%s/%s_all_persistence_histograms.jpg' % (figpath, regname), dpi=300, format='jpg')
    plt.close()


###############################################################################
#
#              DOING THIS NOW FOR THE DIFFERENT MODELS WE HAVE
#              AS WELL AS FOR ERA-INTERIM 
#
############################################################################### 
    
       

#Estimating probs
def name_to_dataset(name):

    if name == 'ECE_lowres':
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_ECEARTH31_base_T255_3ens_DJF_EAT_1979-2008_4pcs.txt' % datapath)
        dataname = 'ECE_Low'
    
    elif name == 'ECE_hires':    
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_ECEARTH31_base_T511_3ens_DJF_EAT_1979-2008_4pcs.txt' % datapath)
        dataname = 'ECE_Hi'

    elif name == 'HADGEM_lowres':    
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_HadGEM3-GA3_base_N216_3ens_DJF_EAT_1986-2011_4pcs.txt' % datapath)
        dataname = 'HADGEM_Low'
        
    elif name == 'HADGEM_hires':    
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_HadGEM3-GA3_base_N512_3ens_DJF_EAT_1986-2011_4pcs.txt' % datapath)
        dataname = 'HADGEM_Hi' 
        
    elif name == 'MRI_lowres':
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_MRIAGCM32_base_TL95_3ens_DJF_EAT_1979-2010_4pcs.txt' % datapath)
        dataname = 'MRI_Low'
    
    elif name == 'MRI_hires':
        data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_MRIAGCM32_base_TL319_3ens_DJF_EAT_1979-2010_4pcs.txt' % datapath)
        dataname = 'MRI_Hi'
    
    return data, dataname
     
     
def estimate_probs(names):

    for name in names:
        print name
        data, dataname = name_to_dataset(name)
        transition_probs(data, dataname, 'NAO+', logodds=False, skip=True, save=True)
        transition_probs(data, dataname, 'NAO-', logodds=False, skip=True, save=True)
        transition_probs(data, dataname, 'ATLR', logodds=False, skip=True, save=True)
        transition_probs(data, dataname, 'BLK', logodds=False, skip=True, save=True)


def estimate_probs_era(firstyear, lastyear):

    data = np.loadtxt('%s/indclORD_4clus_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_%s-%s_4pcs.txt' % (datapath, firstyear, lastyear))
    dataname = 'ERAInterim_%s-%s' % (firstyear, lastyear)
    transition_probs(data, dataname, 'NAO+', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'NAO-', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'ATLR', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'BLK', logodds=False, skip=True, save=True)    


def estimate_probs_ncep():

    data = np.loadtxt('%s/indclORDasREF_4clus_zg500_day_NCEPNCAR_obs_144x73_1ens_DJF_EAT_1979-2010_4pcs.txt' % (datapath))
    dataname = 'NCEP'
    transition_probs(data, dataname, 'NAO+', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'NAO-', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'ATLR', logodds=False, skip=True, save=True)
    transition_probs(data, dataname, 'BLK', logodds=False, skip=True, save=True)   


#Plotting
def plot_pdfs():

    for name in ['ECE', 'HADGEM', 'MRI']:
        plot_lognormal(name, 'NAO+')
        plot_lognormal(name, 'NAO-')
        plot_lognormal(name, 'ATLR')
        plot_lognormal(name, 'BLK')





###############################################################################
#
#                           RUNNING AS MAIN
#
###############################################################################


#estimate_probs_era(1979, 2008)
#estimate_probs_era(1979, 2010)
#estimate_probs_era(1986, 2011)
#estimate_probs_ncep()
 
   
#plot_pdfs()
#plot_lognormal_regime('NAO+')
#plot_lognormal_regime('NAO-')
#plot_lognormal_regime('ATLR')
plot_lognormal_regime('BLK')

#plot_histogram('BLK') 
 
 
        
        
