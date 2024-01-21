import numpy as np
import pandas as pd
import os,glob
import bagpipes as pipes
import pygtc
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from astropy.cosmology import z_at_value
from astropy.table import Table,Column
from astropy.io import fits
from scipy.interpolate import interp1d

import os
import sys
from IPython import display

import corner
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.io.fits as fits
import random

from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
from scipy import integrate
import matplotlib
import copy
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D   
import matplotlib.patches as mpatches

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

nscat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/legac_team_jun22.cat",delim_whitespace=True,comment='#')
nscat['idm'] = nscat['ID'].astype(str) + '+' + nscat['MASK'].astype(str)
pcat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/UVISTA_final_v4.1.cat", delim_whitespace=True)


def get_prospec_sfh(idm):
    
    legac_id = scat['ID_LEGAC'][scat['idm']==idm].iloc[0]
    z = scat['Z_SPEC'][scat['idm']==idm].iloc[0]
    age_of_universe = cosmo.age(z).value   #in Gyrs
    delta_t = 1

    plt.rcParams['xtick.labelsize'] = '20'
    plt.rcParams['ytick.labelsize'] = '20'
    plt.rcParams["figure.figsize"] = [15.0,8.0]

    f = fits.open('/Users/yashakaushal/Documents/legac_19july/angelos_aug22/legac_dr3_final_imf_1/' + str(legac_id) + '_SFH.fits')

    data = f[1].data
    header = f[1].header
    data.columns
    x = np.array(data.lookback_time)
    x = np.insert(x, 0, 0, axis=0)
    y = data.SFR
    yl = y - (data.SFR_m) # 50 -16 
    yu = (data.SFR_p) + y  # 84 - 50

    sfh_x,sfh_y, sfh_yl, sfh_yu = [[] for i in range(4)]

    for i in range(len(x)-1):                                        #make x axis array elements

        globals()['x%s' %i] = np.linspace(x[i],x[i+1],100)
        sfh_x.extend(np.linspace(x[i],x[i+1],100))

    for i in range(len(x)-1):  

        globals()['y%s' %(i)] = y[i]*np.ones(100)
        sfh_y.extend(y[i]*np.ones(100))

        globals()['y%sl' %(i)] = yl[i]*np.ones(100)
        sfh_yl.extend(yl[i]*np.ones(100))

        globals()['y%su' %(i)] = yu[i]*np.ones(100)
        sfh_yu.extend(yu[i]*np.ones(100))

    sfh_x = np.concatenate((x0,x1,x2,x3,x4,x5,x6,x7),axis=0)
    sfh_y = np.concatenate((y0,y1,y2,y3,y4,y5,y6,y7),axis=0)
    sfh_yl = np.concatenate((y0l,y1l,y2l,y3l,y4l,y5l,y6l,y7l),axis=0)
    sfh_yu = np.concatenate((y0u,y1u,y2u,y3u,y4u,y5u,y6u,y7u),axis=0)
    
    return(sfh_x/1000,sfh_y,sfh_yl,sfh_yu)


def get_bp_sfh_spectra(idm):
    
    path = "/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/s_final_fulldata_bj/bp_only_medians/pipes/posterior/" 
    os.chdir(path + '../../')
#     from Age_Z_functions import * 

    scat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/legac_team_mar21.cat",delim_whitespace=True,skiprows=[x for x in range(1,188)])
    scat['idm'] = scat['ID'].astype(str) + '+' + scat['MASK'].astype(str)
    pcat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/UVISTA_final_v4.1.cat", delim_whitespace=True)
    photo_filtlist=np.loadtxt("/Users/yashakaushal/Documents/bagpipes/catalogs/bvrizyj_filt_list.txt",dtype='str')

    ID = idm.split('+')[0]
    mask = idm.split('+')[1]
    z = scat[(scat['ID']== int(ID)) & (scat['MASK']== int(mask))]['Z_SPEC'].iloc[0]
    galaxy = pipes.galaxy(idm,load_both,filt_list=photo_filtlist)
    fit = pipes.fit(galaxy, get_fit_instructions_photo())
    fit.fit(verbose=True,n_live=1000)
    fit.posterior.get_advanced_quantities()
    print(list(fit.posterior.samples))
    
    age_of_universe = cosmo.age(z).value
    ages = list(age_of_universe - (fit.posterior.sfh.ages*10**-9))
    
    sfh50 = list(np.percentile(fit.posterior.samples["sfh"], 50, axis=0))
    sfh16 = list(np.percentile(fit.posterior.samples["sfh"], 16, axis=0))
    sfh84 = list(np.percentile(fit.posterior.samples["sfh"], 84, axis=0))
     
    spec_post = np.percentile(fit.posterior.samples["spectrum_full"],(16, 50,84), axis=0)
    spec_post = spec_post.astype(float)
    spec_50 = (spec_post[1,:])
    waves = fit.posterior.model_galaxy.wavelengths
    
    # plotting photometry
    photometry = np.copy(fit.galaxy.photometry)
    pmask = (photometry[:, 1] > 0.)
    ymax = 1.05*np.max((photometry[:, 1]+photometry[:, 2])[pmask])
    m16 = np.percentile(fit.posterior.samples['photometry'],16,axis=0)
    m50 = np.percentile(fit.posterior.samples['photometry'],50,axis=0)
    m84 = np.percentile(fit.posterior.samples['photometry'],84,axis=0)
    
    return(ages,sfh16,sfh50,sfh84,waves,spec_50,photometry,m16,m50,m84)


def get_bp_posteriors(idm):
    
    path = "/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/s_final_fulldata_bj/bp_only_medians/pipes/posterior/" 
    os.chdir(path + '../../')
#     from Age_Z_functions import * 

    scat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/legac_team_mar21.cat",delim_whitespace=True,skiprows=[x for x in range(1,188)])
    scat['idm'] = scat['ID'].astype(str) + '+' + scat['MASK'].astype(str)
    pcat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/UVISTA_final_v4.1.cat", delim_whitespace=True)
    photo_filtlist=np.loadtxt("/Users/yashakaushal/Documents/bagpipes/catalogs/bvrizyj_filt_list.txt",dtype='str')

    ID = idm.split('+')[0]
    mask = idm.split('+')[1]
    z = scat[(scat['ID']== int(ID)) & (scat['MASK']== int(mask))]['Z_SPEC'].iloc[0]
    galaxy = pipes.galaxy(idm,load_both,filt_list=photo_filtlist)
    fit = pipes.fit(galaxy, get_fit_instructions_photo())
    fit.fit(verbose=True,n_live=1000)
    fit.posterior.get_advanced_quantities()
    print(list(fit.posterior.samples))
    
    mass = list(fit.posterior.samples["formed_mass"])
    sfr = list(fit.posterior.samples["sfr"])
    av = list(fit.posterior.samples["dust:Av"])
    met = list(fit.posterior.samples["dblplaw:metallicity"])
    age = list(fit.posterior.samples["mass_weighted_age"])
    
    return (mass,sfr,av,met,age)

def get_prosp_posteriors(idm):
    
    legac_id = scat['ID_LEGAC'][scat['idm']==idm].iloc[0]
    n_samples=200
    
    path = '/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/prospec_final_fulldata/post_phot_full_final/'
    f = fits.open(path + str(legac_id) + '_post_ages.fits')
    data = f[1].data
    header = f[1].header
    p_mwage = data['mw_age'] 
    p_lwage = data['lw_age_rband']
    
    path = '/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/prospec_final_fulldata/post_phot_full_final/'
    f = fits.open(path + str(legac_id) + '_post_params.fits')
    data = f[1].data
    header = f[1].header

    p_av = random.sample(list(data['dust2'][0]),500)
    p_Z = random.sample(list(data['logzsol'][0]),500)
    p_dust_n = random.sample(list(data['dust_index'][0]),500)
    p_smass = random.sample(list(data['logmass'][0]),500)
    p_sfr = random.sample(list(data['logmass'][0]),500)
    

    return (p_smass,p_sfr,p_av,p_Z,p_mwage)

def get_prospec_model(idm):
    
    legac_id = scat['ID_LEGAC'][scat['idm']==idm].iloc[0]
    z = scat['Z_SPEC'][scat['idm']==idm].iloc[0]

    path = '/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/prospec_final_fulldata/legac_dr3_final_imf/'
    hdul = fits.open(path + str(legac_id) +'_best_model.fits')
    data_a = hdul[1].data
    header_a = hdul[1].header
    hdul.close()
    wl = data_a['wavelength'] * 1e4 # microns to angstroms
    wl_um = data_a['wavelength']
    c_angstrom = 2.998e18  # Angstrom/s
    factor = (np.square(wl) * (1+z)**2 / c_angstrom) * 1e26 # 10**-19 erg s-1 cm-2 Angstrom-1 to mJy ---- observed wl
    spec_flux = np.array(data_a['Flux_density']) / factor
    
    phot_file = '/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/prospec_final_fulldata/legac_dr3_final_imf_1_model_photometry.fits'

    table = Table.read(phot_file)
    df = table.to_pandas()
    df['ID_LEGAC'] = (df['ID_LEGAC']).astype(int)
    
    phot_wls = [ 4427.13421337,  5454.8498548,  6248.95357062,  7645.8890470, 9010.9118943,  10203.27590588, 12499.31100549]
    phot_factor = (np.square(phot_wls) * (1+z)**2 / c_angstrom) * 1e29 # 10**-19 erg s-1 cm-2 Angstrom-1 to mJy ---- observed wl
    phot_flux = np.array(df[df['ID_LEGAC']==legac_id].iloc[0][:-1]) / phot_factor 

    return(wl,spec_flux,phot_flux)

def get_bp_model(idm):
    
    path = "/Users/yashakaushal/Documents/my_papers/paper1_sfh/data/s_final_fulldata_bj/bp_only_medians/pipes/posterior/" 
    os.chdir(path + '../../')
#     from Age_Z_functions import * 

    scat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/legac_team_mar21.cat",delim_whitespace=True,skiprows=[x for x in range(1,188)])
    scat['idm'] = scat['ID'].astype(str) + '+' + scat['MASK'].astype(str)
    pcat = pd.read_csv("/Users/yashakaushal/Documents/bagpipes/catalogs/UVISTA_final_v4.1.cat", delim_whitespace=True)
    photo_filtlist=np.loadtxt("/Users/yashakaushal/Documents/bagpipes/catalogs/bvrizyj_filt_list.txt",dtype='str')

    ID = idm.split('+')[0]
    mask = idm.split('+')[1]
    z = scat[(scat['ID']== int(ID)) & (scat['MASK']== int(mask))]['Z_SPEC'].iloc[0]
    galaxy = pipes.galaxy(idm,load_both,filt_list=photo_filtlist)
    fit = pipes.fit(galaxy, get_fit_instructions_photo())
    fit.fit(verbose=True,n_live=1000)
    fit.posterior.get_advanced_quantities()
    print(list(fit.posterior.samples))

    spec_post = np.percentile(fit.posterior.samples["spectrum_full"],(16, 50,84), axis=0)
    spec_post = spec_post.astype(float)
    spec_50 = (spec_post[1,:])
    waves = fit.posterior.model_galaxy.wavelengths
    
    # plotting photometry
    photometry = np.copy(fit.galaxy.photometry)
    pmask = (photometry[:, 1] > 0.)
    ymax = 1.05*np.max((photometry[:, 1]+photometry[:, 2])[pmask])
    m16 = np.percentile(fit.posterior.samples['photometry'],16,axis=0)
    m50 = np.percentile(fit.posterior.samples['photometry'],50,axis=0)
    m84 = np.percentile(fit.posterior.samples['photometry'],84,axis=0)

    return(waves,spec_50,photometry,m16,m50,m84)

def get_redshift(IDM):
    ID = int(IDM.split("+")[0])
    mask = int(IDM.split("+")[1])
    z = nscat[(nscat["ID"]==ID) & (nscat["MASK"]==mask)]["Z_SPEC"].iloc[0]
    return(z)

def get_legac_spec(IDM):

    ID = IDM.split('+')[0]
    mask = IDM.split('+')[1]
    path = '/Users/yashakaushal/Documents/bagpipes/spectra_v3.11/'

    sname = fits.open(path + 'spec1d_v3.11/' + 'M' + str(mask) + '/skysub' +  "/legac_M" + str(mask) + "_v3.11_spec1d_" + str(ID) + ".fits")
    sheader = sname[0].header
    sdata = sname[0].data
    sname.close()
    wl101 = sheader["CRVAL1"] + np.arange(0,sheader["CD1_1"]*sheader["NAXIS1"],sheader["CD1_1"])
    ivname = fits.open(path + 'wht1d_v3.11/' + 'M' + str(mask) + '/skysub' + "/legac_M" + str(mask) + "_v3.11_wht1d_" + str(ID) + ".fits")
    iv = ivname[0].data
    ivname.close()
    err101 = (iv)**-0.5
#     print(list(err101))
    Mask = (err101!=0) & (~np.isnan(sdata)) & (sdata>0) & (np.isfinite(err101))
    print(np.sum(Mask))
    wl101 = wl101[Mask]
    spectrum101 = sdata[Mask]*1e-19
    err101 = err101[Mask]*1e-19

    return(np.c_[wl101,spectrum101,err101].T)
