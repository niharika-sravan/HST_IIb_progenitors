import os,sys
import traceback
import numpy as np
import json
import pandas as pd
import math
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import emcee, corner
from itertools import product
from scipy.optimize import minimize
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

from astropy import units as u
from synphot import SourceSpectrum
from synphot import ReddeningLaw
from dust_extinction.parameter_averages import CCM89
from synphot.models import BlackBodyNorm1D
from synphot import Observation
import stsynphot as stsyn

from IIb_progenitor import defs

class Supernova():
  def __init__(self,ID:str,distance:float,distance_uncer:float,mags:list,mag_uncers:list,obsmode:list,
              excl_comp_sed:bool=True,obsmode_comp:list=[],
              ntst_ext:list=[],mw_ext:list=[],
              Z:float=None,mt_eff:float=None,
              envM_const:list=[],R_const:list=[],coreM_const:list=[],Mdot_vwind_const:list=[],
              outdir:str='output'
              ):
    """
    Class that records all properties including output posteriors

    Arguments
    ---------
    ID: name of supernova
    distance: distance to object in Mpc
    distance_uncer: distance uncertainity in Mpc
              propogated to abs mag calculated
    mags: list of HST mags, Vega
    mag_uncers: list of associated uncertainities for HST mags
    obsmode: list of corresponding obsmode strings to generate theoretical mags
              https://stsynphot.readthedocs.io/en/latest/stsynphot/appendixb.html#stsynphot-appendixb

    Keyword Arguments
    -----------------
    excl_comp_sed: exclude companion SED when computing theoretical mags. True by default
    obsmode_comp: list of obsmode strings for theoretical companion properties. Skipped by default
              Useful for planning post SN identification.
    ntst_ext: interstellar+host [Rv,E(B-V)], CCM89 model. Ignored if not provided
    mw_ext: Milky-Way [Rv,E(B-V)], CCM89 model. Ignored if not provided
    Z: set progenitor metallicity, options are 0.005 or 0.02 (solar). Uses both if not specified
    mt_eff: set efficiency of mass transfer in binary RLO, options are 0.5 and 0.1. Uses both if not specified

    envM_const: impose constraints on envelope mass in Msun. None by default
              typically would come from SN light curve fitting.
    R_const: impose constraints on progenitor radius in Rsun. None by default
              typically would come from SN light curve fitting.
    coreM_const: impose constraints on core mass in Msun. None by default
              can come from light curve modeling or ejecta calcuations.
    Mdot_vwind_const: impose constraints on Mdot/vwind in Msun yr-1/km s-1. None by default
              typically comes from X-ray and radio measurements.

      for all four above
        if three values in list, interpretted as 1-sigma uncertainties
            In the case of normal distrbution, repeat the upper and lower limits [mean,sigma,sigma]
            otherwise [mean,sigma_minus,sigma_plus], uses split normal funct
        if two values in list, interpretted as range
            [lower,upper]. Skip one value to use as upper or lower limit [None,upper]
        ignored if None or list empty

    outdir: path to folder where outputs need to be written

    Usage:
    ------
    Define the progenitor properties.
    Minimal intialization:
    >>> from main import Supernova
    >>> SN=Supernova('2016gkg',26.4,5.3,
              [23.41,23.08,23.29],
              [0.47,0.33,0.59],
              ['wfpc2,2,f450w,a2d7,cont#50180',
              'wfpc2,2,f606w,a2d7,cont#50180',
              'wfpc2,2,f814w,a2d7,cont#50180'])

    Detailed specification:
    >>> SN=Supernova('2013df',16.6,0.4,
              [24.535,23.144],
              [0.071,0.055],
              ['wfpc2,3,f555w,a2d7,cont#50180',
              'wfpc2,3,f814w,a2d7,cont#50180'],
              excl_comp_sed=False,obsmode_comp=['wfc3,2,'f218w'],
              ntst_ext=[3.1,0.25],mw_ext=[3.1,0.053],
              Z=0.02,mt_eff=0.1,
              envM_const=[0.05,0.09],R_const=[64,169],coreM_const=[2,3.6],Mdot_vwind_const=[10.5e-6,3e-6,3e-6])

    Get theoretical model properties (takes an hour):
    >>> SN.get_mod_props()
    Need to run just once for a choice of extinction and obsmode. If already run once (have *_mod_props.csv
    in outdir) turn on 'saved' flag.

    Run posterior calculation:
    Outputs posterior plots in outdir using preset values for nsteps and burnin
    >>> SN.run_posterior()
    If autocorrelation time is shorter than nsteps over 50 an emcee.autocorr.AutocorrError will be displayed
    You can increase chain length accordingly.

    To ensure using sensible values for burin set diag to True to view chains (written to outdir) and call
    again with better values
    >>> SN.run_posterior(diag=True)
    >>> SN.run_posterior(nsteps=10000,burnin=100)

    """
    self.ID=ID
    self.distance=distance
    self.distance_uncer=distance_uncer
    self.mags=mags
    self.mag_uncers=mag_uncers
    self.obsmode=obsmode
    self.excl_comp_sed=excl_comp_sed
    self.obsmode_comp=obsmode_comp
    self.ntst_ext=ntst_ext
    self.mw_ext=mw_ext
    self.Z=Z
    self.mt_eff=mt_eff
    self.envM_const=envM_const
    self.R_const=R_const
    self.coreM_const=coreM_const
    self.Mdot_vwind_const=Mdot_vwind_const
    self.outdir=outdir

  def get_mod_props(self,saved=False) -> 'Supernova':
    '''
    creates attribute mod_props that contains theoretical mags
    can take an hour for the full set of theoretical models

    Keyword Arguments:
    ------------------
    saved: Read model properties from a file saved in a previous call.
          Disabled by default because each set depends on Supernova object values for obsmode

    '''
    def get_mags(logT_1,logT_2,R_1,R_2,bp):
      spect1=SourceSpectrum(BlackBodyNorm1D,temperature=10.**logT_1)*(np.pi*R_1**2.)
      spect2=SourceSpectrum(BlackBodyNorm1D,temperature=10.**logT_2)*(np.pi*R_2**2.)
      if self.excl_comp_sed: spect=spect1
      else: spect=spect1+spect2
      wav=np.arange(0.1,3,0.001)*u.micron
      if self.ntst_ext: spect=spect*ReddeningLaw(CCM89(Rv=self.ntst_ext[0])).extinction_curve(
                                                self.ntst_ext[1],wavelengths=wav)
      if self.mw_ext: spect=spect*ReddeningLaw(CCM89(Rv=self.mw_ext[0])).extinction_curve(
                                                self.mw_ext[1],wavelengths=wav)
      mag=Observation(spect,bp,binset=bp.binset,force='taper').effstim('vegamag',vegaspec=stsyn.Vega)
      return mag

    if not saved:
      self.mod_props=pd.DataFrame()
      for Z_ in defs.Z_list:
        for eps_ in defs.eps_list:
          data=pd.read_csv(defs.DATA_PATH+'IIb_props_'+str(Z_)+'_'+str(eps_)+'.csv')
          for bp_str in self.obsmode:
            data[bp_str]=''
            data['Z']=Z_
            data['eps']=eps_
            bp=stsyn.band(bp_str)
            for idx,row in data.iterrows():
              data.at[idx,bp_str]=get_mags(row['logT_1'],row['logT_2'],row['R_1(Rsun)'],row['R_1(Rsun)'],bp)
          self.mod_props=pd.concat([self.mod_props,data])
      self.mod_props=self.mod_props.reset_index(drop=True)
      self.mod_props.to_csv(self.outdir+'/'+self.ID+'_mod_props.csv',index=False)
    else:
      self.mod_props=pd.read_csv(self.outdir+'/'+self.ID+'_mod_props.csv')
    return self

  @staticmethod
  def convert_app_abs(mags:list,mag_uncers:list,distance:float,distance_uncer:float) -> tuple[list,list]:
    for i,mag in enumerate(mags):
      mags[i]=mag-5.*np.log10(distance*1e3)
      sig_d=5.*distance_uncer/(distance*np.log(10))
      mag_uncers[i]=np.sqrt(mag_uncers[i]**2.+sig_d**2)
    return mags,mag_uncers

  def setup_interpolants(self,props:pd.DataFrame):
    '''
    returns a list of linear interpolant for theoretical mag values in each band
    enables fast sampling

    Arguments:
    ----------
    props: Dataframe with grid of mag values per Z and mt_eff

    '''
    for bp_str in self.obsmode:
      props[bp_str]=pd.to_numeric(props[bp_str].astype(str).apply(lambda x: x.split('VEGAMAG')[0]))
    M_arr=np.round(np.arange(props['log10(M_1i)(Msun)'].min()-2.*defs.dellogM1,
                            props['log10(M_1i)(Msun)'].max()+2.*defs.dellogM1,
                            defs.dellogM1),2)
    q_arr=np.round(np.arange(props['q_i(M_2i/M_1i)'].min()-2*defs.delq,
                            props['q_i(M_2i/M_1i)'].max()+2*defs.delq,
                            defs.delq),3)
    if props['Z'].unique()==[0.005]:
      lowP_arr=np.round(np.arange(props['log10(P_i)(days)'].min()-2.*defs.dellogP_low,
                                  2.7+defs.dellogP_low, #2.7 included
                                  defs.dellogP_low),2)
      highP_arr=np.round(np.arange(2.7,
                                  props['log10(P_i)(days)'].max()+2.*defs.dellogP_high,
                                  defs.dellogP_high),2)
      space_lowP=pd.DataFrame(list(product(M_arr,q_arr,lowP_arr)),columns=defs.params)
      space_highP=pd.DataFrame(list(product(M_arr,q_arr,highP_arr)),columns=defs.params)
      space=space_highP.merge(space_lowP,how='outer')
      space=space.merge(props,on=defs.params,how='left',indicator='is_IIb').replace(np.nan,-np.inf)
      space['is_IIb']=space['is_IIb'].replace({'left_only':'no','both':'yes'})

      ntr_lowP=[]
      ntr_highP=[]
      for om in self.obsmode:
        ntr_lowP.append(LinearNDInterpolator(space[space['log10(P_i)(days)']<=2.7][['log10(M_1i)(Msun)','q_i(M_2i/M_1i)','log10(P_i)(days)']],
                                            space[space['log10(P_i)(days)']<=2.7][om]))
        ntr_highP.append(LinearNDInterpolator(space[space['log10(P_i)(days)']>=2.7][['log10(M_1i)(Msun)','q_i(M_2i/M_1i)','log10(P_i)(days)']],
                                            space[space['log10(P_i)(days)']>=2.7][om]))
      return [ntr_lowP,ntr_highP]
    else:
      P_arr=np.round(np.arange(props['log10(P_i)(days)'].min()-2*defs.dellogP_high,
                              props['log10(P_i)(days)'].max()+2*defs.dellogP_high,
                              defs.dellogP_high),2)
      space=pd.DataFrame(list(product(M_arr,q_arr,P_arr)),columns=defs.params)
      space=space.merge(props,on=defs.params,how='left',indicator='is_IIb').replace(np.nan,-np.inf)
      space['is_IIb']=space['is_IIb'].replace({'left_only':'no','both':'yes'})

      ntr=[]
      for om in self.obsmode:
        ntr.append(LinearNDInterpolator(space[['log10(M_1i)(Msun)','q_i(M_2i/M_1i)','log10(P_i)(days)']],
                                        space[om]))
      return ntr

  def run_posterior(self,nsteps:int=15000,burnin:int=1000,diag=False):
    '''
    Main function that runs sampling and creates plots for posteriors and summary statistics

    Keyword Arguments:
    ------------------
    nsteps: Steps in MCMC. This should be 50x greater than estimated autocorrelation time.
    burnin: Number of intial steps to drop in posterior estimates.
            This is typically few times autocorr length

    '''
    def log_prior(theta):
      logM,q,logP=theta
      #barring unphysical values or grid extremes
      if logM<0.9 or q<=0.1 or logP<0.8 or logM>1.7 or q>=1.0 or logP>3.8:
        return -np.inf
      return (-1*defs.alpha*np.log(10.**logM)+defs.beta*np.log(q)+defs.gamma*np.log(logP))

    def log_likelihood(theta):
      logM,q,logP=theta
      mod_mags=np.zeros(len(mags))
      if Z==0.005:
        if logP<2.7:
          for i in range(len(mags)):
            mod_mags[i]=ntr[0][i](theta)
            if math.isnan(mod_mags[i]): return -np.inf
          else:
            for i in range(len(mags)):
              mod_mags[i]=ntr[0][i](theta)
              if math.isnan(mod_mags[i]): return -np.inf
      else:
        for i in range(len(mags)):
          mod_mags[i]=ntr[i](theta)
          if math.isnan(mod_mags[i]): return -np.inf
      return -0.5*np.sum((np.array(mags)-mod_mags)**2/np.array(mag_uncers)**2)

    def log_probability(theta):
      logM,q,logP=theta
      return log_prior(theta)+log_likelihood(theta)

    mags,mag_uncers=self.convert_app_abs(self.mags,self.mag_uncers,self.distance,self.distance_uncer)
    for Z_eps,props in self.mod_props.groupby(['Z','eps']):
      Z,eps=Z_eps[0],Z_eps[1]
      if self.Z:
        if self.Z!=Z: continue
      if self.mt_eff:
        if self.mt_eff!=eps: continue
      ntr=self.setup_interpolants(props)
      nll=lambda *args: -log_likelihood(*args)
      guess=np.array([props['log10(M_1i)(Msun)'].mean(),
                      props['q_i(M_2i/M_1i)'].mean(),
                      props['log10(P_i)(days)'].mean()])
      soln=minimize(nll,guess)
      pos=soln.x+defs.disturb*soln.x*np.random.randn(defs.nwalkers,defs.ndim)
      smplr=emcee.EnsembleSampler(defs.nwalkers,defs.ndim,log_probability)
      smplr.run_mcmc(pos,nsteps,progress=True)
      flat_samples=smplr.get_chain(discard=burnin,thin=defs.thin,flat=True)
      fig=corner.corner(flat_samples,labels=defs.params)
      for i in range(defs.ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "$\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}$"
        txt = txt.format(mcmc[1],q[0],q[1],defs.params[i])
        plt.figtext(0.7,0.9-i*0.05,txt)
      fig.savefig(self.outdir+'/'+self.ID+'_progenitor_dist_'+str(Z_eps)+'.png')
      print('chain length',nsteps)
      try: print('autocorrelation time',max(smplr.get_autocorr_time()))
      except Exception: print(traceback.format_exc())
      if diag:
        fig,axes=plt.subplots(defs.ndim, figsize=(10, 7), sharex=True)
        samples=smplr.get_chain()
        for i in range(defs.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0,len(samples))
            ax.set_ylabel(defs.params[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
        fig.savefig(self.outdir+'/'+self.ID+'_chains_'+str(Z_eps)+'.png')
    return
