import pkg_resources
DATA_PATH = pkg_resources.resource_filename('IIb_progenitor', 'data/')

Z_list=[0.005,0.02]
eps_list=[0.1,0.5]

alpha=2.3
beta=-1.
gamma=-0.22

#grid properties
#Note: grid M and P are base 10 log
params=['log10(M_1i)(Msun)','q_i(M_2i/M_1i)','log10(P_i)(days)']
dellogM=0.0005
dellogM1=0.02
delq=0.05
dellogP_low=0.1
dellogP_high=0.02

#emcee parameters; nsteps and burnin can be passed to run_posterior
nwalkers=32
disturb=0.1
ndim=len(params)
thin=200 #this is around half autocorr length
