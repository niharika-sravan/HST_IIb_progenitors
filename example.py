from main import Supernova

SN=Supernova('2016gkg',26.4,5.3,
              [23.41,23.08,23.29],
              [0.47,0.33,0.59],
              ['wfpc2,2,f450w,a2d7,cont#50180',
              'wfpc2,2,f606w,a2d7,cont#50180',
              'wfpc2,2,f814w,a2d7,cont#50180'],
              mw_ext=[3.1,0.017])

SN=Supernova('2013df',16.6,0.4,
               [24.535,23.144],
               [0.071,0.055],
               ['wfpc2,2,f555w,a2d7,cont#50180',
               'wfpc2,2,f814w,a2d7,cont#50180'],
               mw_ext=[3.1,0.017],ntst_ext=[3.1,0.08])

#Detector not sure, and other obsmode kws

SN=Supernova('2011dh',8.58,0.1,
               [23.39,22.36,21.83,21.28,21.20],
               [0.25,0.02,0.04,0.04,0.03],
               ['wfpc2,2,f336w,a2d7,cont#50180',
               'acs,wfc1,f555w',
               'acs,wfc1,f658n',
               'acs,wfc1,f814w'],
               mw_ext=[3.1,0.07])

#wfc1 or 2? Need mjd or aperture?


SN.get_mod_props(saved=True).run_posterior(diag=True)
