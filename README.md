# HST IIb progenitors
CAUTION: In some cases, the sampler will experience label switching and the chains will not converge. A fix is in the works...

## Installation

In the directory where you want to install package (not necessarily where you want to run inference):

```
git clone https://github.com/niharika-sravan/HST_IIb_progenitors.git
cd HST_IIb_progenitors
conda env create -f environment.yml
conda activate IIb_prog
pip install .
```

Download HST Files:
https://stsynphot.readthedocs.io/en/latest/stsynphot/data_hst.html

```
tar xzf your/downloads/directory/synphot1.tar.gz -C IIb_progenitor
mkdir IIb_progenitor/grp/redcat/trds/calspec
cp IIb_progenitor/data/*.fits grp/redcat/trds/calspec/.
```
The last moves necessary spectrum files in the right location.

Finally, either add this to your shell start up or repeat this for every new shell session:

```
export PYSYN_CDBS=/install/directory/IIb_progenitor/grp/redcat/trds
```

## Simple example execution

In the folder you want the results in:
```
>>> from IIb_progenitor import main
>>> SN = main.Supernova('2016gkg',26.4,5.3,
                     [23.41,23.08,23.29],
                     [0.47,0.33,0.59],
                     ['wfpc2,2,f450w,a2d7,cont#50180',
                     'wfpc2,2,f606w,a2d7,cont#50180',
                     'wfpc2,2,f814w,a2d7,cont#50180'],
                     mw_ext=[3.1,0.017])
>>> SN.get_mod_props()
>>> SN.run_posterior(diag=True)
```

See main.py for more examples.
