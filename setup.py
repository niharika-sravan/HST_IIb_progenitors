from setuptools import setup, find_packages

setup(
    name='IIb_progenitors',
    version='0.0.1',
    packages=find_packages(),
    author='Niharika Sravan',
    author_email='niharika.sravan@gmail.com',
    install_requires=['matplotlib',
                      'pandas',
                      'astropy',
                      'emcee',
                      'stsynphot',
                      'corner',
                      'dust_extinction',
                      'tdqm'
                     ]
)
