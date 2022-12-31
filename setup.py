from setuptools import setup, find_packages

setup(
    name='IIb_progenitor',
    version='0.0.1',
    packages=find_packages(),
    package_data={'IIb_progenitor': ['data/*']},
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
