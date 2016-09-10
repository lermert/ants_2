#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name = 'ants_2',
    version = '0.0.0a0',
    description = 'Package to download, preprocess and correlate ambient noise \
    data in a simple way.',
    #long_description =
    # url = 
    author = 'ETH CSE Noise',
    author_email  = 'lermert@student.ethz.ch',
    # license
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Topic :: Seismology',
        'Programming Language :: Python :: 2',
    ],
    keywords = 'Ambient seismic noise',
    packages = find_packages(),
    #package_data = ,
    install_requires = [
        "obspy>=1.0.1",
        "geographiclib",
        "click",
        "h5py",
        "mpi4py>=2.0.0"],
    
    # ToDo: Add entry points for test suite
    entry_points = {
        'console_scripts': [
            'ants = ants_2.main:run'            
        ]
    },
)