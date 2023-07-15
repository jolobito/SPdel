#!/usr/bin/env python3

import sys
import platform
import setuptools
from distutils.core import setup


if platform.architecture()[0] != '64bit':
    sys.stderr.write('Architecture requires 64bit')
    sys.stderr.flush()
    exit()


myos = sys.platform

if myos == 'darwin':
    bins = ['./ext_bin/abgd_darwin','./ext_bin/asap_darwin','./ext_bin/mptp_darwin']

elif myos == 'linux' or myos == "linux2":
    bins = ['./ext_bin/abgd_linux','./ext_bin/asap_linux','./ext_bin/mptp_linux']

elif myos == 'win32':
    bins = ['./ext_bin/abgd.exe','./ext_bin/asap.exe','./ext_bin/mptp.exe']
    
else:
    sys.stderr.write('Package does not work with %s operative system'  % myos)
    sys.stderr.flush()
    exit()


dependencies = [
                "biopython",
                "scipy",
                "matplotlib",
                "pandas",
                "plotly",
                "numpy",
                "ipykernel",
                "nbformat"
                ]

scripts = [
           "SPdel.py",
           "spdelib/summary.py",
           "spdelib/PTP.py",
           "spdelib/bPTP.py",
           "spdelib/GMYC.py",
           "spdelib/Diagnoser.py"           
           ]

setup(
    name = "spdel",
    author='Jorge L. Ramirez',
    version  = "1.0",
    packages = ["spdelib"],
    package_dir={"spdelib": "spdelib"},
    package_data = {'spdelib': bins},
    data_files = [ ('bin', bins) ],
    include_package_data=True,
    zip_safe = False,
    scripts = scripts,
    install_requires = dependencies,
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
