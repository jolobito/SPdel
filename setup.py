#!/usr/bin/env python3

import setuptools
from distutils.core import setup

dependencies = [
                "biopython",
                "scipy",
                "matplotlib",
                "pandas",
                "plotly",
                "numpy",
                "python-kaleido"
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
    scripts = scripts,
    install_requires = dependencies,
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
