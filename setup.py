#!/usr/bin/env python3

import setuptools
from distutils.core import setup

dependencies = [
                "biopython",
                "ete3",
                "scipy",
                "matplotlib",
                "seaborn",
                "nexus",
                "pandas",
                "scipy",
                "plotly",
                "numpy",
                "toytree",
                "toyplot",
                "PyQt5"
                ]

scripts = [
           "SPdel.py",
           "spdelib/summary.py",
           "spdelib/PTP.py",
           "spdelib/bPTP.py",
           "spdelib/GMYC.py"
           ]

setup(
    name = "spdel",
    author='Jorge L. Ramirez',
    version  = "1.0",
    packages = ["spdelib"],
    package_dir={"spdelib": "spdelib"},
    # entry_points = {
    #     'console_scripts': [
    #         'SPdel  = fishlifescript.joinexonfiles:main'
    #         ]
    # },
    scripts = scripts,
    install_requires = dependencies,
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )