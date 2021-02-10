#!/usr/bin/env python3

import setuptools
from distutils.core import setup

no_pip_deps = [
  "ete3 @ https://github.com/etetoolkit/ete/archive/3.1.2.zip",
  "toytree @ https://github.com/eaton-lab/toytree/archive/v1.0.0.zip"
]

pip_deps = [
  "biopython",
  "scipy",
  "matplotlib",
  "seaborn",
  "pandas",
  "plotly",
  "numpy",
  "six",
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
    scripts = scripts,
    install_requires = pip_deps + no_pip_deps,
    zip_safe = False,
    python_requires='>3,<3.9', # <3.9 to avoid WarningSytanx messages
    classifiers = [
        'Programming Language :: Python :: 3'
        ]
    )
