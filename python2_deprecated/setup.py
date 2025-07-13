# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

#Change directory to setup.py directory
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

setup(name='cobrafunctions',
      version='0.1',
      python_requires='<3',
      description='cobrafunctions',
      author='Carles Foguet',
      author_email='',
      url='',
      install_requires=["openpyxl==2.3.3","cobra==0.9.0", "future==0.18.2", "mpmath==1.2.1", "numpy==1.16.6", "optlang==1.5.2", "pandas==0.24.2", "python-dateutil==2.8.2", "pytz==2021.3", "ruamel.ordereddict==0.4.15", "ruamel.yaml==0.14.12", "six==1.16.0", "swiglpk==5.0.3", "sympy==1.5.1", "tabulate==0.8.9","scipy==0.19.0","cobra==0.9.0","cobrafunctions"],
      scripts=["gim3e_and_sampling.py","run_qMTA.py","map_expression_to_reactions.py"],
      packages=find_packages(),
      package_dir = {'cobrafunctions': 'cobrafunctions'},
     )

