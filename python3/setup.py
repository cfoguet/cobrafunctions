# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

#Change directory to setup.py directory
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

setup(name='cobrafunctions',
      version='1.2',
      description='cobrafunctions',
      author='Carles Foguet',
      author_email='',
      url='',
      install_requires=["openpyxl==3.0.9","cobra==0.29.1","scipy==1.14.1","cobrafunctions"],
      scripts=["gim3e_and_sampling.py","run_qMTA_individual_samples.py","map_expression_to_reactions.py"],
      packages=find_packages(),
      package_dir = {'cobrafunctions': 'cobrafunctions'},
     )

