# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

#Change directory to setup.py directory
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

setup(name='cobrafunctions',
      version='0.1',
      description='cobrafunctions',
      author='Carles Foguet',
      author_email='',
      url='',
      install_requires=["openpyxl==3.0.9","cobra","scipy","cobrafunctions"],
      scripts=["gim3e_and_sampling.py","run_qMTA.py","map_expression_to_reactions.py"],
      packages=find_packages(),
      package_dir = {'cobrafunctions': 'cobrafunctions'},
     )

