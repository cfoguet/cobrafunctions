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
      install_requires=["anyio==3.4.0", "appdirs==1.4.4", "async-generator==1.10", "charset-normalizer==2.0.9", "colorama==0.4.4", "commonmark==0.9.1", "contextvars==2.4", "dataclasses==0.8", "depinfo==1.7.0", "diskcache==5.3.0", "et-xmlfile==1.1.0", "future==0.18.2", "h11==0.12.0", "httpcore==0.14.3", "httpx==0.21.1", "idna==3.3", "immutables==0.16", "importlib-metadata==4.8.2", "importlib-resources==5.4.0", "mpmath==1.2.1", "numpy==1.19.5", "openpyxl==3.0.9", "optlang==1.5.2", "pandas==1.1.5", "pydantic==1.8.2", "Pygments==2.10.0", "python-dateutil==2.8.2", "python-libsbml==5.19.0", "pytz==2021.3", "rfc3986==1.5.0", "rich==10.15.2", "ruamel.yaml==0.17.17", "ruamel.yaml.clib==0.2.6", "six==1.16.0", "sniffio==1.2.0", "swiglpk==5.0.3", "sympy==1.9", "typing_extensions==4.0.1", "zipp==3.6.0","scipy==1.5.4","cobra==0.22.1","cobrafunctions"],
      scripts=["gim3e_and_sampling.py","run_qMTA.py","map_expression_to_reactions.py"],
      packages=find_packages(),
      package_dir = {'cobrafunctions': 'cobrafunctions'},
     )

#["create_and_solve_iso2flux_model.py"]
#      python_requires='<3',  # Your supported Python ranges

