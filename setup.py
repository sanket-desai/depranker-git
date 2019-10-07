#!/usr/bin/env python

from distutils.core import setup

setup(name='DepRanker',
      version='0.1.0',
      description='A gene impact score calculator (prioritization method) for RNAi / CRISPR screen results',
      author='Sanket Desai',
      author_email='desai,sanket12@gmail.com',
      url='https://github.com/sanket-desai/depranker-git',
      license = "MIT",
      long_description=open('README.md').read()
      #nstall_requires=['numpy']
     )
