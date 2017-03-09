#!/usr/bin/python
# setup.py
# written on 3/7/2017 by Xiangchen Li

import os
from setuptools import setup

setup(
    name='codonPY',
    version='v0.2.0',
    scripts=[os.path.join('bin', 'codonPY.py')],
    packages=['src'],
    package_data={'src': ['data/*.fna']},
    url='https://github.com/cvn001/codonPY',
    license='MIT',
    author='Xiangchen Li',
    author_email='lixiangchenxy@outlook.com',
    description='A Python package for codon usage bias research'
)
