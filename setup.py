#!/usr/bin/python
# setup.py
# written on 3/7/2017 by Xiangchen Li

import os
import re
from setuptools import setup


current_version = ''
with open('src/__init__.py') as f:
    for line in f.readlines():
        m = re.search(r"^__version__ = '(.*)'$", line.strip())
        if m:
            current_version = m.group(1)
            break

setup(
    name='codonPY',
    version=current_version,
    scripts=[os.path.join('bin', 'codonPY.py')],
    packages=['src'],
    package_data={'src': ['data/*.fna']},
    url='https://github.com/cvn001/codonPY',
    license='MIT',
    author='Xiangchen Li',
    author_email='lixiangchenxy@outlook.com',
    description='A Python package for codon usage bias research'
)
