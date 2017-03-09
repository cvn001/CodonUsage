#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: build setup script
# Created by galaxy on 2017/3/9 0009 17:19

# try using setuptools or distutils.
import os
import re
import sys


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# parse version from package/module without importing or evaluating the code
current_version = ''
with open('src/__init__.py') as f:
    for line in f.readlines():
        m = re.search(r"^__version__ = '(.*)'$", line.strip())
        if m:
            current_version = m.group(1)
            break

if sys.version_info <= (2, 7):
    sys.stderr.write("ERROR: codonPY requires Python 2.7, or Python 3.3 or later. " +
                     "Python {0} detected. Exiting.".format(sys.version_info[:2]))
    sys.exit(1)

setup(
    name="codonPY",
    version=current_version,
    author="Xiangchen Li",
    author_email="lixiangchenxy@outlook.com",
    description=''.join(["src provides a package and script for " +
                         "calculation of codon usage bias."]),
    license="MIT",
    keywords="genome bioinformatics sequence",
    platforms="Posix; MacOS X",
    url="http://cvn001.github.io/codonPY/",  # project home page
    download_url="https://github.com/cvn001/codonPY/releases",
    scripts=[os.path.join('bin', 'codonPY.py')],
    packages=['codonPY'],
    package_data={'codonPY': ['data/']},
    include_package_date=True,
    install_requires=['biopython',
                      'scipy',
                      'numpy'],
    classifiers=[
        'Development Status :: Dev',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )
