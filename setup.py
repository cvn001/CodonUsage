#!/usr/bin/python
# setup.py
# written on 3/7/2017 by Xiangchen Li

import os
import re
import sys
from setuptools import setup


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
    name='codonPY',
    version=current_version,
    scripts=[os.path.join('bin', 'codonPY.py')],
    packages=['src'],
    package_data={'src': ['data/*.fna']},
    url='https://github.com/cvn001/codonPY',
    license='MIT',
    platforms="Posix; MacOS X",
    author='Xiangchen Li',
    author_email='lixiangchenxy@outlook.com',
    description='A Python package for codon usage bias research',
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
