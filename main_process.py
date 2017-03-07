#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: 本程序用于
# Created by galaxy on 17-3-7 下午10:08

import os
import sys
import time
import shutil
import subprocess
import logging.handlers
import traceback
from multiprocessing import Pool, cpu_count
from argparse import ArgumentParser


def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="main_process.py")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
                        help="How many threads will be used? [default all]")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                        help="Give verbose output")
    return parser.parse_args()


