#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction:
# Created by galaxy on 17-3-7 10:08 pm

import os
import sys
import time
import shutil
import subprocess
import logging.handlers
import traceback
from multiprocessing import Pool, cpu_count
from argparse import ArgumentParser
from get_RSCU_from_fasta import get_rscu
from global_items import genetic_code
from calculate_relative_adaptiveness import calculate_wi
from Bio import SeqIO
from get_CAI_from_fasta import calculate_cai


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


def each_step(step, input_file, fname):
    if step == 1:
        logger.info('Step 1: get relative synonymous codon usage (RSCU) from fasta')
        step_dir = os.path.join(args.outdirname, 'RSCU')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        output_file = os.path.join(step_dir, fname + '.rscu')
        get_rscu(input_file, output_file, args.interesting)
        logger.info('Step 1 done')
    elif step == 2:
        logger.info('Step 2: calculate relative adaptiveness (W) from fasta')
        step_dir = os.path.join(args.outdirname, 'R_A')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        output_file = os.path.join(step_dir, fname + '.w')
        calculate_wi(input_file, output_file)
    elif step == 3:
        logger.info('Step 3: get Codon Adaptation Index (CAI) from fasta')
        step_dir = os.path.join(args.outdirname, 'CAI')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        # calculate_cai(input_file, output_file, )


def run_in_order(step):
    for root, dirs, files in os.walk(args.outdirname):
        for each_file in files:
            fname = os.path.splitext(each_file)[0]
            input_file = os.path.join(args.outdirname, each_file)
            try:
                handle = open(input_file)
                SeqIO.parse(handle, 'fasta')
            except TypeError:
                logger.error('Not correct FASTA format, please check.')
                sys.exit(1)
            if step == 1:
                pass
            elif step == 2:
                pass


if __name__ == '__main__':
    # Run as script
    # Parse command-line
    args = parse_cmdline()
    # Set up logging
    logger = logging.getLogger('get_ENC_GC3s_from_fasta.py: %s' % time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            log_stream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(log_stream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.info("Could not open %s for logging", args.logfile)
            logger.error(last_exception())
            sys.exit(1)
    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    if not os.path.exists(args.outdirname):
        os.makedirs(args.outdirname)
    logger.info("Output directory: %s", args.outdirname)
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))


