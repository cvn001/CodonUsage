#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to calculate ENC (Nc) and GC3s from each genome.
# Created by galaxy on 2017/3/2 0002 21:38

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


def run_codonw(each_type, fasta_file, out_file):
    devnull = open(os.devnull, 'w')
    cmd = '{0} {1} {2} -noblk -nomenu -nowarn -silent -human -fasta -enc'.format(each_type, fasta_file, out_file)
    try:
        subprocess.call(cmd, shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to run codonw but failed, please check the dependency.'.format(args.part))
        logger.error(last_exception())
        sys.exit(1)


def multi_process(each_type):
    run_list = []
    type_dir = os.path.join(args.outdirname, each_type)
    for root, dirs, files in os.walk(args.indirname):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            fasta_file = os.path.join(args.indirname, each_file)
            out_file = os.path.join(type_dir, strain_id + '.out')
            run_list.append([fasta_file, out_file])
    p = Pool(args.threads)
    for each_item in run_list:
        p.apply_async(run_codonw, args=(each_type, each_item[0], each_item[1]))
    p.close()
    p.join()


def main():
    type_list = ['gc3s']
    for each_type in type_list:
        logger.info('Processing {0} Calculating'.join(each_type.upper()))
        multi_process(each_type)


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
            logger.error("Could not open %s for logging", args.logfile)
            sys.exit(1)
    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    if not os.path.exists(args.outdirname):
        os.makedirs(args.outdirname)
    else:
        shutil.rmtree(args.outdirname)
        os.makedirs(args.outdirname)
    logger.info("Output directory: %s", args.outdirname)
    main()
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
