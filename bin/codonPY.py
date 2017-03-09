#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: codonPY main script
# Created by galaxy on 17-3-7 10:08 pm

import logging.handlers
import os
import sys
import time
import traceback
from argparse import ArgumentParser
from Bio import SeqIO
from src.calculate_relative_adaptiveness import calculate_wi
from src.concatenate_fasta import concatenate
from src.get_CAI_from_fasta import calculate_cai
from src.get_ENC_GC3s_from_fasta import get_enc_gc3s
from src.get_RSCU_from_fasta import get_rscu
from src.integrate_results import integrate_enc


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
    parser = ArgumentParser(prog="src.py")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None,
                        type=str, help="Input directory name")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None,
                        type=str, help="Output directory")
    parser.add_argument("-p", "--step", type=int, dest="step", default=0,
                        help="Which part will be run? [0|1|2|3|4]")
    # parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
    #                     help="How many threads will be used? [default all]")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None,
                        type=str, help="Logfile location")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument('-c', '--col', required=False, default=3, dest='col', type=int,
                        help='The column number where the wi values are in the relative adaptiveness table. '
                             'Codons are assumed to be in the first column. Default = 3rd column.')
    parser.add_argument('-e', '--exclude_cpg_codons', required=False, default=None, dest='exclude_cpg', type=str,
                        help='Argument for whether to exclude CpG codons in the calculation of CAI. '
                             'Default is to include all codons. '
                             'Change to "yes" to exclude CpG codons from the calculation. '
                             'Change to "only" to base the CAI measures only on CpG codons.')
    parser.add_argument('-n', '--interest', required=False, default='None', dest='interest', type=str,
                        help='Enter a codon or amino acid (three letter notation) you are interested '
                             'in to view its count data individually in standard output.')
    parser.add_argument('-s', '--sub', required=False, default="None", dest='sub', type=str,
                        help='A subset of sequence names to pull and concatenate')
    return parser.parse_args()


def each_step(step, input_file, fname):
    if step == 1:
        logger.info('Step 1: get relative synonymous codon usage (RSCU) from {0}'.format(input_file))
        step_dir = os.path.join(args.outdirname, 'RSCU')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        output_file = os.path.join(step_dir, fname + '.rscu')
        get_rscu(input_file, output_file, args.interest)
        logger.info('Step 1 done')
    elif step == 2:
        logger.info('Step 2: calculate relative adaptiveness (W) from {0}'.format(input_file))
        step_dir = os.path.join(args.outdirname, 'Relative_Adaptiveness')
        rscu_dir = os.path.join(args.outdirname, 'RSCU')
        rscu_file = os.path.join(rscu_dir, fname + '.rscu')
        if not os.path.exists(rscu_file):
            logger.warning('The essential wi values files are not exist, so go back step 2.')
            each_step(1, input_file, fname)
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        output_file = os.path.join(step_dir, fname + '.wi')
        calculate_wi(rscu_file, output_file)
    elif step == 3:
        logger.info('Step 3: get Codon Adaptation Index (CAI) from {0}'.format(input_file))
        step_dir = os.path.join(args.outdirname, 'CAI')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        wi_dir = os.path.join(args.outdirname, 'Relative_Adaptiveness')
        wi_file = os.path.join(wi_dir, fname + '.wi')
        if not os.path.exists(wi_file):
            logger.warning('The essential wi values files are not exist, so go back step 2.')
            each_step(2, input_file, fname)
        output_file = os.path.join(step_dir, fname + '.cai')
        calculate_cai(input_file, output_file, wi_file, args.col, args.exclude_cpg)
    elif step == 4:
        logger.info('Step 4: get expected ENC and GC3s from {0}'.format(input_file))
        step_dir = os.path.join(args.outdirname, 'ENC_GC3s')
        if not os.path.exists(step_dir):
            os.makedirs(step_dir)
        output_file = os.path.join(step_dir, fname + '.enc_gc3s')
        get_enc_gc3s(input_file, output_file)
        message = integrate_enc(args.outdirname)
        logger.info(message)


def separate_run(step, input_dir):
    logger.info('Run in step: {0}'.format(str(step)))
    for rt, ds, fs in os.walk(input_dir):
        for f in fs:
            fname = os.path.splitext(f)[0]
            input_file = os.path.join(input_dir, f)
            each_step(step, input_file, fname)


def auto_run(input_dir):
    logger.info('Run all steps automatically.')
    step_list = [1, 2, 3, 4]
    for i in step_list:
        separate_run(i, input_dir)


# def multi_process():
#     p = Pool(args.threads)
#     out_put = os.path.join()
#     return


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
    if not os.path.exists(args.indirname) or not os.path.getsize(args.indirname):
        logger.info('Input directory {0} is missing or empty. Please check.'.format(args.indirname))
        logger.error(last_exception())
        sys.exit(1)
    if not os.path.exists(args.outdirname):
        os.makedirs(args.outdirname)
    logger.info("Output directory: %s", args.outdirname)
    concatenate_dir = os.path.join(args.outdirname, 'concatenated_fasta')
    if not os.path.exists(concatenate_dir):
        os.makedirs(concatenate_dir)
    if args.sub == "None":
        logger.info("\nNo subset file detected. Concatenating the entire fasta")
    for root, dirs, files in os.walk(args.indirname):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            in_file = os.path.join(args.indirname, each_file)
            try:
                handle = open(in_file)
                SeqIO.parse(handle, 'fasta')
            except TypeError:
                logger.info('Not correct FASTA format, please check.')
                logger.error(last_exception())
                sys.exit(1)
            out_file = os.path.join(concatenate_dir, strain_id + '.fasta')
            if not os.path.exists(out_file):
                concatenate(in_file, out_file, args.sub)
    if args.step == 0:
        auto_run(concatenate_dir)
    elif 1 <= args.step <= 4:
        separate_run(args.step, concatenate_dir)
    else:
        logger.info('You picked a wrong step [0|1|2|3|4]. Stop.')
        logger.error(last_exception())
        sys.exit(1)
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
