#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line.upper()
        else:
            seq.append(line.upper())

    if len(seq) == 2:
        yield seq
    fp.close()


def deal_with_pfam(files, describe="SusC"):

    seqids = set()

    for file in files:
        for seqid,seq in read_fasta(file):
            if seqid in seqids:
                continue
            seqids.add(seqid)
            seqid = "%s#%s" % (seqid.replace("/","|"), describe)
            print(">%s\n%s" % (seqid, seq))

    return 0


def add_hlep_args(parser):

    parser.add_argument('input', nargs='+', metavar='FILE', type=str,
        help='Input protein sequence file, format(fasta, fa.gz')
    parser.add_argument('-d', '--describe', metavar='STR', type=str, required=True,
        help='Input gene family name')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    deal_with_pfam: Process protein sequence and add it to pfam database.
attention:
    deal_with_pfam PF00953_full.txt -d SusC >SusC.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    deal_with_pfam(args.input, args.describe)


if __name__ == "__main__":

    main()
