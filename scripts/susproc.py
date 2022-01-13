#!/usr/bin/env python
# -*- coding: utf-8 -*-



import re
import sys
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def suslproc(file):

    print("#qseqid\tsseqid\tGene family")
    for line in read_tsv(file, "\t"):
        sseqid, family = line[1].split("#", 1)
        print("%s\t%s\t%s" % (line[0], sseqid, family))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar='FILE', type=str,
        help="Input annotation result file(sus.out).")

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
    suslproc.py: Annotate statistical SUS analysis results

attention:
    suslproc.py  sus.out >sus.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()
    suslproc(args.input)


if __name__ == "__main__":

    main()
