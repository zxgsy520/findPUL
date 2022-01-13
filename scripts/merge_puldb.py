#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []



def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_sus(file):

    r = {}

    for line in read_tsv(file, "\t"):
        r[line[0]] = line[-1]

    return r


def read_cazy(file):

    r = {}

    for line in read_tsv(file, "\t"):
        type = list(set(line[3].split(";")))
        if ("GH" not in type) and ("CE" not in type):
            continue
        r[line[0]] = ";".join(type)

    return r


def read_pul(file):

    r = {}

    for line in read_tsv(file, "\t"):
        r[line[0]] = [line[2], line[4]]

    return r


def sorted_seqid(seqids):

    r = {}

    for i in set(seqids):
        i = i.split(".")
        seqid = ".".join(i[0:-1])
        if seqid not in r:
            r[seqid] = []
        r[seqid].append(int(i[-1]))

    ids = []
    for line in sorted(r.items(), key=lambda x:x[0]):
        for j in sorted(line[1]):
            ids.append("%s.%s" % (line[0], j))

    return ids
     

def merge_puldb(susfile, cazyfile, pulfile):

    dsus =  read_sus(susfile)
    dcazy = read_cazy(cazyfile)
    dpul = read_pul(pulfile)

    seqids = dsus.keys() + dcazy.keys() + dpul.keys()
    
    print("#Seqid\tGene family\tPulid\tDegradation/Biosynthesis")
    for i in sorted_seqid(seqids):
        family = "."
        if i in dsus:
            family = dsus[i]
        if i in dcazy:
            if "." != family:
                family += ";%s" % dcazy[i]
            else:
                family = dcazy[i]

        pulid = "."
        des = ""
        if i in dpul:
            pulid, des = dpul[i]
        print("%s\t%s\t%s\t%s" % (i, family.strip(";"), pulid, des))

    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input pul annotation result file(pul.tsv).')
    parser.add_argument('--cazy', metavar='FILE', type=str, required=True,
        help="Input cazy annotation result file(cazy.tsv)")
    parser.add_argument("--sus", metavar='FILE', type=str, required=True,
        help="Input gene annotation result file(sus.tsv)")

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
    merge_puldb.py --Merged polysaccharide utilization site database

attention:

    merge_puldb.py pul.tsv --cazy cazy.tsv --sus sus.tsv >merge_puldb.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_puldb(args.sus, args.cazy, args.input)


if __name__ == "__main__":
    main()
