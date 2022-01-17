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

__version__ = "1.2.0"
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


def process_geneid(file):

    data = OrderedDict()
    r = {}

    for line in read_tsv(file, "\t"):
        geneid = line[0]
        sample, contig, position = geneid.split(".")
        did = "%s.%s" % (sample, contig)
        if did not in data:
            data[did] = []
        data[did].append(int(position))
        gene = "PUL"
        if line[1] != ".":
            gene = line[1]
        r[geneid] = [gene, line[-1]]

    return data, r


def find_gene_family(glist, gaplen=2, minegene=2, gaps=1):

    r = []
    temp = []
    gs = 0

    for i in sorted(glist):
        if not temp:
            temp.append(i)
            continue
        if i <= (max(temp)+1):
            temp.append(i)
        elif i <= (max(temp)+gaplen+1):
            gs += 1
            if gs >= (gaps+1):
                if len(temp) >= minegene:
                    r.append(temp)
                gs = 0
                temp = [i]
            else:
                temp.append(i)
        else:
            if len(temp) >= minegene:
                r.append(temp)
            gs = 0
            temp = [i]
    if len(temp) >= minegene:
        r.append(temp)

    return r


def pul_annotation(seqid, start, end, gene_dict):

    strus = []
    funcs = []

    for i in range(start, end+1, 1):
        geneid = "%s.%s" % (seqid, i)
        if geneid not in gene_dict:
            strus.append("_")
            continue
        stru, func = gene_dict[geneid]
        strus.append(stru)
        if func not in ["degradation", "biosynthesis"]:
            continue
        funcs.append(func)
    if funcs:
        funcs = max(funcs, key=funcs.count)
    else:
        funcs = ""

    return "|".join(strus), funcs



def find_pul(file, gaplen=1, minegene=2, gaps=1, anrate=75):

    data, gene_dict = process_geneid(file)

    print("#Seq_id\tPULs\tStart\tEnd\tGene Number\tAnnotation Genes\tAnnotation rate(%)\tPul structure\tFunction")
    n = 0
    for seqid in data:
        for temp in find_gene_family(data[seqid], gaplen, minegene, gaps):
            start, end = min(temp), max(temp)
            n += 1
            stru, func = pul_annotation(seqid, start, end, gene_dict)
            gene_anrate = len(temp)*100.0/(end-start+1)
            if gene_anrate <= anrate:
                continue
            if ("GH" not in stru) and ("CE" not in stru) and ("Sus" not in stru):
                continue
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(seqid, n, start, end,
                end-start+1, len(temp), gene_anrate, stru, func
                )
            )

    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input pul annotation result file(pul.tsv).')
    parser.add_argument('-gl', '--gaplen', metavar='INT', type=int, default=4,
        help="Maximum allowable gap size, default=4")
    parser.add_argument("-mg", '--minegene', metavar='INT', type=int, default=2,
        help="Minimum number of genes allowed, default=2")
    parser.add_argument('-g', '--gaps', metavar='INT', type=int, default=4,
        help="Minimum number of gaps allowed, default=4")
    parser.add_argument('-ar', '--anrate', metavar='FLOAT', type=float, default=60.0,
        help="Minimum genome annotation rate to be met, default=60.0")

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
    find_pul.py -- Polysaccharide utilization site prediction

attention:

    find_pul.py pul.tsv >stat_pul.tsv
optional arguments:
  --gaplen INT             允许插入其他基因的最大数目, default=4
  -mg INT, --minegene INT
                        允许最小的中靶基因组数, default=2
  --gaps INT            允许最小的开口数目, default=4
  -ar FLOAT, --anrate FLOAT
                        注释上的基因占比, default=60.0

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    find_pul(args.input, args.gaplen, args.minegene, args.gaps, args.anrate)


if __name__ == "__main__":
    main()
