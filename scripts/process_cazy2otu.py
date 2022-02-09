#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import math
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def generate_list(length):

    r = []

    for i in range(length):
        r.append(0)

    return r


def read_stat_cazy(files):

    samples = []
    data = {}
    genes = []
    n = 0

    for file in files:
        n += 1
        sample = file.split("/")[-1].split(".")[0]
        samples.append(sample)
        
        gene = 0
        for line in read_tsv(file):
            gene += int(line[1])
            if line[0] not in data:
                data[line[0]] = generate_list(n-1)
            data[line[0]].append(int(line[1]))
        genes.append(gene)

    return samples, genes, data


def process_unite(genes, data):
    
    r = {}
    meangene = sum(genes)*1.0/len(genes)

    for i in data:
        temp = []
        k = 0
        for j in data[i]:
            temp.append(j*1.0*meangene/genes[k])
            k += 1
        r[i] = temp

    return r
       

def process_zscore(otulist):

    total = sum(otulist)
    meanotu = total*1.0/len(otulist)
    sd = 0

    for i in otulist:
        sd += math.pow((i-meanotu), 2)
    sd = math.sqrt(sd/len(otulist))

    r = []
    for i in otulist:
        r.append("{:.6f}".format((i-meanotu)/sd))

    return r


def process_cazy2otu(files, maxrow=20):

    samples, genes, data = read_stat_cazy(files)
    #print(genes)
    data = process_unite(genes, data)
    n = 0

    print("#ID\t%s" % "\t".join(samples))
    for otuid, line in sorted(data.items(), key=lambda x:min(x[1]), reverse=True):
        n += 1
        line = process_zscore(line)
        print("%s\t%s" % (otuid, "\t".join(line)))
        if n >= maxrow:
            break

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", nargs='+', metavar="FILE", type=str,
        help="Input CAZy database classification statistics, stat_cazy.tsv")
    parser.add_argument("-maxr", "--maxrow", metavar="INT", type=int, default=20,
        help="Output the top abundance number,  default=20.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
For exmple:
        python process_cazy2otu.py *.stat_cazy.tsv >cazy_otu.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    process_cazy2otu(args.input, args.maxrow)


if __name__ == "__main__":

    main()
