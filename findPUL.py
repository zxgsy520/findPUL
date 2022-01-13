#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "1.0.0"
__all__ = []


from dagflow import Task, ParallelTask, DAG, do_dag
from common import check_path, mkdir
from seqkit.split import seq_split
QUEUE = "-q all.q,s01"
ROOT = "/Work/user/zhangxg/pipeline/findPUL"
SCRIPTS = os.path.join(ROOT, "scripts")
DIAMOND_BIN = "/Work/pipeline/software/Base/diamond/v2.0.3/"
BLAST_BIN = "/Work/pipeline/software/Base/blast+/bin/"

CAZY = "/Work/database/CAZy/v10/"
CAZY_DB = os.path.join(CAZY, "CAZy.dmnd")
CAZY_ACTIV = os.path.join(CAZY, "CAZy.activities.txt")
CAZY_SUBFAM = os.path.join(CAZY, "CAZy.subfam.txt")

PHI_DB = "/Work/database/phi/v4-12/phi"
PUL_DB = "/Work/database/dbCAN-PUL/v202010/PUL"
GDB = "/Work/user/zhangxg/pipeline/findPUL/database/SUS"


def create_cazy_task(proteins, prefix, evalue, coverage, threads, job_type,
                     work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "cazy"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={diamond}:$PATH
time diamond blastp --query {{proteins}} --db {db} \\
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \\
--max-target-seqs 5 --evalue 1e-05 --threads {threads} --out {{prefixs}}.CAZy.m6
""".format(diamond=DIAMOND_BIN,
           db=CAZY_DB,
           evalue=evalue,
           coverage=coverage,
           threads=threads),
        proteins=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_CAZy",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.CAZy.m6 > {prefix}.CAZy.m6
time {script}/blast_filter.py {prefix}.CAZy.m6 \\
  --outfmt std qlen slen stitle --out qseqid sseqid qstart qend evalue bitscore stitle \\
  --min_qcov {coverage} --min_scov 0 --evalue {evalue} --best >{prefix}.CAZy.out
{script}/cazyproc.py {prefix}.CAZy.out --activ {activ} \\
  --subfam {subfam} -o {prefix}.cazy_classify.tsv >{prefix}.cazy.tsv
{script}/plot_cazy.py {prefix}.cazy_classify.tsv -p {prefix}
cp {prefix}.CAZy.m6 {prefix}.CAZy.out {prefix}.cazy.tsv {out_dir}
cp {prefix}.cazy_classify.tsv {prefix}.cazy.png {prefix}.cazy.pdf {out_dir}
""".format(id=id,
           prefix=prefix,
           activ=CAZY_ACTIV,
           subfam=CAZY_SUBFAM,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, os.path.join(work_dir, "%s.cazy.tsv" % prefix)


def create_pul_task(proteins, prefix, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "pul"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={blast}:$PATH
time blastp -query {{proteins}} -db {db} \\
-outfmt '6 std qlen slen stitle' \\
-max_target_seqs 5 -evalue 1e-05 -num_threads {threads} -out {{prefixs}}.pul.m6 \
""".format(blast=BLAST_BIN,
           db=PUL_DB,
           evalue=evalue,
           coverage=coverage,
           threads=threads),
        proteins=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_pul",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.pul.m6 > {prefix}.pul.m6
time {script}/blast_filter.py {prefix}.pul.m6 \\
  --outfmt std qlen slen stitle --out qseqid sseqid qstart qend stitle evalue bitscore \\
  --min_qcov {coverage} --min_scov 0 --evalue {evalue} --best >{prefix}.pul.out
{script}/pulproc.py {prefix}.pul.out \\
  -d {db}.txt >{prefix}.pul.tsv 2>{prefix}.stat_pul.tsv
cp {prefix}.pul.m6 {prefix}.pul.out {prefix}.pul.tsv {prefix}.stat_pul.tsv {out_dir}
""".format(id=id,
           db=PUL_DB,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def create_sus_task(proteins, prefix, evalue, coverage, threads, job_type,
                    work_dir="", out_dir=""):

    prefixs = [os.path.basename(i) for i in proteins]
    id = "sus"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={blast}:$PATH
time blastp -query {{proteins}} -db {db} \\
-outfmt '6 std qlen slen stitle' \\
-max_target_seqs 5 -evalue 1e-05 -num_threads {threads} -out {{prefixs}}.sus.m6 \
""".format(blast=BLAST_BIN,
           db=GDB,
           evalue=evalue,
           coverage=coverage,
           threads=threads),
        proteins=proteins,
        prefixs=prefixs,
    )

    join_task = Task(
        id="merge_pul",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/*.pul.m6 > {prefix}.sus.m6
time {script}/blast_filter.py {prefix}.sus.m6 \\
  --outfmt std qlen slen stitle --out qseqid sseqid qstart qend stitle evalue bitscore \\
  --min_qcov {coverage} --min_scov 0 --evalue {evalue} --best >{prefix}.sus.out
{script}/pulproc.py {prefix}.sus.out \\
  -d {db}.txt >{prefix}.pul.tsv 2>{prefix}.stat_pul.tsv
cp {prefix}.pul.m6 {prefix}.pul.out {prefix}.pul.tsv {prefix}.stat_pul.tsv {out_dir}
""".format(id=id,
           db=PUL_DB,
           prefix=prefix,
           script=SCRIPTS,
           evalue=evalue,
           coverage=coverage,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def run_findpul(protein, prefix, evalue, coverage, threads,
                 job_type, concurrent, refresh, work_dir="", out_dir=""):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    protein = check_path(protein)

    work_dict = {
        "split": "00_data",
        "cazy": "01_CAZy",
        "pul": "02_PUL",
        "sus": "03_SUS",
    }
    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    proteins = seq_split([protein], mode="length", num=5000000, output_dir=work_dict["split"])

    dag = DAG("run_findpul")
    cazy_tasks, cazy_join, cazy_stat = create_cazy_task(
        proteins=proteins,
        prefix=prefix,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["cazy"],
        out_dir=out_dir)
    dag.add_task(*cazy_tasks)
    dag.add_task(cazy_join)

    pul_tasks, pul_join = create_pul_task(
        proteins=proteins,
        prefix=prefix,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["pul"],
        out_dir=out_dir
    )
    dag.add_task(*pul_tasks)
    dag.add_task(pul_join)

    sus_tasks, sus_join = create_sus_task(
        proteins=proteins,
        prefix=prefix,
        evalue=evalue,
        coverage=coverage,
        threads=threads,
        job_type=job_type,
        work_dir=work_dict["sus"],
        out_dir=out_dir
    )
    dag.add_task(*sus_tasks)
    dag.add_task(sus_join)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return 0


def add_hlep_args(parser):

    parser.add_argument("protein", metavar='FILE', type=str,
        help="Input protein sequence.")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="out",
        help="Input sample name.")
    parser.add_argument("-e", "--evalue", metavar="NUM", type=float, default=1e-05,
        help="Evalue cutoff of blastp for Refseq and KEGG, default=1e-05.")
    parser.add_argument("-c", "--coverage", metavar="NUM", type=float, default=30,
        help="Coverage cutoff of blastp for Refseq and KEGG,default=30.")
    parser.add_argument("-t", "--thread", metavar='INT', type=int, default=4,
        help="Set the running thread, default=4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

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
    run_findpul.py :Polysaccharide Utilization Loci Prediction
attention:
    run_findpul.py protein.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    run_findpul(args.protein, args.prefix, args.evalue, args.coverage,
                 args.thread, args.job_type, args.concurrent,
                 args.refresh, args.work_dir, args.out_dir)


if __name__ == "__main__":

    main()
