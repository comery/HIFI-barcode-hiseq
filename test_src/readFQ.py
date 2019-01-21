#!/usr/bin/env python
# --*-- coding:utf-8 --*--

import os
import sys
import time
import gzip
import subprocess
import argparse

t = time.time()

try:
    import Bio
except:
    sys.exit("package biopython not found! Please install it!")
else:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna


###############################################################################
#####------------------------- parameters --------------------------------#####

## common group  ##
common_parser = argparse.ArgumentParser(add_help=False)

common_group = common_parser.add_argument_group("common arguments")

common_group.add_argument(
    "-outpre",
    metavar="<STR>",
    required=True,
    help="prefix for output files",
)
## index group ##
index_parser = argparse.ArgumentParser(add_help=False)

index_group = index_parser.add_argument_group("index arguments")

index_group.add_argument(
    "-index",
    metavar="INT",
    type=int,
    required=True,
    help="the length of tag sequence in the ends of primers",
)

# Software group
soft_parser = argparse.ArgumentParser(add_help=False)

soft_group = soft_parser.add_argument_group("software path")

soft_group.add_argument(
    "-vsearch",
    metavar="<STR>",
    help="vsearch path" + "(only needed if vsearch is not in $PATH)",
)

soft_group.add_argument(
    "-threads", metavar="<INT>", default=2, help="threads for vsearch, default=2"
)

soft_group.add_argument(
    "-cid",
    metavar="FLOAT",
    type=float,
    default=0.98,
    dest="cluster_identity",
    help="identity for clustering, default=0.98",
)

## filter group  ##
filter_parser = argparse.ArgumentParser(
    add_help=False,
    description="Use the whole raw" + "dataset (Only adapters should be removed)!",
)

filter_group = filter_parser.add_argument_group("filter arguments")

filter_group.add_argument(
    "-q1",
    metavar="<STR>",
    required=True,
    dest="fastq1",
    help="input raw paired-end left fastq file, and\n"
    + "only adapters should be removed;\nsupposed on\n"
    + "Phred score system (Hiseq platform)",
)

filter_group.add_argument(
    "-q2",
    metavar="<STR>",
    required=True,
    dest="fastq2",
    help="input raw paired-end right fastq file, and\n"
    + "only adapters should be removed;\nsupposed on\n"
    + "Phred score system (Hiseq platform)",
)
filter_group.add_argument(
    "-phred",
    metavar="<INT>",
    type=int,
    dest="phred",
    default=33,
    help="Phred score system, default=64",
)
filter_group.add_argument(
    "-e",
    metavar="<INT>",
    type=int,
    dest="expected_err",
    help="expected error threshod, default=10\n"
    + "see more: http://drive5.com/usearch/manual/exp_errs.html",
)

filter_group.add_argument(
    "-q",
    metavar="<INT>",
    type=int,
    dest="quality",
    nargs=2,
    help="filter by base quality; for example: '20 5' means\n"
    + "dropping read which contains more than 5 percent of \n"
    + "quality score < 20 bases.",
)

filter_group.add_argument(
    "-trim",
    dest="trim",
    action="store_true",
    help="whether to trim 5' end of read, it adapts to -e mode\n"
    + "or -q mode",
)
filter_group.add_argument(
    "-n",
    metavar="<INT>",
    type=int,
    default=1,
    help="remove reads containing [INT] Ns, default=1",
)

# ------------------------------------------------------------------------------

## assign group ##
assign_parser = argparse.ArgumentParser(
    add_help=False,
    description="assing clean reads to +\
                                        samples by unique tag sequence with +\
                                        100% similarity",
)

assign_group = assign_parser.add_argument_group("assign arguments")

assign_group.add_argument(
    "-primer",
    metavar="<STR>",
    required=True,
    help="taged-primer list, on following format:\n"
    + "Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA\n"
    + "For001   AAGCGGTCAACAAATCATAAAGATATTGG\n"
    + "...\n"
    + "this format is necessary!",
)

assign_group.add_argument(
    "-outdir",
    metavar="<STR>",
    default="assigned",
    help="output directory for assignment," + 'default="assigned"',
)

assign_group.add_argument(
    "-tmis",
    metavar="<INT>",
    type=int,
    dest="tag_mismatch",
    default=0,
    help="mismatch number in tag when demultiplexing, default=0",
)

assign_group.add_argument(
    "-pmis",
    metavar="<INT>",
    type=int,
    dest="primer_mismatch",
    default=1,
    help="mismatch number in primer when demultiplexing, default=1",
)



## only assign need
only_assign_parser = argparse.ArgumentParser(add_help=False)
only_assign_group = only_assign_parser.add_argument_group(
    "when only run assign arguments"
)

only_assign_group.add_argument(
    "-fq", metavar="<STR>", required=True, help="cleaned fastq file "
)

# ------------------------------------------------------------------------------

## assembly group ##
assembly_parser = argparse.ArgumentParser(
    description="Due to connect HIFI barcode sequence by overlaping two"
    + "consensus sequences which generated from clustering method,"
    + "in this program I use VSEARCH.  You can define the length of overlap,"
    + "how many reads used to make clusters, and whether to check codon"
    + "translation for PCG.",
    add_help=False,
)

assembly_group = assembly_parser.add_argument_group("assembly arguments")

assembly_group.add_argument(
    "-min",
    metavar="INT",
    type=int,
    default=80,
    dest="min_overlap",
    help="minimun length of overlap, default=80",
)

assembly_group.add_argument(
    "-max",
    metavar="INT",
    type=int,
    default=90,
    dest="max_overlap",
    help="maximum length of overlap, default=90",
)

assembly_group.add_argument(
    "-oid",
    metavar="FLOAT",
    type=float,
    default=0.95,
    dest="overlap_identity",
    help="minimun similarity of overlap region, default=0.95",
)

assembly_group.add_argument(
    "-tp",
    metavar="INT",
    type=int,
    dest="cluster_number_needKeep",
    help="how many clusters will be used in" + "assembly, recommendation=2",
)

assembly_group.add_argument(
    "-ab",
    metavar="INT",
    type=int,
    dest="abundance_threshod",
    help="keep clusters to assembly if its abundance >=INT ",
)

assembly_group.add_argument(
    "-seqs_lim",
    metavar="INT",
    type=int,
    default=0,
    help="reads number limitation. by default,\n" + "no limitation for input reads",
)

assembly_group.add_argument(
    "-len",
    metavar="INT",
    type=int,
    default=400,
    dest="standard_length",
    help="standard read length, default=400",
)

assembly_group.add_argument(
    "-ds",
    dest="drop_short_read",
    action="store_true",
    help="drop short reads away before assembly",
)

assembly_group.add_argument(
    "-mode",
    metavar="INT",
    type=int,
    choices=[1, 2],
    default=1,
    help="1 or 2; modle 1 is to cluster and keep\n"
    + "most [-tp] abundance clusters, or clusters\n"
    + "abundance more than [-ab], and then make a \n"
    + "consensus sequence for each cluster.\n"
    + "modle 2 is directly to make only one consensus\n"
    + "sequence without clustering. default=1",
)

assembly_group.add_argument(
    "-rc",
    dest="reads_check",
    action="store_true",
    help="whether to check amino acid translation\n" + "for reads, default not",
)

# translation need
trans_parser = argparse.ArgumentParser(add_help=False)
trans_group = trans_parser.add_argument_group(
    "translation arguments(when set -rc or -cc)"
)

trans_group.add_argument(
    "-codon",
    metavar="INT",
    type=int,
    dest="codon_table",
    default=5,
    help="codon usage table used to check" + "translation, default=5",
)

trans_group.add_argument(
    "-frame",
    metavar="INT",
    type=int,
    choices=[0, 1, 2],
    default=1,
    help="start codon shift for amino acid" + "translation, default=1",
)


## only assembly need
only_assembly_parser = argparse.ArgumentParser(add_help=False)
only_assembly_group = only_assembly_parser.add_argument_group(
    "only run assembly arguments(not all)"
)

only_assembly_group.add_argument(
    "-list",
    metavar="FILE",
    type=str,
    required=True,
    help="input file, fastq file list. [required]",
)

polish_parser = argparse.ArgumentParser(
    description="polish all assemblies, \n"
    + "to make a confident COI barcode"
    + "reference.",
    add_help=False,
)
polish_group = polish_parser.add_argument_group(
    "polish arguments"
)

polish_group.add_argument(
    "-i",
    metavar="STR",
    type=str,
    dest="coi_input",
    required=True,
    help="COI barcode assemblies",
)

polish_group.add_argument(
    "-cc",
    dest="coi_check",
    action="store_false",
    help="whether to check final COI contig's\n"
    + "amino acid translation, default yes",
)

polish_group.add_argument(
    "-cov",
    metavar="INT",
    type=int,
    dest="min_coverage",
    default=5,
    help="minimun coverage of 5' or 3' end allowed, default=5",
)

polish_group.add_argument(
    "-l",
    metavar="INT",
    type=int,
    dest="coi_length",
    default=650,
    help="minimun length of COI barcode allowed, default=650",
)

# ------------------------------------------------------------------------------------------------

###############################################################################
#####----------------------- main subcommand parsers --------------------######

description = """

Description

    An automatic pipeline for HIFI-SE400 project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences).

Version

    1.0.2 2018-12-10  Add "-trim" function in filter;
        accept mismatches in tag or primer sequence,
        when demultiplexing; accept uneven reads to
        assembly; add "-ds" to drop short reads before
        assembly.
    1.0.1 2018-12-2  Add "polish" function
    1.0.0 2018-11-22 formated as PEP8 style
    0.0.1 2018-11-3

Author
    yangchentao at genomics.cn, BGI.
    mengguanliang at genomics.cn, BGI.
"""

parser = argparse.ArgumentParser(
    prog="HIFI-SE",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-v", "--version",
    action="version",
    version="%(prog)s 1.0.2"
)

subparsers = parser.add_subparsers(dest="command")

########## subcommmands ###########

## all subcommand
parser_all = subparsers.add_parser(
    "all",
    parents=[
        common_parser,
        index_parser,
        soft_parser,
        filter_parser,
        assign_parser,
        assembly_parser,
        trans_parser,
    ],
    formatter_class=argparse.RawTextHelpFormatter,
    help="run filter, assign and assembly",
)

## filter subcommand
parser_filter = subparsers.add_parser(
    "filter",
    parents=[common_parser, filter_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="filter raw reads",
)

## assign subcommand
parser_assign = subparsers.add_parser(
    "assign",
    parents=[common_parser,
             index_parser,
             only_assign_parser,
             assign_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="assign reads to samples",
)

## assembly subcommand
parser_assembly = subparsers.add_parser(
    "assembly",
    parents=[common_parser,
             index_parser,
             only_assembly_parser,
             soft_parser,
             assembly_parser,
             trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do assembly from input fastq\n" + "reads, output HIFI barcodes.",
)

## polish subcommand
parser_polish = subparsers.add_parser(
    "polish",
    parents=[polish_parser,
            index_parser,
            trans_parser,],
    formatter_class=argparse.RawTextHelpFormatter,
    help="polish COI barcode assemblies,\n"
    + "output confident barcodes."
)

## BOLD_identification
parser_bold = subparsers.add_parser(
    "bold_identification",
    parents=[],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do taxa identification\n" + "on BOLD system,\n",
)

###############################################################################
#####---------------------- program execution start ----------------------#####

args = parser.parse_args()
# -----------------------BOLD identification----------------------#
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

if sys.argv[1] == "bold_identification":
    # if args.command == 'bold_identification':
    sys.argv = sys.argv[1:]
    sys.exit(bold_identification())


# -----------------------arguments checking-----------------------#
## softwares and databases
def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return 0
    else:
        print(cmd + " not found!", file=sys.stderr)
        return 1


def files_exist_0_or_1(filelist):
    '''
    check files involed whether are existing!
    '''
    num = 0
    for file in filelist:
        if os.path.exists(file):
            num += 1
        else:
            print("%s doesn't exist!" % file, file=sys.stderr)
    if len(filelist) == num:
        return 0
    else:
        return 1


def check_and_open(file):
    if os.path.exists(file):
        print("WARRNING: " + file + " exists! now overwriting")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

# ----------------------------------------------------------------
## file existing check
errors_found = 0
if args.command == "all":
    errors_found += files_exist_0_or_1(
        [args.fastq2, args.fastq2, args.primer]
    )
elif args.command == "filter":
    errors_found += files_exist_0_or_1([args.fastq1, args.fastq2])
elif args.command == "assign":
    errors_found += files_exist_0_or_1([args.primer])
elif args.command == "assembly":
    errors_found += files_exist_0_or_1([args.list])
elif args.command == "polish":
    errors_found += files_exist_0_or_1([args.coi_input])
else:
    parser.print_help()
    parser.exit()

if args.command in ["all", "assembly"]:
    vsearch = "vsearch"
    if hasattr(args, "vsearch"):
        if args.vsearch:
            vsearch = args.vsearch
    errors_found += check_program_involed(vsearch)

if errors_found > 0:
    parser.exit("Errors found! Exit!")


if hasattr(args, "outpre") and args.outpre.endswith("/"):
    print("outpre is in bad format! no \"/\"")
    exit()

# -----------------------functions for filtering------------------#

def deleteDuplicatedElementFromList3(listA):
        #return list(set(listA))
        return sorted(set(listA), key = listA.index)

def parse_pe_fastq(fq1, fq2, phred=64):
    from numpy import fromstring,byte
    while True:
        name1 = fq1.readline().split(" ")[0]
        name2 = fq2.readline().split(" ")[0]
        if not name1 or not name2:
            break
        read1, nothing1, qual1 = fq1.readline()[:-1], fq1.readline(), fq1.readline()[:-1]
        read2, nothing2, qual2 = fq2.readline()[:-1], fq2.readline(), fq2.readline()[:-1]
        qual1 = fromstring(qual1, dtype=byte) - phred
        print(qual1)
        qual2 = fromstring(qual2, dtype=byte) - phred
        assert name1 == name2, 'fastq1, fastq2 is not paired-end'
        yield name1, read1, read2, qual1, qual2


def exp_e(q, phred):
    '''
    expected error number(E*) = sum(P), where P is
    the probability that the base call is incorrect.
    '''
    exp = 0
    tmp = list(q)
    ascill = [ord(n) - phred for n in tmp]

    for i in ascill:
        exp += 10 ** (-i / 10)

    return exp

def lowquality_rate(qual_str, cut_off, phred):
    '''
    calculate the rate of low quality base in read.
    '''
    low_base = 0
    tmp = list(qual_str)
    ascill = [ord(n) - phred for n in tmp]

    for i in ascill:
        if i < cut_off:
            low_base += 1

    low_rate = low_base / len(qual_str)
    return low_rate


# ----------------------functions for assigning-------------------#
def complementation(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper() ## [a bug fixed], reported by Wu Ping 20181129
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence

def comp_rev(sequence):
    # make a sequence complement and reversed #
    sequence = complementation(sequence)
    return sequence[::-1]


def detect_mis(f, r, dict):
    tag_mis = 0
    primer_mis = 0
    strs = [f, r]
    while(strs):
        s1 = strs.pop()
        for s2 in dict.keys():
            if len(s1) == len(s2):
                tag_mis = 0
                primer_mis = 0
                for base in range(len(s1)):
                    if s1[base] is not s2[base]:
                        if base < args.index:
                            tag_mis += 1
                        else:
                            primer_mis += 1
                if (tag_mis <= args.tag_mismatch
                    and primer_mis <= args.primer_mismatch):
                    goal = dict[s2]
                    break
                else:
                    goal = ''
    
    if len(goal) > 0:
        return goal
    else:
        return False


def dis_barcode(barcode_list):
    # ----count matched bases of two barcodes----#
    dis = []
    while(barcode_list):
        b1 = barcode_list.pop()
        for b2 in barcode_list:
            mismatch = 0
            for base in range(len(b1)):
                if b1[base] is not b2[base]:
                    mismatch += 1
            dis.append(mismatch)
    min_dis = min(dis)
    max_dis = max(dis)
    return (min_dis, max_dis)


# --------------------------------------------------------
if args.command in ["all", "filter"]:

    filtered_outfile1 = args.outpre + "_filter_highqual_1.fastq"
    filtered_outfile2 = args.outpre + "_filter_highqual_2.fastq"
    out1 = check_and_open(filtered_outfile1)
    out2 = check_and_open(filtered_outfile2)

    # ini state
    total = 0
    clean = 0
    nn = 0
    low = 0

    log = check_and_open(args.outpre + "_filter_log.txt")

    if args.expected_err and args.quality:
        print(
            "Bad arguments:\n\t"
            + " -e argument is confilicting with -q,"
            + " can not using in the same time"
        )
        exit()
    elif args.quality:
        high_qual = args.quality[0]
        low_qual_cont = args.quality[1] / 100
        filter_type = 2
        log.write("Filtering by quality score: {}, content:"
                  + " {}%".format(args.quality, args.quality[1])
                  + "\n")
    else:
        filter_type = 1
        if not args.expected_err:
            args.expected_err = 10
        log.write("Filtering by expected_err: {}".format(args.expected_err) + "\n")

    print_time("Filtering start:")

    if args.fastq1.endswith("gz"):
        fq1 = gzip.open(args.fastq1, 'rt')
    else:
        fq1 = open(args.fastq1, 'r')
    if args.fastq2.endswith("gz"):
        fq2 = gzip.open(args.fastq2, 'rt')
    else:
        fq2 = open(args.fastq2, 'r')

    for i in parse_pe_fastq(fq1, fq2, args.phred):
        name, seq1, seq2, qua1, qua2 = i

        if not name[0].startswith("@"):
            print("ERROR: input fastq is not a correct fastq format")
            exit()

        total += 1
        N_count1 = seq1.count("N")
        N_count2 = seq2.count("N")

        if N_count1 < args.n or N_count2 < args.n:
            if filter_type == 1:
                if (exp_e(qual1, args.phred) <= args.expected_err
                    and exp_e(qual2, args.phred) <= args.expected_err):
                    out1.write(name + "-1\n" + seq1 + "\n" + "+\n" + qual1 + "\n")
                    out2.write(name + "-2\n" + seq2 + "\n" + "+\n" + qual2 + "\n")
                    clean += 1
                else:
                    low += 1
            else:
                # filter_type == 2
                if (lowquality_rate(qual1, high_qual, args.phred) > low_qual_cont
                    and lowquality_rate(qual2, high_qual, args.phred) > low_qual_cont):
                    out1.write(name + "-1\n" + seq1 + "\n" + "+\n" + qual1 + "\n")
                    out2.write(name + "-2\n" + seq2 + "\n" + "+\n" + qual2 + "\n")

                    clean += 1
                else:
                    low += 1
        else:
            nn += 1

    log.write("total reads:\t{}".format(total) + "\n")
    log.write("clean reads:\t{}".format(clean) + "\n")
    log.write("low quality reads:\t{}".format(low) + "\n")
    log.write("containing N reads:\t{}".format(nn) + "\n")

    log.close()
    fq1.close()
    fq2.close()
    out1.close()
    out2.close()

    print_time("Filtering done:")
