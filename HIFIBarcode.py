#!/usr/bin/env python
# --*-- coding:utf-8 --*--

import os
import sys
import time
import gzip
import subprocess
import argparse
from icecream import ic

from bold_identification.BOLD_identification import (
        main as bold_identification,
)

t = time.time()

try:
    import Bio
except:
    sys.exit("package biopython not found! Please install it!")
else:
    from Bio.Seq import Seq



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
    "-threads",
    metavar="<INT>",
    default=2,
    help="threads for vsearch, default=2",
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
    help="input raw paired-end left fastq file, and only adapters \n"
    + "should be removed; supposed on Phred score system (Hiseq platform)",
)

filter_group.add_argument(
    "-q2",
    metavar="<STR>",
    required=True,
    dest="fastq2",
    help="input raw paired-end right fastq file"
)
filter_group.add_argument(
    "-phred",
    metavar="<INT>",
    type=int,
    dest="phred",
    default=33,
    help="Phred score system, default=33",
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
    help="filter by base quality; for example: '20 5' means dropping\n"
    + "read which contains more than 5 percent of quality score<20\n"
    + "bases.",
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
    description="assing clean reads to" +\
    "samples by unique tag sequence with" +\
    "100% similarity",
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
    "-tmis",
    metavar="<INT>",
    type=int,
    dest="tag_mismatch",
    default=1,
    help="mismatch number in tag when demultiplexing, default=1",
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
    "-fq1", metavar="<STR>", required=True, help="cleaned fastq1 file "
)
only_assign_group.add_argument(
    "-fq2", metavar="<STR>", required=True, help="cleaned fastq2 file "
)

# ------------------------------------------------------------------------------

## buildend group -------------------------------------------------------
buildend_parser = argparse.ArgumentParser(
    description="Due to connect HIFI barcode sequence by overlaping two"
    + "consensus sequences which generated from clustering method,"
    + "in this program I use VSEARCH.  You can define the length of overlap,"
    + "how many reads used to make clusters, and whether to check codon"
    + "translation for PCG.",
    add_help=False,
)

buildend_group = buildend_parser.add_argument_group("buildend arguments")

buildend_group.add_argument(
    "-min",
    metavar="INT",
    type=int,
    default=10,
    dest="min_overlap",
    help="minimun length of overlap, default=10",
)

buildend_group.add_argument(
    "-max",
    metavar="INT",
    type=int,
    default=150,
    dest="max_overlap",
    help="maximum length of overlap, default=150",
)

buildend_group.add_argument(
    "-tp",
    metavar="INT",
    type=int,
    dest="cluster_number_needKeep",
    help="how many clusters will be used in" + "buildend, recommendation=2",
)

buildend_group.add_argument(
    "-ab",
    metavar="INT",
    type=int,
    dest="abundance_threshod",
    help="keep clusters to buildend if its abundance >=INT ",
)

buildend_group.add_argument(
    "-seqs_lim",
    metavar="INT",
    type=int,
    default=0,
    help="reads number limitation. by default, no limitation",
)

buildend_group.add_argument(
    "-len",
    metavar="INT",
    type=int,
    default=150,
    dest="standard_length",
    help="standard read length, default=150",
)

buildend_group.add_argument(
    "-ds",
    dest="drop_short_read",
    action="store_true",
    help="drop short reads away before buildend",
)

buildend_group.add_argument(
    "-rc",
    dest="reads_check",
    action="store_true",
    help="whether to check amino acid translation\n" + "for reads, default not",
)

buildend_group.add_argument(
    "-bom",
    metavar="INT",
    type=int,
    dest="buildend_mismatch",
    default=1,
    help="mismatches allow in overlap",
)

## only buildend need 
only_buildend_parser = argparse.ArgumentParser(add_help=False)
only_buildend_group = only_buildend_parser.add_argument_group(
    "only run buildend arguments(not all)"
)

only_buildend_group.add_argument(
    "-list",
    metavar="FILE",
    type=str,
    required=True,
    help="input file, fastq file list. [required]",
)

# chain group ---------------------------------------------
chain_parser = argparse.ArgumentParser(add_help=False)
chain_group = chain_parser.add_argument_group(
    "chain arguments")

chain_group.add_argument(
    "-mi",
    metavar="INT",
    type=int,
    dest="min_insertsize",
    default=180,
    help="minimun length of connected fragment",
)

chain_group.add_argument(
    "-com",
    metavar="INT",
    type=int,
    dest="chain_mismatch",
    default=1,
    help="mismatches allow in overlap",
)

## only chain need
only_chain_parser = argparse.ArgumentParser(add_help=False)
only_chain_group = only_chain_parser.add_argument_group(
    "when only use chain command"
)

only_chain_group.add_argument(
    "-ms",
    metavar="STR",
    type=str,
    dest="middle_input",
    required=True,
    help="middle fasta in ssam format",
)


## ----------------gapfill for all mode --------------------
gapfill_parser = argparse.ArgumentParser(add_help=False)
gapfill_group = gapfill_parser.add_argument_group(
    "gapfill arguments when run all command"
)
gapfill_group.add_argument(
    "-gmax",
    metavar="INT",
    type=int,
    dest="gmax",
    default=700,
    help="maximal length of sequence assembly(700)"
)

gapfill_group.add_argument(
    "-gmin",
    metavar="INT",
    type=int,
    dest="gmin",
    default=450,
    help="minimal length of sequence assembly(450)"
)

gapfill_group.add_argument(
    "-kl",
    metavar="INT",
    type=int,
    dest="lowkmer",
    default=55,
    help="kmer lower limit(55)"
)

gapfill_group.add_argument(
    "-i",
    metavar="INT",
    type=int,
    dest="kmer_interval",
    default=10,
    help="kmer interval(10)"
)

gapfill_group.add_argument(
    "-k",
    metavar="INT",
    type=int,
    dest="kmer",
    default=127,
    help="length of kmer(121),permit 127mer"
)

gapfill_group.add_argument(
    "-S",
    metavar="INT",
    type=int,
    dest="samp_num",
    default=10,
    help="sample number in a subfile"
)

gapfill_group.add_argument(
    "-t",
    metavar="INT",
    type=int,
    dest="cpu",
    default=4,
    help="CPU number(4)"
)


## -------------only gapfill nedd---------------------------
only_gapfill_parser = argparse.ArgumentParser(add_help=False)
only_gapfill_group = only_gapfill_parser.add_argument_group(
    "arguments when only run gapfill"
)
only_gapfill_group.add_argument(
    "-ends",
    metavar="STR",
    type=str,
    required=True,
    dest="gapfill_ends",
    help="ends fasta for gap filling"
)

only_gapfill_group.add_argument(
    "-middle",
    metavar="STR",
    type=str,
    required=True,
    dest="gapfill_mid",
    help="the middle fasta's path\n"
)

mkout_parser = argparse.ArgumentParser(add_help=False)
mkout_group = mkout_parser.add_argument_group(
    "rename gap filling result to raw barcodes"
)

mkout_group.add_argument(
    "-d",
    metavar="STR",
    type=str,
    required=True,
    dest="contigDir",
    help="the outdir of gapfill step",
)

# ------------------translation need -------------------------
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

## -----------------------polish -----------------------------
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


###############################################################################
#####----------------------- main subcommand parsers --------------------######

description = """

Description

    An automatic pipeline for HIFI-hiseq project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences).

Version

    2.1.0 2021-09-14 The first version.

Author

    yangchentao at genomics.cn, BGI.
    zhouchengran at genomics.cn, BGI.
    liushanlin at genomics.cn, BGI.

"""

parser = argparse.ArgumentParser(
    prog="HIFI-hiseq",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-v", "--version",
    action="version",
    version="%(prog)s 2.1.0"
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
        buildend_parser,
        chain_parser,
        gapfill_parser,
        trans_parser,
    ],
    formatter_class=argparse.RawTextHelpFormatter,
    help="run filter, assign, buildend, chain, and gapfill",
)

## filter subcommand
parser_filter = subparsers.add_parser(
    "filter",
    parents=[common_parser, filter_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="filter raw reads by quality or expected_err",
)

## assign subcommand
parser_assign = subparsers.add_parser(
    "assign",
    parents=[common_parser,
             index_parser,
             only_assign_parser,
             assign_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="assign clean reads to samples",
)

## buildend subcommand
parser_buildend = subparsers.add_parser(
    "buildend",
    parents=[common_parser,
             only_buildend_parser,
             soft_parser,
             buildend_parser,
             trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="buildends for each sample, output DNA fragment with tag",
)

## chain subcommand
parser_chain = subparsers.add_parser(
    "chain",
    parents=[common_parser,
             chain_parser,
             only_chain_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="connect middle meta-paried reads to longer contigs",
)

## gapfill subcommand
parser_gapfill = subparsers.add_parser(
    "gapfill",
    parents=[common_parser,
            gapfill_parser,
            only_gapfill_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="gap filling to generate raw contigs",
)

## mkout subcommand
parser_mkout = subparsers.add_parser(
    "mkout",
    parents=[common_parser,
             mkout_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="rename raw contigs to final COI barcodes"
)

## polish subcommand
parser_polish = subparsers.add_parser(
    "polish",
    parents=[polish_parser,
            index_parser,
            trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="polish COI barcode assemblies, output confident barcodes"
)

## BOLD_identification, taxonomy
parser_bold = subparsers.add_parser(
    "taxonomy",
    parents=[],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do taxa identification on BOLD system",
)

###############################################################################
#####---------------------- program execution start ----------------------#####

# -----------------------BOLD identification----------------------#
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

if sys.argv[1] == "taxonomy":
    sys.argv = sys.argv[1:]
    sys.exit(bold_identification())

args = parser.parse_args()
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
        return True
    else:
        return False

def files_exist_0_or_1(filelist):
    '''
    check files involed whether are existing!
    '''
    num = 0
    for file in filelist:
        if os.path.exists(file):
            num += 1
        else:
            print("[Error]: %s doesn't exist!" % file, file=sys.stderr)
    if len(filelist) == num:
        return 0
    else:
        return 1

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("[WARNING]: " + file + " exists! now overwriting")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out


## input file existing check
errors_found = 0
if args.command == "all":
    errors_found += files_exist_0_or_1(
        [args.fastq2, args.fastq2, args.primer]
    )
elif args.command == "filter":
    errors_found += files_exist_0_or_1([args.fastq1, args.fastq2])
elif args.command == "assign":
    errors_found += files_exist_0_or_1([args.primer])
elif args.command == "buildend":
    errors_found += files_exist_0_or_1([args.list])
elif args.command == "chain":
    errors_found += files_exist_0_or_1([args.middle_input])
elif args.command == "gapfill":
    errors_found += files_exist_0_or_1([args.gapfill_ends, args.gapfill_mid])
elif args.command == "mkout":
    errors_found += files_exist_0_or_1([args.contigDir])
elif args.command == "polish":
    errors_found += files_exist_0_or_1([args.coi_input])
else:
    parser.print_help()
    parser.exit()

if args.command in ["all", "buildend"]:
    vsearch = "vsearch"
    if hasattr(args, "vsearch"):
        if args.vsearch:
            vsearch = args.vsearch
            if check_program_involed(vsearch):
                print("[INFO]: find vesearch in {}".format(vsearch))
            else:
                print("[ERROR]: can not find vesearch in your $PATH")
                exit()

if args.command in ["all", "gapfill"]:
    # find cmd in the same folder of main script by default.
    bin_path = os.path.dirname(sys.argv[0])
    barcode_exe = os.path.join(bin_path, 'barcode')

    if check_program_involed(barcode_exe):
        print("[INFO]: find barcode in {}".format(barcode_exe) )
    else:
        print("[ERROR]: can not find soapbarcode program in your environmental $PATH ")
        exit()


if errors_found > 0:
    parser.exit("Errors found, Exit!")

if hasattr(args, "outpre") and args.outpre.endswith("/"):
    print("[ERROR]: outpre is in bad format! no \"/\"")
    exit()


# -----------------------functions for common---------------------#
def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def checkDirMkdir(fold):
    outdir = os.path.abspath(args.outpre + "_assign")
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    return outdir

def print_time(info):
    print(info + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def run_time(info):
    print("[INFO]: {0} total run time: {1:.2f}".format(info,
                                                       time.time() - t) + "s")

# -----------------------functions for filtering------------------#

def deleteDuplicatedElementFromList3(listA):
        #return list(set(listA))
        return sorted(set(listA), key = listA.index)

def parse_pe_fastq(fq1, fq2):
    while True:
        name1 = fq1.readline().split(" ")[0]
        name2 = fq2.readline().split(" ")[0]
        if not name1 or not name2:
            break
        read1, nothing1, qual1 = fq1.readline()[:-1], fq1.readline(), fq1.readline()[:-1]
        read2, nothing2, qual2 = fq2.readline()[:-1], fq2.readline(), fq2.readline()[:-1]
        assert name1 == name2, 'fastq1, fastq2 is not paired-end'
        yield name1, read1, read2, qual1, qual2

def exp_e(q_list, phred):
    '''
    expected error number(E*) = sum(P), where P is
    the probability that the base call is incorrect.
    '''
    exp = 0
    quals = [ord(n) - phred for n in q_list]
    for i in quals:
        exp += 10 ** (-i / 10)

    return exp

def lowquality_rate(q_list, cut_off, phred):
    '''
    calculate the rate of low quality base in read.
    '''
    low_base = 0
    quals = [ord(n) - phred for n in q_list]
    for i in quals:
        if i < cut_off:
            low_base += 1

    low_rate = low_base / len(quals)
    return low_rate


# ----------------------functions for assigning-------------------#

def complementation(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    sequence = sequence.translate(transtable)
    return sequence # [BUG20210914]


def comp_rev(sequence):
    # make a sequence complement and reversed #
    sequence = complementation(sequence)
    return sequence[::-1]

def decode(head_flen, head_rlen, dict):

    if head_flen in dict.keys():
        target = dict[head_flen]
    elif head_rlen in dict.keys():
        target = dict[head_rlen]
    else:
        target = False

    return target

def decode_mismatch(f, r, dict):
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
                    if s1[base] != s2[base]:
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
                if b1[base] != b2[base]:
                    mismatch += 1
            dis.append(mismatch)
    min_dis = min(dis)
    max_dis = max(dis)
    return (min_dis, max_dis)

# ----------------------functions for buildend------------------#

def parse_ssam(ssam):
    while True:
        tmp = ssam.readline().strip().split()
        name = tmp[0]
        if not name:
            break
        read1, read2 = tmp[1], tmp[3]
        qual1, qual2 = tmp[2], tmp[4]
        if len(read1) == len(read2) and len(qual1) == len(qual2):
            yield name, read1, read2, qual1, qual2

def comp_rev_list(reads_list):
    # make a list of sequences reverse and complement #
    new_reads_list = []
    for read in reads_list:
        read = comp_rev(read)
        new_reads_list.append(read)

    return new_reads_list

def mismatch(str1, str2):
    # ----count matched bases of two sequences----#
    mismatches = 0
    for base in range(len(str1)):
        if str1[base] != str2[base]:
            mismatches += 1

    return mismatches

def translate_dnaseq(seq, codon):
    # ---------translate_dnaseq------------#
    l_dna = len(seq)
    if l_dna % 3 != 0:
        seq = seq[: -(l_dna % 3)]
        # print("your sequence lenght is not tripple" + \
        # "but no worries, I have trimmed well format")
    #coding_dna = Seq(seq, generic_dna)
    #protein = coding_dna.translate(table=codon)
    protein = Seq(seq).translate()
    if "*" in protein:
        return False
    else:
        return True

def connectMetapairReads(record,
                         min_overlap=10,
                         max_overlap=150,
                         overlap_mismatch=1,
                         standard_length=150):

    records = record.strip().split()
    name = records[0]
    seq1, seq2 = records[1], comp_rev(records[3]) # get complementary and reverse seq2
    forward_qual, reverse_qual = records[2], records[4][::-1] # reverse seq2's quality
    singal = 0
    overlaps = {}
    for s in range(min_overlap, max_overlap + 1):
        l0 = seq1[-s:]
        l1 = seq2[0:s]
        tmp_mismatch = mismatch(l0, l1)
        if tmp_mismatch == 0:
            overlaps[s] = 1
            # find best result, so exit loop #
            break

        elif tmp_mismatch <= overlap_mismatch:
            tmp_identity = 1 - (tmp_mismatch / standard_length)
            overlaps[s] = tmp_identity

    # find best overlaping result in all potenial positions
    # candidates = sorted(overlaps.items(),
    # lambda x, y: cmp(x[1], y[1]), reverse=True)
    candidates = sorted(
        overlaps, key=overlaps.__getitem__, reverse=True
    )

    if len(candidates) > 0:

        potenial = candidates[0]  # overlap similarity top 1
        s0 = seq1[-potenial:]
        s1 = seq2[0:potenial]

        corrected = ""

        # compare each base from forward and reverse to keep one
        # forward == reverse
        # forward ne reverse, quality(forward) > quality(reverse)
        # forward ne reverse, quality(forward) < quality(reverse)

        for p in range(len(s0)):
            # site is changed, be careful!#
            tmp_loca0 = standard_length - potenial + p

            if s0[p] == s1[p]:
                corrected += s0[p]
                info = s0[p] + "=" + s1[p]
            else:
                for_quality = forward_qual[tmp_loca0]
                rev_quality = reverse_qual[p]
                if for_quality >= rev_quality:
                    corrected += s0[p]
                else:
                    corrected += s1[p]

        makeup_consensus = (
            seq1[: standard_length - potenial]
            + corrected
            + seq2[potenial - standard_length :]
        )

        len_makeup_consensus = len(makeup_consensus)
        this_oid = overlaps[potenial] * 100
        this_oid = str("%.2f" % this_oid)

        return makeup_consensus

def sortLengthSizeTrim(seqs, ori=1):
    """
    sort list by count of item
    """
    # when connect ends, I reverse and complement all 2 end,
    # so for reverse primer end, must return very begining
    # state.
    if ori == 2:
        seqs = comp_rev_list(seqs)
    # trim low coverage bases from tail
    seqlens = [len(n) for n in seqs]
    max_length = max(seqlens)
    effect_length = max_length
    """from tail to head, check the coverage of each base"""
    for i in range(max_length, 0, -1):
        coverage_this_base = 0
        for item in seqs:
            if len(item) >= i:
                coverage_this_base += 1
        if coverage_this_base >= 5:
            break
        else:
            effect_length -= 1
    new_seqs = []
    for s in seqs:
        new_seqs.append(s[:effect_length])

    # sort trimed seq list by abundance
    count = {}
    for s in new_seqs:
        if s in count.keys():
            count[s] += 1
        else:
            count[s] = 1
    #just sort by abundance.
    #seqs = sorted(new_seqs, key=lambda k : count[k], reverse=True)

    #sort by length and abundance
    seqs = sorted(new_seqs, key=lambda k : (len(k), count[k]), reverse=True)

    if ori == 1:
        return seqs
    else:
        return comp_rev_list(seqs)

# --------------------function for gapfill---------------------#
def parse_fasta(fa_fh):
    while True:
        name = fa_fh.readline().strip()
        if not name or name[0] != ">":
            break
        name = name.replace(">", "")
        seq = fa_fh.readline().strip()
        yield name, seq

def format_ends(ends, subsam, split_ends_dir):
    endfile = []
    seqsf = {}
    seqsr = {}
    abundance = {}
    all_samples = []
    with open(ends, 'r') as fh:
        for i in parse_fasta(fh):
            head, seq = i
            tmp = head.split("_")
            tag1 = int(tmp[0][-3:])
            if tag1 not in all_samples:
                all_samples.append(tag1)
            ab = tmp[2]
            if "For" in head:
                if tag1 not in seqsf.keys():
                    seqsf[tag1] = [seq]
                else:
                    seqsf[tag1].append(seq)
                abundance[seq] = ab
            elif "Rev" in head:
                if tag1 not in seqsr.keys():
                    seqsr[tag1] = [seq]
                else:
                    seqsr[tag1].append(seq)
                abundance[seq] = ab
            else:
                print("[ERROR]: wrong format in head of ends fasta")
                exit()

    #sub_len = int(len(seqsf.keys()) / subsam)
    sub_len = int(subsam)
    current_len = 0
    current_cont = ""
    file_mark = 1
    sorted_all_samples = sorted(all_samples)
    for s in sorted_all_samples:
        sample_id = format(s, '03d')
        current_len += 1
        order_sam = 0
        for f in range(len(seqsf[s])):
            fs = f + 1
            abf = abundance[seqsf[s][f]]
            for r in range(len(seqsr[s])):
                rs = r + 1
                order_sam += 1
                abr = abundance[seqsr[s][r]]
                end = ">{0}_{1}-{2};{5};{3}-{4}\n".format(
                    sample_id, fs, rs, abf, abr, order_sam)
                end += "{0}\n{1}".format(seqsf[s][f], seqsr[s][r])
                if len(current_cont) == 0:
                    current_cont = end
                else:
                    current_cont += "\n" + end

        if current_len >= sub_len:
            outfile = split_ends_dir + "/ends." + str(file_mark)
            endfile.append(outfile)
            with open(outfile, 'w') as endf:
                endf.write(current_cont)
                current_len = 0
                current_cont = ""
                file_mark += 1
    ##make the last file
    if current_cont:
        outfile = split_ends_dir + "/ends." + str(file_mark)
        endfile.append(outfile)
        with open(outfile, 'w') as endf:
            endf.write(current_cont)

    return endfile

# -------------------functions for translation-----------------#
def coi_check(contig, codon):
    # ---------------coi_check------------#
    for_trim = args.index + 25 + 1
    rev_trim = args.index + 26
    contig = contig[for_trim:]
    # contig = contig[0:for_trim] #bug!
    contig = contig[:-rev_trim]
    if translate_dnaseq(contig, codon):
        return True
    else:
        return False

#-------------------------run suncommand------------------------#
if args.command in ["all", "filter"]:
    t = time.time()
    begining_all = time.time()
    print_time("[INFO]: Filtering starts:")
    filter_outdir = checkDirMkdir(args.outpre + "_filter")
    filtered_outfile1 = filter_outdir + "/filter_highqual_1.fastq"
    filtered_outfile2 = filter_outdir + "/filter_highqual_2.fastq"
    out1 = check_and_open_outhandle(filtered_outfile1)
    out2 = check_and_open_outhandle(filtered_outfile2)

    # ini state
    total = 0
    clean = 0
    nn = 0
    low = 0

    log = check_and_open_outhandle(filter_outdir + "/filter_log.txt")

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

    fq1 = open_input(args.fastq1)
    fq2 = open_input(args.fastq2)

    for i in parse_pe_fastq(fq1, fq2):
        name, seq1, seq2, qual1, qual2 = i

        if not name[0].startswith("@"):
            print("[ERROR]: input fastq is not a correct fastq format")
            exit()

        total += 1
        N_count1 = seq1.count("N")
        N_count2 = seq2.count("N")

        if N_count1 < args.n or N_count2 < args.n:
            if filter_type == 1:
                if (exp_e(qual1, args.phred) <= args.expected_err
                    and exp_e(qual2, args.phred) <= args.expected_err):
                    out1.write(name + " 1\n" + seq1 + "\n" + "+\n" + qual1 + "\n")
                    out2.write(name + " 2\n" + seq2 + "\n" + "+\n" + qual2 + "\n")
                    clean += 1
                else:
                    low += 1
            else:
                # filter_type == 2
                if (lowquality_rate(qual1, high_qual, args.phred) > low_qual_cont
                    and lowquality_rate(qual2, high_qual, args.phred) > low_qual_cont):
                    out1.write(name + " 1\n" + seq1 + "\n" + "+\n" + qual_str1 + "\n")
                    out2.write(name + " 2\n" + seq2 + "\n" + "+\n" + qual_str2 + "\n")

                    clean += 1
                else:
                    low += 1
        else:
            nn += 1

    log.write("total reads:\t{}".format(total) + "\n")
    log.write("clean reads:\t{}".format(clean) + "\n")
    log.write("low quality reads:\t{}".format(low) + "\n")
    log.write("containing N (> {}) reads:\t{}".format(args.n, nn) + "\n")

    log.close()
    fq1.close()
    fq2.close()
    out1.close()
    out2.close()

    print_time("[INFO]: Filtering done:")
    run_time("Filtering")

#----------------assign----------------------#

if args.command in ["all", "assign"]:
    t = time.time()
    print_time("[INFO]: Assigning starts:")

    if args.command == "all":
        args.fq1 = filtered_outfile1
        args.fq2 = filtered_outfile2

    assign_outdir = checkDirMkdir(args.outpre + "_assign")

    ErrFile = check_and_open_outhandle(assign_outdir + "/assign_err.fasta")
    middle_assigned = assign_outdir + "/assign_middle.ssam"
    Middle_ssam = check_and_open_outhandle(middle_assigned)

    indexlen = args.index

    pris = {}
    indp = {}
    FH = {}
    barcodes = []

    with open(args.primer, "r") as p:
        primer_lines = 0
        for i in p.readlines():
            primer_lines += 1
            i = i.strip()
            arr = i.split()
            if len(arr) != 2:
                print("primer set is not well-formated")
                exit()
            sam = arr[0]
            ipr = arr[1]
            FH[ipr] = arr[0]
            if sam in indp.keys():
                print(arr[0] + "show twice in primer set")
            else:
                ori = sam[0:3]
                num = sam[-3:]

                pris[num] = {}
                pris[num][ori] = ipr
                indp[sam] = ipr

            if ori == "For":
                plenf = len(ipr)
                primerF = ipr[indexlen:]
                # save barcodes in a array, analyze later
                barcodes.append(ipr[0:indexlen])

            if ori == "Rev":
                plenr = len(ipr)
                primerR = ipr[indexlen:]
    # check primer lines
    if primer_lines % 2 != 0:
        print("primer lines ({}) is not even number".format(primer_lines))
        print("the primer file need to have each forward and reverse primer")
        exit()

    # analysis barcodes and demultiplex argument(mismatch)
    (min_dis, max_dis) = dis_barcode(barcodes)
    print("[INFO]: min distance among barcodes is {}".format(min_dis))
    print("[INFO]: max distance among barcodes is {}".format(max_dis))
    print("[INFO]: mismatches allowed in barcodes is {}".format(args.tag_mismatch))

    if args.tag_mismatch and args.tag_mismatch > (min_dis - 1):
        print("mismatch you set is too large to demultiplex, it must be smaller"
             + " than {},".format(min_dis) + " because min distance among barcodes"
             + " is {}".format(min_dis)
             )
        exit()
    # ini dict of count_total and count_assigned for statistic.
    count_assigned = {}
    for fh in FH.values():
        count_assigned[fh] = 0
    count_assigned['middle'] = 0

    neg_priF = comp_rev(primerF)
    neg_priR = comp_rev(primerR)
    assigned_list =  assign_outdir + "/assign.list"

    with open(assigned_list, "w") as ls:
        sorted_sample = sorted(pris.keys())
        for s in sorted_sample:
            ls.write(assign_outdir
                     + "/For" + s + ".ssam"+ "\n"
                     + assign_outdir
                     + "/Rev" + s + ".ssam" + "\n")
    # open all assigned files
    filehandle = {}
    for sam in indp.keys():
        filehandle[sam] = open(assign_outdir + "/" + sam + ".ssam", "w")

    err = 0
    seqnum = 0
    assigned = 0

    # open clean paired fastq
    fh1 = open_input(args.fq1)
    fh2 = open_input(args.fq2)

    for i in parse_pe_fastq(fh1, fh2):
        name, seq1, seq2, qual1, qual2 = i
        seqnum += 1
        fmatch = 0
        rmatch = 0
        head1f = seq1[0:plenf] # head length(forward primer) of read1
        head1r = seq1[0:plenr] # head length(reverse primer) of read1
        head2f = seq2[0:plenf] # head length(forward primer) of read2
        head2r = seq2[0:plenr] # head length(reverse primer) of read2
        # cut head (max of for and rev) to make a tmp sequence
        len_head_cut = max(plenf, plenr)
        tmp1 = seq1[len_head_cut:]
        tmp2 = seq2[len_head_cut:]
        # if primer in the wrong position, remove this reads.
        if (primerF in tmp1 or primerR in tmp1
            or neg_priF in tmp1 or neg_priR in tmp1):
            ErrFile.write(">" + name + "_priErr_1\n" + seq1 + "\n")
            ErrFile.write(">" + name + "_2\n" + seq2 + "\n")
            err += 1
            continue
        if (primerF in tmp2 or primerR in tmp2
            or neg_priF in tmp2 or neg_priR in tmp2):
            ErrFile.write(">" + name + "_1\n" + seq1 + "\n")
            ErrFile.write(">" + name + "_priErr_2\n" + seq2 + "\n")
            err += 1
            continue
        # deal forward read
        fmatch = decode(head1f, head1r, FH)
        fmatch_withMis = decode_mismatch(head1f, head1r, FH)
        rmatch = decode(head2f, head2r, FH)
        rmatch_withMis = decode_mismatch(head2f, head2r, FH)

        if (fmatch or fmatch_withMis) and (rmatch == False
                                           and rmatch_withMis == False):
            assigned += 1
            targetf = fmatch or fmatch_withMis
            count_assigned[targetf] += 1
            if targetf[:3] == "For":
                filehandle[targetf].write(
                    targetf + "_" + str(seqnum) + " "
                    + seq1 + " " + qual1 + " ")
                filehandle[targetf].write(
                    seq2 + " " + qual2 + "\n")
            else:
                filehandle[targetf].write(
                    targetf + "_" + str(seqnum) + " "
                    + seq2 + " " + qual2 + " ")
                filehandle[targetf].write(
                    seq1 + " " + qual1 + "\n")

        elif (rmatch or rmatch_withMis) and (fmatch == False
                                             and fmatch_withMis == False):
            assigned += 1
            targetr = rmatch or rmatch_withMis
            count_assigned[targetr] += 1
            if targetr[:3] == "For":
                filehandle[targetr].write(
                    targetr + "_" + str(seqnum) + " "
                    + seq2 + " " + qual2 + " ")
                filehandle[targetr].write(
                    seq1 + " " + qual1 + "\n")
            else:
                filehandle[targetr].write(
                    targetr + "_" + str(seqnum) + " "
                    + seq1 + " " + qual1 + " ")
                filehandle[targetr].write(
                    seq2 + " " + qual2 + "\n")
        elif (fmatch == False
              and fmatch_withMis ==False
              and rmatch == False
              and rmatch_withMis == False):
            assigned += 1
            count_assigned['middle'] += 1
            Middle_ssam.write(
                ">" + str(seqnum) + "_1 " + seq1 + " " + qual1
                + " " + seq2 + " " + qual2 + "\n"
            )

        else:
            err += 1

    ErrFile.close()
    # close all assigned files
    for fh in filehandle.values():
        fh.close()

    # report assignment information
    with open(assign_outdir + "/assign.log", "w") as log:
        log.write("total reads:\t{}\n".format(seqnum))
        log.write("err reads:\t{}\n".format(err))
        log.write("assigned:\t{}\n".format(assigned))
        for i in sorted(count_assigned.keys()):
            log.write(
                i + "\t" + str(count_assigned[i]) + "\n"
            )

    print_time("[INFO]: Assigning done:")
    run_time("Assigning")

#-------------------buildend----------------------------#
if args.command in ["all", "buildend"]:
    t = time.time()
    print_time("[INFO]: Building ends starts:")

    if args.min_overlap < 10 or args.min_overlap > 120:
        print("[ERROR]: "
            + "For COI barcodes, by and large, overlaped length is 30~120 bp, so\
              {} is not proper!".format(args.min_overlap)
        )
        exit()
    if args.max_overlap < args.min_overlap:
        print("[ERROR]: maximum overlap length must be large than minimun")
        exit()

    if args.command == "all":
        args.list = assigned_list  # list generated from assign step

    buildend_outdir = checkDirMkdir(args.outpre + "_buildend")
    buildends_result = buildend_outdir + "/buildends.fasta"
    fh_out = check_and_open_outhandle(buildends_result)
    fh_log = check_and_open_outhandle(buildend_outdir + "/buildends.log")

    fh_log.write("## assigned reads list file = " + args.list + "\n")

    if args.seqs_lim:
        fh_log.write("## reads input limitation: " + str(args.seqs_lim) + "\n")
    else:
        fh_log.write("## Using all reads to make consensus, no limitation\n")

    fh_log.write("## input read length = " + str(args.standard_length) + "\n")
    fh_log.write("## vsearch path = " + vsearch + "\n")

    if args.reads_check:
        fh_log.write("## check codon translation = yes\n")
    else:
        fh_log.write("## check codon translation = no\n")

    fh_log.write("## clustering identity = " + str(args.cluster_identity) + "\n")

    fh_log.write("## overlaping mismatch allow = " + str(args.buildend_mismatch) + "\n")

    fh_log.write("## min overlap = " + str(args.min_overlap) + "\n")
    fh_log.write("## max overlap = " + str(args.max_overlap) + "\n")

    # --------------main-----------------------#

    try:
        with open(args.list) as fh_list:
            lines = fh_list.readlines()
    except FileNotFoundError:
        print("[ERROR]: can not find " + args.list)
        exit(0)
    total_pairs = 0
    connected_pairs = 0
    for line in lines:
        line = line.rstrip()
        name = os.path.basename(line).split(".")[0]
        short_outname = name[:6]
        if "For" in name:
            ori = 1
        else:
            ori = 2
        fh_log.write("//processing " + name + " done\n")
        if os.path.exists(line) == False:
            print("[ERROR]: can not find " + line)
            exit()

        # if file is empty, continue
        if os.path.getsize(line) == 0:
            fh_log.write("! file is empty!")
            fh_log.write("\n")
            continue
        success_connected = []
        with open(line,'r') as fh:
            records = fh.readlines()
            for r in records:
                total_pairs += 1
                if (len(r.strip().split())) != 5:
                    continue
                makeup_consensus = connectMetapairReads(
                    r,
                    min_overlap=args.min_overlap,
                    max_overlap=args.max_overlap,
                    overlap_mismatch=args.buildend_mismatch,
                    standard_length=args.standard_length,
                    )

                if makeup_consensus:
                    # if check result, and ok so write into output_checked file #
                    success_connected.append(makeup_consensus)
                    connected_pairs += 1
        # sort success_connected by length and size, meanwhile trim end by
        # (coverage < 5)
        if success_connected:
            success_connected = sortLengthSizeTrim(success_connected, ori)

            # output successfully connected ends to a temp file
            pid = os.getpid()
            temp_fasta = "temp.fa" + "." + str(pid)
            temp_uc = "temp.uc" + "." + str(pid)
            with open(temp_fasta, "w") as TM:
                for i in range(len(success_connected)):
                    TM.write(">" + str(i) + "\n" + success_connected[i] + "\n")

            # cluster these ends
            vsearch_cmd = (
                vsearch
                + " --cluster_smallmem "
                + temp_fasta
                + " -usersort "
                + " --threads "
                + str(args.threads)
                + " --quiet "
                + " --uc "
                + temp_uc
                + " --id "
                + str(args.cluster_identity)
            )
            subprocess.call(vsearch_cmd, shell=True)
            # to store clusters, clusters[represent ID] = [clustered IDs]
            clusters = {}
            # to store each cluster's abundance, for sorting clusters
            # e.g. count[represent ID] = 100
            count = {}

            # open temp.uc and statistic each cluster's abundance.
            with open(temp_uc, "r") as uc:
                # there are "H","S","C" in the head of line
                for line in uc.readlines():
                    if line[0] != "H":
                        continue
                    array = line.split()
                    if array[9] in clusters.keys():
                        clusters[array[9]].append(array[8])
                        count[array[9]] += 1
                    else:
                        clusters[array[9]] = []
                        clusters[array[9]].append(array[9])
                        count[array[9]] = 1

            # sorting clusters by abundance.#
            # sorted_clusters = sorted(count, key=count.__getitem__,reverse=True)
            sorted_clusters = sorted(count, key=lambda k: (count[k], k), reverse=True)

            if args.cluster_number_needKeep:
                # if set "-tp", keep top N clusters to buildend
                keep = args.cluster_number_needKeep
                sorted_clusters = sorted_clusters[0:keep]
            elif args.abundance_threshod:
                # if set "-ab", keep all clusters of abundance > ab
                while sorted_clusters:
                    item = sorted_clusters.pop()
                    if count[item] < args.abundance_threshod:
                        pass
                    else:
                        sorted_clusters.append(item)
                        break
            elif len(sorted_clusters) > 1:
                # if set nothing, I will set it to -tp 2
                # if second most abundant sequence less than 1/10 of first,
                # remove it!# of course it is just for when tp==2
                sorted_clusters = sorted_clusters[0:2]
                if count[sorted_clusters[1]] < count[sorted_clusters[0]] / 10:
                    sorted_clusters.pop()
            order = 0
            for k in sorted_clusters:
                order += 1
                fh_out.write(">"
                            + short_outname
                            + "_" + str(order)
                            + "_" + str(count[k]) + "\n"
                            + success_connected[int(k)] + "\n"
                            )
            fh_log.write("{0} {1} {2:.3f}".format(total_pairs,
                                            connected_pairs,
                                            connected_pairs/total_pairs))
        else:
            print(f"[WARNING]: No read pair can be connected in {short_outname}! please check your data!")

    fh_out.close()
    fh_log.close()
    if success_connected:
        rm_tmp_cmd = "rm temp.fa.* temp.uc.*"
        os.system(rm_tmp_cmd)
    print_time("[INFO]: Building ends done:")
    run_time("Building")

# -------------------chain ---------------------#
if args.command in ["all", "chain"]:
    t = time.time()
    print_time("[INFO]: Middle Chaining starts:")
    if args.command == "all":
        inputMiddle = middle_assigned
        if args.min_insertsize < args.standard_length:
            print("[ERROR]: min_insertsize can not be less than " + args.standard_length)
            exit()
    else:
        inputMiddle = args.middle_input
        if args.min_insertsize < 150:
            print("[ERROR]: min_insertsize can not be less than 150 bp")
            exit()

    chain_outdir = checkDirMkdir(args.outpre + "_chain")
    chain_out = chain_outdir + "/chain.fasta"
    chain_log = chain_outdir + "/chain.log"
    cout = open(chain_out, 'w')
    clog = open(chain_log, 'w')
    clog.write("total_pairs	connected_pairs	connect_ratio(%)\n")
    total_pairs = 0
    connected_pairs = 0
    with open_input(inputMiddle) as fm:
        for r in fm:
            total_pairs += 1
            tmp = r.strip().split()
            if len(tmp) != 5:
                continue
            else:
                name = tmp[0].split("_")[0]
                name = name.replace(">", "")
            middle_connected = connectMetapairReads(r,
                                                    overlap_mismatch=args.chain_mismatch)

            if middle_connected and len(middle_connected) >= args.min_insertsize:
                cout.write(">cop" + name + "\n"
                           + middle_connected + "\n")
                connected_pairs += 1

    clog.write("{0} {1} {2:.3f}".format(total_pairs,
                                        connected_pairs,
                                        connected_pairs/total_pairs))
    cout.close()
    clog.close()

    print_time("[INFO]: Middle Chaining done:")
    run_time("Middle Chaining")

# -------------------gap filling---------------------------#
if args.command in ["all", "gapfill"]:
    t = time.time()

    print_time("[INFO]: Gap filling starts:")

    gapfill_outdir = checkDirMkdir(args.outpre + "_gapfill")
    result_outdir = checkDirMkdir(gapfill_outdir + "/result")
    split_ends_dir = checkDirMkdir(gapfill_outdir + "/ends")
    shell_outdir = checkDirMkdir(gapfill_outdir + "/shell")
    middle_lis = gapfill_outdir + "/middle.lis"

    if args.command == "all":
        ends = buildends_result
    else:
        ends = args.gapfill_ends
        chain_out = os.path.abspath(args.gapfill_mid)

    with open(middle_lis, 'w') as mf:
        mf.write(">\n"
                 +"f=" + chain_out)


    # format ends fasta
    # split ends fasta to several subfiles
    if args.samp_num < 1:
        print("samp_num can not be 0")
        exit()
    split_ends = format_ends(ends, args.samp_num, split_ends_dir)
    subfile = 0
    for s in split_ends:
        s = os.path.abspath(s)
        subfile += 1
        shell = shell_outdir + "/barcodes." + str(subfile) + ".sh"
        barcode_out = result_outdir + "/barcodes." + str(subfile)
        barcode_out = os.path.abspath(barcode_out)
        with open(shell, 'w') as sh:
            sh.write(barcode_exe
                     + " -e " + s
                     + " -r " + middle_lis
                     + " -o " + barcode_out
                     + " -x " + str(args.gmax)
                     + " -n " + str(args.gmin)
                     + " -l " + str(args.lowkmer)
                     + " -v " + str(args.kmer_interval)
                     + " -k " + str(args.kmer)
                     + " -t " + str(args.cpu)
                    )

    print("[RUN]: step 1: run shell files in " + shell_outdir)

    print_time("[INFO]: Making gap-filling shell files done:")
    run_time("Gap-filling")

# ---------------------make output -------------------------#
if args.command == "mkout":
    t = time.time()
    print_time("[INFO]: Generating result starts:")

    if os.path.exists(args.contigDir) == False:
        print("[ERROR]: can not find contigDir")
        exit()
    mkout_outdir = checkDirMkdir(args.outpre + "_mkout")
    for e in os.listdir(args.contigDir + "/ends/"):
        ef = args.contigDir + "/ends/" + e
        file_maker = ef.split(".")[-1]
        count = 0
        order = {}
        with open(ef, 'r') as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1
                    order[count] = line.strip().replace(">", "")
        cooresponding_result = (args.contigDir
                                +"/result/barcodes." + file_maker + ".contig")
        if os.path.exists(cooresponding_result) == False:
            print("[ERROR]: can not find " + cooresponding_result)
            continue
        new_cooresponding_result = cooresponding_result + ".add"
        new_out = check_and_open_outhandle(new_cooresponding_result)
        with open(cooresponding_result) as fh:
            for i in parse_fasta(fh):
                head, seq = i
                tmp = head.split()[0].split("_")
                kmer_info = tmp[0].replace(">", "")
                sample_order = int(tmp[1])
                haplotype = int(tmp[2])
                # just output the first haplotype
                if haplotype == 1:
                    head_items = [order[sample_order], "k="+ kmer_info]
                    new_head = ";".join(head_items)
                    new_out.write(">" + new_head + "\n" + seq + "\n")
        new_out.close()
    # cat results to one file
    generate_coi = "cat " + args.contigDir + "/result/*.add > " +\
            mkout_outdir + "/all_barcodes.fasta"
    subprocess.call(generate_coi, shell=True)
    print_time("[INFO]: Generating result done:")
    run_time("Generating result")

# ----------------------polish assemblies----------------------------#
if args.command == "polish":
    t = time.time()
    print_time("[INFO]: Polishing starts:")

    polish_outdir = checkDirMkdir(args.outpre + "_polish")
    polish_outfile = polish_outdir + "/all.barcodes.polished"
    coiout = open(polish_outfile,'w')

    total_barcodes = 0
    good_barcodes = 0
    with open(args.coi_input, "r") as handle:
        kmer_info = {}
        seq_abu = {}
        sample_seq = {}
        for i in parse_fasta(handle):
            name, seq = i
            total_barcodes += 1
            items = name.split(";")
            kmer = items[-1]
            sample_tag = items[0].split("_")
            sample = sample_tag[0]
            coverages = items[2].split("-")
            for_coverage = int(coverages[0])
            rev_coverage = int(coverages[1])
            abu = int((for_coverage + rev_coverage) / 2)
            length = len(seq)
            # remove low record with low coverage or short length
            if (for_coverage < args.min_coverage
                or rev_coverage < args.min_coverage
                or length < args.coi_length):
                continue

            elif args.coi_check and coi_check(str(record.seq), args.codon_table) == False:
                print(str(sample_tag) + " translation failed")
                continue

            else:
                if sample in sample_seq.keys():
                    sample_seq[sample].append(seq)
                else:
                    sample_seq[sample] = [seq]
                seq_abu[seq] = abu
                kmer_info[seq] = kmer

    # pick up most confident COI barcode.
    sorted_samples = sorted(sample_seq.keys())
    polished_count = len(sorted_samples)
    for s in sorted_samples:
        sorted_seq = sorted(sample_seq[s], key=lambda k: seq_abu[k], reverse=True)
        seq = sorted_seq[0]
        top_abu = seq_abu[seq]
        kmer = kmer_info[seq]
        seqlen = len(seq)
        head = ">{0};size={1};{2};l={3}".format(s, top_abu, kmer, seqlen)
        coiout.write(head + "\n" + seq + "\n")
        good_barcodes += 1
    coiout.close()
    print("All barocdes generated: {}".format(total_barcodes))
    print("All barocdes polished: {}".format(good_barcodes))
    print_time("[INFO]: Polishing done:")
    run_time("Polishing")

if args.command == "all":
    total_time = time.time() - begining_all
    print("[INFO]: Total run time: {0:.2f}".format(total_time + "s"))
