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
    default=30,
    dest="min_overlap",
    help="minimun length of overlap, default=30",
)

buildend_group.add_argument(
    "-max",
    metavar="INT",
    type=int,
    default=120,
    dest="max_overlap",
    help="maximum length of overlap, default=120",
)

buildend_group.add_argument(
    "-oid",
    metavar="FLOAT",
    type=float,
    default=0.95,
    dest="overlap_identity",
    help="minimun similarity of overlap region, default=0.95",
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
    "-mode",
    metavar="INT",
    type=int,
    choices=[1, 2],
    default=1,
    help="1 or 2; modle 1 is to cluster and keep most [-tp] abundance\n"
    + "clusters, or clusters abundance more than [-ab], and then make\n"
    + "a consensus sequence for each cluster. modle 2 is directly to \n"
    + "make only one consensus sequence without clustering. default=1\n",
)

buildend_group.add_argument(
    "-rc",
    dest="reads_check",
    action="store_true",
    help="whether to check amino acid translation\n" + "for reads, default not",
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
    "-cmr",
    metavar="STR",
    type=str,
    dest="cmr",
    help="the path of cmr software",
)

## only chain need
only_chain_parser = argparse.ArgumentParser(add_help=False)
only_chain_group = only_chain_parser.add_argument_group(
    "when only use chain command"
)

only_chain_group.add_argument(
    "-1",
    metavar="STR",
    type=str,
    dest="middle_R1",
    required=True,
    help="read 1 of middle assignment",
)

only_chain_group.add_argument(
    "-2",
    metavar="STR",
    type=str,
    dest="middle_R2",
    required=True,
    help="read 2 of middle assignment",
)

## gapfill for all mode ---------------------------------------
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
    default=8,
    help="CPU number(8)"
)


## only gapfill nedd
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

# translation need ---------------------------------------
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

## polish --------------------------------------------------
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

    An automatic pipeline for HIFI-hiseq project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences).

Version

    0.0.1 2018-12-13 The first version.

Author

    yangchentao at genomics.cn, BGI.
    zhouchengran at genomics.cn, BGI.
    liushanlin at genomics.cn, BGI.
    mengguanliang at genomics.cn, BGI.

"""

parser = argparse.ArgumentParser(
    prog="HIFIBarcode",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-v", "--version",
    action="version",
    version="%(prog)s 0.0.1"
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

## buildend subcommand
parser_buildend = subparsers.add_parser(
    "buildend",
    parents=[common_parser,
             only_buildend_parser,
             soft_parser,
             buildend_parser,
             trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do buildend from input fastq\n" + "reads, output HIFI barcodes.",
)

## chain subcommand
parser_chain = subparsers.add_parser(
    "chain",
    parents=[common_parser,
             chain_parser,
             only_chain_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="connect middle meta-paried reads to longer contigs.",
)

## gapfill subcommand
parser_gapfill = subparsers.add_parser(
    "gapfill",
    parents=[common_parser,
            gapfill_parser,
            only_gapfill_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="gap filling to make complete COI barcodes.",
)

## polish subcommand
parser_polish = subparsers.add_parser(
    "polish",
    parents=[polish_parser,
            index_parser,
            trans_parser],
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

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("WARRNING: " + file + " exists! now overwriting")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

# ----------------------------------------------------------------
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
    errors_found += files_exist_0_or_1([args.middle_R1, args.middle_R2])
elif args.command == "gapfill":
    errors_found += files_exist_0_or_1([args.gapfill_ends, args.gapfill_mid])
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
    errors_found += check_program_involed(vsearch)

if args.command in ["all", "chain"]:
    cmr = "cmr"
    if hasattr(args, "cmr"):
        if args.cmr:
            cmr = args.cmr
    errors_found += check_program_involed(cmr)

if errors_found > 0:
    parser.exit("Errors found! Exit!")

if hasattr(args, "outpre") and args.outpre.endswith("/"):
    print("outpre is in bad format! no \"/\"")
    exit()

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
        qual1 = list(qual1)
        qual2 = list(qual2)
        #qual1 = fromstring(qual1, dtype=byte) - phred
        #qual2 = fromstring(qual2, dtype=byte) - phred
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
    sequence = sequence.upper() ## [a bug fixed], reported by Wu Ping 20181129
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence

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

# ----------------------functions for buildend------------------#

def comp_rev_list(reads_list):
    # make a list of sequences reverse and complement #
    new_reads_list = []
    for read in reads_list:
        read = comp_rev(read)
        new_reads_list.append(read)

    return new_reads_list

def match(str1, str2):
    # ----count matched bases of two sequences----#
    matched = 0
    for base in range(len(str1)):
        if str1[base] == str2[base]:
            matched += 1
    identity = matched / len(str1)
    return identity

def translate_dnaseq(seq, codon):
    # ---------translate_dnaseq------------#
    l_dna = len(seq)
    if l_dna % 3 is not 0:
        seq = seq[: -(l_dna % 3)]
        # print("your sequence lenght is not tripple" + \
        # "but no worries, I have trimmed well format")
    coding_dna = Seq(seq, generic_dna)
    protein = coding_dna.translate(table=codon)
    if "*" in protein:
        return False
    else:
        return True


# --------------------function for gapfill---------------------#
def parse_fasta(fa_fh):
    while True:
        name = fa_fh.readline().strip()
        if len(name) == 0:
            break
        seq = fa_fh.readline().strip()
        yield name, seq


def strpy(w1=96, w2=2):
    x = []
    for i in range(w1):
        arr = []
        for j in range(w2):
            arr.append(0)
        x.append(arr)
    return x

def format_ends(ends, subsam, split_ends_dir):
    endfile = []
    seqsf = {}
    seqsr = {}
    abundance = {}
    with open(ends, 'r') as fh:
        for i in parse_fasta(fh):
            head, seq = i
            tmp = head.split("_")
            tag1 = int(tmp[0][-3:]) - 1
            tag2 = int(tmp[1]) - 1
            ab = tmp[2]
            if "For" in head:
                if tag1 not in seqsf.keys():
                    seqsf[tag1] = []
                seqsf[tag1].append(seq)
                abundance[seq] = ab
            elif "Rev" in head:
                if tag1 not in seqsr.keys():
                    seqsr[tag1] = []
                seqsr[tag1].append(seq)
                abundance[seq] = ab
            else:
                print("[ERROR]: wrong format in head of ends fasta")
                exit()
    sub_len = int(len(seqsf) / subsam)
    current_len = 0
    current_cont = ""
    file_mark = 1
    for s in range(96):
        #sample_id = "{0:.3d}".format(s+1)
        sample_id = format(s+1, '03d')
        current_len += 1
        for f in range(len(seqsf[s])):
            fs = f + 1
            abf = abundance[seqsf[s][f]]
            for r in range(len(seqsr[s])):
                rs = r + 1
                abr = abundance[seqsr[s][r]]
                end = ">{0}_{1}-{2};{3}-{4}\n".format(
                    sample_id, fs, rs, abf, abr)
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

#--------------------------------------------------------
if args.command in ["all", "filter"]:

    filtered_outfile1 = args.outpre + "_filter_highqual_1.fastq"
    filtered_outfile2 = args.outpre + "_filter_highqual_2.fastq"
    out1 = check_and_open_outhandle(filtered_outfile1)
    out2 = check_and_open_outhandle(filtered_outfile2)

    # ini state
    total = 0
    clean = 0
    nn = 0
    low = 0

    log = check_and_open_outhandle(args.outpre + "_filter_log.txt")

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

    print_time("[INFO]: Filtering start:")

    if args.fastq1.endswith("gz"):
        fq1 = gzip.open(args.fastq1, 'rt')
    else:
        fq1 = open(args.fastq1, 'r')
    if args.fastq2.endswith("gz"):
        fq2 = gzip.open(args.fastq2, 'rt')
    else:
        fq2 = open(args.fastq2, 'r')

    for i in parse_pe_fastq(fq1, fq2):
        name, seq1, seq2, qual1, qual2 = i

        if not name[0].startswith("@"):
            print("[ERROR]: input fastq is not a correct fastq format")
            exit()

        total += 1
        N_count1 = seq1.count("N")
        N_count2 = seq2.count("N")
        qual_str1 = "".join(qual1)
        qual_str2 = "".join(qual2)

        if N_count1 < args.n or N_count2 < args.n:
            if filter_type == 1:
                if (exp_e(qual1, args.phred) <= args.expected_err
                    and exp_e(qual2, args.phred) <= args.expected_err):
                    out1.write(name + " 1\n" + seq1 + "\n" + "+\n" + qual_str1 + "\n")
                    out2.write(name + " 2\n" + seq2 + "\n" + "+\n" + qual_str2 + "\n")
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

#----------------assign------------------------------------------------------------

if args.command in ["all", "assign"]:
    print_time("[INFO]: Assigning start:")

    if args.command == "all":
        args.fq1 = filtered_outfile1
        args.fq2 = filtered_outfile2
        assigned_outdir = os.path.abspath(args.outpre + "_assign")
    elif args.command == "assign":
        assigned_outdir = os.path.abspath(args.outdir)

    ErrFile = check_and_open_outhandle(assigned_outdir + "_err.fasta")
    Middle_F = check_and_open_outhandle(assigned_outdir + "_assign_midF.fasta")
    Middle_R = check_and_open_outhandle(assigned_outdir + "_assign_midR.fasta")

    indexlen = args.index

    if os.path.exists(assigned_outdir) == False:
        os.mkdir(assigned_outdir)

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
    print("min distance among barcodes is {}".format(min_dis))
    print("max distance among barcodes is {}".format(max_dis))
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
    assigned_list =  args.outpre + "_assign.list"

    with open(assigned_list, "w") as ls:
        sorted_sample = sorted(pris.keys())
        for s in sorted_sample:
            ls.write(
                assigned_outdir
                + "/For"
                + s
                + ".ssam"
                + "\n"
                + assigned_outdir
                + "/Rev"
                + s
                + ".ssam"
                + "\n"
            )
    # open all assigned files
    filehandle = {}
    for sam in indp.keys():
        filehandle[sam] = open(assigned_outdir + "/" + sam + ".ssam", "w")

    seqnum = 0
    err = 0
    assigned = 0

    # open clean paired fastq
    if args.fq1.endswith(".gz"):
        fh1 = gzip.open(args.fq1, "rt")
    else:
        fh1 = open(args.fq1, "r")

    if args.fq2.endswith(".gz"):
        fh2 = gzip.open(args.fq2, "rt")
    else:
        fh2 = open(args.fq2, "r")

    for i in parse_pe_fastq(fh1, fh2):
        name, seq1, seq2, qual1, qual2 = i
        qual_str1 = "".join(qual1)
        qual_str2 = "".join(qual2)
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
                    + seq1 + " " + qual_str1 + " ")
                filehandle[targetf].write(
                    seq2 + " " + qual_str2 + "\n")
            else:
                filehandle[targetf].write(
                    targetf + "_" + str(seqnum) + " "
                    + seq2 + " " + qual_str2 + " ")
                filehandle[targetf].write(
                    seq1 + " " + qual_str1 + "\n")

        elif (rmatch or rmatch_withMis) and (fmatch == False
                                             and fmatch_withMis == False):
            assigned += 1
            targetr = rmatch or rmatch_withMis
            count_assigned[targetr] += 1
            if targetr[:3] == "For":
                filehandle[targetr].write(
                    targetr + "_" + str(seqnum) + " "
                    + seq2 + " " + qual_str2 + " ")
                filehandle[targetr].write(
                    seq1 + " " + qual_str1 + "\n")
            else:
                filehandle[targetr].write(
                    targetr + "_" + str(seqnum) + " "
                    + seq1 + " " + qual_str1 + " ")
                filehandle[targetr].write(
                    seq2 + " " + qual_str2 + "\n")
        elif (fmatch == False
              and fmatch_withMis ==False
              and rmatch == False
              and rmatch_withMis == False):
            assigned += 1
            count_assigned['middle'] += 1
            Middle_F.write(
                ">" + str(seqnum) + "_1\n" + seq1 + "\n"
            )
            Middle_R.write(
                ">" + str(seqnum) + "_2\n" + seq2 + "\n"
            )
        else:
            err += 1

    ErrFile.close()
    # close all assigned files
    for fh in filehandle.values():
        fh.close()

    # report assignment information
    with open(args.outpre + ".assign.log", "w") as log:
        log.write("total reads:\t{}\n".format(seqnum))
        log.write("err reads:\t{}\n".format(err))
        log.write("assigned:\t{}\n".format(assigned))
        for i in sorted(count_assigned.keys()):
            log.write(
                i + "\t" + str(count_assigned[i]) + "\n"
            )

    print_time("[INFO]: Assigning done:")

#-----------------------------------------------------------------------------
if args.command in ["all", "buildend"]:
    print_time("[INFO]: Building ends start:")

    if args.min_overlap > 120:
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

    buildends_result = args.outpre + "_buildends.fasta"
    fh_out = check_and_open_outhandle(buildends_result)
    fh_log = check_and_open_outhandle(args.outpre + "_buildends.log")

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

    fh_log.write("## consensus mode = " + str(args.mode) + "\n")

    if args.mode == 1:
        fh_log.write("## clustering identity = "
                     + str(args.cluster_identity)
                     + "\n")

    fh_log.write("## overlaping identity = "
                 + str(args.overlap_identity)
                 + "\n")
    fh_log.write("## min overlap = " + str(args.min_overlap) + "\n")
    fh_log.write("## max overlap = " + str(args.max_overlap) + "\n")

    # --------------main-----------------------#

    try:
        with open(args.list) as fh_list:
            lines = fh_list.readlines()
    except FileNotFoundError:
        print("[ERROR]: can not find " + args.list)
        exit(0)

    for line in lines:
        line = line.rstrip()
        name = os.path.basename(line).split(".")[0]
        short_outname = name[:6]
        fh_log.write("//processing " + name + " done\n")
        # if file is empty, continue
        if os.path.getsize(line) == 0:
            fh_log.write("! file is empty!")
            fh_log.write("\n")
            continue
        success_connected = []
        with open(line,'r') as fh:
            records = fh.readlines()
            for r in records:
                both_ends = r.split()
                if len(both_ends) != 5:
                    continue
                name = both_ends[0]
                forward_read = both_ends[1]
                forward_qual = both_ends[2]
                reverse_read = both_ends[3]
                reverse_qual = both_ends[4]
                reverse_read = comp_rev(reverse_read)
                reverse_qual = reverse_qual[::-1]

                ##-----------anchoring overlap site--------#

                read0 = forward_read[-args.max_overlap :]
                read1 = reverse_read[0 : args.max_overlap]

                singal = 0
                overlaps = {}
                for s in range(args.min_overlap, args.max_overlap + 1):
                    l0 = read0[-s:]
                    l1 = read1[0:s]
                    tmp_identity = match(l0, l1)
                    if tmp_identity == 1:
                        overlaps[s] = 1
                        # find best result, so exit loop #
                        break

                    elif tmp_identity >= args.overlap_identity:
                        overlaps[s] = tmp_identity

                # find best overlaping result in all potenial positions
                # candidates = sorted(overlaps.items(),
                # lambda x, y: cmp(x[1], y[1]), reverse=True)
                candidates = sorted(
                    overlaps, key=overlaps.__getitem__, reverse=True
                )

                if len(candidates) > 0:

                    potenial = candidates[0]  # overlap similarity top 1
                    s0 = read0[-potenial:]
                    s1 = read1[0:potenial]

                    corrected = ""

                    # compare each base from forward and reverse to keep one
                    # forward == reverse
                    # forward ne reverse, quality(forward) > quality(reverse)
                    # forward ne reverse, quality(forward) < quality(reverse)

                    for p in range(len(s0)):
                        # site is changed, be careful!#
                        tmp_loca0 = args.standard_length - potenial + p

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
                        forward_read[: args.standard_length - potenial]
                        + corrected
                        + reverse_read[potenial - args.standard_length :]
                    )

                    len_makeup_consensus = len(makeup_consensus)
                    this_oid = overlaps[potenial] * 100
                    this_oid = str("%.2f" % this_oid)

                    # if check result, and ok so write into output_checked file #
                    success_connected.append(makeup_consensus)
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
            + " --cluster_fast "
            + temp_fasta
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
                if line[0] is not "H":
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

    fh_out.close()
    fh_log.close()
    rm_tmp_cmd = "rm temp.fa.* temp.uc.*"
    os.system(rm_tmp_cmd)

    print_time("[INFO]: Building ends done:")

# connect middle fasta
if args.command in ["all", "chain"]:
    print_time("[INFO]: Middle Chaining start:")
    if args.command == "all":
        read1 = Middle_F
        read2 = Middle_R
    else:
        read1 = args.middle_R1
        read2 = args.middle_R2
    failed1 = args.outpre + "_failed1.fa"
    failed2 = args.outpre + "_failed2.fa"
    chain_out = args.outpre + "_chainout.fa"
    chain_out = os.path.abspath(chain_out)
    chain_cmr = cmr +\
                " -a " + read1 +\
                " -b " + read2 +\
                " -2 " + failed1 +\
                " -3 " + failed2 +\
                " -o " + chain_out +\
                " -m 0 >cec.log 2>cec.error"
    print("Run:" + chain_cmr)
    subprocess.call(chain_cmr, shell=True)
    print_time("[INFO]: Middle Chaining done:")

# gap filling
if args.command in ["all", "gapfill"]:
    if check_program_involed("barcode"):
        print("can not find soapbarcode program in $BIN path")
        exit()
    print_time("[INFO]: Gap filling start:")
    if args.command == "all":
        ends = buildends_result
        middle_lis = os.path.abspath(args.outpre + "_middle_lis")
    else:
        ends = args.gapfill_ends
        chain_out = os.path.abspath(args.gapfill_mid)
        middle_lis = os.path.abspath(args.outpre + "_middleFasta.lis")
    with open(middle_lis, 'w') as mf:
        mf.write(">\n"
                 +"f=" + chain_out)


    # outdir
    barcode_outdir = args.outpre + "_barcode"
    result_outdir = barcode_outdir + "/result"
    split_ends_dir = barcode_outdir + "/ends"
    shell_outdir = barcode_outdir + "/shell"
    for d in [barcode_outdir, result_outdir, split_ends_dir, shell_outdir]:
        if os.path.exists(d) == False:
            os.mkdir(d)

    # format ends fasta
    # split ends fasta to several subfiles
    split_ends = format_ends(ends, args.samp_num, split_ends_dir)
    subfile = 0
    for s in split_ends:
        s = os.path.abspath(s)
        subfile += 1
        shell = shell_outdir + "/barcodes." + str(subfile) + ".sh"
        barcode_out = result_outdir + "/barcodes." + str(subfile) + ".fa"
        barcode_out = os.path.abspath(barcode_out)
        with open(shell, 'w') as sh:
            sh.write("barcode"
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
    print("Now you can run shell files in " + shell_outdir)
    if args.samp_num > 1:
        print("And, finally run: cat " + result_outdir + "/*.fa >all.barcodes.fa")
    else:
        print("the final barcodes in " + result_outdir)

