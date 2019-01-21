#!/usr/bin/env python
# --*-- coding:utf-8 --*--

import os
import sys
import time
import gzip

def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def parse_pe_fasta(fa1, fa2):
    while True:
        name1 = fa1.readline().split("_")[0]
        name2 = fa2.readline().split("_")[0]

        read1 = fa1.readline()[:-1]
        read2 = fa2.readline()[:-1]
        if not name1 or name2:
            break
        assert name1 == name2, 'fastq1, fastq2 is not paired-end'
        yield name1, read1, read2

def reverseComplement(s):
    complement = { 'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A',
                  'a' : 'T', 'g' : 'C', 'c' : 'G', 't': 'A',
                  'N' : 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def match(str1, str2):
    # ----count matched bases of two sequences----#
    matched = 0
    for base in range(len(str1)):
        if str1[base] == str2[base]:
            matched += 1
    identity = matched / len(str1)
    return identity

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# -----------------------------------------------------------


if len(sys.argv) < 4:
    usage = "Usage:\n\tpython3 {0} <fasta1> <fasta2> <output>\n".format("cmr.py")
    print(usage)
    exit()
else:
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    min_overlap = 30
    max_overlap = 120
    overlap_identity = 0.98
    standard_length = 150

print_time("[INFO]: start running:")


with open_input(fasta1) as fa1:

    with open_input(fasta2) as fa2:

        output = open(sys.argv[3] ,'w')
        for i in parse_pe_fasta(fa1, fa2):
            name, seq1, seq2 = i

            if len(seq1) == 0 or len(seq2) == 0:
                print("[warnning]: " + name + " with wrong length!")
                continue
            seq2 = reverseComplement(seq2)

            singal = 0
            overlaps = {}
            for s in range(min_overlap, max_overlap + 1):
                l0 = seq1[-s:]
                l1 = seq2[0:s]
                tmp_identity = match(l0, l1)
                if tmp_identity == 1:
                    overlaps[s] = 1
                    # find best result, so exit loop #
                    break

                elif tmp_identity >= overlap_identity:
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
                        #for_quality = forward_qual[tmp_loca0]
                        #rev_quality = reverse_qual[p]
                        #if for_quality >= rev_quality:
                        corrected += s0[p]
                        #else:
                        #    corrected += s1[p]

                makeup_consensus = (
                    seq1[: standard_length - potenial]
                    + corrected
                    + seq2[potenial - standard_length :]
                )

                len_makeup_consensus = len(makeup_consensus)
                this_oid = overlaps[potenial] * 100
                this_oid = str("%.2f" % this_oid)

                print(name + "\n" + makeup_consensus, file=output)

        output.close()
        print_time("[INFO]: finish running:")
