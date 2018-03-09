#!/usr/bin/env python

import pysam
import sys

def print_fastq(read, out_r1, out_r2):
    out_fh = out_r1 if read.is_read1 else out_r2
    out_fh.write("@{} {}:N:0\n".format(read.query_name, '1' if read.is_read1 else '2' ))
    out_fh.write("{}\n".format(read.query_sequence))
    out_fh.write("+\n")
    qual = "".join([ chr(c+33) for c in read.query_qualities])
    out_fh.write("{}\n".format(qual))


samples = ["Sample_75641", "Sample_79162"]
regions = [("chr38", 1000000, 1030000),
           ("chr37", 1000000, 1030000),
           ("chr36", 1000000, 1030000),
           ("chrX",  1000000, 1030000)]

for sample in samples:
    for region in regions:
        print("Processing {} {}".format(sample, region))
        samfile = pysam.AlignmentFile("{}.bam".format(sample), "rb")

        out_r1="{}_{}_{}_{}_R1.fq".format(sample, *region)
        out_r2="{}_{}_{}_{}_R2.fq".format(sample, *region)


        fh1 = open(out_r1, 'w')
        fh2 = open(out_r2, 'w')

        iter = samfile.fetch(*region)

        reads = dict()

        for x in iter:
            read = '1' if x.is_read1 else '2'
            if x.query_name not in reads:
                reads[ x.query_name ] = dict()
            reads[ x.query_name ][read] = x

        search_for = []

        for (r,v) in reads.items():
            ms = list(v.keys())
            if len(ms) < 2:
                search_for.append(list(v.values())[0])
                continue
            for read in v.values():
                print_fastq(read, fh1, fh2)

        for s in search_for:
            try:
                m = samfile.mate(s)
                print_fastq(s, fh1, fh2)
                print_fastq(m, fh1, fh2)
            except ValueError as e:
                print("  Exception: {}".format(e))

        fh1.close()
        fh2.close()
