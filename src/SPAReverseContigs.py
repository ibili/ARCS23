#!/usr/bin/python

import sys

def complement(read):
    ret = ""
    for c in read:
        if c == "A":
            ret += "T"
        elif c == "C":
            ret += "G"
        elif c == "G":
            ret += "C"
        elif c == "T":
            ret += "A"
        else:
            ret += "N"
    return ret[::-1]

def read_contig(f):
    line = f.readline().strip()  
    state = 0
    contig = {}
    label = {}
    seq = "" 
    fout = open(sys.argv[2], 'w')
    while line:
        cols = line.split("_")
        if len(cols) > 2:
            if len(seq) > 0:
                contig[contigId] = seq
            state = 0
            seq = ""
        if state == 0 and len(cols) > 2:
            contigId = int(cols[1])
            label[contigId] = line
            state = 1
        if state == 1 and len(cols) == 1:
            seq = seq + line
        line = f.readline().strip()
    if contigId not in contig:
        contig[contigId] = seq
    print "contig number",len(contig)
    for key in contig:
        fout.write("%s\n" % label[key])
        fout.write("%s\n" % contig[key])
    for key in contig: 
        newLabel = '_'.join(label[key].split('_')[2:])
        fout.write(">NODE_%d_%s\n" % (key + len(contig), newLabel))
        fout.write("%s\n" % complement(contig[key]))

if len( sys.argv ) != 3:
    print "Usage: " + sys.argv[0] + " contigs.fasta reContigs.fasta"
    exit(0)

with file( sys.argv[1] ) as f:
    read_contig(f)
