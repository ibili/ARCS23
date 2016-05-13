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
def getReverseFile(f):
    line = f.readline()
    fout = open( sys.argv[2], "w" )
    i = 0
    while line:
        fout.write(">read_%d" % i)
        line = f.readline()
        fout.write("_" + str(len(line.strip())) + "\n")
        fout.write(line.strip() + "\n")
        line = f.readline().strip()
        i += 1
    f.close()
    f = open( sys.argv[1], "r" )
    line = f.readline()
    #i = 1
    while line:
        fout.write(">read_%d" % i)
        line = f.readline() 
        fout.write("_" + str(len(line.strip())) + "\n")
        fout.write(complement(line.strip()) + "\n")
        line = f.readline().strip();
        i += 1
if len( sys.argv ) != 3:
    print "Usage: " + sys.argv[0] + " pacbio_read.fasta reverse_pabio_read.fasta"
    exit(0)

with file( sys.argv[1] ) as f:
    getReverseFile(f)
