#!/usr/bin/python

import sys

#read contigs.fasta , write uniqueContig.fasta, contig_len, component_0 

def read_contig(f,threshold):
    uniOut = open(sys.argv[2], "w")
    lenOut = open(sys.argv[3], "w")
    comOut = open(sys.argv[4], "w")
    line = f.readline().strip()  
    state = 0
    contig = {}
    label = {}
    coverage = {}
    seq = "" 
    while line:
        #print line
        cols = line.split("_")
        if len(cols) > 2:
            if len(seq) > 0:
                #print "test:",contigId
                contig[contigId] = seq
                #print len(contig)
            state = 0
            seq = ""
        #print "aa",len(cols)
        if state == 0 and len(cols) > 2:
            #print cols
            contigId = int(cols[1])
            label[contigId] = line
            coverage[contigId] = float(cols[5]) 
            state = 1
        if state == 1 and len(cols) == 1:
            seq = seq + line
        line = f.readline().strip()
    if contigId not in contig:
        contig[contigId] = seq
    #print "final contig", contig 
    #print "final coverage", coverage
    num = 1
    print "contig number",len(contig)
    for key in contig:
        print key,len(contig[key])
        lenOut.write("%d %d\n" % (key,len(contig[key])))
        # select unique strategy ????
        #if ( len(contig[key])>200 ):
        if (coverage[key] <= threshold and len(contig[key])>200) or len(contig[key]) > 10000:
        #if (coverage[key] <= threshold and len(contig[key])>200):
        #if (coverage[key]/len(contig[key]) <= threshold/1000 and len(contig[key])>200):
            uniOut.write("%s\n" % label[key])
            uniOut.write("%s\n" % contig[key])
            #comOut.write(">component  %d\n" % num)
        
            if key <= len(contig)/2:
                comOut.write(">component  %d\n" % num)
            elif key == len(contig)/2 + 1:
                num = 1
                comOut.write(">component  %d\n" % -num)
            else:
                comOut.write(">component  %d\n" % -num)
            
            comOut.write("%d\n\n" % key)
            num = num + 1
if len( sys.argv) != 5:
    print "Usage: " + sys.argv[0] + " contigs.fasta, uniqueContig.fasta, contig_len, component_0"
    exit(0)

#get average,median
with file(sys.argv[1]) as f:
     line = f.readline().strip()
     coverage = []
     while line:
         cols = line.split("_") 
         #print len(cols)
         if len(cols) > 2:
             coverage.append(float(cols[5]))
             #print float(cols[5]),len(coverage)
         line = f.readline().strip()
average = sum(coverage)/len(coverage)
median = sorted(coverage)[len(coverage)/2]

print "average", average
print "median", median

#simple think: copy number less than median is unique, can be changed
with file(sys.argv[1]) as f:
    read_contig(f,average)


#./blasr contig.fasta contig.fasta -out bb -m 0 -nproc 8
