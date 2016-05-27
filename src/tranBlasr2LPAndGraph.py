#!/usr/bin/python

import sys
EDGE_LENGTH_CUTOFF = 150
CORVAGE_PERCENT_CUTOFF = 0.8
#COVER_NUM_CUTOFF = 10
COVER_NUM_CUTOFF = 5
SIM_CUTOFF = 85
SCORE_CUTOFF = -600
#SIM_UPPER = 99


###pay attention contigsId from 0 or 1
def read_contigLen(f):
    contigLen = []
    line = f.readline().strip()
    while line:
        contigLen.append( int(line) );
        line = f.readline().strip()
    return contigLen


def read_component(f):
    with file(sys.argv[1] ) as contigLenFile:
        contigLen = read_contigLen(contigLenFile)
    fout = open(sys.argv[6], "w")
    #print "contigId contigLength"
    #print contigLen
    contigInComp = {}
    line = f.readline().strip()
    componentNum = 0
    while line:
        componentId = int(line.split()[1].strip())
        componentNum = componentNum + 1
        contigs = f.readline().strip().split()
        gaps    = f.readline().strip().split()
        contigInComp[ int(contigs[0]) ] = (componentId, 0)
        last = 0
        for i in range(1, len(contigs)):
            cur = last + contigLen[ int(contigs[i-1]) ] + int(gaps[i-1])
            contigInComp[ int(contigs[i]) ] = (componentId, cur)
            last = cur
        
        #contigId is from 1 or from 0
        #fout.write(str(last + contigLen[ int(contigs[ len(contigs)-1 ]) -1 ]) + "\n")
        fout.write(str(last + contigLen[ int(contigs[ len(contigs)-1 ])]) + "\n")
        line = f.readline().strip()
    fout = open(sys.argv[7], "w")
    fout.write("EDGE_CLUSTER_NUM=%d\n" % (componentNum))
    return contigLen, contigInComp, componentNum
    
def read_blasr(f):
    with file( sys.argv[2] ) as componentFile:
        contigLen, contigInComp, componentNum = read_component(componentFile)
    graph = {}
    read2Contigs = {}
    line = f.readline().strip()
    print "readId , contigId, contig in read, readStart, readEnd, readLength, contigStart, contigEnd, contigLength, sim, score"
    while line:
        fileEnd = False
        while line and ( not line.strip().startswith("%sim") ):
            line = f.readline()
            if not line:
                fileEnd = True
                break
        if fileEnd:
            break
        sim = line.strip().split(" ")[1]
        line = f.readline();
        score = int(line.strip().split(" ")[1])
        line = f.readline().strip()
        readId = int(line.split("/")[0].split("d")[1].strip())
        line = f.readline().strip()
        contigId = int(line.split("_")[1].strip())
        for i in range(0, 3):
            f.readline()
        line = f.readline().strip()
        readDirection = int(line.split(":")[1].strip()) 
        line = f.readline().strip()
        contigDirection = int(line.split(":")[1].strip())
        if (readDirection !=0 or contigDirection != 0):
            line = f.readline();
            continue;
        line = f.readline().strip()
        readStart = int(line.split()[1].strip())
        readEnd   = int(line.split()[3].strip())
        readLength = int(line.split()[5].strip())
        line = f.readline().strip()
        contigStart = int(line.split()[1].strip())
        contigEnd   = int(line.split()[3].strip())
        contigLength= int(line.split()[5].strip())
        #print "sim >="
        #print SIM_CUTOFF  
        if (score > SCORE_CUTOFF or contigLength < EDGE_LENGTH_CUTOFF or sim < SIM_CUTOFF):
            line = f.readline()
            continue
        #if (readStart > 5 and readLength-readEnd > 5) and ((contigEnd - contigStart) < contigLength * CORVAGE_PERCENT_CUTOFF or contigLength < EDGE_LENGTH_CUTOFF or sim < SIM_CUTOFF):
        #    line = f.readline()
        #    continue

        #print "after filter"
        if not read2Contigs.has_key(readId):
            read2Contigs[ readId ] = [];
        if readId == 100:
            print readId,"\t", contigId,"\t", readStart -contigStart,"\t",readStart, "\t", readEnd,"\t", readLength, " ", contigStart, " ", contigEnd, " ", contigLength, " ", sim, " ", score 
        read2Contigs[ readId ].append(( contigId, readStart - contigStart))
        line = f.readline()
    
    #print read2Contigs
    #print "readId , contigId, contig in read"
    for idx in read2Contigs:
        read2Contigs[idx].sort(key=lambda x:x[1])
        for i in range(len( read2Contigs[idx] )):
            contig_i, readPos_i = read2Contigs[idx][i]
            if not contigInComp.has_key( contig_i ):
                continue
            com_i, comPos_i = contigInComp[ contig_i ]
            for j in range(i+1, len( read2Contigs[idx] )):
                contig_j, readPos_j = read2Contigs[idx][j]
                if not contigInComp.has_key( contig_j ):
                    continue
                com_j, comPos_j = contigInComp[ contig_j ]
                if com_j == com_i:
                    continue
                #if readPos_i == readPos_j and contigLen[contig_i] == contigLen[contig_j]:  // read same position, length equal, must be forward backward
                    #continue
                #if contigLen[ contig_i ] + readPos_i < readPos_j:
                    #continue
                
                dis = readPos_j - readPos_i + comPos_i - comPos_j
                #if com_i < com_j:
                key = '%d_%d' % (com_i, com_j)
                if not graph.has_key(key):
                    graph[key] = (dis, 1)
                    if (idx ==  100):
                        print "com_i,com_j,readId,dis,cov:",com_i, com_j, idx, graph[key], graph[key][0]/graph[key][1]
                else:
                    x, y = graph[key]
                    graph[key] = (x + dis, y + 1)
                    if (idx ==  100):
                        print "com_i,com_j,readId,dis,cov:",com_i, com_j, idx, graph[key], graph[key][0]/graph[key][1]
    eleForDelete = []
    #print graph["2_515"]
    for key in graph:
        components = key.split("_")
        com_i = int(components[0])
        com_j = int(components[1])
        dis, cov = graph[key]
        if cov < COVER_NUM_CUTOFF:
            eleForDelete.append(key)
    for key in eleForDelete:
        del graph[key]
    #print graph["2_515"]
    return graph, componentNum

if len( sys.argv) != 8:
    print "Usage: " + sys.argv[0] + " contig_len component_n blasr.out contig_arc_graph_after_remove_ambigous_arcs_n position_lp_n.math edge_cluster_len_n scaffold_parameter_n"
    exit(0)

with file(sys.argv[3]) as f:
    graph, componentNum = read_blasr(f)
    print "graph"
    for key in graph:
        print key,graph[key]
    
    '''
    #filter degree too larger
    #see newARCS3
    '''
#write LP
fout = open(sys.argv[4], "w")
for key in graph:
    components = key.split("_")
    com_i = components[0]
    com_j = components[1]
    dis, cov = graph[key]
    fout.write(com_i + "\t" + com_j + "\t" + str(int(dis/cov)) + "\t" + str(cov) + "\t0.0\n")
fout = open(sys.argv[5], "w")
for i in range(componentNum):
    fout.write("var x_%d;\n" % i)
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write("var e_%d_%d;\n" % (com_i, com_j))
    fout.write("var E_%d_%d;\n" % (com_i, com_j))
fout.write("minimize z: \n")
count = 0
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write(" E_%d_%d + " % (com_i, com_j))
    count += 1
    if count % 10 == 0:
        fout.write("\n")
fout.write("0;\n")

index = 1
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    dis, cov = graph[key]
    fout.write("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;\n" % (index , com_j , com_i , com_i , com_j , int(dis/cov)))
    index += 1

for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write("s.t. con%d : E_%d_%d + e_%d_%d >= 0;\n" % (index , com_i , com_j, com_i , com_j))
    index += 1
    fout.write("s.t. con%d : E_%d_%d - e_%d_%d >= 0;\n" % (index , com_i , com_j, com_i , com_j))
    index += 1
fout.write("end;\n")



