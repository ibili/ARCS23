#!/usr/bin/python

import sys
EDGE_LENGTH_CUTOFF = 100
COVRAGE_PERCENT_CUTOFF = 0.9 # to ful, no need change
#COVER_NUM_CUTOFF = 10
COVER_NUM_CUTOFF = 3 # to ful, no need change
SIM_CUTOFF = 85
SCORE_CUTOFF = 0
#SIM_UPPER = 99


###pay attention contigsId from 0 or 1
def read_contigLen(f):
    contigLen = []
    line = f.readline().strip()
    while line:
        contigLen.append( int(line) );
        line = f.readline().strip()
    return contigLen

def read_component2Contig(f):
    with file(sys.argv[1] ) as contigFile:
        contigLen = read_contigLen(contigFile)
    components = {}
    line = f.readline().strip()
    while line:
        componentId = int(line.split()[1].strip())
        contigs = f.readline().strip().split()
        gaps    = f.readline().strip().split()
        if not components.has_key(componentId):
            components[ componentId ] = [];
        components[ componentId ].append( (int(contigs[0]), 0) )
        last = 0
        for i in range(1, len(contigs)):
            cur = last + contigLen[ int(contigs[i-1]) - 1] + gaps[i-1] #in contigs.fasta contigId from 1
            components[ componentId ].append( (int(contigs[i]), cur) )
            last = cur
        line = f.readline().strip()
    return components

def read_component(f):
    with file(sys.argv[1] ) as contigFile:
        contigLen = read_contigLen(contigFile)
    fout = open(sys.argv[6], "w")
    #print contigLen
    contigInComp = {}
    line = f.readline().strip()
    componentNum = 0
    while line:
        componentId = int(line.split()[1].strip())
        componentNum = componentId
        contigs = f.readline().strip().split()
        gaps    = f.readline().strip().split()
        contigInComp[ int(contigs[0]) ] = (componentId, 0)
        last = 0
        for i in range(1, len(contigs)):
            cur = last + contigLen[ int(contigs[i-1]) - 1 ] + int(gaps[i-1])
            contigInComp[ int(contigs[i]) ] = (componentId, cur)
            last = cur
        #print "contigs size, contigLen size",len(contigs), len(contigLen) 
        fout.write(str(last + contigLen[ int(contigs[ len(contigs)-1 ]) - 1 ]) + "\n")
        line = f.readline().strip()
    fout = open(sys.argv[7], "w")
    fout.write("EDGE_CLUSTER_NUM=%d\n" % (componentNum+1))
    return contigLen, contigInComp, componentNum+1
    
def read_blasr(f):
    print "yes1"
    with file( sys.argv[2] ) as componentFile:
        contigLen, contigInComp, componentNum = read_component(componentFile)
    print contigInComp
    print componentNum
    graph = {}
    read2Contigs = {}
    line = f.readline().strip()
    
    #print "contigId, readId , mapLength, contigStart, contigEnd, contigLength, readStart, readEnd, sim, score"
    while line:
        contigId = int(line.split("\t")[0].split("_")[1].strip())
        readId = int(line.split("\t")[1].split("_")[1].strip())
        
        sim = float(line.split("\t")[2].strip())         
        mapLength = int(line.split("\t")[3].strip())
        score = float(line.strip().split("\t")[-1])
        
        readStart = int(line.split("\t")[8].strip())
        readEnd   = int(line.split("\t")[9].strip())
        readLength = int(line.split("\t")[1].split("_")[2].strip())


        contigStart = int(line.split("\t")[6].strip())
        contigEnd   = int(line.split("\t")[7].strip())
        contigLength = int(line.split("\t")[0].split("_")[3].strip())
        
        readDirection = readStart < readEnd 
        contigDirection = contigStart < contigEnd
        
        if (readDirection != 1 or contigDirection != 1):
            line = f.readline();
            continue; 
        
        if (score < SCORE_CUTOFF or mapLength < EDGE_LENGTH_CUTOFF or sim < SIM_CUTOFF):
            line = f.readline()
            continue
        shouldLength = mapLength + min(contigStart, readStart) + min(contigLength-contigEnd, readLength-readEnd) 
        # actual map / should map >  COVRAGE_PERCENT_CUTOFF
        if (mapLength < shouldLength * COVRAGE_PERCENT_CUTOFF):
            line = f.readline()
            continue
        
        if not read2Contigs.has_key(readId):
            read2Contigs[ readId ] = [];
        read2Contigs[ readId ].append(( contigId, readStart - contigStart))
        line = f.readline()
    
    for idx in read2Contigs:
        read2Contigs[idx].sort(key=lambda x:x[1])
        if len(read2Contigs[idx]) >= 2:
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
                    dis = readPos_j - readPos_i + comPos_i - comPos_j
                    key = '%d_%d' % (com_i, com_j)
                    if not graph.has_key(key):
                        graph[key] = (dis, 1) 
                        print "com_i,com_j,readId,dis,cov:",com_i, com_j, idx, graph[key],graph[key][0]/graph[key][1]
                    else:
                        x, y = graph[key]
                        graph[key] = (x + dis, y + 1)
                        print "com_i,com_j,readId,dis,cov:",com_i, com_j, idx, graph[key],graph[key][0]/graph[key][1]
    eleForDelete = []
    print "before coverage filter",graph 

    for key in graph:
        components = key.split("_")
        com_i = int(components[0])
        com_j = int(components[1])
        dis, cov = graph[key]
        if cov < COVER_NUM_CUTOFF:
            eleForDelete.append(key)
    for key in eleForDelete:
        del graph[key]
    print "SPA add edge",len(graph)
    
    fout = open("subgraph","w")
    fout.write("digraph {\n")
    for key in graph:
        fout.write("%s->%s[label=\"%s\"]\n" % (key.split("_")[0] , key.split("_")[1] , graph[key]))
        print key,graph[key]
    fout.write("}")



    return graph, componentNum

if len( sys.argv) != 8:
    print "Usage: " + sys.argv[0] + " contig_len component_n blasr.out contig_arc_graph_after_remove_ambigous_arcs_n position_lp_n.math edge_cluster_len_n scaffold_parameter_n"
    exit(0)

'''
with file(sys.argv[1]) as f:
    contigLen = read_contigLen(f)
print contigLen
'''
'''
with file(sys.argv[2]) as f:
    contigInComp, componentNum = read_component(f)
print contigInComp
print componentNum
'''
'''
with file(sys.argv[3]) as f:
    graph, componentNum = read_blasr(f)
print graph
print componentNum
'''
with file(sys.argv[3]) as f:
    graph, componentNum = read_blasr(f)
    '''
    #filter degree too larger
    sortedGraphKey = sorted(graph.keys(), key=lambda t:int(t.split("_")[0]))
    #print sortedGraphKey
    comFirst = {}
    for key in sortedGraphKey:
        keyFirst = key.split("_")[0]
        if keyFirst in comFirst:
            comFirst[keyFirst] += 1
        else:
            comFirst[keyFirst] = 1
    print comFirst
    comFirstFilter = []
    for ii in comFirst:
        if comFirst[ii] >= 4:
            comFirstFilter.append(ii)
    print comFirstFilter

    comSecond = {}
    for key in graph.keys():
        keySecond = key.split("_")[1]
        if keySecond in comSecond:
            comSecond[keySecond] += 1
        else:
            comSecond[keySecond] = 1
    print comSecond
    comSecondFilter = []
    for ii in comSecond:
        if comSecond[ii] >= 5:
            comSecondFilter.append(ii)
    print comSecondFilter

    for key in graph.keys():
        if (key.split("_")[0] in comFirstFilter) and key.split("_")[0] in comSecondFilter:
            del graph[key]
            continue
        if (key.split("_")[1] in comFirstFilter) and key.split("_")[1] in comSecondFilter:
            del graph[key]
    #print graph        
    fout = open("copy_num","w")
    for key in graph.keys():
        fout.write(str(graph[key][1]) + "\n")
    '''
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



