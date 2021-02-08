################################################################################
#required packages:
#   clustal omega
#   paup
#   python: sys, skbio, Bio, scipy, re, math, os
#
#ARGS:
#1: input/output dir
#2: homolog sequence file
#3: target sequence fasta identification
#4: protein or dna; default dna
#
#Example:
#julia evolutionAnalysis.jl FULL_PATH/example APP_refseq_protein.fasta NP_000475.1 protein
#
################################################################################
import sys
import re
import math
import os
from skbio import TreeNode
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
import numpy as np
from scipy import stats
import random


dir = sys.argv[1]
seqType = sys.argv[4]
accession_full = sys.argv[3]
print(seqType)

'''
replace the amino acid with X if it is not among the 20 amino acid codes
'''
def replaceAA(sequence):
    seq = ''
    for char in sequence:
        if char in 'ARNDCQEGHILKMFPSTWYV':
            seq = seq + char
        else:
            seq = seq + 'X'
    return seq

'''
write a sequence to fasta file
'''
def write_fasta(file, accession, sequence):
    file.write('>{}\n'.format(accession))
    if len(sequence) <= 70:
        file.write('>{}\n'.format(sequence))
    else:
        line_num_exact = len(sequence) / 70
        line_num = math.floor(line_num_exact)
        for i in range(line_num):
            start = i * 70
            stop = i * 70 + 70
            file.write('{}\n'.format(sequence[start:stop]))

        start = line_num * 70
        file.write('{}\n'.format(sequence[start:]))

################################################################################
#
#process homologous DNA/protein sequence file
#
################################################################################
in_file = open('{}/{}'.format(dir, sys.argv[2]), 'r')
print('reading input fasta file {}'.format(sys.argv[2]))
accession = ''
sequence = ''
fastaDict = {}
for line in in_file:
    line = line.rstrip()
    if len(line) == 0:
        continue

    if '>' in line:
        #if the first accession, initiate
        if accession == '':
            accession = re.search('>([A-Za-z0-9_]+.[0-9]+)', line)
            accession = accession.group(0)[1:]
            continue

        #store accession-sequence pair in fastaDict
        if seqType == 'protein':
            sequence = replaceAA(sequence)
        fastaDict[accession] = sequence
        sequence = ''
        accession = re.search('>([A-Za-z0-9_]+.[0-9]+)', line)
        accession = accession.group(0)[1:]

    else:
        sequence = sequence + line

if seqType == 'protein':
    sequence = replaceAA(sequence)
fastaDict[accession] = sequence
in_file.close()

################################################################################
#
#         convert the homologous sequence file into fasta format
#
################################################################################
file = open('{}/{}.fasta'.format(dir, accession_full), 'w')
for key in fastaDict.keys():
    write_fasta(file, key, fastaDict[key])
file.close()
fastaDict.clear()
################################################################################
#
#                      run cluster omega on the fasta file
#
################################################################################
os.system('clustalo -i {}/{}.fasta -o {}/{}_aligned.fasta \
    --auto -v --force'.format(dir, accession_full, dir, accession_full))

################################################################################
#
#  read the aligned fasta file, and clean the identifier
#
################################################################################
file = open('{}/{}_aligned.fasta'.format(dir, accession_full), 'r')
seqDict = {}

accession = ''
seq = ''
for line in file:

    line = line.rstrip()
    if len(line) == 0:
        continue

    if '>' in line:
        if len(accession) > 0:
            seqDict[accession] = seq
            seq = ''

        accession = line.replace('_', '')
        accession = accession.replace('.', '')
        accession = accession.replace('>', '')

    else:
        seq = seq + line

seqDict[accession] = seq
# nchar = len(seq)
file.close()

################################################################################
#
#                      convert fasta to nexus file
#
################################################################################
def writeNexus(fileName, seqDict, sequencetype):
    nexus_file = open(fileName, 'w')
    print('generating nexus file {}'.format(nexus_file))
    nexus_file.write('#NEXUS\n\n')
    nexus_file.write('BEGIN TAXA;\n')
    nexus_file.write('\tTITLE Taxa;\n')
    nexus_file.write('\tDIMENSIONS NTAX={};\n'.format(len(seqDict)))

    taxa = ' '.join(sorted(seqDict.keys()))
    dimension = 0
    for key in seqDict.keys():
        dimension = len(seqDict[key])
        break

    nexus_file.write('\tTAXLABELS\n')
    nexus_file.write('\t\t{}\n\t;\n\n'.format(taxa))
    nexus_file.write('END;\n\n')

    nexus_file.write('BEGIN CHARACTERS;\n')
    nexus_file.write('\tTITLE Character_Matrix;\n')
    nexus_file.write('\tDIMENSIONS NCHAR={};\n'.format(dimension))
    nexus_file.write('\tFORMAT DATATYPE = {} GAP = - MISSING = ?;\n'.format(sequencetype))
    nexus_file.write('\tMATRIX\n')

    for key in sorted(seqDict):
        nexus_file.write('\t{} {}\n'.format(key, seqDict[key]))

    nexus_file.write('\n;\nEND;')
    nexus_file.close()

nexusFile = '{}/{}.nex'.format(dir, accession_full)
writeNexus(nexusFile, seqDict, seqType)

################################################################################
#
#                                   run paup
#
################################################################################
def getTreeCMD(nexusFileName, outputFile, dir):
    #TREE SEARCH METHOD
    #tree search method for simultaneous analysis tree & support tests
    treeSearchMethod = 'hsearch nreps=1000 swap=tbr multrees=no'

    #BOOTSTRAP TREE SEARCH METHOD
    #   specify method for bootstrap tree searhc
    #   bootstrapMethod = "search=bandb"
    bootstrapMethod = 'search=heuristic nreps=50'

    createTreeCmdFile = open(outputFile, 'w')
    createTreeCmdFile.write('#NEXUS\n\n')
    createTreeCmdFile.write('set warnReset = no;\n')
    createTreeCmdFile.write('set increase = auto;\n')
    createTreeCmdFile.write('set datastorage=full;\n')
    createTreeCmdFile.write('set criterion=parsimony;\n')
    createTreeCmdFile.write('execute {};\n'.format(nexusFileName))
    createTreeCmdFile.write('{};\n'.format(treeSearchMethod))
    createTreeCmdFile.write('filter best;\n')

    treeFile = outputFile.replace('_getTree.cmd', '_trees.nex')

    createTreeCmdFile.write('savetrees file={} format=nexus replace=yes root=yes;\n'.format(treeFile))
    createTreeCmdFile.write('quit warnTsave=no;\n')
    createTreeCmdFile.close()

getTreeCMD(nexusFile, '{}/{}_getTree.cmd'.format(dir, accession_full), dir)
os.system('paup {}/{}_getTree.cmd'.format(dir, accession_full))
################################################################################
#
#                              Generate variation number
#
################################################################################
with open('{}/{}_aligned.fasta'.format(dir,accession_full)) as f:
    alist = [line.rstrip() for line in f]
seqDict = {}

accession = ''
seq = ''
for line in alist:
    if '>' in line:
        if accession != '':
            seqDict[accession] = seq
        accession = line.replace('>', '')
        accession = accession.replace('_', '')
        accession = accession.replace('.', '')
        seq = ''
    else:
        seq = seq + line
seqDict[accession] = seq
#print(len(seqDict))

homoAccession = accession_full
homoAccession = homoAccession.replace('_', '')
homoAccession = homoAccession.replace('.', '')
homoSeq = seqDict[homoAccession]
seqLength = len(homoSeq)

################################################################################
#
#                       convert phylogenetic tree
#
################################################################################
Phylo.convert('{}/{}_trees.nex'.format(dir, accession_full), 'nexus', '{}/{}_tree.tree'.format(dir, accession_full), 'newick')

f = open('{}/{}_tree.tree'.format(dir, accession_full), 'r')
tree = f.readline()
#print(tree)

tree = re.sub(r':\d.\d+', '', tree)
#print(tree)
tree = TreeNode.read(StringIO(tree))

################################################################################
#
#     find set of characters given a list of taxa l and position number i
#
################################################################################
def findSet(l, i, seqDict):
    ll = []
    for j in l:
        ll.append(str(seqDict[j][i]))
    ll = list(set(ll))
    return ll

def updateVN(node, child, variation_number, seqDict, length):

    allClades = re.findall(r'[a-zA-Z0-9]+', str(node))
    #print("all")
    #print(allClades)
    for c in child:
        cClades = re.findall(r'[a-zA-Z0-9]+', str(c))
        subClades = list(set(allClades) - set(cClades))
        #print("sub")
        #print(subClades)
        for i in range(length):
            cSet = findSet(cClades, i, seqDict)
            subSet = findSet(subClades, i, seqDict)
            for item in cSet:
                if item not in subSet:
                    variation_number[i] = variation_number[i] + 1
    return variation_number

def generateVN(tree, seqDict, seqLength):
    variation_number = np.zeros((seqLength,), dtype=int)
    queue = []
    queue.append(tree)

    while len(queue) > 0:
        node = queue.pop()
        if node.is_tip():
            continue
        else:
            child = node.children
            variation_number = updateVN(node, child, variation_number, seqDict, seqLength)
            #print('variation')
            #print(variation_number)
            for c in child:
                queue.append(c)
    return variation_number

variation_number = generateVN(tree, seqDict, seqLength)

homoIndexList = []
f_vn = []
for i in range(len(homoSeq)):
    if str(homoSeq[i]) != '-':
        homoIndexList.append(i)
        f_vn.append(variation_number[i])
#print(len(homoIndexList))

outputFile = open("{}/vn_{}.txt".format(dir,accession_full), 'w')

vn_max = max(f_vn)
vn_min = min(f_vn)
#print(vn_min, vn_max)
#exit()
for i in range(len(homoIndexList)):
    j = i + 1
    vn = variation_number[homoIndexList[i]]
    vn = (vn - vn_min) / (vn_max - vn_min)
    outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))

#outputFile.write("\n")
outputFile.close()
