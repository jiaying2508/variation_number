#run XGBoost Classifier to predict SAV disease pathogenicity
'''
python vn.py gene outputDir sequenceType refseqID 

python version:
    python3.7.4 or higher

args:
    [1] gene
    [2] outputDir fullPath
    [3] sequence type: protein or nucleotide
    [4] refseqID

required python packages:
    Bio
    skbio
    BeautifulSoup
    numpy

required external packages:
    clustal omega
    paup
    
'''

__author__ = "Jiaying Lai"
__date__ = "August, 2021"
__maintainer__ = "Jiaying Lai"
__email__ = "jiaying_lai@brown.edu"

import argparse
import sys
import re
import math
import os
from skbio import TreeNode
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
import numpy as np
from Bio import Entrez
import json
from bs4 import BeautifulSoup
import requests
import time

'''
replace the amino acid with X if it is not among the 20 amino acid codes
'''
def replaceAA(sequence, seqType):
    if seqType == 'nucleotide':
        return sequence
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

def writeNexus(fileName, seqDict, sequencetype):
    nexus_file = open(fileName, 'w')
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

'''
functions to calculate variation number
'''
def findSet(l, i, seqDict):
    ll = []
    for j in l:
        ll.append(str(seqDict[j][i]))
    ll = list(set(ll))
    return ll

def updateVN(node, child, variation_number, seqDict, length):

    allClades = re.findall(r'[a-zA-Z0-9]+', str(node))
    for c in child:
        cClades = re.findall(r'[a-zA-Z0-9]+', str(c))
        subClades = list(set(allClades) - set(cClades))
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

def processVN(file, outputDir, accession_full, seqType):

    ################################################################################
    #process homologous protein sequence file
    ################################################################################
    try:
        in_file = open('{}/{}'.format(outputDir, file, 'r'))

        print('# Generating Variation Number for {}'.format(file))

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
                sequence = replaceAA(sequence, seqType)
                fastaDict[accession] = sequence
                sequence = ''
                accession = re.search('>([A-Za-z0-9_]+.[0-9]+)', line)
                accession = accession.group(0)[1:]

            else:
                sequence = sequence + line

        sequence = replaceAA(sequence, seqType)
        fastaDict[accession] = sequence
        in_file.close()

        ################################################################################
        #         convert the homologous sequence file into fasta format
        ################################################################################
        file1 = open('{}/{}.fasta'.format(outputDir, accession_full), 'w')
        for key in fastaDict.keys():
            write_fasta(file1, key, fastaDict[key])
        file1.close()
        fastaDict.clear()

        ################################################################################
        #                      run cluster omega on the fasta file
        ################################################################################
        os.system('clustalo -i {}/{}.fasta -o {}/{}_aligned.fasta \
            --auto -v --force >/dev/null'.format(outputDir, accession_full, outputDir, accession_full))
            
        ################################################################################
        #  read the aligned fasta file, and clean the identifier
        ################################################################################
        file1 = open('{}/{}_aligned.fasta'.format(outputDir, accession_full), 'r')
        seqDict = {}
        accession = ''
        seq = ''
        for line in file1:

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
        file1.close()

        ################################################################################
        #                      convert fasta to nexus file
        ################################################################################
        nexusFile = '{}/{}.nex'.format(outputDir, accession_full)
        writeNexus(nexusFile, seqDict, seqType)

        ################################################################################
        #                                   run paup
        ################################################################################
        getTreeCMD(nexusFile, '{}/{}_getTree.cmd'.format(outputDir, accession_full), outputDir)
        os.system('paup {}/{}_getTree.cmd >/dev/null'.format(outputDir, accession_full))

        ################################################################################
        #                              Generate variation number
        ################################################################################
        with open('{}/{}_aligned.fasta'.format(outputDir,accession_full)) as f:
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

        homoAccession = accession_full
        homoAccession = homoAccession.replace('_', '')
        homoAccession = homoAccession.replace('.', '')
        homoSeq = seqDict[homoAccession]
        seqLength = len(homoSeq)

        ################################################################################
        #                       convert phylogenetic tree
        ################################################################################
        Phylo.convert('{}/{}_trees.nex'.format(outputDir, accession_full), 'nexus', '{}/{}_tree.tree'.format(outputDir, accession_full), 'newick')

        f = open('{}/{}_tree.tree'.format(outputDir, accession_full), 'r')
        tree = f.readline()

        tree = re.sub(r':\d.\d+', '', tree)
        tree = TreeNode.read(StringIO(tree))

        variation_number = generateVN(tree, seqDict, seqLength)

        homoIndexList = []
        f_vn = []
        for i in range(len(homoSeq)):
            if str(homoSeq[i]) != '-':
                homoIndexList.append(i)
                f_vn.append(variation_number[i])

        outputFile = open("{}/vn_{}.txt".format(outputDir,accession_full), 'w')

        vn_max = max(f_vn)
        vn_min = min(f_vn)

        for i in range(len(homoIndexList)):
            j = i + 1
            vn = variation_number[homoIndexList[i]]
            vn = (vn - vn_min) / (vn_max - vn_min)
            outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))

        outputFile.close()
    except:
        pass

def getFasta(geneName, outputDir, seqType, refseqID):
    print('\n# process {} refseq'.format(geneName))
    page = requests.get('https://www.ncbi.nlm.nih.gov/protein/?term={}'.format(geneName))
    geneID = ''
    homo_acc = ''
    if page.status_code == 200:
        soup = BeautifulSoup(page.text, 'html.parser')
        soup.prettify()

        for line in soup:
            geneID1 = re.search('<li>\s?Gene\s?ID\s?:\s?(\d+)\s?<\/li>', str(line))
            if geneID1:
                geneID = geneID1.group(1)
                break
        time.sleep(0.34)
        page = requests.get('https://www.ncbi.nlm.nih.gov/gene/{}/ortholog/?term={}'.format(geneID,geneName))
        if page.status_code == 200:
            refseqList = []
            soup = BeautifulSoup(page.text, 'html.parser')
            soup.prettify()
            # c = 0
            for line in soup:
                if 'appData' in str(line):
                    match = re.search('appData\.genes\s?=\s?([^;]+);', str(line))
                    if match:
                        genes = match.group(1)
                        try:
                            res = json.loads(genes)
                        except:
                            match = re.search('appData\.genes\s?=\s?(.+}]);', str(line))
                            genes = match.group(1)
                            res = json.loads(genes)
                        count = 0
                        for gene in res:
                            count += 1
                            if gene['tax_id'] == 9606:
                                if refseqID != None:
                                    homo_acc = refseqID
                                    refseqList.append(refseqID)
                                    continue
                                else:
                                    for item in gene['refseq_accessions']:
                                        try:
                                            if seqType == 'nucleotide':
                                                homo_acc = item['transcript_acc']
                                                refseqList.append(homo_acc)
                                                break
                                            else:
                                                homo_acc = item['protein_acc']
                                                refseqList.append(homo_acc)
                                                break
                                        except:
                                            continue
                            else:
                                try:
                                    if seqType == 'nucleotide':
                                        refseqList.append(gene['refseq_accessions'][0]['transcript_acc'])
                                    else:
                                        refseqList.append(gene['refseq_accessions'][0]['protein_acc'])
                                except:
                                    pass

            outputFile = open('{}/{}'.format(outputDir, geneName), 'w')
            print('# downloading orthologs files')
            for gene in refseqList:
                time.sleep(0.34)
                if seqType == 'nucleotide':
                    handle = Entrez.efetch(db="nuccore", id=gene, rettype = 'fasta', retmode='text')
                else:
                    handle = Entrez.efetch(db="protein", id=gene, rettype = 'fasta', retmode='text')
                txt = handle.read()
                handle.close()
                seq = ''
                lCount = 0
                for l in txt.split('\n'):
                    if lCount == 0:
                        lCount = 1
                        continue
                    seq = seq + l
                outputFile.write('>{}\n{}\n'.format(gene, seq))
            outputFile.close()
            return homo_acc
    return ''

parser = argparse.ArgumentParser(description='Calculation for the Variation Number')
parser.add_argument('-g', '--gene', type=str, metavar='', required=True, help='Gene Name')
parser.add_argument('-o', '--outputDir', type=str, metavar='', help='Output Directory')
parser.add_argument('-s', '--seqType', type=str, metavar='', required=True, help='sequence type, protein or nucleotide')
parser.add_argument('-a', '--acc', type=str, metavar='', help='refseq accession')

args = parser.parse_args()

if __name__ == '__main__':
    if args.outputDir == None:
        args.outputDir = os.getcwd()
    else:
        try:
            os.mkdir(args.outputDir)
        except:
            pass
    seqType = args.seqType.lower()
    if seqType != 'protein' and seqType != 'nucleotide':
        print('Please specify is the sequence type is protein or nucleotide')
        exit()
    ################################################################################
    '''
    download orthologs file
    '''
    ################################################################################
    Entrez.email = ''

    acc = getFasta(args.gene, args.outputDir, seqType, args.acc)
    if len(acc) == 0:
        print('Orthologs information for {} is not available'.format(args.gene))
    else:
        processVN(args.gene, args.outputDir, acc, seqType)