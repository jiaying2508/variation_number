"""
A module for calculating variation number
"""

__author__ = "Jiaying Lai"
__date__ = "August, 2021"
__maintainer__ = "Jiaying Lai"
__email__ = "jiaying_lai@brown.edu"

import argparse
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
    """

    Parameters
    ----------
    sequence : input sequence
        
    seqType : nucleotide or protein
        

    Returns
    -------
    if the input sequence is protein sequence, replace letter not in 'ARNDCQEGHILKMFPSTWYV' with 'X'

    """
    if seqType == 'nucleotide':
        return sequence
    seq = ''
    for char in sequence:
        if char in 'ARNDCQEGHILKMFPSTWYV-':
            seq = seq + char
        else:
            seq = seq + 'X'
    return seq

'''
write a sequence to fasta file
'''
def write_fasta(file, accession, sequence):
    """

    Parameters
    ----------
    file : output file 
        
    accession : gene/protein accession number
        
    sequence : output sequence
        

    Returns
    -------

    """
    file.write('>{}\n'.format(accession))
    if len(sequence) <= 70:
        file.write('{}\n'.format(sequence))
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
    """

    Parameters
    ----------
    fileName : output file
        
    seqDict : dictionary where key is accession and value is sequence
        
    sequencetype : protein or nucleotide
        

    Returns
    -------

    """
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

def getTreeCMD(nexusFileName, outputFile):
    """

    Parameters
    ----------
    nexusFileName : nexus file name
        
    outputFile : output file name
        

    Returns
    -------

    """
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
    """

    Parameters
    ----------
    node : tree node
        
    child : all the child clade of the current node
        
    variation_number : variation number
        
    seqDict : sequence dictionary
        
    length : sequence length
        

    Returns
    -------
    variation_number(np.array)
    """

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
    """

    Parameters
    ----------
    tree : input tree
        
    seqDict : dictionary where key is accession and value is sequence
        
    seqLength : sequence length
        

    Returns
    -------
    calculated variation number
    """
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

def processVN(file, outputDir, accession_full=None, seqType='protein', aligned=False, alignTool='clustal'):
    """

    Parameters
    ----------
    file : input file in fasta format (full path)
        
    outputDir : output directory
        
    accession_full : accession number; if provided, will output variation number for this sequence without gaps
        
    seqType : nucleotide or protein

    aligned: whether the input file is aligned    

    alignTool: the alignment program to use. clustal: clustalo; mafft: mafft

    Returns
    -------
    vn.txt
    vn_scaled.txt: scaled using min max normalization
    vn_{$accession_full}.txt
    vn_{$accession_full}_scaled.txt: scaled using min max normalization

    """

    if accession_full == None:
        accession_full = '_accession'
    ################################################################################
    #process homologous protein sequence file
    ################################################################################
    in_file = open(file)

    accession = ''
    accession1 = ''
    sequence = ''
    fastaDict = {}
    accessionCount = 1
    accessionDict = {}
    homoAccession = ''
    for line in in_file:
        line = line.rstrip()
        if len(line) == 0:
            continue

        if '>' in line:
            #if the first accession, initiate
            if accession != '':
                #store accession-sequence pair in fastaDict
                sequence = replaceAA(sequence, seqType)
                fastaDict[accession1] = sequence
                sequence = ''

            accession = line[1:]
            accession1 = 's' + str(accessionCount)
            accessionDict[accession1] = accession

            if accession == accession_full:
                homoAccession = accession1
            accessionCount += 1

        else:
            sequence = sequence + line

    sequence = replaceAA(sequence, seqType)
    fastaDict[accession1] = sequence
    in_file.close()

    outputFile = open('{}/accession.txt'.format(outputDir), 'w')
    for key in accessionDict.keys():
        outputFile.write('{}\t{}\n'.format(key, accessionDict[key]))
    outputFile.close()

    ################################################################################
    #         convert the homologous sequence file into fasta format
    ################################################################################
    if aligned:
        file1 = open('{}/_temp.fasta'.format(outputDir), 'w')
        for key in fastaDict.keys():
            write_fasta(file1, key, fastaDict[key])
        file1.close()

        os.system('mv {}/_temp.fasta {}/{}_aligned.fasta'.format(outputDir, outputDir, accession_full))
        fastaDict.clear()
    else:
        file1 = open('{}/{}.fasta'.format(outputDir, accession_full), 'w')
        for key in fastaDict.keys():
            write_fasta(file1, key, fastaDict[key])
        file1.close()
        fastaDict.clear()

        ################################################################################
        #                      run cluster omega on the fasta file
        ################################################################################
        # print('# Performing multiple sequence alignment\n')
        # print('clustalo -i {}/{}.fasta -o {}/{}_aligned.fasta \
        #     --auto -v --force >/dev/null'.format(outputDir, accession_full, outputDir, accession_full))
        
        if alignTool=='clustal':
            print('# Performing multiple sequence alignment using CLUSTALO\n')
            os.system('clustalo -i {}/{}.fasta -o {}/{}_aligned.fasta \
                --auto -v --force >/dev/null'.format(outputDir, accession_full, outputDir, accession_full))
        else:
            print('# Performing multiple sequence alignment using MAFFT\n')
            os.system('mafft --auto --quiet {}/{}.fasta > {}/{}_aligned.fasta'.format(outputDir, accession_full, outputDir, accession_full))
    
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
            accession = line[1:]
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
    getTreeCMD(nexusFile, '{}/{}_getTree.cmd'.format(outputDir, accession_full))
    os.system('paup {}/{}_getTree.cmd >/dev/null'.format(outputDir, accession_full))

    ################################################################################
    #                              Generate variation number
    ################################################################################
    print('# Generating Variation Number for {}'.format(file))

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
    seqLength = len(seq)

    if accession_full == '_accession':
        homoAccession = None
    else:
        homoSeq = seqDict[homoAccession]
        

    ################################################################################
    #                       convert phylogenetic tree
    ################################################################################
    Phylo.convert('{}/{}_trees.nex'.format(outputDir, accession_full), 'nexus', '{}/{}_tree.tree'.format(outputDir, accession_full), 'newick')

    f = open('{}/{}_tree.tree'.format(outputDir, accession_full), 'r')
    tree = f.readline()

    tree = re.sub(r':\d.\d+', '', tree)
    tree = TreeNode.read(StringIO(tree))

    variation_number = generateVN(tree, seqDict, seqLength)

    outputFile = open("{}/vn.txt".format(outputDir), 'w')
    for i in range(len(variation_number)):
        j = i + 1
        vn = variation_number[i]
        outputFile.write('{}\t{}\n'.format(str(j), str(vn)))
    outputFile.close()

    outputFile = open("{}/vn_scaled.txt".format(outputDir), 'w')
    vn_max = max(variation_number)
    vn_min = min(variation_number)
    for i in range(len(variation_number)):
        j = i + 1
        vn = variation_number[i]
        vn = (vn - vn_min) / (vn_max - vn_min)
        outputFile.write('{}\t{}\n'.format(str(j), str(vn)))
    outputFile.close()

    if homoAccession != None:
        homoIndexList = []
        f_vn = []
        for i in range(len(homoSeq)):
            if str(homoSeq[i]) != '-':
                homoIndexList.append(i)
                f_vn.append(variation_number[i])

        outputFile = open("{}/vn_{}.txt".format(outputDir,accession_full), 'w')
        for i in range(len(homoIndexList)):
            j = i + 1
            vn = variation_number[homoIndexList[i]]
            outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))
        outputFile.close()
    
        outputFile = open("{}/vn_{}_scaled.txt".format(outputDir,accession_full), 'w')
        vn_max = max(f_vn)
        vn_min = min(f_vn)
        for i in range(len(homoIndexList)):
            j = i + 1
            # vn = variation_number[homoIndexList[i]]
            vn = f_vn[i]
            vn = (vn - vn_min) / (vn_max - vn_min)
            outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))
        outputFile.close()

def getFasta(geneName, outputDir, seqType, refseqID=None, email=''):
    """

    Parameters
    ----------
    geneName : gene name
        
    outputDir : output directory
        
    seqType : protein or nucleotide
        
    refseqID : accession number
        

    Returns
    -------
    accession number for homo sapiens
    """
    Entrez.email = email
    print('\n# process {} refseq'.format(geneName))
    page = requests.get('https://www.ncbi.nlm.nih.gov/protein/?term={}'.format(geneName))
    geneID = ''
    homo_acc = ''
    taxidDict = {}
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
                                    taxidDict[homo_acc] = 9606
                                    continue
                                else:
                                    for item in gene['refseq_accessions']:
                                        try:
                                            if seqType == 'nucleotide':
                                                homo_acc = item['transcript_acc']
                                                refseqList.append(homo_acc)
                                                taxidDict[homo_acc] = 9606
                                                break
                                            else:
                                                homo_acc = item['protein_acc']
                                                refseqList.append(homo_acc)
                                                taxidDict[homo_acc] = 9606
                                                break
                                        except:
                                            continue
                            else:
                                try:
                                    seqAcc = gene['tax_id']
                                    if seqType == 'nucleotide':
                                        taxidDict[gene['refseq_accessions'][0]['transcript_acc']] = seqAcc
                                        refseqList.append(gene['refseq_accessions'][0]['transcript_acc'])
                                    else:
                                        taxidDict[gene['refseq_accessions'][0]['protein_acc']] = seqAcc
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

            outputFile = open('{}/taxid.txt'.format(outputDir), 'w')
            for key in taxidDict.keys():
                outputFile.write('{}\t{}\n'.format(key, taxidDict[key]))
            outputFile.close()
            return homo_acc
    return ''
