# Variation Number

A package for calculating the variation number of nucleotide/protein sequence using sequence orthologs.

Characteristic Attribute Organization System (CAOS) discovers rules associated with a given phylogenetic tree. A pure (Pu) rule or character attribute (CA) is a state that exists in all elements of a clade but not the alternate clade; a private (Pr) CA is present in some members of a clade but absent in the alternate clade. A variation number (VN) is defined as the number of occurrences of a position as a CA in all the tree clades.

The method is described in the publication:  
Lai, J., & Sarkar, I. N. (2021). A Phylogenetic Approach to Analyze the Conservativeness of BRCA1 and BRCA2 Mutations. AMIA ... Annual Symposium proceedings. AMIA Symposium, 2020, 677–686. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8075528/

## Features

- Download orthologs
- Build phylogenetic trees
- Generate variation numbers

## Required python packages
Python packages (most of which can be installed using pip) needed to run LYRUS include:
- skbio: http://scikit-bio.org
- numpy: https://numpy.org/install/
- Bio: https://biopython.org/wiki/Download
- BeautifulSoup: https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup

## Required external packages
In order to run **vn.py**, please install **command line** version for:
1. Clustal Omega: http://www.clustal.org/omega/
2. Mafft: https://mafft.cbrc.jp/alignment/software/
3. PAUP: http://phylosolutions.com/paup-test/

## Running instructions for installation using pip
variation_number is published on [PyPI](https://pypi.org/). Use the following command to install variation_number using pip:
```console
$ pip install variation-number
```

### Usage
```
import variation_number as vn
import os
gene = 'BRCA1'
seqtype =' protein'
outputDir = '{}/output'.format(os.getcwd())

# Download orthologs from NCBI orthologs database (optional; can use user provided sequence file)
acc = vn.getFasta(gene, outputDir,seqtype,refseqID=None)

# Calculate variation number using clustal omega
vn.processVN(file=gene, outputDir, accession_full=acc, seqType=seqtype, aligned=False, alignTool='clustal')
```