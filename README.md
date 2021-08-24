# Variation Number
Characteristic Attribute Organization System (CAOS) discovers rules associated with a given phylogenetic tree. A pure (Pu) rule or character attribute (CA) is a state that exists in all elements of a clade but not the alternate clade; a private (Pr) CA is present in some members of a clade but absent in the alternate clade. A variation number (VN) is defined as the number of occurrences of a position as a CA in all the tree clades.

The method is described in the publication:  
Lai, J., & Sarkar, I. N. (2021). A Phylogenetic Approach to Analyze the Conservativeness of BRCA1 and BRCA2 Mutations. AMIA ... Annual Symposium proceedings. AMIA Symposium, 2020, 677–686. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8075528/

## Required python packages
Python packages (most of which can be installed using pip) needed to run LYRUS include:
- skbio: http://scikit-bio.org
- numpy: https://numpy.org/install/
- Bio: https://biopython.org/wiki/Download
- BeautifulSoup: https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup

## Required external packages
In order to run **vn.py**, please install **command line** version for:
1. Clustal Omega: http://www.clustal.org/omega/
2. PAUP: http://phylosolutions.com/paup-test/

## Running Instructions
Download **vn.py**, and run use python version 3.7.4 or higher
```console
$ python vn.py -g <geneName> -s <sequenceType: nucleotide or protein>
```
Example:
```console
$ python -g BRCA1 -s nucleotide
```

## Optional Parameters
```
-o output directory: need to be full path
-a RefSeq accession
```
