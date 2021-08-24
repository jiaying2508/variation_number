# Variation Number
Characteristic Attribute Organization System (CAOS) discovers rules associated with a given phylogenetic tree. A pure (Pu) rule or character attribute (CA) is a state that exists in all elements of a clade but not the alternate clade; a private (Pr) CA is present in some members of a clade but absent in the alternate clade. A variation number (VN) is defined as the number of occurrences of a position as a CA in all the tree clades.

## Required external packages
In order to run **vn.py**, please install **command line** version for:
1. Clustal Omega: http://www.clustal.org/omega/
2. PAUP: http://phylosolutions.com/paup-test/

## Running Instructions
Download **vn.py**, and run use python version 3.7.4 or higher
```console
$ python vn.py -g <geneName> -s <sequenceType>
```
Example:
```console
$ python -g NAT2 -s nucleotide
```
Please note that the sequence type need to be either **nucleotide** or **protein**

## Optional Parameters
```
-o
-a
```
