# tree_monophyly

This script will calculate the proportion of gene trees for which a species is monophyletic. For a gene tree to be considered for a given species, two or more individuals from that species must be sampled. Two output files are generated: 1) **results_taxa_monophyly.txt** contains all sampled species, the proportion of gene trees for which the species were monophyletic, and the number of tested gene trees that met sampling criteria of two or more individuals, and 2) **results_loci_nonmonphyly.txt** contains each locus, the proportion of sampled species that were *not* monophyletic, and the species that were *not* monophyletic.  

You need to specify a traits file (assigning individuals to species; see below) and provide the path to your folder of individual tree files. Tree files are expected to have a file ending starting with ".tr*".

usage:  
```python
    python monophylyTest.py traits.file /path/to/trees/
```

***
The traits file is tab-delimited and assigns individuals to species:

ind1&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind2&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind3&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind4&nbsp;&nbsp;&nbsp;&nbsp;sp2  
ind5&nbsp;&nbsp;&nbsp;&nbsp;sp2  
ind6&nbsp;&nbsp;&nbsp;&nbsp;sp3  
ind7&nbsp;&nbsp;&nbsp;&nbsp;sp3  
ind8&nbsp;&nbsp;&nbsp;&nbsp;sp3  
...
***
References:

Sukumaran, J. and Mark T. Holder. 2010. DendroPy: A Python library for phylogenetic computing. Bioinformatics 26: 1569-1571.
