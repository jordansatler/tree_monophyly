#!/usr/bin/env python

"""
Get proportion of gene trees for which a given species is monophyletic.
Requires at least two individuals sampled per species for a gene
tree to be considered, and tests whether the species is monophyletic. 
Outputs proportion of trees for which they are monophyletic, as well as
total number of trees for which they were present with minimum sampling.
Also reports proportion of species not monophyletic for each locus.

author: J. Satler
date: 21 Nov 2019
version: 1

usage: python monophylyTest.py traits.file /path/to/trees/

"""
import sys
import glob
import dendropy

def traits(infile):
    """read traits file and assign individuals to OTUs"""
    with open(infile, 'r') as t:
        return {i.split()[0]:i.split()[1] for i in t}

def tree_files(files):
    """get list of tree files"""
    return glob.glob(files + "/*.tr*")

def get_tree(infile):
    """get trees and return as list"""
    return dendropy.TreeList.get(path=infile,
                                 schema="newick",
                                 rooting="force-unrooted",
                                 preserve_underscores=True)

def tree_sampling(traits, tree):
    """map individuals to species in tree"""
    species = {sp:[] for sp in set(traits.values())}
    for i in tree.taxon_namespace:
        species[traits[i.label]].append(i.label)
    return species

def test_for_monophyly(sp_map, trees):
    """test if species is monophyletic"""
    res = {}
    for sp, inds in sp_map.items():
        # make sure there are at least two individuals for a species
        if len(inds) > 1:
            res[sp] = trees.frequency_of_bipartition(labels=inds)
    return res, len(res)

def get_non_mono_taxa(mono):
    """get taxa that are not monophyletic for locus"""
    return [k for k, v in mono.items() if v == 0.0]

def proportion_monophyletic(res):
    """get proportion of monophyletic trees"""
    with open("results_taxa_monophyly.txt", "w") as out:
        header = "species\tproportion monophyletic\tN gene trees"
        out.write(header + "\n")
        for sp, mono in sorted(res.items()):
            m = [y for y in mono if y == 1.0]
            # skip if a species was not sampled in any tree
            if len(mono) == 0:
                continue
            out.write("{0}\t{1:0.2f}\t{2}\n".format(sp,
                                               float(len(m))/len(mono),
                                               len(mono)))

def monophyly_of_loci(non_mono):
    """output numbers on loci and monophyly"""
    with open("results_loci_nonmonphyly.txt", "w") as out:
        header = "locus\tproportion Not monophyletic\tspecies"
        out.write(header + "\n")
        for k, v in sorted(non_mono.items(),
                           key = lambda x: x[1][1],
                           reverse = True):
            out.write("{0}\t{1:0.2f}\t{2}\n".format(k,
                                                    v[1],
                                                    ', '.join(sorted(v[0]))))

def main():
    if len(sys.argv) != 3:
        print("python monophylyTest.py traits.file /path/to/trees/")
        sys.exit()

    tr = traits(sys.argv[1])
    treef = tree_files(sys.argv[2])

    results = {sp:[] for sp in set(tr.values())}
    non_monophyly = {}

    for t in treef:
        locus = t.split("/")[-1].split(".")[0]
        tr_in = get_tree(t)
        sp_map = tree_sampling(tr, tr_in)

        # test if species are monophyletic
        r, taxa = test_for_monophyly(sp_map, tr_in)
        non_mono = get_non_mono_taxa(r)
        # record non monophyletic taxa
        if non_mono:
            non_monophyly[locus] = (non_mono, len(non_mono) / float(taxa))

        for k, v in r.items():
            results[k].append(v)

    # write results to file
    proportion_monophyletic(results)
    monophyly_of_loci(non_monophyly)

if __name__ == '__main__':
    main()
