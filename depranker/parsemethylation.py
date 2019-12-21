#!/usr/bin/env python

import numpy as np
import os
import sys

__author__ = "Sanket S Desai"
__license__ = "MIT"
__version__ = "0.2"
__email__ = "desai.sanket12@gmail.com"

class MethylationParser(object):
    def __init__(self, efn):
        if not os.path.exists(efn) and not os.path.isfile(efn):
            print("Methylation beta value file %s not found or unreadable, please check!!" %(enf))
            sys.exit(0)
        efi=open(efn)
        self.header_=efi.readline().strip().replace('\"','').split("\t")
        self.gene_cnv_map_={}
        for line in efi:
            sline=line.strip().split("\t")
            self.gene_cnv_map_[sline[0]]=sline[1]
        efi.close()

    def get_methylation_beta_values(self, gene):
        try:
            return self.gene_cnv_map_[gene]
        except ValueError :
            return ""

    def get_number_of_genes(self):
        return len(self.gene_cnv_map_)

    def get_gene_betavalue_map(self):
        return self.gene_cnv_map_

    def get_gene_betavalue_map(self, genelist): #given a gene list extract the values
        gem={}
        for g in genelist:
            try:
                gem[g]=self.gene_cnv_map_[g]
            except KeyError:
                #print("Gene name \'%s\' not found in the CNV dataset!! Please check!!" %(g))
                pass
        return gem

    #def get_gene_cnv_map(self):
    #    return self.gene_cnv_map_

    def get_gene_methylation_score_map(self):
        gemap=self.get_gene_betavalue_map()
        gexpnp=np.array(gemap.values()).astype(float)
        gexpnp_std= (gexpnp - gexpnp.min(axis=0)) / (gexpnp.max(axis=0) - gexpnp.min(axis=0))
        gexpnp_scaled = gexpnp_std * (10 - 0) + 0
        #gexpnp_scaled=10-avglfcnp_scaled #inverse is true
        grs_map={}
        ind=0
        for g in gemap.keys():
            grs_map[g]=float(gexpnp_scaled[ind])
            ind=ind+1
        return grs_map
'''
def main():
    cnvp=CopyNumberVariationParser(sys.argv[1])
    print(cnvp.get_gene_cnv_score_map())

if __name__=="__main__":
    main()
'''
