#!/usr/bin/env python

import numpy as np
import os
import sys

__author__ = "Sanket S Desai"
__license__ = "MIT"
__version__ = "0.2"
__email__ = "desai.sanket12@gmail.com"

class DemeterParser(object):
    def __init__(self, efn):
        if not os.path.exists(efn) and not os.path.isfile(efn):
            print("DEMETER score file %s not found or unreadable, please check!!" %(enf))
            sys.exit(0)
        efi=open(efn)
        self.header_=efi.readline().strip().replace('\"','').split("\t")
        self.gene_demeter_map_={}
        for line in efi:
            sline=line.strip().split("\t")
            self.gene_demeter_map_[sline[0]]=sline[1]
        efi.close()

    def get_demeter_score(self, gene):
        try:
            return self.gene_demeter_map_[gene]
        except ValueError :
            return ""

    def get_number_of_genes(self):
        return len(self.gene_demeter_map_)

    def get_gene_demeterscore_map(self):
        return self.gene_demeter_map_

    def get_gene_demeterscore_map(self, genelist): #given a gene list extract the values
        gem={}
        for g in genelist:
            try:
                gem[g]=self.gene_demeter_map_[g]
            except KeyError:
                #print("Gene name \'%s\' not found in the expression dataset!! Please check!!" %(g))
                pass
        return gem

    def get_gene_demeter_score_map(self):
        gemap=self.get_gene_demeterscore_map()
        gexpnp=np.array(list(gemap.values())).astype(float)
        gexpnp_std= (gexpnp - gexpnp.min(axis=0)) / (gexpnp.max(axis=0) - gexpnp.min(axis=0))
        gexpnp_scaled = gexpnp_std * (10 - 0) + 0
        gexpnp_scaled=10-avglfcnp_scaled #inverse is true
        grs_map={}
        ind=0
        for g in gemap.keys():
            grs_map[g]=float(gexpnp_scaled[ind])
            ind=ind+1
        return grs_map
