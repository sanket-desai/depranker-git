#!/usr/bin/env python

import numpy as np
import os
import sys

__author__ = "Sanket S Desai"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "desai.sanket12@outlook.com"
'''
RoastParser - class to parse the ROAST (DOI: 10.1093/bioinformatics/btq401) output file (a bioconductor package)
Format description: Example record below
"NGenes"	"PropDown"	"PropUp"	"Direction"	"PValue"	"FDR"	"PValue.Mixed"	"FDR.Mixed"
"DSTYK"	5	0	0	"Down"	0.118	0.9995	0.903	0.9995
'''
class RoastParser(object):
    def __init__(self, rfn):
        if not os.path.exists(rfn) and not os.path.isfile(rnf):
            print("Roast file %s not found or unreadable, please check!!" %(rnf))
            sys.exit(0)
        rfi=open(rfn)
        records=[]
        self.header_=rfi.readline().strip().replace('\"','').split("\t")
        for line in rfi:
            sline=line.strip().replace('\"','').split("\t")
            records.append(sline)
        self.rnpmatrix_=np.array(records)
        #print(self.rnpmatrix_)
        rfi.close()

    def get_number_of_genes(self):
        return len(self.rnpmatrix_)

    def get_roast_genes(self):
        return list(self.rnpmatrix_[range(len(self.rnpmatrix_)),0])

    def get_ngenes(self):
        return list(self.get_column_as_nparray(1))

    def get_number_of_records(self):
        return self.get_number_of_genes()

    def get_header(self):
        return self.header_

    def get_column_as_nparray(self, i):
        return self.rnpmatrix_[range(len(self.rnpmatrix_)),i]

    def get_gene_as_nparray(self):
        return self.get_column_as_nparray(0)

    def get_ngenes_as_nparray(self):
        return self.get_column_as_nparray(1)

    def get_propdown_as_nparray(self):
        return self.get_column_as_nparray(2)

    def get_propup_as_nparray(self):
        return self.get_column_as_nparray(3)

    def get_direction_as_nparray(self):
        return self.get_column_as_nparray(4)

    def get_pvalue_as_nparray(self):
        return self.get_column_as_nparray(5)

    def get_fdr_as_nparray(self):
        return self.get_column_as_nparray(6)

    def get_pvaluemixed_as_nparray(self):
        return self.get_column_as_nparray(7)

    def get_fdrmixed_as_nparray(self):
        return self.get_column_as_nparray(8)

    def get_up_direction_records_as_npmatrix(self):
        return self.rnpmatrix_[self.get_direction_as_nparray()=="Up",]

    def get_down_direction_records_as_npmatrix(self):
        return self.rnpmatrix_[self.get_direction_as_nparray()=="Down",]

    def get_up_direction_records_as_npmatrix(self, npmat):
        return npmat[npmat[range(len(npmat)),4]=="Up",]

    def get_down_direction_records_as_npmatrix(self,npmat):
        return npmat[npmat[range(len(npmat)),4]=="Down",]

    def get_pvalue_sorted_records_as_npmatrix(self):
        pvals=self.get_pvalue_as_nparray().astype(np.float)
        return self.rnpmatrix_[pvals.argsort(),]

    def get_pvalue_sorted_records_as_npmatrix(self, npmat):
        pvals=npmat[range(len(npmat)),5].astype(np.float)
        return npmat[pvals.argsort(),]

    def get_pvalue_filtered_records(self, pfilter): #pfilter less than or equal to pval
        pvals=self.get_pvalue_as_nparray().astype(np.float)
        return self.rnpmatrix_[pvals <= float(pfilter),]

    def get_pvalue_filtered_records(self, npmat, pfilter): #pfilter less than or equal to pval
        pvals=npmat[range(len(npmat)),5].astype(np.float)
        return self.npmat[pvals <= float(pfilter),]

    def get_fdr_sorted_records_as_npmatrix(self):
        fdr=self.get_fdr_as_nparray().astype(np.float)
        return self.rnpmatrix_[fdr.argsort(),]

    def get_fdr_sorted_records_as_npmatrix(self, npmat):
        fdr=npmat[range(len(npmat)),5].astype(np.float)
        return npmat[fdr.argsort(),]

    def get_fdr_filtered_records(self, pfilter): #pfilter less than or equal to pval
        fdr=self.get_fdr_as_nparray().astype(np.float)
        return self.rnpmatrix_[fdr <= float(pfilter),]

    def get_fdr_filtered_records(self, npmat, pfilter): #pfilter less than or equal to pval
        fdr=npmat[range(len(npmat)),5].astype(np.float)
        return self.npmat[fdr <= float(pfilter),]

    def get_records_as_npmatrix_with_minimum_ngenes(self, mingenes):
        ngenes=self.get_ngenes_as_nparray().astype(np.int)
        return self.rnpmatrix_[ngenes >= mingenes,]

#    def get_records_as_npmatrix_with_minimum_ngenes(self, npmat, mingenes):
#        ngenes=npmat[range(len(npmat)),1].astype(np.int)
#        return self.npmat[ngenes >= mingenes,]

    def get_gene_rank_map(self, mingene=2, direction="Down", sorttype="pvalue"):
    #Rank of the gene in a screen is based on NGENE, Direction and pvalue sorted
        filteredmat=self.get_records_as_npmatrix_with_minimum_ngenes(mingene)
        if direction=="Down":
            filteredmat=self.get_down_direction_records_as_npmatrix(filteredmat)
        elif direction=="Up":
            filteredmat=self.get_up_direction_records_as_npmatrix(filteredmat)
        if sorttype=="pvalue":
            filteredmat=self.get_pvalue_sorted_records_as_npmatrix(filteredmat)
        elif sorttype=="fdr":
            filteredmat=self.get_fdr_sorted_records_as_npmatrix(filteredmat)
        gene_rank={}
        rnk=1
        for fg in filteredmat[range(len(filteredmat)),0]:
            gene_rank[str(fg)]=rnk
            rnk=rnk+1
        return gene_rank
    #Return a normalized (scale 0-10) score for individual genes
    def get_gene_roast_score_map(self, mingene=2, direction="Down", sorttype="pvalue"):
        mgene_rank=self.get_gene_rank_map(mingene, direction, sorttype)
        ranksnp=np.array(mgene_rank.values()).astype(float)
        ranksnp_std= (ranksnp - ranksnp.min(axis=0)) / (ranksnp.max(axis=0) - ranksnp.min(axis=0))
        ranksnp_scaled = ranksnp_std * (10 - 0) + 0
        ranksnp_scaled=10-ranksnp_scaled #inverse is true
        grs_map={}
        ind=0
        for g in mgene_rank.keys():
            grs_map[g]=float(ranksnp_scaled[ind])
            ind=ind+1
        return grs_map
