#!/usr/bin/env python

import numpy as np
import os
import sys

__author__ = "Sanket S Desai"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "desai.sanket12@outlook.com"
'''
EdgeRTopTagsParser - class to parse the top tag results as given output by edgeR package. ProcessHairpin Function
Format Description - as shown in the example below (10 columns including the ID repeated twice at every record)
"ID"	"Sequences"	"Gene.symbol"	"Plate"	"logFC"	"logCPM"	"LR"	"PValue"	"FDR"
"RLGH-GN47331_0.103"	"RLGH-GN47331_0.103"	"CCTTTTGTAGGATGCATTCCT"	"PDIK1L"	"247-0001"	-5.20958150765904	10.436427154501	4.51499597689099	0.0335989651347352	0.999923418061804
Column->data
0 ID_rowname
1 ID
2 Sequence
3 Gene.symbol
4 Plate
5 logFC
6 logCPM
7 LR
8 PValue
9 FDR
'''

class EdgeRTopTagsParser(object):
    def __init__(self, toptagsfn):
        if not os.path.exists(toptagsfn) and not os.path.isfile(toptagsfn):
            print("EdgeR top tags file %s not found or unreadable, please check!!" %(rnf))
            sys.exit(0)
        toptagsfi=open(toptagsfn)
        records=[]
        self.header_=toptagsfi.readline().strip().replace('\"','').split("\t")
        for line in toptagsfi:
            sline=line.replace('\"','').split("\t")
            records.append(sline)
        self.toptagsmatrix_=np.array(records)
        toptagsfi.close()

    def get_column_as_nparray(self, i):
        return self.toptagsmatrix_[range(len(self.toptagsmatrix_)),i]

    def get_unique_genes(self):
        return list(np.unique(self.get_column_as_nparray(3)))

    def get_toptag_records_as_npmatrix_for_gene(self, gene):
        res=np.array
        try:
            res=self.toptagsmatrix_[ self.get_column_as_nparray(3) == gene ,]
        except IndexError:
            print("Gene %s not found in the toptag results, please check!!" %(gene))
        return res

    def get_fdr_significant_toptags_as_npmatrix(self, fdrthresh=0.05):
        return self.toptagsmatrix_[ self.get_column_as_nparray(9).astype(np.float) <= fdrthresh ,]

    def get_pvalue_significant_toptags_as_npmatrix(self, pvalthresh=0.05):
        return self.toptagsmatrix_[ self.get_column_as_nparray(8).astype(np.float) <= pvalthresh ,]

    def get_gene_average_logfc(self, gene):
        generecords=self.get_toptag_records_as_npmatrix_for_gene(gene)
        return float(np.mean(generecords[range(len(generecords)),5].astype(np.float)))

    def get_fdr_significant_gene_average_logfc(self, gene, fdrthresh=0.05):
        generecords=self.get_toptag_records_as_npmatrix_for_gene(gene)
        return float(np.mean(generecords[ generecords[range(len(generecords)),9] <= fdrthresh ,5]))

    def get_fdr_significant_gene_average_logfc(self, gene, pvalthresh=0.05):
        generecords=self.get_toptag_records_as_npmatrix_for_gene(gene)
        return float(np.mean(generecords[ generecords[range(len(generecords)),8] <= pvalthresh ,5]))

    def get_ids_as_nparray(self):
        return self.get_column_as_nparray(0)

    def get_sequence_by_id_as_nparray(self, id_):
        return self.toptagsmatrix_[ self.get_column_as_nparray(0)==id_ , 2 ]

    def get_gene_average_logfc_map(self):
        genes=self.get_unique_genes()
        gene_avlfc={}
        for g in genes:
            av=self.get_gene_average_logfc(g)
            gene_avlfc[g]=av
        return gene_avlfc

    def get_gene_average_logfc_score_map(self):
        gene_avlfc=self.get_gene_average_logfc_map()
        avglfcnp=np.array(gene_avlfc.values()).astype(float)
        avglfcnp_std= (avglfcnp - avglfcnp.min(axis=0)) / (avglfcnp.max(axis=0) - avglfcnp.min(axis=0))
        avglfcnp_scaled = avglfcnp_std * (10 - 0) + 0
        avglfcnp_scaled=10-avglfcnp_scaled #inverse is true
        grs_map={}
        ind=0
        for g in gene_avlfc.keys():
            grs_map[g]=float(avglfcnp_scaled[ind])
        return grs_map
