#!/usr/bin/env python

import numpy as np
from scipy.stats import rankdata
import os
import sys
import argparse
from parseedgertoptags import *
from parseroast import *
from parseexpression import *
from parsecnv import *

__author__ = "Sanket S Desai"
__license__ = "MIT"
__version__ = "0.2"
__email__ = "desai.sanket12@gmail.com"
#Smaller ranks are better
class DepRanker(object):
    def __init__(self, toptagsfn, roastfn, exprsfn, cnvfn):
        toptags=EdgeRTopTagsParser(toptagsfn)
        roast=RoastParser(roastfn)
        self.roast_ranks_= self.get_gene_rank_map(roast.get_gene_rank_map(),-1)
        self.roast_scores_=roast.get_gene_roast_score_map()
        self.toptags_ranks_= self.get_gene_rank_map(toptags.get_gene_average_logfc(), -1)
        self.toptags_scores_=toptags.get_gene_average_logfc_score_map()
        expression=ExpressionParser(exprsfn)
        copynumber=CopyNumberVariationParser(cnvfn)
        topgenes=roast.get_roast_genes()
        topgenesexpmap=expression.get_gene_expression_map(topgenes)
        topgenescnvmap=copynumber.get_gene_cnv_map(topgenes)
        self.expression_ranks_=self.get_gene_rank_map(topgenes, 1)
        self.expression_scores_=self.get_gene_score_map(topgenesexpmap,1)
        self.cnv_ranks_=self.get_gene_rank_map(topgenescnvmap, 1)
        self.cnv_scores_=self.get_gene_score_map(topgenescnvmap,1)
        self.genes_=topgenes
#        if len(self.roast_scores_)!=len(self.expression_scores_) or len(self.toptags_scores_) != len(self.roast_scores_) or len(self.roast_scores_) != len(cnv_scores_) :
#            print("DepRanker error: please contact the author!!")
#            sys.exit()

    #if direction 1, defined as lower is lower, elif direction -1, lower is higher: In general higher ranks are better
    def get_gene_rank_map(self, gene_value_map, direction=1):
        grmap={}
        if direction==1:
            gvals=np.array(gene_value_map.values()).astype(float)
            gvalranks=rankdata(np.array(gvals))
            ind=0
            for g in gene_value_map.keys():
                grmap[g]=float(gvalranks[ind])
                ind=ind+1
        elif direction==-1:
            gvalranks=rankdata(- np.array(gvals))
            ind=0
            for g in gene_value_map.keys():
                grmap[g]=float(gvalranks[ind])
                ind=ind+1
        else:
            print("Direction of ranking required; please mention 1 or -1 !!")
            sys.exit(0)
        return grmap

    def get_gene_score_map(self, gene_value_map, direction): #direction defines if high should be considered high or low (as in case of pooled screen logFC)
        grmap={}
        if direction==1:
            avglfcnp=np.array(gene_value_map.values()).astype(float)
            avglfcnp_std= (avglfcnp - avglfcnp.min(axis=0)) / (avglfcnp.max(axis=0) - avglfcnp.min(axis=0))
            avglfcnp_scaled = avglfcnp_std * (10 - 0) + 0
            ind=0
            for g in gene_value_map.keys():
                grmap[g]=float(avglfcnp_scaled[ind])
                ind=ind+1
        elif direction==-1:
            avglfcnp=np.array(gene_value_map.values()).astype(float)
            avglfcnp_std= (avglfcnp - avglfcnp.min(axis=0)) / (avglfcnp.max(axis=0) - avglfcnp.min(axis=0))
            avglfcnp_scaled = avglfcnp_std * (10 - 0) + 0
            avglfcnp_scaled=10-avglfcnp_scaled #inverse is true
            ind=0
            for g in gene_value_map.keys():
                grmap[g]=float(avglfcnp_scaled[ind])
                ind=ind+1
        else:
            print("Direction of ranking required; please mention +1 or -1 !!")
            sys.exit(0)
        return grmap

    def get_gene_rank_impact_score(self, gene):
        rimpscore=0
        if gene in self.roast_ranks_:
            try:
                rimpscore=self.roast_ranks_[gene]
            except KeyError:
                rimpscore=rimpscore+0
        if gene in self.toptags_ranks_:
            try:
                rimpscore=self.toptags_ranks_[gene]
            except KeyError:
                rimpscore=rimpscore+0
        if gene in self.expression_ranks_:
            try:
                rimpscore=self.expression_ranks_[gene]
            except KeyError:
                rimpscore=rimpscore+0
        if gene in self.cnv_ranks_:
            try:
                rimpscore=self.cnv_ranks_[gene]
            except KeyError:
                rimpscore=rimpscore+0
        return rimpscore

    def get_gene_impact_score(self, gene):
        impscore=0
        #currently only consideres scoring if the gene score is available for all the 4 types
        if gene in self.roast_scores_:
            try:
                impscore=self.roast_scores_[gene]
            except KeyError:
                impscore=impscore+0
        if gene in self.toptags_scores_:
            try:
                impscore=impscore+self.toptags_scores_[gene]
            except KeyError:
                impscore=impscore+0
        if gene in self.expression_scores_:
            try:
                impscore=impscore+self.expression_scores_[gene]
            except KeyError:
                impscore=impscore+0
        if gene in self.cnv_scores_:
            try:
                impscore=impscore+self.cnv_scores_[gene]
            except KeyError:
                impscore=impscore+0
        return impscore

    def write_score_table(self, fn):
        fo=open(fn,'w')
        fo.write("Gene\tRoastScore\tPooledTopTagScore\tExpressionScore\tCNVScore\tImpactScore\n")
        for g in self.genes_:
            try:
                #print([g,self.roast_scores_[g],self.toptags_scores_[g],self.expression_scores_[g],self.cnv_scores_[g],self.get_gene_impact_score(g)])
                fo.write("%s\t%f\t%f\t%f\t%f\t%f\n" %(g,self.roast_scores_[g],self.toptags_scores_[g],self.expression_scores_[g],self.cnv_scores_[g],self.get_gene_impact_score(g)))
            except KeyError as ke:
                print("Impact Score cannot be computed for gene: %s" %(g))
        fo.close()
    def write_score_table(self, fn):
        fo=open(fn,'w')
        fo.write("Gene\tRoastRankScore\tPooledTopTagRankScore\tExpressionRankScore\tCNVRankScore\tRankImpactScore\n")
        for g in self.genes_:
            try:
                #print([g,self.roast_scores_[g],self.toptags_scores_[g],self.expression_scores_[g],self.cnv_scores_[g],self.get_gene_impact_score(g)])
                fo.write("%s\t%f\t%f\t%f\t%f\t%f\n" %(g,self.roast_ranks_[g],self.toptags_ranks_[g],self.expression_ranks_[g],self.cnv_ranks_[g],self.get_gene_rank_impact_score(g)))
            except KeyError as ke:
                print("Impact Rank Score cannot be computed for gene: %s" %(g))
        fo.close()

def main():
    aparser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description= '''DepRanker: A gene impact score calculator (prioritization method) for RNAi / CRISPR screen results
Developed by Dutt lab
Version v0.1.0''')
    aparser.add_argument('-roast', action='store', dest='roast_result_file', help='Roast result file')
    aparser.add_argument('-toptags', action='store', dest='edgeR_toptags_file', help='EdgeR toptags result file')
    aparser.add_argument('-exprs', action='store', dest='expression_file', help='Gene expression file')
    aparser.add_argument('-cnv', action='store', dest='copy_number_variation_file', help='Copy number variation file')
    aparser.add_argument('-out', action='store', dest='output_file', help='Ouput file')
    pargs=aparser.parse_args()
    try:
        depranks=DepRanker(pargs.edgeR_toptags_file, pargs.roast_result_file, pargs.expression_file, pargs.copy_number_variation_file)
        depranks.write_score_table(pargs.output_file)
    except TypeError as t:
        aparser.print_help()
        print(t)

if __name__=="__main__":
    main()
