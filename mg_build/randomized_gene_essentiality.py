#!/usr/bin/python

import networkx as nx
import numpy as np
from collections import Counter 
import openpyxl as xl
import sys
from cascadingfailureanalysis import *


def main(argv):
    net_file,experimental_file,p_step,p_max,target_reaction,replica = tuple(argv)
    p_step = float(p_step)
    p_max = float(p_max)
    g = nx.read_gml(net_file)
    
    experimental_classification = {}
    label = {"Y":True,"N":False}
    wb = xl.load_workbook(experimental_file)
    sh=wb['Sheet1']
    for i in range(1,sh.max_row):
        essential = label[sh['E'+str(i+1)].value]
        mol = sh['C'+str(i+1)].value
        if mol != None:
            experimental_classification[mol] = essential

    percentages = np.arange(0.0,p_max,p_step)

    mol_nodes = [n for n,m in g.nodes(data=True) if m['type']=='m']
    #react_nodes = [n for n,m in g.nodes(data=True) if m['type']=='r' and m['reacttype'] not in ['DNAReplicationReaction']]
    react_nodes = [n for n,m in g.nodes(data=True) if m['type']=='r']
    
    results = []
    for p in percentages:
        if p:
            randomize_bipartite_network(g,p_step,mol_nodes,react_nodes)
        local_gene_is_essential = {m:cascade_failure(g,m,target_reaction=target_reaction)[1] for m in experimental_classification.keys()}
        counts = Counter(local_gene_is_essential.values())
        total_essential = counts[True]
        classifications = []
        for gene in experimental_classification.keys():
            if local_gene_is_essential[gene] == experimental_classification[gene]:
                classifications.append('correct')
            elif experimental_classification[gene]:
                classifications.append('fn')
            else:
                classifications.append('fp')
        counts = Counter(classifications)
        correct = float(counts['correct'])/len(classifications)
        false_positive = float(counts['fp'])/len(classifications)
        false_negative = float(counts['fn'])/len(classifications)
        print(p,replica,total_essential,correct,false_positive,false_negative,sep='\t')

if __name__ == "__main__":
    main(sys.argv[1:])