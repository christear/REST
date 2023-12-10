# countGene
# python=3.6-3.10
# Bin Zhang
# Data: Dec 10, 2023
# version: v1
# ENV: 
#

import re
import sys
import pandas as pd
from utils.PAS_utils import get_gene

# based on pas count and pas annotation to aggregate PAS count for each gene 
def countGene (pas_count,pas_anno,output,select_types):
    print('### aggregate PAS count for each gene ')
    pc_df = pd.read_csv(pas_count,sep = '\t',comment = '#')
    annopas_df = pd.read_csv(pas_anno,header = None,sep = '\t')
    #select_types = 'last_exon,intron_anno,intron_unanno,ups_exon'
    genepas_df = annopas_df[annopas_df.iloc[:,6].isin(select_types.split(','))].copy()
    genepas_df['genes'] = genepas_df.iloc[:,9].apply(get_gene,sep = ':')
    ugenes = genepas_df['genes'].unique()
    gene_count = []
    i = 1
    for each_gene in ugenes:
        if i % 1000 == 0:
            print(f'### proccessed {i} genes')
        pas = genepas_df[genepas_df['genes'] == each_gene].iloc[:,3]
        each_count = pc_df[pc_df.iloc[:,0].isin(pas)].copy()
        dat = each_count.iloc[:,6:each_count.shape[1]]
        #each_out = pd.DataFrame()
        #each_out['gene'] = each_gene
        each_gc = dat.sum(axis = 0)
        gene_count.append(each_gc)
        i += 1
        #
    gc_df = pd.DataFrame({'genes':ugenes})
    gc_df = pd.concat([gc_df,pd.DataFrame(gene_count)],axis = 1)
    if output != None:
        gc_df.to_csv(output,sep = '\t',index = False)
    return gc_df
    
if len(sys.argv) < 3:
    print("Usage:python countGene.py pas_count_txt pas_anno_txt output_gene_count")
    sys.exit(1)
else:
    select_types = 'last_exon,intron_anno,intron_unanno,ups_exon,ncRNA'
    pau_df = countGene(sys.argv[1],sys.argv[2],sys.argv[3],select_types)