# callCluster
# python=3.6-3.10
# Bin Zhang
# Data: Nov 30, 2023
# version: v1
# ENV: 
#

import re
import sys
import pandas as pd
from collections import Counter

# get gene id with multiple filed seprated with ':'
def get_gene (gene_str,sep):
    inf = gene_str.split(sep)
    gene = inf[0] + ':' + inf[1]
    return gene
    
# the counts_txt can be the output file from featureCounts
# with comments line staring with '#' and the first 6 columns are: 'Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length' 
# pas_anno is the output from annotatePAS, bed-like format, the column 7 to 10 should be: 'pas_type', 'feature_start', 'feature_end', 'gene_inf'
def calPAU(counts_txt,pas_anno,pau_txt,select_types):
    print(f'### calculating PAU for {counts_txt}')
    counts_df = pd.read_csv(counts_txt,sep = '\t',comment = '#')
    annopas_df = pd.read_csv(pas_anno,header = None,sep = '\t')
    genepas_df = annopas_df[annopas_df.iloc[:,6].isin(select_types.split(','))].copy()
    #annopas_df = annopas_df[annopas_df.iloc[:,9].notna()].copy()
    #genepas_df = annopas_df[(annopas_df.iloc[:,6] == 'last_exon') | (annopas_df.iloc[:,6] == 'intron_anno') | (annopas_df.iloc[:,6] == 'intron_unanno') | (annopas_df.iloc[:,6] == 'ups_exon')].copy()
    genepas_df['genes'] = genepas_df.iloc[:,9].apply(get_gene,sep = ':')
    gene_pas = genepas_df.drop_duplicates(subset=[genepas_df.columns[3], genepas_df.columns[10]])
    ugenes = gene_pas['genes'].unique()
    #pau_df = pd.DataFrame()
    pau_list = []
    i = 1
    for each_gene in ugenes:
        if i % 1000 == 0:
            print(f'### proccessed {i} genes')
        pas = gene_pas[gene_pas['genes'] == each_gene].iloc[:,3]
        each_count = counts_df[counts_df.iloc[:,0].isin(pas)].copy()
        dat = each_count.iloc[:,6:each_count.shape[1]]
        each_out = pd.DataFrame()
        each_out['PAS_id'] = each_count['Geneid']
        each_out['gene'] = each_gene
        prop = dat.div(dat.sum(axis=0), axis=1)
        each_out = pd.concat([each_out,prop],axis = 1)
        pau_list.append(each_out)
        i += 1
    #
    pau_df = pd.concat(pau_list,axis = 0)
    if pau_txt != None:
        pau_df.to_csv(pau_txt,sep = '\t',index = False)
    return pau_df
    
# length from the last exon start to the PAS 
def le_length(df_row):
    row_list = df_row.values
    if row_list[5] == '+':
        _l = row_list[2] - row_list[7]
    else:
        _l = row_list[8] - row_list[2]
    return _l+1
    
# calculate weighted 3' UTR length index (WULI)    
def calWULI (counts_txt,pas_anno,wuli_txt):
    pau_df = calPAU(counts_txt,pas_anno,None,'last_exon')
    gene_counts = Counter(pau_df['gene'].to_list())
    annopas_df = pd.read_csv(pas_anno,header = None,sep = '\t')
    lepas_dfo = annopas_df[annopas_df.iloc[:,6] == 'last_exon'].copy()
    lepas_len = lepas_dfo.apply(le_length,axis = 1)
    lepas_len[lepas_len < 0] = 1
    anno_lelen = lepas_dfo.iloc[:,8] - lepas_dfo.iloc[:,7]
    genes = lepas_dfo.iloc[:,9].apply(get_gene,sep = ':')
    len_ratio = lepas_len/anno_lelen
    len_ratio[len_ratio > 1] = 1
    lepas_df = pd.DataFrame({'PAS_id':lepas_dfo.iloc[:,3],'gene':genes,'le_len':anno_lelen,'pas_len':lepas_len,'len_ratio':len_ratio})
    mpas_gene = [key for key, value in gene_counts.items() if value > 1]
    i = 1
    wuli_list = []
    print(f'### calculating WULI for {counts_txt}')
    for each_gene in mpas_gene:
        if i % 1000 == 0:
            print(f'### processed {i} genes')
        each_pau = pau_df[pau_df['gene'] == each_gene].copy()
        each_len = lepas_df[lepas_df['gene'] == each_gene].copy()
        #each_len = each_len[each_len.iloc[:,3].isin(each_pau['PAS_id'])].copy()
        if each_pau.iloc[:,0].tolist() == each_len.iloc[:,0].tolist():
            wuli_df = each_pau.iloc[:,2:each_pau.shape[1]].mul(each_len['len_ratio'].tolist(),axis = 0)
            wuli = wuli_df.sum(axis = 0)
        else:
            print(f'### PAS from gene {each_gene} don\'t match')
        wuli_list.append([each_gene,max(each_len['le_len'])] + wuli.tolist())
        i += 1
    wuli_df = pd.DataFrame(wuli_list)
    wuli_df.columns = ['gene','lastexon_len'] + pau_df.columns[2:pau_df.shape[1]].tolist()
    if wuli_txt != None:
        wuli_df.to_csv(wuli_txt,sep = '\t',index = False)
    return wuli_df
    
    
if len(sys.argv) < 3:
    print("Usage:python calPAS.py counts_txt pas_anno pau_txt")
    sys.exit(1)
else:
    select_types = 'last_exon,intron_anno,intron_unanno,ups_exon'
    pau_df = calPAU(sys.argv[1],sys.argv[2],sys.argv[3],select_types)
    #tmp_df = calWULI(counts_txt,pas_anno,None)
    #cd47wuli = tmp_df[tmp_df['gene'].str.contains('CD47')]
    #cd47wuli.values
    

