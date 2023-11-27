# annotate PAS with gene annotation gtf
# python=3.6-3.10
# Bin Zhang
# Data: Nov 26, 2023
# version: 1.0
# ENV:

import os
import sys
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
#
from gtfparse import read_gtf
import portion as P

# extract last exons from each transcript for each gene and merge the overlapped ones 
def collapseExon(gene_id,exon_df):
    gene_df = exon_df[exon_df['gene_id'] == gene_id].copy()
    _chr = gene_df['seqname'].iloc[0]
    _strand = gene_df['strand'].iloc[0]
    _gene_name = gene_df['gene_name'].iloc[0]
    #txn_uid = gene_df['transcript_id'].unique()
    regs = P.empty()
    for i,row in gene_df.iterrows():
        regs = regs | P.closed(row['start'],row['end'])
    # format bed output
    j = 1
    out_bed = []
    if len(regs) == 1:
        each_r = regs[0]
        out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
    else:
        for each_r in regs:
            if j == 1:
                if _strand == '+':
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'first_exon',len(regs)]),j,_strand])
                else:
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
            elif j == len(regs):
                if _strand == '+':
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'last_exon',len(regs)]),j,_strand])
                else:
                    out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'first_exon',len(regs)]),j,_strand])
            else:
                out_bed.append([_chr,each_r.lower,each_r.upper,':'.join(str(_e) for _e in [gene_id,_gene_name,'mid_exon',len(regs)]),j,_strand])
            j += 1
    # 
    out_bed = pd.DataFrame(out_bed)
    return out_bed
        
# split the genome into 7 distinct categories, including: terminal_exon, ups_exon, intron, ncRNA, ups_tss1k, extend_pas10k, intergenic 
def annotatePAS(input_pas,anno_gtf,output,header):
    pas_df = pd.read_csv(input_pas,sep = '\t',header = header)
    #pas_bed = BedTool.from_dataframe(pas_df)
    gtf_df = read_gtf(anno_gtf)
    # split gtf for protein_coding and non-coding gene
    pcg_df = gtf_df[gtf_df['gene_type'] == 'protein_coding'].copy()
    ncg_df = gtf_df[gtf_df['gene_type'] != 'protein_coding'].copy()
    # 
    pcg_exon_df = pcg_df[pcg_df['feature'] == 'exon'].copy()
    # to speed up, split dataframe with chromosome 
    #chrs = pcg_exon_df['seqname'].unique()
    chrs = pas_df['chr']
    for _chr in chrs:
        #chr_pas_df = pas_df[pas_df.iloc[:,0] == _chr].copy()
        chr_pas_bed = BedTool.from_dataframe(pas_df[pas_df.iloc[:,0] == _chr])
        chr_pcg_exon_df = pcg_exon_df[pcg_exon_df['seqname'] == _chr].copy()
        chr_gids = chr_pcg_exon_df['gene_id'].unique()
        bed_list = list(map(lambda _gid: collapseExon(_gid,chr_pcg_exon_df),chr_gids))
        collapse_exon_df = pd.concat(bed_list)
        cexon_bed = BedTool.from_dataframe(collapse_exon_df)
        pas_cexon = chr_pas_bed.intersect(cexon_bed,s = True,wa = True, wb = True)
        pas_cexon_df = pas_cexon.to_dataframe()
        last_exon_pas = pas_cexon_df[pas_cexon_df['blockCount'].str.contains('last_exon')].copy()
        
        
        
        
    
        
    
    