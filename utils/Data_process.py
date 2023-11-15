# functions for basic data processing
# python=3.6-3.10
# Bin Zhang
# Data: Nov 14, 2023
# version: 1.0
# ENV:

import pandas as pd
from utils.PAS_utils import cluster2bed,toUpper,seq2kmer,addLabel2bed
#import pybedtools 

def data_preprocessing(input_file,flank_len,strand,reference,kmer,distance,annotation):
    if_filter_chr = True
    bed_seq = cluster2bed(input_file,flank_len,strand,if_filter_chr)
    print('[INFO] extracting sequence')
    #bed_seq = bed_seq.sequence(fi=reference,s=True,name = True,tab = True)
    # it seems getfasta much faster than sequence
    #bed_seq = bed_seq.getfasta(fi=reference,s=True,name = True,tab = True,fo = 'cluster.seq.txt')
    bed_seq = bed_seq.getfasta(fi=reference,s=True,name = True,tab = True)
    bed_seq_df = pd.read_csv(bed_seq.seqfn,sep = '\t',header = None)
    print('[INFO] finished extracting sequence')
    bed_seq_df.columns = ['info','sequence']
    bed_seq_df['sequence'] = bed_seq_df['sequence'].apply(toUpper)
    bed_seq_df['sequence'] = bed_seq_df['sequence'].apply(seq2kmer,k=kmer)
    bed_seq_df['label'] = 0
    bed_seq_df = bed_seq_df.loc[:,['sequence','label','info']]
    if annotation != None:
        print(f'[INFO] annotation:{annotation}')
        # overlap with annotation to label data
        bed1 = cluster2bed(input_file,distance,strand,if_filter_chr)
        bed1_df = addLabel2bed(bed1,annotation)
        bed_seq_df['label'] = bed1_df['label']
        #bed_seq_df['info'] = bed1_df['name']
    pas_bed = cluster2bed(input_file,1,strand,if_filter_chr)
    pas_bed_df = pas_bed.to_dataframe()
    pas_bed_df['anno'] = bed_seq_df['label']
    return [bed_seq_df,pas_bed_df]
    