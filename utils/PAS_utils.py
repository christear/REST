# relabelPAS
# python=3.6-3.10
# Bin Zhang
# Data: Nov 12, 2023
# version: 1.0
# ENV: 

import re
import sys
import pandas as pd

# connect kmer to sequence: wrote by myself 
def kmer2seq2 (kmer_seq):
    kmer_list = re.split(' ',kmer_seq)
    _cs = ''
    for i in range(0,len(kmer_list) - 1):
        _cs += kmer_list[i][0]
    _cs += kmer_list[len(kmer_list) - 1]
    return _cs

# copied from DNABERT/motif/motif_utils.py 
def seq2kmer(seq, k):
    """
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space

    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers
    
def kmer2seq(kmers):
    """
    Convert kmers to original sequence
    
    Arguments:
    kmers -- str, kmers separated by space.
    
    Returns:
    seq -- str, original sequence.

    """
    kmers_list = kmers.split(" ")
    bases = [kmer[0] for kmer in kmers_list[0:-1]]
    bases.append(kmers_list[-1])
    seq = "".join(bases)
    assert len(seq) == len(kmers_list) + len(kmers_list[0]) - 1
    return seq
    
    
# determine internal primming (IP)
# 7 consecutive A or more than 8 A among 10 nucleotide 
def determineIP (seq):
    _subseq = seq[100:110]
    #_subseq = seq
    if re.match(r'.*A{7}',_subseq):
        _ip = 1
    else:
        _na = 0
        for _n in _subseq:
            if _n == "A":
                _na += 1
        if _na > 7:
            _ip = 1
        else:
            _ip = 0
    return _ip

### determine PAS motif
# motif_list is an array of PAS motifs ranked by frequency, human and mouse may be different from the third  
# k: only check the top k motifs  
def determineMotif (seq,motif_list,k):
    # only using the upstream sequences of a PAS 
    _subseq = seq[0:104]
    _nm = 0
    for _i in range(0,k):
        if motif_list[_i] in _subseq:
            _nm += 1
            #_nm = 1
    return _nm
    
def relabelPred (pred_txt,motif_file,output):
    pred_df = pd.read_csv(pred_txt,sep = '\t')
    cseq_list = pred_df["sequence"].apply(kmer2seq)
    ip_res = cseq_list.apply(determineIP)
    motif_df = pd.read_csv(motif_file,sep = '\t',header = None)
    motif_list = motif_df[0].values.tolist()
    # the top 3 PAS motif
    motif_res = cseq_list.apply(determineMotif,motif_list = motif_list,k = 3)
    # all possible PAS motif 
    motif_res2 = cseq_list.apply(determineMotif,motif_list = motif_list,k = len(motif_list))
    relabel = []
    #n = 0
    for i in range(0,len(motif_res)):
        # relabel PAS in FN without PAS motif and with internal primming as negative
        if pred_df["label"][i] == 1 and pred_df["pred_01"][i] == 0 and motif_res2[i] == 0 and ip_res[i] == 1:
            relabel.append(0)
        # relabel predicted PAS not from ground truth (FP)
        elif pred_df["label"][i] == 0 and pred_df["pred_01"][i] == 1:
            # relabel FP with PAS motif and without internal primming as positive
            if motif_res[i] > 0 and ip_res[i] == 0:
                relabel.append(1)
            # relabel FP without PAS motif or with internal primming as negative
            else:
                relabel.append(0)
        else:
            #n += 1
            relabel.append(pred_df["label"][i])
    out_df = pred_df
    out_df['relabel'] = relabel
    if output != None:
        repred.to_csv(output, index=False, sep = "\t")
    #pred_df['relabel'] = relabel
    return out_df