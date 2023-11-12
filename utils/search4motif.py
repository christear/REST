# search4motif
# python=3.9
# Bin Zhang
# Data: Sep 15, 2021
# version: 1.0
# ENV: 

import re
import sys

# read motif file, each line include one motif 
def readMotif(motif_file):
    motifs = []
    with open(motif_file,'r') as m:
        for _line in m:
            _line = _line.rstrip()
            motifs.append(_line)
    return motifs
            
# search for motifs in the provided fasta file 
def search4motif (input_fa,motif_file,start,end,out_tab):
    motifs = readMotif(motif_file)
    print(motifs)
    sequences = {}
    _start = int(start)
    _end = int(end)
    with open(input_fa,'r') as f:
        #_id = ''
        _id = 0
        _seq = ''
        _n = 0
        for _line in f:
            _line = _line.rstrip()
            if '>' in _line:
                #_id = _line[1:]
                #_n += 1
                if _id != 0:
                    sequences[str(_id)] = _seq
                    _seq = ''
                _id += 1
                    #_seq = ''
            else:
                _subs = _line.upper()
                _seq = _seq + _subs
        sequences[str(_id)] = _seq
    with open(out_tab,'w') as w:
        _sep = '\t'
        for _id in sequences:
            outs = []
            outs.append(_id)
            _seq = sequences[_id]
            #print(_id,_seq)
            _subs = _seq[_start : _end]
            #positions = []
            for _m in motifs:
                _p = -1
                if _m in _subs:
                    _p = _subs.index(_m)
                    #print(_m,_p)
                #positions.append(_p)
                outs.append(_p + 1)
            _out_line = _sep.join(str(e) for e in outs) + '\n'
            w.write(_out_line)
        
#"""
if len(sys.argv) < 3:
    print("Usage:python search4motif.py input_fa motif_file start end output_tab")
    sys.exit(1)
else:
    search4motif(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    
        