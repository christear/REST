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

# split the genome into 7 distinct categories, including: terminal_exon, ups_exon, intron, ncRNA, ups_tss1k, extend_pas10k, intergenic 
def annotatePAS(input_pas,anno_gtf,output,header):
    pas_df = pd.read_csv(input_pas,sep = '\t')
    pas_bed = BedTool.from_dataframe(pas_df)
    
    