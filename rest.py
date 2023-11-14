# REST(RNA 3'End Sequencing data analysis Toolkit)
# python=3.6-3.10
# Bin Zhang
# Data: Nov 14, 2023
# version: 1.0
# ENV:

import os
import sys
import argparse

def call_cluster_main(args):
    print(f'call cluster')
    
def filter_cluster_main(args):
    print(f'filter cluster')


### main function 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "REST:RNA 3\'End Sequencing data analysis Toolkit", add_help = 'REST Availbale tools: call_cluster, filter_cluster')
    subparsers = parser.add_subparsers(help = 'sub-command help')
    
    # subfunction: call_cluster
    call_cluster = subparsers.add_parser('call_cluster', help='call cluster from alignment bam/bed file')
    call_cluster.add_argument('--input_file', required=True, help = 'input file in bam or bed format')
    call_cluster.add_argument('--input_format', default = 'bed', help = 'format of the input file')
    call_cluster.add_argument('--strand', type = int, default = 2, help = 'strand of the sequencing data that generate the input file. 1: forward strand, 2: reverse strand, 0: strandless')
    call_cluster.add_argument('--output', default = 'tmp.cluster', help = 'output file of the identified cluster')
    call_cluster.add_argument('--output_dis', default = 'tmp.dis', help = 'distance of the output cluster')
    call_cluster.set_defaults(func = call_cluster_main)
    
    # subfunction: filter_cluster
    filter_cluster = subparsers.add_parser('filter_cluster', help = 'filter cluster file based on fine-tuned DNABERT models')
    filter_cluster.add_argument('--run', default = 'pred', help = 'pred/train: only predict true/false for fintering, or also train a fine-tuned model')
    filter_cluster.add_argument('--input_file', required = True, help = 'input file of cluster in bed-like format')
    filter_cluster.add_argument('--rounds', type = int, default = 5, help = 'rounds of the training, eg: 5')
    filter_cluster.add_argument('--model', default='model/5-new-12w-0', help='eg: /active_r0_m5-0, default=model/5-new-12w-0')
    filter_cluster.add_argument('--kmer', type = int,default=5,help='eg: 5, default=5')
    filter_cluster.add_argument('--batch_size',type=int,default=100,help='eg: 100, default=100')
    filter_cluster.add_argument('--train_epoch',type=int,default=5,help='number of train epochs, eg: 5, default=5')
    filter_cluster.add_argument('--max_seq_len',type=int,default=200,help='eg: 200, default=200')
    filter_cluster.add_argument('--n_process',type=int,default=100,help='number of processors, eg: 100, default=100')
    filter_cluster.add_argument('--motif_file',default='mouse.PAS.motif',help='path to motif file, default=mouse.PAS.motif')
    filter_cluster.add_argument('--DNABERT_path',default='DNABERT/examples', help='path to DNABERT script, default=DNABERT/examples')
    filter_cluster.add_argument('--out_dir',default='test_out', help='directory for output, default=test_out')
    filter_cluster.add_argument('--strand',default=2, help='strand of the cluster, 1: forward strand, 2; reverse strand, 0: strandless')
    filter_cluster.add_argument('--distance',default=50, help='distance threshold to define the overlap with annotation, default = 50')
    filter_cluster.add_argument('--annotation',default=None, help='PAS annotation to define true and false for training the model')
    filter_cluster.add_argument('--reference',default=None, help='the reference genome used to extract sequence flanking peaks of each cluster')
    filter_cluster.set_defaults(func=filter_cluster_main)
    #
    args = parser.parse_args()
    if 'call_cluster' in args:
        print(call_cluster)
    elif 'filter_cluster' in args:
        print(filter_cluster)
    elif 'count' in args:
        print('counting')
    else:
        sys.exit()
    print(args)
    args.func(args)
