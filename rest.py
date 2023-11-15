# REST(RNA 3'End Sequencing data analysis Toolkit)
# python=3.6-3.10
# Bin Zhang
# Data: Nov 14, 2023
# version: 1.0
# ENV:

import os
import sys
import argparse
from utils.Cluster import callCluster
from utils.Data_process import data_preprocessing

def call_cluster_main(args):
    print(f'[INFO] call cluster based on {args.input_file}')
    callCluster(args.input_file,args.input_format,args.strand,args.output,args.output_dis)
    
def filter_cluster_main(args):
    print(f'[INFO] filter cluster in {args.input_file}')
    run=args.run
    if run=='pred':
        model=args.model
    #file=args.file
    input_file=args.input_file
    rounds=args.rounds
    kmer=args.kmer
    batch_size=args.batch_size
    train_epoch=args.train_epoch
    max_seq_len=args.max_seq_len
    n_process=args.n_process
    motif_file=args.motif_file
    #
    DNABERT_path = args.DNABERT_path
    out_dir = args.out_dir
    strand = args.strand
    distance = args.distance
    annotation = args.annotation
    reference = args.reference
    # checking 
    if reference == None:
        print(f'[INFO] the reference genome is required for both training and prediction')
        sys.exit()
    print(f'[INFO] reference:{reference}')
    #
    flank_len = int(max_seq_len/2)
    bed_seq_df,pas_bed_df = data_preprocessing(input_file,flank_len,strand,reference,kmer,distance,annotation,)
    print('[INFO] finished data processing')
    #
    if run=='train':
        if annotation == None:
            print(f'[INFO] annotation is required for training')
            sys.exit()
        print('[INFO] start training pipeline')
        print(f'[INFO] N rounds: {rounds}')
        print(f'[INFO] kmer: {kmer}')
        #
        os.system(f'mkdir -p {out_dir}')
        print(f'[INFO] creat output folder: {out_dir}')
        for r in range(rounds):
            print(f'[INFO] **start training round {r}**')
            output_path=f'{out_dir}/active_r{r}_m{kmer}-0'
            # it the fine-tuned model exists, skipped the training 
            if os.path.isfile(f'{output_path}/pytorch_model.bin'):
                print(f'[INFO] a pre-trained model of round {r} exists')
                continue
            os.system('mkdir -p {}'.format(output_path))
            print(f'[INFO] save files in this round under: {output_path}')
            if r==0:
                bed_seq_df.to_csv(f'{output_path}/train.tsv',index=False,sep='\t')
                print(f'[INFO] saved {output_path}/train.tsv')      
            elif r!=0:
                print('[INFO] train.tsv should exist, from the relabel of last round')
            print('[INFO] start finetune model')
            ### checking if this requires for 
            init_model = 
            if kmer==3:
                init_model='model/3-new-12w-0'
            elif kmer==4:
                init_model='model/4-new-12w-0'
            elif kmer==5:
                init_model='model/5-new-12w-0'
            elif kmer==6:
                init_model='model/6-new-12w-0'
                
        
    elif run=='pred':
        print('[INFO] start pred pipeline')
        print(f'[INFO] save files to prediction folder {our_dir}')
        pred_path=out_dir
        os.system('mkdir -p {}'.format(pred_path))
        bed_seq_df.to_csv(f'{pred_path}/dev.tsv',index=False,sep='\t')
        print('[INFO] start prediction')
        
    


### main function 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "REST:RNA 3\'End Sequencing data analysis Toolkit", add_help = 'REST Availbale tools: call_cluster, filter_cluster')
    subparsers = parser.add_subparsers(help = 'sub-command help')
    
    # subfunction: call_cluster
    call_cluster = subparsers.add_parser('call_cluster', help='call cluster based on alignment bam/bed file')
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
    filter_cluster.add_argument('--strand',default=1, help='strand of the cluster, 1: forward strand, 2; reverse strand, 0: strandless')
    filter_cluster.add_argument('--distance',default=50, help='distance threshold to define the overlap with annotation, default = 50')
    filter_cluster.add_argument('--annotation',default=None, help='PAS annotation to define true and false for training the model')
    filter_cluster.add_argument('--reference',default=None, help='the reference genome used to extract sequence flanking peaks of each cluster')
    filter_cluster.set_defaults(func=filter_cluster_main)
    #
    if len(sys.argv) < 2:
        sys.exit()
    args = parser.parse_args()
    #print(args)
    args.func(args)
