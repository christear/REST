# REST(RNA 3'End Sequencing data analysis Toolkit)
# python=3.6-3.10
# Bin Zhang
# Data: Nov 14, 2023
# Last modidication: Nov 28, 2023
# version: 1.0
# ENV:

import os
import sys
import argparse
import numpy as np 
import pandas as pd
from utils.Cluster import callCluster, mergeCluster
from utils.PAS_utils import relabelPred, annotatePAS
from utils.Data_process import data_preprocessing
from utils.PAS_utils import calPAU, calWULI
from utils.PAS_utils import get_gene

def call_cluster_main(args):
    print(f'[INFO] call cluster based on {args.input_file}')
    callCluster(args.input_file,args.input_format,args.strand,args.output,args.output_dis)

def merge_cluster_main(args):
    print(f'[INFO] merge cluster based on {args.file_list}')
    if args.with_header != True:
        header = None
    else:
        header = 'infer'
    mergeCluster(args.file_list,args.output,args.read,args.sam_num,args.distance,header,args.sam_label)

def annotate_PAS_main(args):
    print(f'[INFO] annotate PAS for {args.input_pas}')
    if args.with_header != None:
        header = 'infer'
    else:
        header = args.with_header
    annotatePAS(args.input_pas,args.gtf,args.output,header,args.distance)

def calcu_PAS_main(args):
    counts_txt = args.pas_count
    pas_anno = args.pas_anno
    output = args.output
    method = args.method
    pas_type = args.pas_type
    pos = args.pos
    if 'PAU' in method or 'pau' in method:
        print('[INFO] calculate PAS usage (PAU)')
        pau_df = calPAU(counts_txt,pas_anno,output,pas_type)
    elif 'WULI' in method or 'wuli' in method:
        print('[INFO] calculate weighted 3\' UTR length index')
        wuli_df = calWULI(counts_txt,pas_anno,output,pos = pos)
    else:
        print(f'[Error] {method} has not be defined')
        
from utils.PAS_utils import get_gene

# based on pas count and pas annotation to aggregate PAS count for each gene 
# def countGene (pas_count,pas_anno,output):
def count_gene_main (args):
    print('[INFO] aggregate PAS count for each gene ')
    pas_count = args.pas_count
    pas_anno = args.pas_anno
    output = args.output
    select_types = args.pas_type
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
    
def count_PAS_main(args):
    bam = args.bam
    data_type = args.data_type
    pas = args.pas
    win = args.win
    output = args.output
    temp = args.temp
    strand = args.strand
    end = args.end
    thread = args.thread
    #multiple = args.multiple
    if data_type == 'bulk':
        print('[INFO] counting for bulk RNA-seq data')
        if args.with_header != None:
            header = 'infer'
        else:
            header = args.with_header
        #
        pas_df = pd.read_csv(args.pas,sep = '\t',header = header)
        saf_df = pas_df.iloc[:,[3,0,1,2,5]].copy()
        saf_df.columns = ['GeneID','Chr','Start','End','Strand']
        saf_df['Start'] = saf_df['Start'] - win + 1
        saf_df['End'] = saf_df['End'] + win
        saf_df.to_csv(args.temp,index = False, sep = '\t')
        if args.multiple == True:
            print('### counting with mutlipe-mapped reads')
            if end == 5 or end == 3:
                print(f'### counting based on read {end} prime end')
                os.system(f'featureCounts -a {temp} -F SAF -o {output} -s {strand} --read2pos {end} -T {thread} -M {bam}')
                #os.system(f'rm -f {temp}')
            else:
                print(f'### counting based on whole read rather than end')
                os.system(f'featureCounts -a {temp} -F SAF -o {output} -s {strand} -T {thread} -M {bam}')
        else:
            print('### counting without mutlipe-mapped reads')
            if end == 5 or end == 3:
                print(f'### counting based on read {end} prime end')
                os.system(f'featureCounts -a {temp} -F SAF -o {output} -s {strand} --read2pos {end} -T {thread} {bam}')
                #os.system(f'rm -f {temp}')
            else:
                print(f'### counting based on whole read rather than end')
                os.system(f'featureCounts -a {temp} -F SAF -o {output} -s {strand} -T {thread} {bam}')
    elif 'sc' in data_type:
        print('[INFO] counting for scRNA-seq data')
    else:
        print(f'[INFO] unknown data type{data_type}')
    print('### counting Done')

def filter_cluster_main(args):
    print(f'[INFO] filter cluster in {args.input_file}')
    run=args.run
    if run=='pred':
        model=args.model
    else:
        init_model = args.model
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
    if DNABERT_path is None:
        #repo_root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        DNABERT_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'DNABERT/examples')
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
    bed_seq_df,pas_bed_df = data_preprocessing(input_file,flank_len,strand,reference,kmer,distance,annotation)
    #print([input_file,flank_len,strand,reference,kmer,distance,annotation])
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
            if str(kmer) not in init_model:
                print(f'[INFO] error: kmer{kmer} does not match with the model{init_model}')
            os.system(f'python {DNABERT_path}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{kmer} \
                --model_name_or_path {init_model} \
                --task_name dnaprom \
                --do_train \
                --data_dir {output_path} \
                --max_seq_length {max_seq_len} \
                --per_gpu_eval_batch_size={batch_size}   \
                --per_gpu_train_batch_size={batch_size}   \
                --learning_rate 2e-4 \
                --num_train_epochs {train_epoch} \
                --output_dir {output_path} \
                --save_steps 4000 \
                --warmup_percent 0.1 \
                --hidden_dropout_prob 0.1 \
                --overwrite_output \
                --weight_decay 0.01 \
                --n_process {n_process}')
            print(f'[INFO] finish finetune model under: {output_path}')
            print('[INFO] start prediction with finetuned model')
            pred_path=f'{output_path}/pred'
            os.system(f'mkdir -p {pred_path}')
            os.system(f'cp {output_path}/train.tsv {pred_path}/dev.tsv')
            os.system(f'python {DNABERT_path}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{kmer} \
                --model_name_or_path {output_path} \
                --task_name dnaprom \
                --do_predict \
                --data_dir {pred_path}  \
                --max_seq_length {max_seq_len} \
                --per_gpu_pred_batch_size={batch_size}   \
                --output_dir {output_path} \
                --predict_dir {pred_path} \
                --n_process {n_process}')
            print(f'[INFO] finished prediction under: {pred_path}')
            print('[INFO] start process prediction result')
            data=pd.read_csv(f'{pred_path}/dev.tsv',sep='\t')
            _pred=np.load(f'{pred_path}/pred_results.npy')
            threshold=0.5
            _pred01=np.int8(_pred>=threshold)
            data['pred_value']=_pred
            data['pred_01']=_pred01
            data.to_csv(f'{pred_path}/pred.tsv',index=False,sep='\t')
            print('[INFO] finished process prediction result')

            print('[INFO] start relabel')
            relabel_data = relabelPred(f'{pred_path}/pred.tsv',motif_file,f'{pred_path}/pred.relabeled.tsv')
            print('[INFO] finished relabel')
            print('[INFO] start process relabeled results to new train data for next round')
            output_path='{}/active_r{}_m{}-0'.format(out_dir,r+1,kmer)
            os.system(f'mkdir -p {output_path}')
            relabel_data=relabel_data.loc[:,['sequence','relabel','info']]
            relabel_data.columns=['sequence', 'label', 'info']
            relabel_data.to_csv(f'{output_path}/train.tsv',sep='\t',index=False)
            print('[INFO] finish process')
            print('[INFO] **finished training round {}**'.format(r))
        # summarize multiple rounds predictions 
        print('[INFO] summarize predictions')
        for r in range(rounds):
            round_df = pd.read_csv(f'{out_dir}/active_r{r}_m{kmer}-0/pred/pred.relabeled.tsv',sep = '\t')
            pas_bed_df[f'inputlabel_r{r}'] = round_df['label']
            pas_bed_df[f'predvalue_r{r}'] = round_df['pred_value']
            pas_bed_df[f'pred01_r{r}'] = round_df['pred_01']
            pas_bed_df[f'outputlabel_r{r}'] = round_df['relabel']
        #
        pas_bed_df.to_csv(f'{out_dir}/train.summary',sep = '\t',index = False)
    elif run=='pred':
        print('[INFO] start pred pipeline')
        print(f'[INFO] save files to prediction folder {out_dir}')
        pred_path=out_dir
        os.system('mkdir -p {}'.format(pred_path))
        bed_seq_df.to_csv(f'{pred_path}/dev.tsv',index=False,sep='\t')
        #ap = bed_seq_df['label'].sum()
        #print(f'{ap} pas overlapped with annotation')
        print('[INFO] start prediction')
        os.system(f'python {DNABERT_path}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{kmer} \
                --model_name_or_path {model} \
                --task_name dnaprom \
                --do_predict \
                --data_dir {pred_path}  \
                --max_seq_length {max_seq_len} \
                --per_gpu_pred_batch_size={batch_size}   \
                --output_dir {model} \
                --predict_dir {pred_path} \
                --n_process {n_process}')
        print(f'[INFO] finished prediction under: {pred_path}')
        print('[INFO] start process prediction result')
        _pred=np.load(f'{pred_path}/pred_results.npy')
        threshold=0.5
        _pred01=np.int8(_pred>=threshold)
        pas_bed_df['pred_value'] = _pred
        pas_bed_df['pred_01'] = _pred01
        pas_bed_df.to_csv(f'{pred_path}/prediction.tsv',index = False, sep = '\t')
        print(f'[INFO] prediction saved to {out_dir}/prediction.tsv')
        print('[INFO] finished process prediction result')
        #
    elif run=='retrain':
        # make prediction first and then finetune the model by relabeling the prediction
        print('[INFO] start retrain pipeline')
        print(f'[INFO] save files to retrain folder {out_dir}')
        #retrain_path=out_dir
        os.system(f'mkdir -p {out_dir}')
        initial_pred=f'{out_dir}/ini_pred'
        if os.path.isfile(f'{initial_pred}/dev.tsv'):
            print(f'[INFO] using the existing prediction from {initial_pred}')
        else:
            os.system(f'mkdir -p {initial_pred}')
            bed_seq_df.to_csv(f'{initial_pred}/dev.tsv',index=False,sep='\t')
            print('[INFO] start prediction')
            os.system(f'python {DNABERT_path}/run_finetune.py \
                    --model_type dna \
                    --tokenizer_name=dna{kmer} \
                    --model_name_or_path {init_model} \
                    --task_name dnaprom \
                    --do_predict \
                    --data_dir {initial_pred}  \
                    --max_seq_length {max_seq_len} \
                    --per_gpu_pred_batch_size={batch_size}   \
                    --output_dir {init_model} \
                    --predict_dir {initial_pred} \
                    --n_process {n_process}')
            print(f'[INFO] finished initial prediction under: {initial_pred}')
            print('[INFO] start process prediction result')
            _pred=np.load(f'{initial_pred}/pred_results.npy')
            threshold=0.5
            _pred01=np.int8(_pred>=threshold)
            pas_bed_df['pred_value'] = _pred
            pas_bed_df['pred_01'] = _pred01
            pas_bed_df.to_csv(f'{initial_pred}/prediction.tsv',index = False, sep = '\t')
            print(f'[INFO] prediction saved to {initial_pred}/prediction.tsv')
            print('[INFO] finished process prediction result')
        #
        print('[INFO] start relabel')
        data=pd.read_csv(f'{initial_pred}/dev.tsv',sep='\t')
        pred_df = pd.read_csv(f'{initial_pred}/prediction.tsv',sep = '\t')
        
        data['pred_value']=pred_df['pred_value']
        data['pred_01']=pred_df['pred_01']
        data.to_csv(f'{initial_pred}/pred.tsv',index=False,sep='\t')
        relabel_data = relabelPred(f'{initial_pred}/pred.tsv',motif_file,f'{initial_pred}/pred.relabeled.tsv')
        print('[INFO] finished relabel')
        ft_path = f'{out_dir}/finetune_model'
        os.system(f'mkdir -p {ft_path}')
        relabel_data=relabel_data.loc[:,['sequence','relabel','info']]
        relabel_data.columns=['sequence', 'label', 'info']
        relabel_data.to_csv(f'{ft_path}/train.tsv',sep='\t',index=False)
        #
        print('[INFO] finetune the model')
        os.system(f'python {DNABERT_path}/run_finetune.py \
            --model_type dna \
            --tokenizer_name=dna{kmer} \
            --model_name_or_path {init_model} \
            --task_name dnaprom \
            --do_train \
            --data_dir {ft_path} \
            --max_seq_length {max_seq_len} \
            --per_gpu_eval_batch_size={batch_size}   \
            --per_gpu_train_batch_size={batch_size}   \
            --learning_rate 2e-4 \
            --num_train_epochs {train_epoch} \
            --output_dir {ft_path} \
            --save_steps 4000 \
            --warmup_percent 0.1 \
            --hidden_dropout_prob 0.1 \
            --overwrite_output \
            --weight_decay 0.01 \
            --n_process {n_process}')
        print(f'[INFO] finish finetune model under: {ft_path}')
        

### main function 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "REST:RNA 3\'End Sequencing data analysis Toolkit", add_help = 'REST Availbale tools: call_cluster, filter_cluster')
    subparsers = parser.add_subparsers(help = 'sub-command help')
    
    # subfunction: call_cluster
    call_cluster = subparsers.add_parser('call_cluster', help='call cluster based on the alignment bam/bed file')
    call_cluster.add_argument('--input_file', required=True, help = 'input file in bam or bed format')
    call_cluster.add_argument('--input_format', default = 'bed', help = 'format of the input file, default: bed')
    call_cluster.add_argument('--strand', type = int, default = 2, help = 'strand of the sequencing data that generate the input file. 1: forward strand, 2: reverse strand, 0: strandless')
    call_cluster.add_argument('--output', default = 'tmp.cluster', help = 'output file of the identified cluster')
    call_cluster.add_argument('--output_dis', default = 'tmp.dis', help = 'distance of the output cluster')
    call_cluster.set_defaults(func = call_cluster_main)
    
    # subfunction: filter_cluster
    filter_cluster = subparsers.add_parser('filter_cluster', help = 'filter cluster and fine-tuned DNABERT models for true/false PAS')
    filter_cluster.add_argument('--run', default = 'pred', help = 'pred/train/retrain: only predict true/false for fintering, or also train a fine-tuned model')
    filter_cluster.add_argument('--input_file', required = True, help = 'input file of cluster in bed-like format')
    filter_cluster.add_argument('--rounds', type = int, default = 5, help = 'rounds of the training, eg: 5')
    filter_cluster.add_argument('--model', default='model/5-new-12w-0', help='eg: /active_r0_m5-0, default=model/5-new-12w-0')
    filter_cluster.add_argument('--kmer', type = int,default=5,help='eg: 5, default=5')
    filter_cluster.add_argument('--batch_size',type=int,default=100,help='eg: 100, default=100')
    filter_cluster.add_argument('--train_epoch',type=int,default=5,help='number of train epochs, eg: 5, default=5')
    filter_cluster.add_argument('--max_seq_len',type=int,default=200,help='eg: 200, default=200')
    filter_cluster.add_argument('--n_process',type=int,default=100,help='number of processors, eg: 100, default=100')
    filter_cluster.add_argument('--motif_file',default='mouse.PAS.motif',help='path to motif file, default=mouse.PAS.motif')
    filter_cluster.add_argument('--DNABERT_path',default=None, help='path to DNABERT script, default=DNABERT/examples')
    filter_cluster.add_argument('--out_dir',default='test_out', help='directory for output, default=test_out')
    filter_cluster.add_argument('--strand',type = int,default=1, help='strand of the cluster, 1: forward strand, 2; reverse strand, 0: strandless')
    filter_cluster.add_argument('--distance',default=50, help='distance threshold to define the overlap with annotation, default = 50')
    filter_cluster.add_argument('--annotation',default=None, help='PAS annotation to define true and false for training the model')
    filter_cluster.add_argument('--reference',default=None, help='the reference genome used to extract sequence flanking peaks of each cluster')
    filter_cluster.set_defaults(func=filter_cluster_main)
    
    # subfunction: merge_cluster
    merge_cluster = subparsers.add_parser('merge_cluster', help = 'merge the cluster passed PASBERT filtering from multiple samples')
    merge_cluster.add_argument('--file_list', required=True, help = 'input file list seprated by , and each file should be in bed-like format')
    merge_cluster.add_argument('--read', type = int, default = 5, help = 'threshold of supporting read for a eligiable cluster, this information should be stored as score in bed format file')
    merge_cluster.add_argument('--sam_num', type = int, default = 2, help = 'mininal number of detected samples required for keeping the cluster')
    merge_cluster.add_argument('--distance',type = int, default = 25, help = 'distance threshold for merging cluster')
    merge_cluster.add_argument('--with_header', default = True, help = 'whether the bed file has header')
    merge_cluster.add_argument('--output', default = 'test.merged.cluster', help = 'the output merged cluster')
    merge_cluster.add_argument('--sam_label',default = None, help = 'short label of each sample/file')
    merge_cluster.set_defaults(func = merge_cluster_main)
    
    # subfunction: annotate_PAS
    annotate_PAS = subparsers.add_parser('annotate_PAS', help = 'annotate PAS based on gene annotation gtf')
    annotate_PAS.add_argument('--input_pas', required = True, help = 'input file of putative PAS in bed-like format')
    annotate_PAS.add_argument('--gtf', required = True, help = 'gene annotation in gtf format, required gene_id, gene_name, gene_type attributes')
    annotate_PAS.add_argument('--with_header', default = None, help = 'whether the input file has header')
    annotate_PAS.add_argument('--output', default = 'test.anno', help = 'the output PAS with annotation in bed-like format')
    annotate_PAS.add_argument('--distance',default = 24, type = int, help = 'distance to define overlap')
    annotate_PAS.set_defaults(func = annotate_PAS_main)
    
    # subfunction: count_PAS: count the number of supporting read for each PAS from 3'end-seq data, more options need to rewrite the code 
    count_PAS = subparsers.add_parser('count_PAS', help = 'count reads supporting each PAS from bam file')
    count_PAS.add_argument('--bam', required = True, help = 'bam file of 3\'-end sequencing data or single-cell sequencing data')
    count_PAS.add_argument('--data_type', default = 'bulk', help = 'the type of bam data, either bulk or sc(scRNA-seq)')
    count_PAS.add_argument('--pas', required = True, help = 'file of PAS in bed format')
    count_PAS.add_argument('--with_header', default = None, help = 'whether the PAS file has header')
    count_PAS.add_argument('--win', type = int, default = 24, help = 'length of the window flaking the pas for counting')
    count_PAS.add_argument('--output', default = 'pas.count', help = 'the output PAS with read count')
    count_PAS.add_argument('--temp',default = 'temp.saf', help = 'temporary file for counting')
    count_PAS.add_argument('--strand',default = 2, help = 'strand of the sequencing data. 1: forward strand, 2: reverse strand, 3: strandless')
    count_PAS.add_argument('--end',default = 3,type = int, help = 'count based on read 3 or 5 end')
    count_PAS.add_argument('--multiple',default = False, help = 'count multi-mapping reads or not')
    count_PAS.add_argument('--thread',default = 1, help = 'number of the threads')
    count_PAS.set_defaults(func = count_PAS_main)
    
    # subfunction: count_gene: aggerate count of PAS for each gene  
    count_gene = subparsers.add_parser('count_gene', help = 'count reads for each gene based on PAS count and annotation')
    count_gene.add_argument('--pas_count', required = True, help = 'file of PAS count from count_PAS function')
    count_gene.add_argument('--pas_anno', required = True, help = 'file of PAS annotation from annotate_PAS function')
    count_gene.add_argument('--pas_type', default = 'last_exon,intron_anno,intron_unanno,ups_exon,ncRNA', help = 'types of PAS used for calculating')
    count_gene.add_argument('--output', default = 'pas.count', help = 'the output PAS with read count')
    count_gene.set_defaults(func = count_gene_main)
    
    #subfunction: calcu_PAS: calculate PAS usage (PAU) and weighted 3'UTR length index (WULI)
    calcu_PAS = subparsers.add_parser('calcu_PAS', help = 'calculate PAS usage (PAU) or the weighted 3\' UTR length index')
    calcu_PAS.add_argument('--pas_count', required = True, help = 'txt file with read counts of each PAS output from featureCounts')
    calcu_PAS.add_argument('--method', default = 'PAU', help = 'method for the calculation, either PAU or WULI')
    calcu_PAS.add_argument('--pos', default = 'start', help = 'position of the last exon to calculate WULI, either start or mid')
    calcu_PAS.add_argument('--pas_anno', required = True, help = 'file of PAS annotation in bed-like format')
    calcu_PAS.add_argument('--pas_type', default = 'last_exon,intron_anno,intron_unanno,ups_exon', help = 'types of PAS used for calculating')
    calcu_PAS.add_argument('--output', default = 'pau.txt', help = 'the output of PAU or WULI')
    calcu_PAS.set_defaults(func = calcu_PAS_main)
    
    #
    if len(sys.argv) < 2:
        sys.exit()
    args = parser.parse_args()
    #print(args)
    args.func(args)
