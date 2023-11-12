# python=3.6
# Modified from Juexiao Zhou's paspert.py 
# 
# Data: Nov 12, 2023
# Usage: 
# @train python pasbert.py --run train --file data/H1.cluster.data.txt --rounds 5 --kmer 5 --batch_size 100 --train_epoch 5
# @pred python pasbert.py --run pred --file data/H1.cluster.data.txt --model data/H1.cluster.data.txt.tmp/active_r0_m5-0 --kmer 5 --batch_size 100

import os
import argparse
#from DNABERT.motif.motif_utils import seq2kmer
from utils.PAS_utils import seq2kmer
from utils.PAS_utils import relabelPred
import pandas as pd
import numpy as np

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='pasbert')
    parser.add_argument('--run',help='train/pred')
    parser.add_argument('--file',help='eg: data/H1.cluster.data.txt')
    parser.add_argument('--rounds',type=int,default=5,help='eg: 5, default=5')
    parser.add_argument('--model',default='model/5-new-12w-0',help='eg: /active_r0_m5-0, default=model/5-new-12w-0')
    parser.add_argument('--kmer',type=int,default=5,help='eg: 5, default=5')
    parser.add_argument('--batch_size',type=int,default=100,help='eg: 100, default=100')
    parser.add_argument('--train_epoch',type=int,default=5,help='number of train epochs, eg: 5, default=5')
    parser.add_argument('--max_seq_len',type=int,default=200,help='eg: 200, default=200')
    parser.add_argument('--n_process',type=int,default=100,help='number of processors, eg: 100, default=100')
    parser.add_argument('--motif_file',default='mouse.PAS.motif',help='path to motif file, default=mouse.PAS.motif')
    parser.add_argument('--DNABERT_path',default='DNABERT/examples', help='path to DNABERT script, default=DNABERT/examples')
    parser.add_argument('--out_dir',default='test_out', help='directory for output, default=test_out')
    
    args = parser.parse_args()
    run=args.run
    if run=='pred':
        model=args.model
    file=args.file
    rounds=args.rounds
    kmer=args.kmer
    batch_size=args.batch_size
    train_epoch=args.train_epoch
    max_seq_len=args.max_seq_len
    n_process=args.n_process
    motif_file=args.motif_file
    #
    DNABERT_path=args.DNABERT_path
    out_dir=args.out_dir
    
    if run=='train':

        #file='data/H1.cluster.data.txt.10'
        #rounds=2
        #kmer=5
        print('[INFO] start training pipeline')
        print('[INFO] file to train: {}'.format(file))
        print('[INFO] N rounds: {}'.format(rounds))
        print('[INFO] kmer: {}'.format(kmer))
        print('[INFO] start pipeline')

        os.system(f'mkdir -p {out_dir}')
        print(f'[INFO] creat output folder: {out_dir}')
        #os.system('mkdir -p {}.tmp'.format(file))
        #print('[INFO] creat tmp folder: {}.tmp'.format(file))

        for r in range(rounds):
            print('[INFO] **start training round {}**'.format(r))
            #output_path='{}.tmp/active_r{}_m{}-0'.format(file,r,kmer)
            output_path=f'{out_dir}/active_r{r}_m{kmer}-0'
            os.system('mkdir -p {}'.format(output_path))
            print('[INFO] save files in this round under: {}'.format(output_path))
            if r==0:
                os.system('cp {} {}/data.txt'.format(file, output_path))
                with open('{}/train.tsv'.format(output_path),'w') as w:
                    w.write('sequence\tlabel\tinfo\n')
                    with open('{}/data.txt'.format(output_path),'r') as f:
                        for line in f:
                            line=line.rstrip('\n').split('\t')
                            w.write('{}\t{}\t{}\n'.format(seq2kmer(line[1].upper(),kmer),line[2],line[0]))  
                print('[INFO] saved {}/train.tsv'.format(output_path))
                print('[INFO] finished data processing')
            elif r!=0:
                print('[INFO] train.tsv should exist, from the relabel of last round')

            print('[INFO] start finetune model')
            if kmer==3:
                init_model='model/3-new-12w-0'
            elif kmer==4:
                init_model='model/4-new-12w-0'
            elif kmer==5:
                init_model='model/5-new-12w-0'
            elif kmer==6:
                init_model='model/6-new-12w-0'
            os.system('python {}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{} \
                --model_name_or_path {} \
                --task_name dnaprom \
                --do_train \
                --data_dir {} \
                --max_seq_length {} \
                --per_gpu_eval_batch_size={}   \
                --per_gpu_train_batch_size={}   \
                --learning_rate 2e-4 \
                --num_train_epochs {} \
                --output_dir {} \
                --save_steps 4000 \
                --warmup_percent 0.1 \
                --hidden_dropout_prob 0.1 \
                --overwrite_output \
                --weight_decay 0.01 \
                --n_process {}'.format(DNABERT_path,kmer,init_model, output_path, max_seq_len, batch_size, batch_size, train_epoch, output_path, n_process))
            print('[INFO] finish finetune model under: {}'.format(output_path))
            print('[INFO] start prediction with finetuned model')
            pred_path='{}/pred'.format(output_path)
            os.system('mkdir -p {}'.format(pred_path))
            os.system('cp {}/train.tsv {}/dev.tsv'.format(output_path, pred_path))
            os.system('python {}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{} \
                --model_name_or_path {} \
                --task_name dnaprom \
                --do_predict \
                --data_dir {}  \
                --max_seq_length {} \
                --per_gpu_pred_batch_size={}   \
                --output_dir {} \
                --predict_dir {} \
                --n_process {}'.format(DNABERT_path,kmer,output_path,pred_path,max_seq_len ,batch_size, output_path,pred_path, n_process))
            print('[INFO] finished prediction under: {}'.format(pred_path))

            print('[INFO] start process prediction result')
            data=pd.read_csv('{}/dev.tsv'.format(pred_path),sep='\t')
            _pred=np.load('{}/pred_results.npy'.format(pred_path))
            threshold=0.5
            _pred01=np.int8(_pred>=threshold)
            data['pred_value']=_pred
            data['pred_01']=_pred01
            data.to_csv('{}/pred.tsv'.format(pred_path),index=False,sep='\t')
            print('[INFO] finished process prediction result')

            print('[INFO] start relabel')
            #pred_tsv='{}/pred.tsv'.format(pred_path)
            #os.system("grep -v \"^sequence\" " + pred_tsv + " | awk \'{seqs=\"\"; for(i = 1; i < 196; i++){seqs=seqs\"\"substr($i,0,1)} print \">\"NR\"\\n\"seqs$196}' > " + pred_tsv + ".fa")
            #os.system('python search4motif.py {}.fa {} 0 104 {}.pas.motif.pos'.format(pred_tsv,motif_file, pred_tsv))
            #os.system('Rscript relabel.PAS.v2.r {} {}'.format(pred_tsv, motif_file))
            relabel_data = relabelPred(f'{pred_path}/pred.tsv',motif_file,f'{pred_path}/pred.relabeled.tsv')
            print('[INFO] finished relabel')

            print('[INFO] start process relabeled results to new train data for next round')
            output_path='{}.tmp/active_r{}_m{}-0'.format(file,r+1,kmer)
            os.system('mkdir -p {}'.format(output_path))
            #relabel_data=pd.read_csv('{}/pred.relabeled.tsv'.format(pred_path),sep='\t')
            relabel_data=relabel_data.loc[:,['sequence','relabel','info']]
            relabel_data.columns=['sequence', 'label', 'info']
            relabel_data.to_csv('{}/train.tsv'.format(output_path),sep='\t',index=False)
            print('[INFO] finish process')

            print('[INFO] **finished training round {}**'.format(r))

        print('[INFO] start merge all results into')
        hold_all=pd.read_csv(file,sep='\t',header=None)
        hold_all.columns=['info','rawseq','initlabel']
        for r in range(rounds):
            #round_df=pd.read_csv('{}.tmp/active_r{}_m{}-0/pred/pred.relabeled.tsv'.format(file,r,kmer),sep='\t')
            round_df = pd.read_csv(f'{out_dir}/active_r{r}_m{kmer}/pred/pred.relabeled.tsv',sep = '\t')
            round_dt=round_df.loc[:,['info','label','pred_value','pred_01','relabel']]
            round_dt.columns=['info','inputlabel_r{}'.format(r),'predvalue_r{}'.format(r),'pred01_r{}'.format(r),'outputlabel_r{}'.format(r)]
            hold_all=pd.merge(hold_all, round_dt, how='outer',on='info')
        #hold_all.to_csv('{}.summary'.format(file),sep='\t',index=False)
        hold_all.to_csv(f'{out_dir}/train.summary',sep='\t',index=False)
        
    elif run=='pred':
        
        print('[INFO] start pred pipeline')
        #print('[INFO] save files to prediction folder {}.pred.tmp'.format(file))
        print(f'[INFO] save files to prediction folder {our_dir}')
        #pred_path='{}.pred.tmp'.format(file)
        pred_path=out_dir
        os.system('mkdir -p {}'.format(pred_path))
        
        with open('{}/dev.tsv'.format(pred_path),'w') as w:
            w.write('sequence\tlabel\tinfo\n')
            with open('{}'.format(file),'r') as f:
                for line in f:
                    line=line.rstrip('\n').split('\t')
                    w.write('{}\t{}\t{}\n'.format(seq2kmer(line[1].upper(),kmer),line[2],line[0]))  
        
        print('[INFO] start prediction')
        os.system('python {}/run_finetune.py \
                --model_type dna \
                --tokenizer_name=dna{} \
                --model_name_or_path {} \
                --task_name dnaprom \
                --do_predict \
                --data_dir {}  \
                --max_seq_length {} \
                --per_gpu_pred_batch_size={}   \
                --output_dir {} \
                --predict_dir {} \
                --n_process {}'.format(DNABERT_path,kmer,model, pred_path, max_seq_len,batch_size,model, pred_path, n_process))
        print('[INFO] finished prediction under: {}'.format(pred_path))
        
        print('[INFO] start process prediction result')
        hold_all=pd.read_csv(file,sep='\t',header=None)
        columns=list(hold_all.columns).copy()
        columns[0]='info'
        columns[1]='rowseq'
        hold_all.columns=columns
        hold_all=hold_all.loc[:,['info','rowseq']]
        print('[INFO] start process prediction result')
        data=pd.read_csv('{}/dev.tsv'.format(pred_path),sep='\t')
        _pred=np.load('{}/pred_results.npy'.format(pred_path))
        threshold=0.5
        _pred01=np.int8(_pred>=threshold)
        data['pred_value']=_pred
        data['pred_01']=_pred01
        data=data.loc[:,['info','pred_value','pred_01']]
        hold_all=pd.merge(hold_all,data,how='outer',on='info')
        #hold_all.to_csv('{}.pred'.format(file),index=False,sep='\t')
        hold_all.to_csv(f'{pred_path}/prediction.tsv',index = False, sep = '\t')
        #print('[INFO] prediction saved to {}.pred'.format(file))
        print(f'[INFO] prediction saved to {out_dir}/prediction.tsv')
        print('[INFO] finished process prediction result')