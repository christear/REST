# $1: model path
# $2: data path
# $3: kmer

export MODEL_PATH=$1
export ROOT=$2
export KMER=$3
export DATA_PATH=$ROOT
export OUTPUT_PATH=$ROOT

python DNABERT/examples/run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_train \
    --data_dir $DATA_PATH \
    --max_seq_length 200 \
    --per_gpu_eval_batch_size=100   \
    --per_gpu_train_batch_size=100   \
    --learning_rate 2e-4 \
    --num_train_epochs 3 \
    --output_dir $OUTPUT_PATH \
    --save_steps 4000 \
    --warmup_percent 0.1 \
    --hidden_dropout_prob 0.1 \
    --overwrite_output \
    --weight_decay 0.01 \
    --n_process 100
