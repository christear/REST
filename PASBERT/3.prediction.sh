# $1: model path
# $2: data path
# $3: kmer

export MODEL_PATH=$1
export ROOT=$2
export KMER=$3
export DATA_PATH=$ROOT/pred
export PREDICTION_PATH=$ROOT/pred

python DNABERT/examples/run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_predict \
    --data_dir $DATA_PATH  \
    --max_seq_length 200 \
    --per_gpu_pred_batch_size=1000   \
    --output_dir $MODEL_PATH \
    --predict_dir $PREDICTION_PATH \
    --n_process 100
