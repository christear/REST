ROOT=$1
N_ROUND=$2
MER=$3

echo .....start pipeline.....
echo ROOT: $ROOT
echo N_ROUND: $N_ROUND
echo KMER: $MER

for i in `eval echo {0..$N_ROUND}`
do
    echo Round: $i
    echo [*1] start process data from positive/negative.events.txt
    python ./PASBERT/1.data_process.py --root $ROOT --round $i --mer $MER
    echo [1*] finished process train.tsv to active_r$i_m$MER-0
    
    echo [*2] start finetune model
    if [ $i -eq '0' ]
    then
        bash ./PASBERT/2.finetune.sh model/5-new-12w-0 $ROOT/active_r${i}_m$MER-0 $MER
    else
        bash ./PASBERT/2.finetune.sh $ROOT/active_r$((i-1))_m$MER-0 $ROOT/active_r${i}_m$MER-0 $MER
    fi
    echo [2*] finished finetune
    
    echo [*3] start prediction with finetuned model
    mkdir -p $ROOT/active_r${i}_m$MER-0/pred
    cp $ROOT/active_r${i}_m$MER-0/train.tsv $ROOT/active_r${i}_m$MER-0/pred/dev.tsv
    bash ./PASBERT/3.prediction.sh $ROOT/active_r${i}_m$MER-0 $ROOT/active_r${i}_m$MER-0 $MER
    echo [3*] finished predition
    
    echo [*4] start process prediction results
    python ./PASBERT/4.prediction_process.py --root $ROOT --round $i --mer $MER
    echo [4*] finish process prediction results
    
    echo [*5] start relabel
    bash ./PASBERT/5.relabel.sh $ROOT/active_r${i}_m${MER}-0/pred/dev_active_r${i}_m${MER}.tsv
    echo [5*] finish relabel
    
    echo [*6] start relabel to positive/negative.events.txt
    mkdir $ROOT/active_r$((i+1))_m${MER}-0
    cp $ROOT/active_r${i}_m${MER}-0/pred/dev_active_r${i}_m${MER}.relabeled.tsv $ROOT/active_r$((i+1))_m${MER}-0
    python ./PASBERT/6.relabel_posneg.py --root $ROOT --round $((i+1)) --mer $MER
    echo [6*] finish relabel to positive/negative.events.txt
    echo [***] Round: $i finished
done
