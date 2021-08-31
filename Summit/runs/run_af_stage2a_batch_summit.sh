#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage: $0 <target_file> <out_dir>"
    exit 1
fi

target_file=$1
out_dir=$2

base_dir=/gpfs/alpine/bip198/proj-shared/mcoletti/PSP

#seq_file=/gpfs/alpine/proj-shared/bip198/alphafold_test/casp14/seq/target.seq
#fea_dir=/gpfs/alpine/proj-shared/bip198/alphafold_test/casp14/af_reduced_db
#fea_dir=/gpfs/alpine/world-shared/bif135/alphafold_onsummit/alphafold_test/casp14/af_reduced_db
#fea_dir=${base_dir}/casp14/af_reduced_db
fea_dir=/gpfs/alpine/world-shared/bif135/desulfovibrio

preset=reduced_dbs

##########################################################################################

#af_dir=/app/alphafold
#af_dir=/gpfs/alpine/world-shared/bif135/alphafold_onsummit/alphafold
af_dir=${base_dir}/Summit/alphafold
#data_dir=/gpfs/alpine/world-shared/bip198/alphafold_databases
data_dir=${base_dir}/Summit/alphafold_databases
#data_dri=fea_dir=${base_dir}/desulfovibrio/afold_fea

##########################################################################################

## get the names of targets, not the actual sequence because it is already in 'features.pkl'
fasta_path=''
counter=0
while IFS= read -r line
do
  [[ "$line" = "\#*" ]] && continue
  IFS=' '; read -a args <<< "$line"
  target=${args[0]}
  if [ $counter -eq 0 ]; then
    fasta_path=$target.fas
  else
    fasta_path="$fasta_path,$target.fas"
  fi
  counter=$(($counter+1))
done < "$target_file"


##########################################################################################
cd $af_dir  || { echo 'Error: enter $af_dir failed' ; exit 1 ; }
python $af_dir/run_alphafold_stage2a.py --fasta_paths=$fasta_path --preset=$preset \
  --data_dir=$data_dir --output_dir=$out_dir --feature_dir=$fea_dir \
  --model_names=model_1_ptm,model_2_ptm,model_3_ptm,model_4_ptm,model_5_ptm \
  --benchmark

#--model_names=model_1,model_2,model_3,model_4,model_5 \
