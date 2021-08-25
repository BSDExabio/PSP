#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage: $0 <fasta_file> (<template_date>)"
    exit 1
fi
fasta_path=$1

out_dir=$HOME/scratch/afold/output_ptm
if [ ! -z "$2" ]
then
  out_dir=$2
fi

. load_alphafold

af_dir=$HOME/data/tools/alphafold2
data_dir=$HOME/scratch/afold/data
fea_dir=$HOME/scratch/afold/output

cd $af_dir  || { echo 'Error: enter $af_dir failed' ; exit 1 ; }
python $af_dir/run_alphafold_stage2.py --fasta_paths=$fasta_path --preset=casp14 \
  --data_dir=$data_dir --output_dir=$out_dir --feature_dir=$fea_dir \
  --model_names=model_1_ptm,model_2_ptm,model_3_ptm,model_4_ptm,model_5_ptm \
  --benchmark
