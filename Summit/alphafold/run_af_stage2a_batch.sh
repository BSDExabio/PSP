#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage: $0 <fasta_file> <target_lst> <out_dir> (<template_date>)"
    exit 1
fi
fasta_file=$1
target_lst=$2
fea_dir=$3
out_dir=$4


##########################################################################################
. load_alphafold


af_dir=$HOME/data/tools/alphafold2
#DATA_DIR=$HOME/scratch/afold/data

##########################################################################################
### prepare to extract input sequence libraries on a local nvme/ssd disk
HOSTNAME=$(hostname)
USER=$(whoami)
DATA_ROOT=/scratch/$PBS_JOBID/afold
DATA_DIR=$DATA_ROOT/data
if [ ! -d $DATA_ROOT ]; then
  mkdir -p $DATA_ROOT || { echo "Error: Creating $DATA_ROOT failed on ${HOSTNAME}" ; exit 1; }
fi

DATA_TAR_BALL=$HOME/scratch/afold/afold_params_casp14.tar.gz
total_size=3600   ### 3.5GB
if [ ! -f $DATA_TAR_BALL ]; then
  echo "Error: input file $DATA_TAR_BALL not found" ; exit 1;
fi

cd $DATA_ROOT || { echo "Error: Entering $DATA_ROOT failed on ${HOSTNAME}" ; exit 1; }
DONE_FILE=${DATA_DIR}/extraction_done
if [ ! -f $DONE_FILE ]; then
  #rm -rf $DATA_ROOT/*
  echo "Info: copying input data files to $DATA_ROOT on ${HOSTNAME}"

  ram_free_space=$(df -m . |tail -n1 |awk '{print $4}')
  echo "Info: free space in MB: $ram_free_space on ${HOSTNAME}"
  echo "Info: need space in MB: $total_size on ${HOSTNAME}"

  if [ $total_size -lt $ram_free_space ]; then
     #tar -zxf $DATA_TAR_BALL
     tar -I pigz -xf $DATA_TAR_BALL
  else
     echo "Error: No enough space on $DATA_ROOT!"
     exit 1
  fi
  touch $DONE_FILE && echo "Info: input tarball successfully copied to ${HOSTNAME}"
fi
##########################################################################################

cd $PBS_O_WORKDIR || { echo 'Error: enter $PBS_O_WORKDIR failed' ; exit 1 ; }

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
done < "$target_lst"

##########################################################################################
cd $af_dir  || { echo 'Error: enter $af_dir failed' ; exit 1 ; }
python $af_dir/run_alphafold_stage2a.py --fasta_paths=$fasta_path --preset=reduced_dbs \
  --data_dir=$DATA_DIR --output_dir=$out_dir --feature_dir=$fea_dir \
  --model_names=model_1,model_2,model_3,model_4,model_5 \
  --benchmark
