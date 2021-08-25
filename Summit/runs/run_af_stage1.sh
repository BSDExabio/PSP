#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Usage: $0 <fasta_file> (<template_date>)"
    exit 1
fi
fasta_path=$1

max_template_date=2020-05-14
if [ ! -z "$2" ]
then
  max_template_date=$2
fi


. load_alphafold
export HHLIB=$HOME/data/tools/hh-suite/build
export HMMER=/usr/local/pace-apps/spack/packages/0.13/linux-rhel7-cascadelake/intel-19.0.5/hmmer-3.2.1-sngcwm2qjzzxseh42cryf432role4on5

af_dir=$HOME/data/tools/alphafold2
data_dir=$HOME/scratch/afold/data
out_dir=$HOME/scratch/afold/output

#cd $af_dir  || { echo 'Error: enter $af_dir failed' ; exit 1 ; }
python $af_dir/run_alphafold_stage1.py --fasta_paths=$fasta_path --preset=casp14 \
  --data_dir=$data_dir --output_dir=$out_dir      \
  --uniref90_database_path=$data_dir/uniref90/uniref90.fasta \
  --mgnify_database_path=$data_dir/mgnify/mgy_clusters.fa   \
  --uniclust30_database_path=$data_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
  --bfd_database_path=$data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
  --pdb70_database_path=$data_dir/pdb70/pdb70           \
  --template_mmcif_dir=$data_dir/pdb_mmcif/mmcif_files  \
  --max_template_date=2020-05-14                        \
  --obsolete_pdbs_path=$data_dir/pdb_mmcif/obsolete.dat \
  --hhblits_binary_path=$HHLIB/bin/hhblits   \
  --hhsearch_binary_path=$HHLIB/bin/hhsearch \
  --jackhmmer_binary_path=$HMMER/bin/jackhmmer \
  --kalign_binary_path=$HOME/data/tools/kalign_v2/kalign \
  --benchmark
