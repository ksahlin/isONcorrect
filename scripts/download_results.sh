#!/bin/bash
#SBATCH -J dro_rattle
#SBATCH -t 14-00:00:00
#SBATCH --mail-user ksahlin@kth.se
#SBATCH --mail-type=ALL
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem 500G
#SBATCH --exclude=c[8-11,15-16]

BASE="/nfs/brubeck.bx.psu.edu/scratch4/ksahlin/isoncorrect_eval"
OUT="/galaxy/home/ksahlin/prefix/source/isONcorrect/data"
cp $BASE/evaluation_biological/SIRV_full_results_per_read_to_transcriptome.csv  $OUT/isONcorrect_sirv_full.csv
cp $BASE/evaluation_biological/drosophila/999/results_per_read.csv $OUT/isONcorrect_dros_full.csv
cp $BASE/evaluation_simulated/chr6/3500000/uneven/results.csv $OUT/isONcorrect_sim_full.csv

# sim subsample isoncorrect

cp $BASE/evaluation_simulated/chr6/1/uneven/results.csv $OUT/isONcorrect_sim_error_rate.csv
cp $BASE/evaluation_simulated/chr6/1/uneven/abundance.csv $OUT/isONcorrect_sim_overcorrection.csv

# sirv subsample isoncorrect
cp $BASE/sirv_subsampling.csv $OUT/isONcorrect_sirv_subsampling.csv


tar -cvjf $OUT/isoncorrect_eval_files.tar.bz2 $OUT/isONcorrect*
split -b 49M $OUT/isoncorrect_eval_files.tar.bz2 $OUT/isoncorrect_eval_files.tar.bz2.part

# https://www.tecmint.com/split-large-tar-into-multiple-files-of-certain-size/

# cat isoncorrect_eval_files.tar.bz2.part* > isoncorrect_eval_files.tar.bz2.joined

# RATTLE

BASE="/nfs/brubeck.bx.psu.edu/scratch4/ksahlin/"
OUT="/galaxy/home/ksahlin/prefix/source/isONcorrect/data"

# full

cp $BASE/rattle_sirv_full/results_per_read_to_transcriptome.csv $OUT/rattle_sirv_full.csv
cp $BASE/rattle_drosophila_full/results_per_read.csv $OUT/rattle_dros_full.csv

# sim subsample isoncorrect

cp $BASE/rattle_sim_controlled/results.csv $OUT/rattle_sim_error_rate.csv
cp $BASE/rattle_sim_controlled/abundance.csv $OUT/rattle_sim_overcorrection.csv
