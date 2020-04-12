#!/bin/bash
#SBATCH -J rat_eval
#SBATCH -t 14-00:00:00
#SBATCH --mail-user ksahlin@kth.se
#SBATCH --mail-type=ALL
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem 500G
#SBATCH --exclude=c[8-11,15-16]


ROOT_IN="/galaxy/home/ksahlin/prefix/source/isONcorrect/"
BASE_OUT="/nfs/brubeck.bx.psu.edu/scratch4/ksahlin/isoncorrect_eval/"
plot_script_dir="/galaxy/home/ksahlin/prefix/source/isONcorrect/evaluation/"
eval_dir="/galaxy/home/ksahlin/prefix/source/isONcorrect/scripts/"
experiment_dir=$eval_dir/full_transcriptome_exp/
ref="/galaxy/home/ksahlin/prefix/source/isONcorrect/test_data/chr6_nonredundant.fa"


# SIM CONTROLLED READS 12% ERROR
sim_reads=$BASE_OUT"/data/chr6/1/uneven/reads_12.fq"
#python $experiment_dir/simulate_reads.py $ref $sim_reads 1 --uneven --subsampling_experiment --error_level 12

outfolder="/nfs/brubeck.bx.psu.edu/scratch4/ksahlin/isoncorrect_sim_error_rates/12"
isonclust="/galaxy/home/ksahlin/prefix/source/isONclust/"
#python $isonclust/isONclust --t 50 --k 11 --w 20  --q 4.0 --fastq $sim_reads  --outfolder $outfolder"/isonclust/"
#python $isonclust/isONclust write_fastq --clusters $outfolder"/isonclust/final_clusters.tsv" --fastq $sim_reads --outfolder $outfolder"/isonclust/fastq/" --N 1
#python $ROOT_IN/run_isoncorrect --k 8 --fastq_folder $outfolder"/isonclust/fastq/"  --outfolder $outfolder"/isoncorrect/" --set_w_dynamically --t 62
corrected_reads=$outfolder"/isoncorrect/isONcorrect.fq"
#touch $corrected_reads
FILES=$outfolder/isoncorrect/*/corrected_reads.fastq
# for f in $FILES 
# do 
#   cat $f >> $outfolder"/isoncorrect/isONcorrect.fq"
# done

# ERROR RATE ANALYSIS

error_rates_file=$outfolder/results.csv
#python $experiment_dir/get_error_rates.py $ref $sim_reads  $corrected_reads  > $error_rates_file

# OVERCORRECTION ANALYSIS

corrected_reads_aligned=$outfolder"/all_reads_after_correction_aligned.sam"
#/usr/bin/time -v  minimap2 --eqx -t 8 -a -k9 -w1 -f 0.000001  $ref $corrected_reads >  $corrected_reads_aligned

sim_reads_aligned=$outfolder"/sim_reads_aligned.sam" #$BASE_OUT"/data/chr6/1/uneven/reads.sam"
#/usr/bin/time -v  minimap2 --eqx -t 8 -a -k9 -w1 -f 0.000001  $ref $sim_reads >  $sim_reads_aligned
abundance_file=$outfolder/"abundance.csv" 
echo   "read_acc","aligned_to","transcript_abundance","is_tp","read_type","ed_btw_transcripts","ed_read_to_true","ed_read_to_aligned" > $abundance_file
python $experiment_dir/get_abundance.py $corrected_reads_aligned $ref $corrected_reads   'corrected' >> $abundance_file
python $experiment_dir/get_abundance.py $sim_reads_aligned $ref $sim_reads   'original' >> $abundance_file





# SIM CONTROLLED READS 4% ERROR

sim_reads=$BASE_OUT"/data/chr6/1/uneven/reads_4.fq"
#python $experiment_dir/simulate_reads.py $ref $sim_reads 1 --uneven --subsampling_experiment --error_level 4

outfolder="/nfs/brubeck.bx.psu.edu/scratch4/ksahlin/isoncorrect_sim_error_rates/4"
isonclust="/galaxy/home/ksahlin/prefix/source/isONclust/"
#python $isonclust/isONclust --t 50 --k 13 --w 20  --q 4.0 --fastq $sim_reads  --outfolder $outfolder"/isonclust/"
#python $isonclust/isONclust write_fastq --clusters $outfolder"/isonclust/final_clusters.tsv" --fastq $sim_reads --outfolder $outfolder"/isonclust/fastq/" --N 1
#python $ROOT_IN/run_isoncorrect  --fastq_folder $outfolder"/isonclust/fastq/"  --outfolder $outfolder"/isoncorrect/" --set_w_dynamically --t 62
corrected_reads=$outfolder"/isoncorrect/isONcorrect.fq"
#touch $corrected_reads
FILES=$outfolder/isoncorrect/*/corrected_reads.fastq
# for f in $FILES 
# do 
#  cat $f >> $outfolder"/isoncorrect/isONcorrect.fq"
# done

# ERROR RATE ANALYSIS

error_rates_file=$outfolder/results.csv
#python $experiment_dir/get_error_rates.py $ref $sim_reads  $corrected_reads  > $error_rates_file

# OVERCORRECTION ANALYSIS

corrected_reads_aligned=$outfolder"/all_reads_after_correction_aligned.sam"
#/usr/bin/time -v  minimap2 --eqx -t 8 -a -k9 -w1 -f 0.000001  $ref $corrected_reads >  $corrected_reads_aligned

sim_reads_aligned=$outfolder"/sim_reads_aligned.sam" #$BASE_OUT"/data/chr6/1/uneven/reads.sam"
#/usr/bin/time -v  minimap2 --eqx -t 8 -a -k9 -w1 -f 0.000001  $ref $sim_reads >  $sim_reads_aligned
abundance_file=$outfolder/"abundance.csv" 
echo   "read_acc","aligned_to","transcript_abundance","is_tp","read_type","ed_btw_transcripts","ed_read_to_true","ed_read_to_aligned" > $abundance_file
python $experiment_dir/get_abundance.py $corrected_reads_aligned $ref $corrected_reads   'corrected' >> $abundance_file
python $experiment_dir/get_abundance.py $sim_reads_aligned $ref $sim_reads   'original' >> $abundance_file

