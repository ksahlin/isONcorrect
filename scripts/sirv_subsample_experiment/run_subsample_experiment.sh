#!/bin/bash

# RUN scripts e.g. as:   ./run_subsample_experiment.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/sirv_subsample/  /Users/kxs624/Documents/data/ont/lc19_pcs109_subsample_full_length_pychopper2_phmmer.fq

inbase=$1
outbase=$2
original_reads=$3
#cores=$4

experiment_dir=$inbase"/scripts/sirv_subsample_experiment"
eval_dir=$inbase"/scripts/"

mkdir -p $outbase


IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


error_rates_file=$outbase/"error_rates_"$depth".tsv"
corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
original_reads_mappings=$outbase/"original_"$depth".sam"


results_file=$outbase/"results.csv"
summary_file=$outbase/"summary.csv"
plot_file=$outbase/"summary"

echo -n  "id","type","Depth","k","w","q_acc","r_acc","total_errors","error_rate","subs","ins","del"$'\n' > $results_file

# # align original reads with minimap2 
# original_reads_mapped=$outbase/original_reads.sam
# minimap2 -a --eqx -k 14 -w 1 $inbase/test_data/sirv_transcriptome.fasta  $original_reads  > $original_reads_mapped

# # Subsample reads
# echo  python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/fastq
# python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/fastq
# ###############


# Correct reads with isONcorrect
# python $inbase/run_isoncorrect --keep_old --t $cores --fastq_folder $outbase/fastq  --outfolder $outbase/isoncorrect/ --set_w_dynamically  # &> /dev/null            
###############################

for id in $(seq 1 1 10)    
do 
    # mkdir -p $outbase/$id/fastq
    # python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/$id/fastq > /dev/null
    sampled_transcripts=$outbase/$id/fastq/sampled_transcripts.fasta
    for depth in 3 5 10 20 
    do
        echo
        echo $id,$depth
        reads=$outbase/$id/fastq/$depth
        fastq2fasta $reads.fastq $outbase/$id/fastq/$depth.fasta #&> /dev/null


        original_eval_out=$outbase/$id/$depth/original/evaluation
        # mkdir -p $original_eval_out
        # minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $reads.fasta  > $reads.sam 2>/dev/null
        python $experiment_dir/get_error_rates.py $sampled_transcripts  $reads.sam $original_eval_out/results.csv 

        for k_param in 7 8 9
        do
            for window in 0 2 4 
            do
                echo $k_param,$window

                w_param=$(( $k_param + $window ))
                corrected_approximate=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/
                eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/evaluation
                # mkdir -p $eval_out
                # python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_approximate --k $k_param --w $w_param  &> /dev/null            
                # minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_approximate/corrected_reads.fastq  > $corrected_approximate/corrected_reads.sam 2>/dev/null
                python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_approximate/corrected_reads.sam $eval_out/results.csv 
                
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file                
                corrected_exact=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/
                eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/evaluation
                # mkdir -p $eval_out
                # python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_exact --k $k_param --w $w_param --exact  &> /dev/null            
                # minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_exact/corrected_reads.fastq  > $corrected_exact/corrected_reads.sam 2>/dev/null
                python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_exact/corrected_reads.sam  $eval_out/results.csv 
                
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window '{if (NR!=1) {print awk_id",original,"awk_depth","awk_k","awk_w","$0}}'  $original_eval_out/results.csv >> $results_file

            done
        done
    done
done

python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_subs.pdf" subs
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_ind.pdf" ins
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_del.pdf" del

