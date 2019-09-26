#!/bin/bash

# RUN scripts e.g. as:   ./run_subsample_experiment.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/sirv_subsample/ /Users/sahlkris/Documents/data/ONT_SIRV/lc19_pcs109_subsample_full_length.fq 

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

echo -n  "id","type","Depth","mut","q_acc","r_acc","total_errors","error_rate","subs","ins","del","switch","abundance"$'\n' > $results_file

# # align original reads with minimap2 
original_reads_mapped=$outbase/original_reads.sam
# minimap2 -a --eqx -k 14 -w 4 $inbase/test_data/sirv_transcriptome.fasta  $original_reads  > $original_reads_mapped

# # Subsample reads
# echo  python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/fastq
# python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/fastq
# ###############


# Correct reads with isONcorrect
# python $inbase/run_isoncorrect --keep_old --t $cores --fastq_folder $outbase/fastq  --outfolder $outbase/isoncorrect/ --set_w_dynamically  # &> /dev/null            
###############################

for id in $(seq 1 1 1)    
do 
    mkdir -p $outbase/$id/fastq
    # python $experiment_dir/subsample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/$id/fastq
    sampled_transcripts=$outbase/$id/fastq/sampled_transcripts.fasta
    for depth in 3 #5 10 20 50 100 200 500
    do
        fastq_in=$outbase/$id/fastq/$depth.fastq
        # fastq2fasta -i $fastq_in -o $outbase/$id/fastq/$depth.fasta > /dev/null
        fasta_in=$outbase/$id/fastq/$depth.fasta

        original_eval_out=$outbase/$id/$depth/original/evaluation
        minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $fasta_in  > $fasta_in.sam
        python $eval_dir/evaluate_simulated_reads.py $fasta_in $sampled_transcripts  $original_eval_out --deal_with_ties > /dev/null

        for k_param in 7 # 8 9
        do
            for window in 0 #1 2 
            do
                w_param=$(( $k_param + $window ))
                corrected_approximate=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/
                eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/evaluation
                python $inbase/isONcorrect3 --fastq $fastq_in  --outfolder $corrected_approximate --k $k_param --w $w_param  &> /dev/null            
                minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_approximate/corrected_reads.fastq  > $corrected_approximate/corrected_reads.sam
                python $eval_dir/evaluate_simulated_reads.py  $corrected_approximate/corrected_reads.fastq  $sampled_transcripts $eval_out --deal_with_ties > /dev/null
                
                echo -n  $id,approx,$depth,$k_param,$w_param,&& head -n 1 $eval_out/summary.csv 
                echo -n  $id,approx,$depth,$k_param,$w_param, >> $summary_file && head -n 1 $eval_out/summary.csv >> $summary_file
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$w_param  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file

                
                corrected_exact=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/
                eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/evaluation
                python $inbase/isONcorrect3 --fastq $fastq_in  --outfolder $corrected_exact --k $k_param --w $w_param --exact   &> /dev/null            
                minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_exact/corrected_reads.fastq  > $corrected_exact/corrected_reads.sam
                python $eval_dir/evaluate_simulated_reads.py  $corrected_exact/corrected_reads.fastq  $sampled_transcripts  $eval_out --deal_with_ties > /dev/null
                
                echo -n  $id,exact,$depth,$k_param,$w_param,&& head -n 1 $eval_out/summary.csv 
                echo -n  $id,exact,$depth,$k_param,$w_param, >> $summary_file && head -n 1 $eval_out/summary.csv >> $summary_file
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$w_param '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file


                echo -n  $id,original,$depth,$k_param,$w_param,&& head -n 1 $original_eval_out/summary.csv
                echo -n  $id,original,$depth,$k_param,$w_param, >> $summary_file && head -n 1 $original_eval_out/summary.csv  >> $summary_file
                awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$w_param '{if (NR!=1) {print awk_id",original,"awk_depth","awk_k","awk_w","$0}}'  $original_eval_out/results.csv >> $results_file

            done
        done
    done
done


# # Evaluate indiviually corrected clusters    
# for folder in $(find $outbase/$depth/isoncorrect/*  -type d); #c_id in $(seq 0 1 1000) $"${FILES[@]}"
#     do 
#         # echo $folder
#         cl_id=$(basename $folder)
#         echo $cl_id
#         # if (( $cl_id > 50 ));
#         #     then
#         #         break
#         # fi

#         # if [ *$cl_id*  == "evaluation_original isoncorrect"  ]; then
#         #   echo "Strings are equal" $cl_id
#         #   continue
#         # fi

#         if [ -f $folder/corrected_reads.fastq ]; then
#             echo $folder/corrected_reads.fastq
#             python $eval_dir/evaluate_simulated_reads.py --no_ref_sim $folder/corrected_reads.fastq $database_filtered $folder/evaluation  #> /dev/null
#             awk -F "," -v awk_id=$cl_id -v awk_depth=$depth '{if (NR!=1) {print awk_id",isoncorrect,"awk_depth",0,"$0}}'  $folder/evaluation/results.csv >> $results_file
            
#             python $eval_dir/evaluate_simulated_reads.py --no_ref_sim $outbase/$depth/isonclust/fastq/$cl_id.fastq $database_filtered $folder/evaluation_reads  #> /dev/null
#             awk -F "," -v awk_id=$cl_id -v awk_depth=$depth '{if (NR!=1) {print awk_id",original,"awk_depth",0,"$0}}'  $folder/evaluation_reads/results.csv >> $results_file

#             # echo "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq
#             # cat "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq >> $corrected_reads_fastq
#         else
#            echo $folder  "File  does not exist."
#         fi
#     done
# #######################################

# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_abundance_diff.py $results_file $plot_file"_abundance_diff.pdf" 






# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_subs.pdf" Substitutions
# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_ind.pdf" Insertions
# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_del.pdf" Deletions


# # # Evaluate error rates OLD
# error_plot=$outbase/$id/$depth/"error_rates_"$depth".pdf" 
# python $experiment_dir/get_error_rates.py  $database_filtered $outbase/$id/$depth/reads.fq  $corrected_reads_fastq  > $error_rates_file #$outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
# python $experiment_dir/plot_error_rates.py $error_rates_file  $error_plot
# # ######################

# # ## Map reads
# minimap2 -ax splice --eqx $database_filtered -uf -k14 $corrected_reads_fastq > $corrected_reads_mappings
# minimap2 -ax splice --eqx $database_filtered  -uf -k14 $outbase/$id/$depth/reads.fq > $original_reads_mappings
# # ###################

# # # Evaluate chnage in expression levels compared to true expression
# abundance_file=$outbase/$id/$depth/"abundance.tsv"
# > $abundance_file
# echo -n  "id","cov_aln","cov_true","seq","type"$'\n' > $abundance_file
# python $experiment_dir/get_abundance.py  $corrected_reads_mappings "corrected" >> $abundance_file
# python $experiment_dir/get_abundance.py $original_reads_mappings  "original" >>  $abundance_file
# abundance_plot=$outbase/$id/$depth/"abundance.pdf" 
# python $experiment_dir/plot_abundance.py $abundance_file "transcript" $abundance_plot






