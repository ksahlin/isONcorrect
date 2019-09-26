#!/bin/bash

# RUN scripts e.g. as:   ./run_subsample_experiment.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/sirv_subsample/ /Users/sahlkris/Documents/data/ONT_SIRV/lc19_pcs109_subsample_full_length.fq  2

inbase=$1
outbase=$2
original_reads=$3
cores=$4

experiment_dir=$inbase"/scripts/full_transcriptome_exp"
eval_dir=$inbase"/scripts/"

mkdir -p $outbase


IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


error_rates_file=$outbase/"error_rates_"$depth".tsv"
corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
original_reads_mappings=$outbase/"original_"$depth".sam"


results_file=$outbase/"results_"$depth".csv"
plot_file=$outbase/"summary"

echo -n  "id","type","Depth","mut","q_acc","r_acc","total_errors","error_rate","subs","ins","del","switch","abundance"$'\n' > $results_file

# align original reads with minimap2 
original_reads_mapped=$outbase/original_reads.sam
minimap2 -a --eqx -k 14 -w 4 $inbase/test_data/sirv_transcriptome.fasta  $original_reads  > $original_reads_mapped

# Subsample reads
echo  python $experiment_dir/subsample_reads.py $original_reads $original_reads_mapped 
python $experiment_dir/subsample_reads.py $original_reads $original_reads_mapped
###############


# Correct reads with isONcorrect
python $inbase/run_isoncorrect --keep_old --t $cores --fastq_folder $outbase/fastq  --outfolder $outbase/isoncorrect/ --set_w_dynamically  # &> /dev/null            
###############################


# Evaluate indiviually corrected clusters    
for folder in $(find $outbase/$depth/isoncorrect/*  -type d); #c_id in $(seq 0 1 1000) $"${FILES[@]}"
    do 
        # echo $folder
        cl_id=$(basename $folder)
        echo $cl_id
        # if (( $cl_id > 50 ));
        #     then
        #         break
        # fi

        # if [ *$cl_id*  == "evaluation_original isoncorrect"  ]; then
        #   echo "Strings are equal" $cl_id
        #   continue
        # fi

        if [ -f $folder/corrected_reads.fastq ]; then
            echo $folder/corrected_reads.fastq
            python $eval_dir/evaluate_simulated_reads.py --no_ref_sim $folder/corrected_reads.fastq $database_filtered $folder/evaluation  #> /dev/null
            awk -F "," -v awk_id=$cl_id -v awk_depth=$depth '{if (NR!=1) {print awk_id",isoncorrect,"awk_depth",0,"$0}}'  $folder/evaluation/results.csv >> $results_file
            
            python $eval_dir/evaluate_simulated_reads.py --no_ref_sim $outbase/$depth/isonclust/fastq/$cl_id.fastq $database_filtered $folder/evaluation_reads  #> /dev/null
            awk -F "," -v awk_id=$cl_id -v awk_depth=$depth '{if (NR!=1) {print awk_id",original,"awk_depth",0,"$0}}'  $folder/evaluation_reads/results.csv >> $results_file

            # echo "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq
            # cat "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq >> $corrected_reads_fastq
        else
           echo $folder  "File  does not exist."
        fi
    done
#######################################

python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate
python $experiment_dir/plot_abundance_diff.py $results_file $plot_file"_abundance_diff.pdf" 

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






