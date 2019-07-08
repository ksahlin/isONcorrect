#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/full"
corrected_reads_fastq="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/full/isONcorrect.fastq" 
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/full_transcriptome_exp"
isoncorrect_dir="/users/kxs624/Documents/workspace/isONcorrect"
isonclust_dir="/users/kxs624/Documents/workspace/isONclust"

database="/Users/kxs624/Documents/data/simulated/ensemble_transcripts.fa" #/nfs/brubeck.bx.psu.edu/scratch6/ksahlin/data/ENSEMBLE_simulated/ensembl_coding_seqs.fa
mkdir -p $outbase


IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing

depth=$1 
error_rates_file=$outbase/"error_rates_"$depth".tsv"
abundance_file=$outbase/"abundance_"$depth".tsv"
plot_file=$outbase/"results_"$depth
corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
original_reads_mappings=$outbase/"original_"$depth".sam"

# echo -n  "id","type","Depth","mut","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches"$'\n' > $results_file
# echo -n  "id"$'\t'"Depth"$'\t'"mut"$'\t'"transcript_id"$'\t'"abundance_original"$'\t'"abundance_corrected"$'\n' > $results_file2


for id in $(seq 1 1 1)  
do 
#     for depth in 50 100 #200 #20 #20 50 100 #10 20 #50 # 100 200 500 1000 5000 10000
#     do
    
    # # Simulate reads
    # echo  python $experiment_dir/simulate_reads.py $database $outbase/$id/$depth/reads.fq $depth
    # python $experiment_dir/simulate_reads.py $database $outbase/$id/$depth/reads.fq $depth #&> /dev/null
    # ###############
    
    # #  Remove longest transcripts
    # python $experiment_dir/filter_fastq.py $outbase/$id/$depth/reads.fq 20000
    # mv $outbase/$id/$depth/reads.fq_filtered.fq $outbase/$id/$depth/reads.fq
    # ############################

    # # Cluster with isONclust      
    # rm -rf $outbase/$id/$depth/isonclust/sorted.fastq
    # echo python $isonclust_dir/isONclust --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/ --k 13 --w 25 --q 6.0 --t 2
    # python $isonclust_dir/isONclust --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/ --k 13 --w 25 --q 6.0 --t 2 #&> /dev/null            
    # python $isonclust_dir/isONclust write_fastq --clusters $outbase/$id/$depth/isonclust/final_clusters.csv --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/fastq --N 20  #&> /dev/null            
    # ############################

    # # Correct reads with isONcorrect
    # python $isoncorrect_dir/run_isoncorrect --keep_old --t 2 --fastq_folder $outbase/$id/$depth/isonclust/fastq  --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80  # &> /dev/null            
    # ###############################

    # # Merge individually corrected clusters    
    # > $corrected_reads_fastq
    # FILES=("$outbase/$id/$depth/isoncorrect"/*/corrected_reads.fastq)
    # for c_id in $(seq 0 1 300) $"${FILES[@]}"
    #     do 
    #         if [ -f "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq ]; then
    #            echo "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq
    #            cat "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq >> $corrected_reads_fastq
    #         else
    #            echo "File $FILE does not exist."
    #         fi
    #     done
    # #######################################
        
    # # Evaluate error rates
    # python $experiment_dir/get_error_rates.py  $database $outbase/$id/$depth/reads.fq  $corrected_reads_fastq  > $error_rates_file #$outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
    # # python $experiment_dir/plot_error_rates.py  
    # ######################

    # ## Map reads
    # minimap2 -ax map-ont --eqx $database  $corrected_reads_fastq > $corrected_reads_mappings
    # minimap2 -ax map-ont --eqx $database  $outbase/$id/reads.fq > $original_reads_mappings
    # ###################

    # Evaluate chnage in expression levels compared to true expression
    > $abundance_file
    python $experiment_dir/get_abundance.py  $corrected_reads_mappings "corrected" > $abundance_file
    python $experiment_dir/get_abundance.py $original_reads_mappings  "original" >>  $abundance_file

    


    # echo -n  $id,approx,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
    # echo -n  $id,approx,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file
        
    #     for line in $(cat $outbase/$id/$depth/isoncorrect/evaluation/results2.tsv ); do
    #             echo -n  $id$'\t'$depth$'\t'$mut_rate$'\t' >> $results_file2 && echo $line >> $results_file2
    #     done


    #     fastq2fasta $outbase/$id/$depth/reads.fq $outbase/$id/$depth/reads.fa
    #     python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/reads.fa $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation_reads > /dev/null
    #     echo -n  $id,original,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv 
    #     echo -n  $id,original,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv  >> $results_file

    #     if [ "$depth" -gt 21 ];
    #     then 
    #         echo "Depth greater than 20, skipping exact";
    #         continue
    #     else
    #         echo "Depth less or equal to 20";
    #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
    #         python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation > /dev/null
    #         echo -n  $id,exact,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
    #         echo -n  $id,exact,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file
    #     fi;

done


# echo  $experiment_dir/plot_error_rates.py $results_file $plot_file
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" Total
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_subs.pdf" Substitutions
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_ind.pdf" Insertions
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_del.pdf" Deletions
# python $experiment_dir/plot_abundance_diff.py $results_file2 $plot_file"_abundance_diff.pdf" 


