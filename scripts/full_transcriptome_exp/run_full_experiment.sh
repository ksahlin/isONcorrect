#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/chr6"
corrected_reads_fastq="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/chr6/isONcorrect.fastq" 
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/full_transcriptome_exp"
isoncorrect_dir="/users/kxs624/Documents/workspace/isONcorrect"
isonclust_dir="/users/kxs624/Documents/workspace/isONclust"

database="/Users/kxs624/Downloads/chr6_ensemble.fa"
# database="/Users/kxs624/Documents/data/simulated/ensemble_transcripts.fa" #/nfs/brubeck.bx.psu.edu/scratch6/ksahlin/data/ENSEMBLE_simulated/ensembl_coding_seqs.fa
mkdir -p $outbase


IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing

depth=$1 
error_rates_file=$outbase/"error_rates_"$depth".tsv"
# abundance_file_original=$outbase/"abundance_orig_"$depth".tsv"
corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
original_reads_mappings=$outbase/"original_"$depth".sam"




# echo -n  "id","type","Depth","mut","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches"$'\n' > $results_file
# echo -n  "id"$'\t'"Depth"$'\t'"mut"$'\t'"transcript_id"$'\t'"abundance_original"$'\t'"abundance_corrected"$'\n' > $results_file2


for id in $(seq 1 1 1)  
do 
#     for depth in 50 100 #200 #20 #20 50 100 #10 20 #50 # 100 200 500 1000 5000 10000
#     do
    
    # Simulate reads
    echo  python $experiment_dir/simulate_reads.py $database $outbase/$id/$depth/reads.fq $depth
    python $experiment_dir/simulate_reads.py $database $outbase/$id/$depth/reads.fq $depth #&> /dev/null
    ###############
    
    # #  Remove longest transcripts
    # python $experiment_dir/filter_fastq.py $outbase/$id/$depth/reads.fq 20000
    # mv $outbase/$id/$depth/reads.fq_filtered.fq $outbase/$id/$depth/reads.fq
    # ############################

    # Cluster with isONclust      
    rm -rf $outbase/$id/$depth/isonclust/sorted.fastq
    echo python $isonclust_dir/isONclust --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/ --k 13 --w 25 --q 6.0 --t 2
    python $isonclust_dir/isONclust --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/ --k 13 --w 25 --q 6.0 --t 2 #&> /dev/null            
    python $isonclust_dir/isONclust write_fastq --clusters $outbase/$id/$depth/isonclust/final_clusters.csv --fastq $outbase/$id/$depth/reads.fq --outfolder $outbase/$id/$depth/isonclust/fastq --N 1  #&> /dev/null            
    ############################

    # Correct reads with isONcorrect
    python $isoncorrect_dir/run_isoncorrect --keep_old --t 2 --fastq_folder $outbase/$id/$depth/isonclust/fastq  --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80  # &> /dev/null            
    ###############################

    # Merge individually corrected clusters    
    > $corrected_reads_fastq
    FILES=("$outbase/$id/$depth/isoncorrect"/*/corrected_reads.fastq)
    for c_id in $(seq 0 1 1000) $"${FILES[@]}"
        do 
            if [ -f "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq ]; then
               echo "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq
               cat "$outbase/$id/$depth/isoncorrect"/$c_id/corrected_reads.fastq >> $corrected_reads_fastq
            else
               echo "File $FILE does not exist."
            fi
        done
    #######################################
        
    # # Evaluate error rates
    error_plot=$outbase/$id/$depth/"error_rates_"$depth".pdf" 
    python $experiment_dir/get_error_rates.py  $database $outbase/$id/$depth/reads.fq  $corrected_reads_fastq  > $error_rates_file #$outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
    python $experiment_dir/plot_error_rates.py $error_rates_file  $error_plot
    # ######################

    # ## Map reads
    minimap2 -ax splice --eqx $database -uf -k14 $corrected_reads_fastq > $corrected_reads_mappings
    minimap2 -ax splice --eqx $database  -uf -k14 $outbase/$id/$depth/reads.fq > $original_reads_mappings
    # ###################

    # # Evaluate chnage in expression levels compared to true expression
    abundance_file=$outbase/$id/$depth/"abundance.tsv"
    > $abundance_file
    echo -n  "id","cov_aln","cov_true","seq","type"$'\n' > $abundance_file
    python $experiment_dir/get_abundance.py  $corrected_reads_mappings "corrected" >> $abundance_file
    python $experiment_dir/get_abundance.py $original_reads_mappings  "original" >>  $abundance_file
    abundance_plot=$outbase/$id/$depth/"abundance.pdf" 
    python $experiment_dir/plot_abundance.py $abundance_file "transcript" $abundance_plot



done




