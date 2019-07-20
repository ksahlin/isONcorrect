#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/RBMY1A1"
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/family_experiment"
database="/Users/kxs624/Documents/data/NCBI_RNA_database/Y_EXONS/Y_EXONS.fa"
mkdir -p $outbase

exon_file=""
# if [ "$1" == "exact" ]
#     then
#         echo "Using exact mode"
#         results_file=$outbase/"results_exact.tsv"
#         plot_file=$outbase/"results_exact"

#     else
#         echo "Using approximate mode"
#         results_file=$outbase/"results_approximate.tsv"        
#         plot_file=$outbase/"results_approximate"
# fi


IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing

mut_rate=$1  # use 0.01 and 0.001
family_size=$2
abundance=$3  # exp or const
gene_member=$4
results_file=$outbase/"results_"$mut_rate"_"$abundance"_"$gene_member".tsv"
results_file2=$outbase/"abundance_"$mut_rate"_"$abundance"_"$gene_member".tsv"
plot_file=$outbase/"results_"$mut_rate"_"$abundance"_"$gene_member

echo -n  "id","type","Depth","mut","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches"$'\n' > $results_file
echo -n  "id"$'\t'"Depth"$'\t'"mut"$'\t'"transcript_id"$'\t'"abundance_original"$'\t'"abundance_corrected"$'\n' > $results_file2

# python $experiment_dir/get_exons.py $database $outbase

for id in $(seq 1 1 1)  
do 
    # python $experiment_dir/generate_transcripts.py --exon_file $outbase/$gene_member"_exons.fa"  $outbase/$id/biological_material.fa --gene_member $gene_member  --family_size $family_size --isoform_distribution exponential  --mutation_rate $mut_rate  &> /dev/null
    # python $experiment_dir/generate_abundance.py --transcript_file $outbase/$id/biological_material.fa --abundance $abundance $outbase/$id/biological_material_abundance.fa   #&> /dev/null

    for depth in 20 #100 #200 #20 #20 50 100 #10 20 #50 # 100 200 500 1000 5000 10000
    do
        # python $experiment_dir/generate_ont_reads.py $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/reads.fq $depth &> /dev/null

        # python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80  &> /dev/null            
        # python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
        # echo -n  $id,approx,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
        # echo -n  $id,approx,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file
        
        # for line in $(cat $outbase/$id/$depth/isoncorrect/evaluation/results2.tsv ); do
        #         # echo "tester: $line"
        #         echo -n  $id$'\t'$depth$'\t'$mut_rate$'\t' >> $results_file2 && echo $line >> $results_file2
        # done
        # # cat $outbase/$id/$depth/isoncorrect/evaluation/results2.tsv >> $results_file2


        # fastq2fasta $outbase/$id/$depth/reads.fq $outbase/$id/$depth/reads.fa
        # python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/reads.fa $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation_reads > /dev/null
        # echo -n  $id,original,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv 
        # echo -n  $id,original,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv  >> $results_file

        if [ "$depth" -gt 101 ];
        then 
            echo "Depth greater than 100, skipping exact";
            continue
        else
            echo "Depth less or equal to 100";
            python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
            python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation > /dev/null
            echo -n  $id,exact,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
            echo -n  $id,exact,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file
            # cat $outbase/$id/$depth/isoncorrect/evaluation/results2.tsv >> $results_file2
        fi;


        ########################################
        # # Evaluate error rates TMP TEST
        corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
        original_reads_mappings=$outbase/"original_"$depth".sam"
        error_rates_file=$outbase/"error_rates_"$depth".tsv"
        tmp_experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/full_transcriptome_exp"
        corrected_reads_fastq=$outbase/$id/$depth/"isoncorrect"/"corrected_reads.fastq"
        database=$outbase/$id/"biological_material.fa"
        error_plot=$outbase/$id/$depth/"error_rates_"$depth".pdf" 
        python $tmp_experiment_dir/get_error_rates.py  $database $outbase/$id/$depth/reads.fq  $corrected_reads_fastq  > $error_rates_file #$outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
        python $tmp_experiment_dir/plot_error_rates.py $error_rates_file  $error_plot
        # ######################

        # ## Map reads
        minimap2 -ax splice --eqx $database -uf -k14 $corrected_reads_fastq > $corrected_reads_mappings
        minimap2 -ax splice --eqx $database  -uf -k14 $outbase/$id/$depth/reads.fq > $original_reads_mappings
        # ###################

        # # Evaluate chnage in expression levels compared to true expression
        abundance_file=$outbase/$id/$depth/"abundance.tsv"
        > $abundance_file
        echo -n  "id","cov_aln","cov_true","seq","type"$'\n' > $abundance_file
        python $tmp_experiment_dir/get_abundance.py  $corrected_reads_mappings "corrected" >> $abundance_file
        python $tmp_experiment_dir/get_abundance.py $original_reads_mappings  "original" >>  $abundance_file
        abundance_plot=$outbase/$id/$depth/"abundance.pdf" 
        python $tmp_experiment_dir/plot_abundance.py $abundance_file "transcript" $abundance_plot
        ######################################################################


    done
done


echo  $experiment_dir/plot_error_rates.py $results_file $plot_file
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" Total
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_subs.pdf" Substitutions
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_ind.pdf" Insertions
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_del.pdf" Deletions
python $experiment_dir/plot_abundance_diff.py $results_file2 $plot_file"_abundance_diff.pdf" 


