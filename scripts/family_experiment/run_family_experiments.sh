#!/bin/bash

# RUN scripts e.g. as:   ./run_family_experiments.sh 0.01 4 exponential DAZ2 /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/ /Users/kxs624/Documents/data/NCBI_RNA_database/Y_EXONS/Y_EXONS.fa

mut_rate=$1  # use 0.01 and 0.001
family_size=$2
abundance=$3  # exp or const
gene_member=$4
inbase=$5
outroot=$6
database=$7 # "/Users/kxs624/Documents/data/NCBI_RNA_database/Y_EXONS/Y_EXONS.fa"


outbase=$outroot/$gene_member"_"$mut_rate"_"$abundance
experiment_dir=$inbase"/scripts/family_experiment"
eval_dir=$inbase"/scripts/"

mkdir -p $outbase

# if [ "$1" == "exact" ]
#     then
#         echo "Using exact mode"
#         results_file=$outbase/"results_exact.csv"
#         plot_file=$outbase/"results_exact"

#     else
#         echo "Using approximate mode"
#         results_file=$outbase/"results_approximate.csv"        
#         plot_file=$outbase/"results_approximate"
# fi


IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing


summary_file=$outbase/"summary_"$mut_rate"_"$abundance"_"$gene_member".csv"
results_file=$outbase/"results_"$mut_rate"_"$abundance"_"$gene_member".csv"
plot_file=$outbase/"summary_"$mut_rate"_"$abundance"_"$gene_member

echo -n  "id","type","Depth","mut","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions"$'\n' > $summary_file
echo -n  "id","type","Depth","mut","q_acc","r_acc","total_errors","error_rate","subs","ins","del","switch","abundance"$'\n' > $results_file


python $experiment_dir/get_exons.py $database $outbase/"exons" > /dev/null

for id in $(seq 1 1 3)  
do 
    python $experiment_dir/generate_transcripts.py --exon_file $outbase/exons/$gene_member"_exons.fa"  $outbase/$id/biological_material.fa --gene_member $gene_member  --family_size $family_size --isoform_distribution exponential  --mutation_rate $mut_rate  > /dev/null

    for depth in 10 20 #50 #100 #200 #20 #20 50 100 #10 20 #50 # 100 200 500 1000 5000 10000
    do
        python $experiment_dir/generate_abundance.py --transcript_file $outbase/$id/biological_material.fa --abundance $abundance $outbase/$id/$depth/biological_material_abundance.fa   > /dev/null
        python $experiment_dir/generate_ont_reads.py $outbase/$id/$depth/biological_material_abundance.fa $outbase/$id/$depth/reads.fq $depth > /dev/null

        python $inbase/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --T 0.1 --k 7 --w 10 --xmax 80  &> /dev/null            
        python $eval_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/$depth/isoncorrect/evaluation > /dev/null
        echo -n  $id,approx,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/summary.csv 
        echo -n  $id,approx,$depth,$mut_rate, >> $summary_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/summary.csv >> $summary_file
        awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$mut_rate  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/isoncorrect/evaluation/results.csv >> $results_file

        
        # for line in $(cat $outbase/$id/$depth/isoncorrect/evaluation/results2.csv ); do
        #         # echo "tester: $line"
        #         echo -n  $id$'\t'$depth$'\t'$mut_rate$'\t' >> $results_file && echo $line >> $results_file
        # done
        # # cat $outbase/$id/$depth/isoncorrect/evaluation/results2.csv >> $results_file


        fastq2fasta $outbase/$id/$depth/reads.fq $outbase/$id/$depth/reads.fa
        python $eval_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/reads.fa $outbase/$id/biological_material.fa $outbase/$id/$depth/evaluation_reads > /dev/null
        echo -n  $id,original,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/evaluation_reads/summary.csv 
        echo -n  $id,original,$depth,$mut_rate, >> $summary_file && head -n 1 $outbase/$id/$depth/evaluation_reads/summary.csv  >> $summary_file
        awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$mut_rate '{if (NR!=1) {print awk_id",original,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/evaluation_reads/results.csv >> $results_file


        if [ "$depth" -gt 101 ];
        then 
            echo "Depth greater than 100, skipping exact";
            continue
        else
            python $inbase/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect_exact/ --T 0.1 --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
            python $eval_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect_exact/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/$depth/isoncorrect_exact/evaluation > /dev/null
            echo -n  $id,exact,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect_exact/evaluation/summary.csv 
            echo -n  $id,exact,$depth,$mut_rate, >> $summary_file && head -n 1 $outbase/$id/$depth/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$mut_rate '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/isoncorrect_exact/evaluation/results.csv >> $results_file

        fi;


        # ########################################
        # # # Evaluate error rates TMP TEST
        # corrected_reads_mappings=$outbase/"corrected_"$depth".sam"
        # original_reads_mappings=$outbase/"original_"$depth".sam"
        # error_rates_file=$outbase/"error_rates_"$depth".csv"
        # tmp_experiment_dir="$inbase/scripts/full_transcriptome_exp"
        # corrected_reads_fastq=$outbase/$id/$depth/"isoncorrect"/"corrected_reads.fastq"
        # database=$outbase/$id/"biological_material.fa"
        # error_plot=$outbase/$id/$depth/"error_rates_"$depth".pdf" 
        # python $tmp_experiment_dir/get_error_rates.py  $database $outbase/$id/$depth/reads.fq  $corrected_reads_fastq  > $error_rates_file #$outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
        # python $tmp_experiment_dir/plot_error_rates.py $error_rates_file  $error_plot
        # # ######################

        # # ## Map reads
        # minimap2 -ax splice --eqx $database -uf -k14 $corrected_reads_fastq > $corrected_reads_mappings
        # minimap2 -ax splice --eqx $database  -uf -k14 $outbase/$id/$depth/reads.fq > $original_reads_mappings
        # # ###################

        # # # Evaluate chnage in expression levels compared to true expression
        # abundance_file=$outbase/$id/$depth/"abundance.csv"
        # > $abundance_file
        # echo -n  "id","cov_aln","cov_true","seq","type"$'\n' > $abundance_file
        # python $tmp_experiment_dir/get_abundance.py  $corrected_reads_mappings "corrected" >> $abundance_file
        # python $tmp_experiment_dir/get_abundance.py $original_reads_mappings  "original" >>  $abundance_file
        # abundance_plot=$outbase/$id/$depth/"abundance.pdf" 
        # python $tmp_experiment_dir/plot_abundance.py $abundance_file "transcript" $abundance_plot
        # ######################################################################


    done
done


# echo  $experiment_dir/plot_error_rates.py $summary_file $plot_file
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_subs.pdf" Substitutions
# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_ind.pdf" Insertions
# python $experiment_dir/plot_error_rates.py $summary_file $plot_file"_del.pdf" Deletions
python $experiment_dir/plot_abundance_diff.py $results_file $plot_file"_abundance_diff.pdf" 


