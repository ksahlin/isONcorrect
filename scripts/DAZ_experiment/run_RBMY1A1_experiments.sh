#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/RBMY1A1"
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/DAZ_experiment"
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

mut_rate=$1
family_size=$2
abundance=$3
gene_member="RBMY1A1"
results_file=$outbase/"results_"$mut_rate"_"$abundance"_"$gene_member".tsv"
plot_file=$outbase/"results_"$mut_rate"_"$abundance"_"$gene_member

echo -n  "id","type","Depth","p","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches"$'\n' > $results_file

python $experiment_dir/get_exons.py $database $outbase

for id in $(seq 1 1 3)  
do 
    python $experiment_dir/generate_transcripts.py --exon_file $outbase/$gene_member"_exons.fa"  $outbase/$id/biological_material.fa --gene_member $gene_member  --family_size $family_size --isoform_distribution exponential  --mutation_rate $mut_rate  &> /dev/null
    python $experiment_dir/generate_abundance.py --transcript_file $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa --abundance $abundance  &> /dev/null

    for depth in 20 50 100 #10 20 #50 # 100 200 500
    do
        python $experiment_dir/generate_ont_reads.py $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/reads.fq $depth &> /dev/null

        # if [ "$1" == "exact" ]
        #     then
        #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null
        #     else
        #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80   &> /dev/null
        # fi

        python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80  &> /dev/null            
        python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation #&> /dev/null
        echo -n  $id,approx,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
        echo -n  $id,approx,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file


        python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$depth/reads.fq   --outfolder $outbase/$id/$depth/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
        python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/isoncorrect/corrected_reads.fastq  $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation > /dev/null
        echo -n  $id,exact,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv 
        echo -n  $id,exact,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation/results.tsv >> $results_file

        fastq2fasta $outbase/$id/$depth/reads.fq $outbase/$id/$depth/reads.fa
        python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/reads.fa $outbase/$id/biological_material.fa $outbase/$id/biological_material_abundance.fa $outbase/$id/$depth/isoncorrect/evaluation_reads > /dev/null
        echo -n  $id,original,$depth,$mut_rate,&& head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv 
        echo -n  $id,original,$depth,$mut_rate, >> $results_file && head -n 1 $outbase/$id/$depth/isoncorrect/evaluation_reads/results.tsv  >> $results_file
    done
done


echo  $experiment_dir/plot_error_rates.py $results_file $plot_file
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" Total
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_subs.pdf" Substitutions
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_ind.pdf" Insertions
python $experiment_dir/plot_error_rates.py $results_file $plot_file"_del.pdf" Deletions


