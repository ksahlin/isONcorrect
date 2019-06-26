#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/3_bp_exons_isoncorrect3/mutation_experiment"
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/mutation_experiment"
mkdir -p $outbase

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

results_file=$outbase/"results.tsv"
plot_file=$outbase/"results"


# echo -n  "id","type","Depth","p","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches_to_major","switches_to_minor","Minor_Retained","Major_Retained"$'\n' > $results_file
# for id in $(seq 1 1 10)  
# do 
#     for depth in 100 # 20 50 # 50 # 100  #200 500
#     do
#         # which python
#         python $experiment_dir/simulate_reads.py --sim_genome_len 300 --coords 0 100 200 300 --outfolder $outbase/$id/ --probs 1.0 $p 1.0  --nr_reads $depth > /dev/null
#         for p in $(seq 0.1 0.1 0.5) # $(seq 0.1 0.1 1.0)  # $(seq 0.1 0.1 0.2)
#         do
#             python $experiment_dir/simulate_reads.py --isoforms $outbase/$id/isoforms.fa  --coords 0 100 200 300 --outfolder $outbase/$id/$p --probs 1.0 $p 1.0  --nr_reads $depth #> /dev/null
#             # if [ "$1" == "exact" ]
#             #     then
#             #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null
#             #     else
#             #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80   &> /dev/null
#             # fi

#             python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80  &> /dev/null            
#             python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation > /dev/null
#             echo -n  $id,approx,$depth,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv 
#             echo -n  $id,approx,$depth,$p, >> $results_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv >> $results_file


#             python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
#             python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation > /dev/null
#             echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv 
#             echo -n  $id,exact,$depth,$p, >> $results_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv >> $results_file

#             fastq2fasta $outbase/$id/$p/reads.fq $outbase/$id/$p/reads.fa
#             python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$p/reads.fa $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation_reads > /dev/null
#             echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation_reads/results.tsv 
#             echo -n  $id,original,$depth,$p, >> $results_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation_reads/results.tsv  >> $results_file
#         done
#     done
# done


# echo  $experiment_dir/plot_mutation_data.py $results_file $plot_file
# python $experiment_dir/plot_mutation_data.py $results_file $plot_file"_tot.pdf" Total
# python $experiment_dir/plot_mutation_data.py $results_file $plot_file"_subs.pdf" Substitutions
# python $experiment_dir/plot_mutation_data.py $results_file $plot_file"_ind.pdf" Insertions
# python $experiment_dir/plot_mutation_data.py $results_file $plot_file"_del.pdf" Deletions

python $experiment_dir/plot_frac_switches.py $results_file $plot_file"_minor_mut_retained.pdf" $plot_file"_major_mut_retained.pdf"


