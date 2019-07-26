#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/3_bp_exons_isoncorrect3/mutation_experiment"
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/mutation_experiment"
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

summary_file=$outbase/"summary.csv"
all_results_file=$outbase/"results.csv"
plot_file=$outbase/"summary"


echo -n  "id","type","Depth","p","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions"$'\n' > $summary_file
echo -n  "id","type","Depth","p","q_acc","r_acc", "total_errors","identity","subs","ins","del","switch","abundance","mutation_present","minor"$'\n' > $all_results_file

for depth in 20 50 100 # 20 50 # 50 # 100  #200 500 
do 
    for id in $(seq 1 1 10) 
    do
        # python $experiment_dir/simulate_reads.py --sim_genome_len 300 --coords 0 100 200 300 --outfolder $outbase/$depth/$id/ --probs 1.0 $p 1.0  --nr_reads $depth > /dev/null
        for p in $(seq 0.1 0.1 0.5) # $(seq 0.1 0.1 1.0)  # $(seq 0.1 0.1 0.2)
        do
            # python $experiment_dir/simulate_reads.py --isoforms $outbase/$depth/$id/isoforms.fa  --coords 0 100 200 300 --outfolder $outbase/$depth/$id/$p --probs 1.0 $p 1.0  --nr_reads $depth #> /dev/null

            # if [ "$1" == "exact" ]
            #     then
            #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null
            #     else
            #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80   &> /dev/null
            # fi

            # python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80  &> /dev/null            
            # python $experiment_dir/../family_experiment/evaluate_simulated_reads.py  $outbase/$depth/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/isoforms_abundance.fa $outbase/$depth/$id/$p/isoncorrect/evaluation --deal_with_ties > /dev/null
            echo -n  $id,approx,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/isoncorrect/evaluation/summary.csv 
            echo -n  $id,approx,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/isoncorrect/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/isoncorrect/evaluation/results.csv >> $all_results_file

            # python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null            
            python $experiment_dir/../family_experiment/evaluate_simulated_reads.py  $outbase/$depth/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/isoforms_abundance.fa $outbase/$depth/$id/$p/isoncorrect_exact/evaluation --deal_with_ties > /dev/null
            echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/summary.csv 
            echo -n  $id,exact,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/results.csv >> $all_results_file


            # fastq2fasta $outbase/$depth/$id/$p/reads.fq $outbase/$depth/$id/$p/reads.fa
            # python $experiment_dir/../family_experiment/evaluate_simulated_reads.py   $outbase/$depth/$id/$p/reads.fa $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/isoforms_abundance.fa $outbase/$depth/$id/$p/evaluation_reads --deal_with_ties > /dev/null
            echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/evaluation_reads/summary.csv 
            echo -n  $id,original,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/evaluation_reads/summary.csv  >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",original,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/evaluation_reads/results.csv >> $all_results_file

        done
    done
done


# echo  $experiment_dir/plot_mutation_data.py $results_file $plot_file
python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_tot.pdf" Total
python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_subs.pdf" Substitutions
python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_ind.pdf" Insertions
python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_del.pdf" Deletions

python $experiment_dir/plot_frac_switches.py $all_results_file $plot_file"_minor_mut_retained.pdf" $plot_file"_major_mut_retained.pdf"


