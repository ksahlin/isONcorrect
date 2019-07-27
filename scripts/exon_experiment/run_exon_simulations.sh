#!/bin/bash


# RUN scripts e.g. as: ./run_exon_simulations.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/exon_experiment/

inbase=$1/
outroot=$2/


outbase=$outroot/
experiment_dir=$inbase"/scripts/exon_experiment"
eval_dir=$inbase"/scripts/"

mkdir -p $outbase

# if [ "$1" == "exact" ]
#     then
#         echo "Using exact mode"
#         summary_file=$outbase/"results_exact.csv"
#         plot_file=$outbase/"results_exact"

#     else
#         echo "Using approximate mode"
#         summary_file=$outbase/"results_approximate.csv"        
#         plot_file=$outbase/"results_approximate"
# fi

summary_file=$outbase/"summary.csv"
results_file=$outbase/"results.csv"
plot_file=$outbase/"summary"


echo -n  "id","type","Depth","p","tot","err","subs","ins","del","Total","Substitutions","Insertions","Deletions","switches"$'\n' > $summary_file
echo -n  "id","type","Depth","p","q_acc","r_acc","total_errors","error_rate","subs","ins","del","switch","abundance"$'\n' > $results_file


for id in $(seq 1 1 3)  
do 
    for depth in 20 #50 #10 20 50 # 100 200 500
    do
        python $experiment_dir/simulate_reads.py --sim_genome_len 300 --coords 0 100 200 300 --outfolder $outbase/$depth/$id/ --probs 1.0 $p 1.0  --nr_reads $depth > /dev/null
        for p in $(seq 0.1 0.1 0.2) # $(seq 0.1 0.1 1.0)  # $(seq 0.1 0.1 0.2)
        do
            python $experiment_dir/simulate_reads.py --isoforms $outbase/$depth/$id/isoforms.fa  --coords 0 100 200 300 --outfolder $outbase/$depth/$id/$p --probs 1.0 $p 1.0  --nr_reads $depth > /dev/null
            
            # if [ "$1" == "exact" ]
            #     then
            #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null
            #     else
            #         python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80   &> /dev/null
            # fi

            python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80  > /dev/null            
            python $eval_dir/evaluate_simulated_reads.py  $outbase/$depth/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/isoncorrect/evaluation > /dev/null
            echo -n  $id,approx,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/isoncorrect/evaluation/summary.csv 
            echo -n  $id,approx,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/isoncorrect/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/isoncorrect/evaluation/results.csv >> $results_file

            python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$depth/$id/$p/reads.fq   --outfolder $outbase/$depth/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   > /dev/null            
            python $eval_dir/evaluate_simulated_reads.py  $outbase/$depth/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/isoncorrect_exact/evaluation > /dev/null
            echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/summary.csv 
            echo -n  $id,exact,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/isoncorrect_exact/evaluation/results.csv >> $results_file

            fastq2fasta $outbase/$depth/$id/$p/reads.fq $outbase/$depth/$id/$p/reads.fa
            python $eval_dir/evaluate_simulated_reads.py   $outbase/$depth/$id/$p/reads.fa $outbase/$depth/$id/isoforms.fa $outbase/$depth/$id/$p/evaluation_reads  > /dev/null
            echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$depth/$id/$p/evaluation_reads/summary.csv 
            echo -n  $id,original,$depth,$p, >> $summary_file && head -n 1 $outbase/$depth/$id/$p/evaluation_reads/summary.csv  >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",original,"awk_depth","awk_p","$0}}'  $outbase/$depth/$id/$p/evaluation_reads/results.csv >> $results_file

        done
    done
done


# echo  $experiment_dir/plot_exon_data.py $summary_file $plot_file
python $experiment_dir/plot_exon_data.py $results_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_exon_data.py $summary_file $plot_file"_subs.pdf" Substitutions
# python $experiment_dir/plot_exon_data.py $summary_file $plot_file"_ind.pdf" Insertions
# python $experiment_dir/plot_exon_data.py $summary_file $plot_file"_del.pdf" Deletions

python $experiment_dir/plot_frac_switches.py $results_file $plot_file"_overcorrected.pdf" 

