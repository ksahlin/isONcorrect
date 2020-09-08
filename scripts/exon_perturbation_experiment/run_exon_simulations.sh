#!/bin/bash

# RUN scripts e.g. as: ./run_exon_simulations.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/exon_perturbation_experiment/

inbase=$1/
outroot=$2/


outbase=$outroot/
experiment_dir=$inbase"/scripts/exon_perturbation_experiment"
# eval_dir=$inbase"/scripts/"


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
results_file=$outbase/"results.csv"
plot_file=$outbase/"summary"


echo -n  "id","type","Depth","p","tot","err","Total","mut_retained"$'\n' > $summary_file
echo -n  "id","type","Depth","p","q_acc","r_acc","total_errors","error_rate","abundance","preserved","minor"$'\n' > $results_file


RTN4IP1="ATGCTATTGGTCTTCGAATGCCGGAGTGGGAGCTGTTGAACGGTTCTAAGTAGGAGAGTGAAAGACCAGCATGGAGGAGAAATCCGGCCAGTTATTGAACAAACCTTTCCTTTTTCTAAAGTTCCAGAAGCCTTCCTGAAGGTGGAAAGAGGACACGCACGAGGAAAGACTGTAATTAATGTTGTTTAAATAAAAATGCAGTTTAGTG"
# MAP3K4="CTTCCTGGAACGCAGAGTGTAATCCCTGCATTTGGCAAGAGGGCACATTGGGAGATTCAGCTTGCAAGAGTCCTGAATCTGATCTAGAAGACTTCTCCGATGAAACAAATACAGAGAATCTTTATGGTACCTCTCCCCCCAGCACACCTCGACAGATGAAACGCATGTCAACCAAACATCAGAGGAATAATGTGGGGAGGCCAGCCAGTCGGTCTAATTTGAAAG"
for id in $(seq 1 1 10) #  $(seq 1 1 10)
do 
    # sim random mutation
    python $experiment_dir/simulate_reads.py --sim_exon_removal --exon_length 20 --ref $RTN4IP1 --outfolder $outbase/$id/  > /dev/null

    for depth in 10 20 40 60 80 100
    do
        for p in 0.1 0.2 0.3 0.4 0.5
        do
            python $experiment_dir/simulate_reads.py --isoforms $outbase/$id/isoforms.fa --outfolder $outbase/$id/$depth/$p --probs $p  --nr_reads $depth > /dev/null

            # python $inbase/isONcorrect --fastq $outbase/$id/$depth/$p/reads.fq   --outfolder $outbase/$id/$depth/$p/isoncorrect/  &> /dev/null            
            # python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/$p/isoncorrect/corrected_reads.fastq  $outbase/$id/isoforms.fa  $outbase/$id/$depth/$p/isoncorrect/evaluation > /dev/null
            # echo -n  $id,approx,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/isoncorrect/evaluation/summary.csv 
            # echo -n  $id,approx,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/isoncorrect/evaluation/summary.csv >> $summary_file
            # awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/isoncorrect/evaluation/results.csv >> $results_file

            python $inbase/isONcorrect --fastq $outbase/$id/$depth/$p/reads.fq   --outfolder $outbase/$id/$depth/$p/isoncorrect_exact/ &> /dev/null            
            python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/$p/isoncorrect_exact/corrected_reads.fastq  $outbase/$id/isoforms.fa  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation  > /dev/null
            echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv 
            echo -n  $id,exact,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/results.csv >> $results_file


            # fastq2fasta $outbase/$id/$depth/$p/reads.fq $outbase/$id/$depth/$p/reads.fa
            # python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/$p/reads.fa $outbase/$id/isoforms.fa  $outbase/$id/$depth/$p/evaluation_reads > /dev/null
            # echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/evaluation_reads/summary.csv 
            # echo -n  $id,original,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/evaluation_reads/summary.csv  >> $summary_file
            # awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",original,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/evaluation_reads/results.csv >> $results_file

        done
    done
done


#python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_tot.pdf" error_rate
python $experiment_dir/plot_frac_switches.py $results_file $plot_file"_minor_isoform_retained.pdf" 


