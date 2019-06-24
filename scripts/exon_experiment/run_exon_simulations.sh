#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/3_bp_exons_isoncorrect3/exon_experiment"
experiment_dir="/Users/kxs624/Documents/workspace/isONcorrect/scripts/exon_experiment"
mkdir -p $outbase

if [ "$1" == "exact" ]
    then
        echo "Using exact mode"
        results_file=$outbase/"results_exact.tsv"
        plot_file=$outbase/"results_exact.pdf"

    else
        echo "Using approximate mode"
        results_file=$outbase/"results_approximate.tsv"        
        plot_file=$outbase/"results_approximate.pdf"
fi



echo -n  "id","type","d","p","tot","err","subs","ins","del","rate","switches"$'\n' > $results_file
for id in $(seq 1 1 1)  
do 
    for depth in 10 #20 50 100 200 500
    do
        # which python
        python $experiment_dir/simulate_reads.py --sim_genome_len 300 --coords 0 100 200 300 --outfolder $outbase/$id/ --probs 1.0 $p 1.0  --nr_reads $depth > /dev/null
        for p in $(seq 0.1 0.1 0.2) # $(seq 0.1 0.1 1.0)  # $(seq 0.1 0.1 0.2)
        do
            python $experiment_dir/simulate_reads.py --isoforms $outbase/$id/isoforms.fa  --coords 0 100 200 300 --outfolder $outbase/$id/$p --probs 1.0 $p 1.0  --nr_reads $depth #> /dev/null
            if [ "$1" == "exact" ]
                then
                    # echo "Using exact mode"
                    python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80 --exact   &> /dev/null
                else
                    # echo "Using approximate mode"
                    python /users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ --k 7 --w 10 --xmax 80   &> /dev/null
            fi
            # python /Users/kxs624/Documents/workspace/isONcorrect/isONcorrect3 --fastq $outbase/$id/$p/isoncorrect/corrected_reads_parasail_1.fasta   --outfolder $outbase/$id/$p/isoncorrect/ > /dev/null
            python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$p/isoncorrect/corrected_reads.fastq  $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation > /dev/null
            echo -n  $id,corrected,$depth,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv 
            echo -n  $id,corrected,$depth,$p, >> $results_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv >> $results_file

            # python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$p/isoncorrect/spoa_ref.fa  $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation_spoa > /dev/null
            # echo -n  $id,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation_spoa/results.tsv 
            # echo -n  $id,$p, >> $results_spoa_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation_spoa/results.tsv >> $results_spoa_file

            fastq2fasta $outbase/$id/$p/reads.fq $outbase/$id/$p/reads.fa
            python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$p/reads.fa $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation_reads > /dev/null
            echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$id/$p/isoncorrect/evaluation_reads/results.tsv 
            echo -n  $id,original,$depth,$p, >> $results_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation_reads/results.tsv  >> $results_file
            # echo -n  $id,$p, >> $results_spoa_file && head -n 1 $outbase/$id/$p/isoncorrect/evaluation_spoa/results.tsv >> $results_spoa_file

            # head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv 
            # ( head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv  && $id && $p ) | cat
        done
    done
done
python $experiment_dir/plot_sim_data.py $results_file $plot_file
# python $experiment_dir/plot_sim_data.py $results_spoa_file $plot_spoa_file


# python simulate_reads.py --coords 0 50 100 150 200 250 300 --ref /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/ref_sim_6exons.fa --outfile /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/reads_0.5/100_reads.fq --probs 1.0 0.5 1.0 1.0 0.5 1.0 --nr_reads 100
# time ./isONcorrect --fastq /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/reads_0.1/100_reads.fastq  --outfolder /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/
# python evaluate_simulated_reads.py  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/corrected_reads_parasail_second1.fasta  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoforms.fa /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/evaluation