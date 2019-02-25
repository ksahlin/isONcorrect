#!/bin/bash

outbase="/Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/3_bp_exons"
for id in 1 #2 3 4 5 6 7 8 9 10
do
    # which python
    python /Users/kxs624/Documents/workspace/isONcorrect/scripts/simulate_reads.py --sim_genome_len 150 --coords 0 50 100 150 --outfolder $outbase/$id/ --probs 1.0 $p 1.0  --nr_reads 100
    for p in $(seq 0.1 0.1 1.0)
    do
        # echo $id $p
        python /Users/kxs624/Documents/workspace/isONcorrect/scripts/simulate_reads.py --ref $outbase/$id/reference.fa  --coords 0 50 100 150 --outfolder $outbase/$id/$p --probs 1.0 $p 1.0  --nr_reads 100 > /dev/null
        python /Users/kxs624/Documents/workspace/isONcorrect/isONcorrect --fastq $outbase/$id/$p/reads.fq   --outfolder $outbase/$id/$p/isoncorrect/ > /dev/null
        python /Users/kxs624/Documents/workspace/isONcorrect/scripts/evaluate_simulated_reads.py  $outbase/$id/$p/isoncorrect/corrected_reads_parasail_1.fasta  $outbase/$id/isoforms.fa $outbase/$id/$p/isoncorrect/evaluation > /dev/null
        echo -n  $id $p "    " && head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv
        # head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv 
        # ( head -n 1 $outbase/$id/$p/isoncorrect/evaluation/results.tsv  && $id && $p ) | cat
    done
done


# python simulate_reads.py --coords 0 50 100 150 200 250 300 --ref /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/ref_sim_6exons.fa --outfile /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/reads_0.5/100_reads.fq --probs 1.0 0.5 1.0 1.0 0.5 1.0 --nr_reads 100
# time ./isONcorrect --fastq /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/reads_0.1/100_reads.fastq  --outfolder /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/
# python evaluate_simulated_reads.py  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/corrected_reads_parasail_second1.fasta  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoforms.fa /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/6_50_bp_exons/isoncorrect_out/reads_0.1/evaluation