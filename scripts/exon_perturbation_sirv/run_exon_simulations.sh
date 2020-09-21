#!/bin/bash

# RUN scripts e.g. as: ./run_exon_simulations.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/exon_perturbation_sirv/

inbase=$1/
outroot=$2/


outbase=$outroot/
experiment_dir=$inbase"/scripts/exon_perturbation_sirv"
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


# all sirv fastq ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.fq 
# all sirv alignments: ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam
# take only relevant alignnments:
# grep -E "0\tSIRV616" ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam > $outbase/"SIRV616.sam"
# grep -E "0\tSIRV606" ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam > $outbase/"SIRV606.sam"
# grep -E "0\tSIRV511" ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam > $outbase/"SIRV511.sam"
# grep -E "0\tSIRV506" ~/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam > $outbase/"SIRV506.sam"
# samtools fastq $outbase/"SIRV616.sam" > $outbase/"SIRV616.fastq"
# samtools fastq $outbase/"SIRV606.sam" > $outbase/"SIRV606.fastq"
# samtools fastq $outbase/"SIRV511.sam" > $outbase/"SIRV511.fastq"
# samtools fastq $outbase/"SIRV506.sam" > $outbase/"SIRV506.fastq"

# alignments=""
# fastq=""
exon_length=6 # either 6 or 14

if [ $exon_length == 6 ]
then
    >> $outbase/isoforms_6bp.fa
    isoform_506="GCATACTGCAAGACTTGTAGGCCCATGAATTACCCGCTTAGGCAATGGGTAGACCTTTCCTTGCCGCGTGGAACAAACCGCCATGAGTGTTTTACGTTTACGGTCGGCGACCGGAGTCACAGCGCGACCAACGGGGCCAGGCAAAGGCCCATAGGCTTTTCCATCGGCACGCGACCGATTCCATCAAGTATAGCCACGGTCGTAACTTCGACGGCGGCAATAGTTTATGGCGCTAACACCTTTCCATCGCTGAGTCCAAACGTATGACGTGAGGTGTTAGGTGCAATCCCCATTCGGGACGGACACGCGAGTAAAGGTCTAACATGTTTCAACTTTCCAGTTAGTCGCCGTCCTATATTGGCAAGACCGGATGGCCGTCAACCAGTTCCATTCTGACTATTTGCCACATTGATTAAAGAAGTAGGAGCGTTAGTCGGCGCGTGCGCGGGACTATTTAGCGAAGGTCGCCGTATCGCGGGTGGCGGCGTTTCCGCCATAACAACAAGCGTTTAAGGTCAGAGTCCGCGACGAATGAAACGCCGTGCAAGCGGC"  
    isoform_511="GCATACTGCAAGACTTGTAGGCCCATGAATTACCCGCTTAGGCAATGGGTAGACCTTTCCTTGCCGCGTGGAACAAACCGCCATGAGTGTTTTACGTTTACGGTCGGCGACCGGAGTCACAGCGCGACCAACGGGCAAAGGCCCATAGGCTTTTCCATCGGCACGCGACCGATTCCATCAAGTATAGCCACGGTCGTAACTTCGACGGCGGCAATAGTTTATGGCGCTAACACCTTTCCATCGCTGAGTCCAAACGTATGACGTGAGGTGTTAGGTGCAATCCCCATTCGGGACGGACACGCGAGTAAAGGTCTAACATGTTTCAACTTTCCAGTTAGTCGCCGTCCTATATTGGCAAGACCGGATGGCCGTCAACCAGTTCCATTCTGACTATTTGCCACATTGATTAAAGAAGTAGGAGCGTTAGTCGGCGCGTGCGCGGGACTATTTAGCGAAGGTCGCCGTATCGCGGGTGGCGGCGTTTCCGCCATAACAACAAGCGTTTAAGGTCAGAGTCCGCGACGAATGAAACGCCGTGCAAGCGGC"  
    echo ">isoform_511" >> $outbase/isoforms_6bp.fa
    echo $isoform_511 >> $outbase/isoforms_6bp.fa
    echo ">isoform_506" >> $outbase/isoforms_6bp.fa
    echo $isoform_506 >> $outbase/isoforms_6bp.fa
    isoforms=$outbase/isoforms_6bp.fa
    cat $outbase/"SIRV506.sam" $outbase/"SIRV511.sam" > $outbase/"reads_to_isoforms_6bp.sam"
    alignments=$outbase/"reads_to_isoforms_6bp.sam"

fi

if [ $exon_length == 14 ]
then
    >> $outbase/isoforms_14bp.fa
    isoform_606="GAATATTGTTTCAATGGTCGTGGCGCTGGTCAGACTGCTCGCTAAAGTCGTCGAAGAAACAGCTAGACCATCTATTCTCTGAACTTCCTAAGTGGTCGGGTTGGTCGTGGCGCTCATAGTTAGCCCAGAACGTCCAGTTTATGTGGTGGTGATCTTGGTTGCCAGACTCAGTATGGTCAAAATAGTTGTGCTCGCCATTTCAACTAGTGCAGTTGCTAAAGTGCCAACTAAGGTTCGTCTTGAAGCAATCGGTAGTATCTGCAGTCAAGTCGTAGTCGTAGGTATGGCACTTGTCGTGCGTCTAGGACGCGTTGAAGCTATCTTATTGCCCGAAGGAGCCGAAAGTCTGGGACACGTCGTGGGACGTATCGTTAGGACTATGGGACTAGTCGTCTTCGGTCTTACAAATGTTGCAATGCCCTGTGAGCTACTTATGAAACATGAATGGTCGGTCTTGTGGTCGCTTTTGTCGAAATCACCGAAGTGCGTATAATTGATGAACGAGTCGAAGTCGAAGTTGTTCAACAAGAGAGTTAGGTAGTGCC"  
    isoform_616="GAATATTGTTTCAATGGTCGTGGCGCTGGTCAGACTGCTCGCTAAAGTCGTCGAAGAAACAGCTAGACCATCTATTCTCTGAACTTCCTAAGTGGTCGGGTTGGTCGTGGCGCTCATAGTTAGCCCAGAACGTCCAGTTTATGTGGTGGTGATCTTGGTTGCCAGACTCAGTATGGTCAAAATAGTTGTGCTCGCCATTTCAACTAGTGCAGTTGCTAAAGTGCCAACTAAGGTTCGTCTTGAAGCAATCGGTAGTATCTGCAGTCAAGTCGTAGTCGTAGGTATGGCACTTGTCGTGCGTCTAGGACGCGTTGAAGCTATCTTATTGCCCGAAGGAGCCGAAAGTCTGGGACACGTCGTGGGACGTATCGTTAGGACTATGGGACTAGTCGTCTTCGGTCTTACAAATTGAGCTACTTATGAAACATGAATGGTCGGTCTTGTGGTCGCTTTTGTCGAAATCACCGAAGTGCGTATAATTGATGAACGAGTCGAAGTCGAAGTTGTTCAACAAGAGAGTTAGGTAGTGCC"  
    echo ">isoform_606" >> $outbase/isoforms_14bp.fa
    echo $isoform_606 >> $outbase/isoforms_14bp.fa
    echo ">isoform_616" >> $outbase/isoforms_14bp.fa
    echo $isoform_616 >> $outbase/isoforms_14bp.fa
    isoforms=$outbase/isoforms_14bp.fa
    cat $outbase/"SIRV606.sam" $outbase/"SIRV616.sam" > $outbase/"reads_to_isoforms_14bp.sam"
    alignments=$outbase/"reads_to_isoforms_14bp.sam"

fi



for id in $(seq 1 1 10) #  $(seq 1 1 10)
do 

    for depth in 10 20 #40 60 80 100
    do
        for p in 0.1 #0.2 0.3 0.4 0.5
        do
            echo $depth
            python $experiment_dir/sirv_sample_reads.py $isoforms $alignments  $p $depth $outbase/$id/$depth/$p
            # python $experiment_dir/simulate_reads.py --isoforms $outbase/$id/isoforms.fa --outfolder $outbase/$id/$depth/$p --probs $p  --nr_reads $depth > /dev/null


            # python $inbase/isONcorrect --fastq $outbase/$id/$depth/$p/reads.fq   --outfolder $outbase/$id/$depth/$p/isoncorrect_exact/ &> /dev/null            
            # python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/$p/isoncorrect_exact/corrected_reads.fastq  $outbase/$id/isoforms.fa  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation  > /dev/null
            # echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv 
            # echo -n  $id,exact,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            # awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/results.csv >> $results_file


        done
    done
done


#python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_frac_switches.py $results_file $plot_file"_minor_isoform_retained.pdf" 


