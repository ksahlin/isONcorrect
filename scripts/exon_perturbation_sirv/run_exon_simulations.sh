#!/bin/bash

# RUN scripts e.g. as: ./run_exon_simulations.sh /Users/kxs624/Documents/workspace/isONcorrect/  /Users/kxs624/tmp/ISONCORRECT/SIMULATED_DATA/exon_perturbation_sirv/

inbase=$1/
outroot=$2


outbase=$outroot
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




exon_length=14 # either 6 or 14

if [ $exon_length == 6 ]
then
    > $outbase/isoforms_6bp.fa
    SIRV506="GCATACTGCAAGACTTGTAGGCCCATGAATTACCCGCTTAGGCAATGGGTAGACCTTTCCTTGCCGCGTGGAACAAACCGCCATGAGTGTTTTACGTTTACGGTCGGCGACCGGAGTCACAGCGCGACCAACGGGGCCAGGCAAAGGCCCATAGGCTTTTCCATCGGCACGCGACCGATTCCATCAAGTATAGCCACGGTCGTAACTTCGACGGCGGCAATAGTTTATGGCGCTAACACCTTTCCATCGCTGAGTCCAAACGTATGACGTGAGGTGTTAGGTGCAATCCCCATTCGGGACGGACACGCGAGTAAAGGTCTAACATGTTTCAACTTTCCAGTTAGTCGCCGTCCTATATTGGCAAGACCGGATGGCCGTCAACCAGTTCCATTCTGACTATTTGCCACATTGATTAAAGAAGTAGGAGCGTTAGTCGGCGCGTGCGCGGGACTATTTAGCGAAGGTCGCCGTATCGCGGGTGGCGGCGTTTCCGCCATAACAACAAGCGTTTAAGGTCAGAGTCCGCGACGAATGAAACGCCGTGCAAGCGGC"  
    SIRV511="GCATACTGCAAGACTTGTAGGCCCATGAATTACCCGCTTAGGCAATGGGTAGACCTTTCCTTGCCGCGTGGAACAAACCGCCATGAGTGTTTTACGTTTACGGTCGGCGACCGGAGTCACAGCGCGACCAACGGGCAAAGGCCCATAGGCTTTTCCATCGGCACGCGACCGATTCCATCAAGTATAGCCACGGTCGTAACTTCGACGGCGGCAATAGTTTATGGCGCTAACACCTTTCCATCGCTGAGTCCAAACGTATGACGTGAGGTGTTAGGTGCAATCCCCATTCGGGACGGACACGCGAGTAAAGGTCTAACATGTTTCAACTTTCCAGTTAGTCGCCGTCCTATATTGGCAAGACCGGATGGCCGTCAACCAGTTCCATTCTGACTATTTGCCACATTGATTAAAGAAGTAGGAGCGTTAGTCGGCGCGTGCGCGGGACTATTTAGCGAAGGTCGCCGTATCGCGGGTGGCGGCGTTTCCGCCATAACAACAAGCGTTTAAGGTCAGAGTCCGCGACGAATGAAACGCCGTGCAAGCGGC"  
    echo ">SIRV511" >> $outbase/isoforms_6bp.fa
    echo $SIRV511 >> $outbase/isoforms_6bp.fa
    echo ">SIRV506" >> $outbase/isoforms_6bp.fa
    echo $SIRV506 >> $outbase/isoforms_6bp.fa
    isoforms=$outbase/isoforms_6bp.fa

    # filter so that all reads have strictly better fit to one of the isoforms
    python $experiment_dir/get_supporting_reads.py $outbase/"SIRV506.sam" $isoforms 6 $outbase/"SIRV506_filtered.sam"
    python $experiment_dir/get_supporting_reads.py $outbase/"SIRV511.sam" $isoforms 6 $outbase/"SIRV511_filtered.sam"

    cat $outbase/"SIRV506_filtered.sam" $outbase/"SIRV511_filtered.sam" > $outbase/"reads_to_isoforms_6bp.sam"
    alignments=$outbase/"reads_to_isoforms_6bp.sam"

fi

if [ $exon_length == 14 ]
then
    > $outbase/isoforms_14bp.fa
    SIRV606="GAATATTGTTTCAATGGTCGTGGCGCTGGTCAGACTGCTCGCTAAAGTCGTCGAAGAAACAGCTAGACCATCTATTCTCTGAACTTCCTAAGTGGTCGGGTTGGTCGTGGCGCTCATAGTTAGCCCAGAACGTCCAGTTTATGTGGTGGTGATCTTGGTTGCCAGACTCAGTATGGTCAAAATAGTTGTGCTCGCCATTTCAACTAGTGCAGTTGCTAAAGTGCCAACTAAGGTTCGTCTTGAAGCAATCGGTAGTATCTGCAGTCAAGTCGTAGTCGTAGGTATGGCACTTGTCGTGCGTCTAGGACGCGTTGAAGCTATCTTATTGCCCGAAGGAGCCGAAAGTCTGGGACACGTCGTGGGACGTATCGTTAGGACTATGGGACTAGTCGTCTTCGGTCTTACAAATGTTGCAATGCCCTGTGAGCTACTTATGAAACATGAATGGTCGGTCTTGTGGTCGCTTTTGTCGAAATCACCGAAGTGCGTATAATTGATGAACGAGTCGAAGTCGAAGTTGTTCAACAAGAGAGTTAGGTAGTGCC"  
    SIRV616="GAATATTGTTTCAATGGTCGTGGCGCTGGTCAGACTGCTCGCTAAAGTCGTCGAAGAAACAGCTAGACCATCTATTCTCTGAACTTCCTAAGTGGTCGGGTTGGTCGTGGCGCTCATAGTTAGCCCAGAACGTCCAGTTTATGTGGTGGTGATCTTGGTTGCCAGACTCAGTATGGTCAAAATAGTTGTGCTCGCCATTTCAACTAGTGCAGTTGCTAAAGTGCCAACTAAGGTTCGTCTTGAAGCAATCGGTAGTATCTGCAGTCAAGTCGTAGTCGTAGGTATGGCACTTGTCGTGCGTCTAGGACGCGTTGAAGCTATCTTATTGCCCGAAGGAGCCGAAAGTCTGGGACACGTCGTGGGACGTATCGTTAGGACTATGGGACTAGTCGTCTTCGGTCTTACAAATTGAGCTACTTATGAAACATGAATGGTCGGTCTTGTGGTCGCTTTTGTCGAAATCACCGAAGTGCGTATAATTGATGAACGAGTCGAAGTCGAAGTTGTTCAACAAGAGAGTTAGGTAGTGCC"  
    echo ">SIRV606" >> $outbase/isoforms_14bp.fa
    echo $SIRV606 >> $outbase/isoforms_14bp.fa
    echo ">SIRV616" >> $outbase/isoforms_14bp.fa
    echo $SIRV616 >> $outbase/isoforms_14bp.fa
    isoforms=$outbase/isoforms_14bp.fa

    # filter so that all reads have strictly better fit to one of the isoforms
    python $experiment_dir/get_supporting_reads.py $outbase/"SIRV606.sam" $isoforms 7 $outbase/"SIRV606_filtered.sam"
    python $experiment_dir/get_supporting_reads.py $outbase/"SIRV616.sam" $isoforms 7 $outbase/"SIRV616_filtered.sam"

    cat $outbase/"SIRV606_filtered.sam" $outbase/"SIRV616_filtered.sam" > $outbase/"reads_to_isoforms_14bp.sam"
    alignments=$outbase/"reads_to_isoforms_14bp.sam"

fi




for id in $(seq 1 1 10) #  $(seq 1 1 10)
do 

    for depth in 10 20 40 60 80 100
    do
        for p in 0.1 0.2 0.3 0.4 0.5 
        do

            python $experiment_dir/sirv_sample_reads.py $isoforms $alignments  $p $depth $outbase/$id/$depth/$p/"reads.fq"

            python $inbase/isONcorrect --fastq $outbase/$id/$depth/$p/reads.fq   --outfolder $outbase/$id/$depth/$p/isoncorrect_exact/  &> /dev/null            
            python $experiment_dir/evaluate_simulated_reads.py  $outbase/$id/$depth/$p/isoncorrect_exact/corrected_reads.fastq  $isoforms  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation  #> /dev/null
            echo -n  $id,exact,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv 
            echo -n  $id,exact,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/summary.csv >> $summary_file
            awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/isoncorrect_exact/evaluation/results.csv >> $results_file

            # python $experiment_dir/evaluate_simulated_reads.py   $outbase/$id/$depth/$p/reads.fq $isoforms $outbase/$id/$depth/$p/evaluation_reads #> /dev/null
            # echo -n  $id,original,$depth,$p,&& head -n 1 $outbase/$id/$depth/$p/evaluation_reads/summary.csv 
            # echo -n  $id,original,$depth,$p, >> $summary_file && head -n 1 $outbase/$id/$depth/$p/evaluation_reads/summary.csv  >> $summary_file
            # awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_p=$p '{if (NR!=1) {print awk_id",original,"awk_depth","awk_p","$0}}'  $outbase/$id/$depth/$p/evaluation_reads/results.csv >> $results_file

        done
    done
done


# python $experiment_dir/plot_mutation_data.py $summary_file $plot_file"_tot.pdf" error_rate
python $experiment_dir/plot_frac_switches.py $results_file $plot_file"_minor_isoform_retained.pdf" 


