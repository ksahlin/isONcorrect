#!/bin/bash
set -e
if [ $# -ne 3 ]; then
        echo "Usage: `basename $0`  <raw_reads.fq>  <outfolder>  <num_cores> "
        exit -1
fi

raw_reads=$1
outfolder=$2
num_cores=$3
echo "Running: " `basename $0` $raw_reads $outfolder $num_cores
mkdir -p $outfolder

echo
echo "Will run pychopper (cdna_classifier.py), isONclust, and isONcorrect. Make sure you have these tools installed." 
echo "For installation see: https://github.com/ksahlin/isONcorrect#installation "
echo 

echo
echo "Running pychopper"
echo 

cdna_classifier.py  $raw_reads $outfolder/full_length.fq -t $num_cores 

echo
echo "Finished pychopper"
echo 



echo
echo "Running isONclust"
echo 

isONclust  --t $num_cores  --ont --fastq $outfolder/full_length.fq \
             --outfolder $outfolder/clustering

isONclust write_fastq --N 1 --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $outfolder/full_length.fq --outfolder  $outfolder/clustering/fastq_files
echo
echo "Finished isONclust"
echo 



echo
echo "Running isONcorrect"
echo 

run_isoncorrect --t $num_cores  --fastq_folder $outfolder/clustering/fastq_files  --outfolder $outfolder/correction/ 

echo
echo "Finished isONcorrect"
echo 


echo
echo "Merging reads back to single file. Corrected reads per cluster are still stored in: " $outfolder/correction/
echo 

# OPTIONAL BELOW TO MERGE ALL CORRECTED READS INTO ONE FILE
touch $outfolder/all_corrected_reads.fq
OUTFILES=$outfolder"/correction/"*"/corrected_reads.fastq"
for f in $OUTFILES
do 
  echo $f
  cat $f >> $outfolder/all_corrected_reads.fq
done

echo
echo "Finished with pipeline and wrote corrected reads to: " $outfolder/all_corrected_reads.fq
echo 

