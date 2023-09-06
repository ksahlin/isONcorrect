isONcorrect
===========

isONcorrect is a tool for error-correcting Oxford Nanopore cDNA reads. It is designed to handle highly variable coverage and exon variation within reads and achieves about a 0.5-1% median error rate after correction. It leverages regions shared between reads from different isoforms achieve low error rates even for low abundant transcripts. See [paper](https://www.nature.com/articles/s41467-020-20340-8) for details. 

**Update:** Since v0.0.8, isONcorrect uses different default parameters compared to what was used in the [paper](https://www.nature.com/articles/s41467-020-20340-8). The new parameters make isONcorrect 2-3 times faster and use 3-8 times less memory with only a small cost of increased median post-correction error rate. With the new parameter setting the correction accuracy is 98.5-99.3% instead of 98.9–99.6% on the data used in the paper. Current default uses `--k 9 --w 20 --max_seqs 2000`. To invoke settings used in paper, set parameters `--k 9 --w 10 --max_seqs 1000`.

Processing and error correction of full-length ONT cDNA reads is achieved by the pipeline of running [pychopper](https://github.com/nanoporetech/pychopper) --> [isONclust](https://github.com/ksahlin/isONclust) --> [isONcorrect](https://github.com/ksahlin/isONcorrect). All these steps can be run in one go with [this script](https://github.com/ksahlin/isONcorrect/blob/master/scripts/correction_pipeline.sh). See below for installation and usage. 


isONcorrect is distributed as a python package supported on Linux / OSX with python v>=3.4. [![Build Status](https://travis-ci.org/ksahlin/isONcorrect.svg?branch=master)](https://travis-ci.com/ksahlin/isONcorrect).

Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Using pip](#Using-pip)
    * [Downloading source from GitHub](#Downloading-source-from-github)
    * [Dependencies](#Dependencies)
    * [Testing installation](#testing-installation)
  * [USAGE](#USAGE)
    * [Running](#Running)
    * [Output](#Output)
    * [Parallelization across nodes](#Parallelization-across-nodes)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
=================

Typical install time on a desktop computer is about 5 minutes with conda for this software.

## Using conda
Conda is the preferred way to install isONcorrect.

1. Create and activate a new environment called isoncorrect

```
conda create -n isoncorrect python=3.9 pip 
conda activate isoncorrect
```

2. Install isONcorrect and its dependency `spoa`.

```
pip install isONcorrect
conda install -c bioconda spoa
```
3. You should now have 'isONcorrect' installed; try it:
```
isONcorrect --help
```

Upon start/login to your server/computer you need to activate the conda environment "isonclust" to run isONcorrect as:
```
conda activate isoncorrect
```

4. You probably want to install `pychopper` and `isONclust` in the isoncorrect environmment as well to run the complete correction pipeline if you haven't already. This can be done with:

```
pip install isONclust
conda install -c bioconda "hmmer>=3.0"
conda install -c bioconda "pychopper>=2.0"
```

You are now set to run the [correction_pipeline](https://github.com/ksahlin/isONcorrect/blob/master/scripts/correction_pipeline.sh). See [USAGE](https://github.com/ksahlin/isONcorrect#usage).


### Dependencies

isONcorrect has the following dependencies (the three first are automatically installed with `pip`)
* [edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python)
* [NumPy](https://numpy.org/) 
* [parasail](https://github.com/jeffdaily/parasail-python)
*  [spoa](https://github.com/rvaser/spoa) 


## Testing installation

You can verify successul installation by running isONcorrect on [these](https://github.com/ksahlin/isONcorrect/tree/master/test_data/isoncorrect/) two small datasets of 100 reads. Download the two datasets and put in a folder `test_data` run, e.g,

```
isONcorrect --fastq test_data/0.fastq \
            --outfolder [output path]
```
Expected runtime for this test data is about 15 seconds. The output will be found in `[output path]/corrected_reads.fastq` where the 100 reads have the same headers as in the original file, but with corrected sequence. Testing the paralleized version (by separate clusters) of isONcorrect can be done by running

```
./run_isoncorrect --t 3 --fastq_folder test_data/ \
                  --outfolder [output path]
```
 
This will perform correction on `0.fastq` and `1.fastq` in parallel. Expected runtime for this test data is about 15 seconds. The output will be found in `[output path]/0/corrected_reads.fastq` and `[output path]/1/corrected_reads.fastq` where the 100 reads in each separate cluster have the same headers as in the respective original files, but with corrected sequence. 


USAGE
=================
 
## Running

### Using correction_pipeline.sh

You can simply run `./correction_pipeline.sh <raw_reads.fq>  <outfolder>  <num_cores> ` which will perform the steps 1-5 below for you. The `correction_pipeline.sh` script is available in this repository [here](https://github.com/ksahlin/isONcorrect/blob/master/scripts/correction_pipeline.sh). Simply download the reposotory or the individual [correction_pipeline.sh file](https://github.com/ksahlin/isONcorrect/blob/master/scripts/correction_pipeline.sh). 

For a fastq file with raw ONT cDNA reads, the following pipeline is recommended:
1.  Produce full-length reads (with [pychopper](https://github.com/nanoporetech/pychopper) (a.k.a. `cdna_classifier`))
2.  Cluster the full length reads into genes/gene-families ([isONclust](https://github.com/ksahlin/isONclust))
3.  Make fastq files of each cluster (`isONclust write_fastq` command)
4.  Correct individual clusters ([isONcorrect](https://github.com/ksahlin/isONcorrect))
5.  Join reads back to a single fastq file (This is of course optional)


### Manually

The contents of the `correction_pipeline.sh` is (roughly) provided below. If you want more individual control over the steps than what the `correction_pipeline.sh` can do for you (such as different parameters in each step), you can modify/remove arguments as needed in `correction_pipeline.sh` or in the below script.  

```
#!/bin/bash

# Pipeline to get high-quality full-length reads from ONT cDNA sequencing

# Set path to output and number of cores
root_out="outfolder"
cores=20

mkdir -p $root_out

cdna_classifier.py  raw_reads.fq $root_out/full_length.fq -t $cores 

isONclust  --t $cores  --ont --fastq $root_out/full_length.fq \
             --outfolder $root_out/clustering

isONclust write_fastq --N 1 --clusters $root_out/clustering/final_clusters.tsv \
                      --fastq $root_out/full_length.fq --outfolder  $root_out/clustering/fastq_files 

run_isoncorrect --t $cores  --fastq_folder $root_out/clustering/fastq_files  --outfolder $root_out/correction/ 

# OPTIONAL BELOW TO MERGE ALL CORRECTED READS INTO ONE FILE
touch $root_out/all_corrected_reads.fq
OUTFILES=$root_out"/correction/"*"/corrected_reads.fastq"
for f in $OUTFILES
do 
  cat $f >> $outfolder/all_corrected_reads.fq
done
```

isONcorrect does not need ONT reads to be full-length (i.e., produced by `pychopper`), but unless you have specific other goals, it is advised to run pychopper for any kind of downstream analysis to guarantee full-length reads. 

## Output

The output of `run_isoncorrect` is one file per cluster with identical headers to the original reads.

### Few large clusters

For some datasets, e.g. targeted data, `isONclust` can produce highly uneven clusters, i.e., a few very large clusters and some/many small ones. In such cases, runtime can be reduced if the argument `--split_wrt_batches` is specified to `run_isoncorrect`.


## Parallelization across nodes

isONcorrect currently supports parallelization across cores on a node (parameter `--t`), but not across several nodes. There is a way to overcome this limitation if you have access to multiple nodes as follows. The `run_isoncorrect` step can be parallilized across n nodes by (in bash or other environment, e.g., snakemake) parallelizing the following commands

```
run_isoncorrect --fastq_folder outfolder/clustering/fastq_files  --outfolder /outfolder/correction/ --split_mod n --residual 0
run_isoncorrect --fastq_folder outfolder/clustering/fastq_files  --outfolder /outfolder/correction/ --split_mod n --residual 1
run_isoncorrect --fastq_folder outfolder/clustering/fastq_files  --outfolder /outfolder/correction/ --split_mod n --residual 2
...
run_isoncorrect --fastq_folder outfolder/clustering/fastq_files  --outfolder /outfolder/correction/ --split_mod n --residual n-1
```
Which tells isONcorrect to only work with distinct cluster IDs.

CREDITS
----------------

Please cite [1] when using isONcorrect.

1. Sahlin, K., Medvedev, P. Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. Nat Commun 12, 2 (2021). https://doi.org/10.1038/s41467-020-20340-8  [Link](https://www.nature.com/articles/s41467-020-20340-8).

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/isONcorect/blob/master/LICENCE.txt).

