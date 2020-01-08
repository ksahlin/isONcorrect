isONcorrect
===========

isONcorrect is a tool for error-correcting  Oxford Nanopore cDNA reads. It is designed to handle highly variable coverage and exon variation within reads and achieves about a 0.5-1% median error rate after correction (see [preprint]() for details). It leverages regions shared between reads from different isoforms achieve low error rates even for low abundant transcripts. See [preprint]() for details.  

Processing and error correction of full-length ONT cDNA reads is acheved by the pipeline of running [pychopper](https://github.com/nanoporetech/pychopper) --> [isONclust](https://github.com/ksahlin/isONclust) --> [isONcorrect](https://github.com/ksahlin/isONcorrect) 


isONcorrect is distributed as a python package supported on Linux / OSX with python v>=3.4. [![Build Status](https://travis-ci.org/ksahlin/isONcorrect.svg?branch=master)](https://travis-ci.org/ksahlin/isONcorrect).

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
    * [Parameters](#Parameters)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
----------------

### Using conda
Conda is the preferred way to install isONcorrect.

1. Create and activate a new environment called isoncorrect

```
conda create -n isoncorrect python=3 pip 
source activate isoncorrect
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
source activate isoncorrect
```

### Using pip 

`pip` is pythons official package installer and is included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 

To install isONcorrect, run:
```
pip install isONcorrect
```
Then install [spoa](https://github.com/rvaser/spoa).


### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as isONcorrect has not been tested with earlier versions of these libraries. However, isONcorrect may also work with earliear versions of these libaries.
* [spoa](https://github.com/rvaser/spoa) (1.1.5)
* [edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python) (1.1.2)
* [NumPy](https://numpy.org/) (1.16.2)

In addition, please make sure you use python version >=3.

With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/isONcorrect.git
cd isONcorrect
./isONcorrect
```

### Testing installation

You can verify successul installation by running isONcorrect on this [small dataset](https://github.com/ksahlin/isONcorrect/tree/master/test/sample_alz_2k.fastq). Simply download the test dataset and run:

```
isONcorrect --fastq [path to isONcorrect basefolder]/test_data/isoncorrect/0.fastq --outfolder [output path]
```


USAGE
-------
 
### Running

For a file with raw ONT cDNA reads the following pipeline is recommended (bash script provided below)
1.  Get full-length ONT cDNA sequences produced by [pychopper](https://github.com/nanoporetech/pychopper) (a.k.a. `cdna_classifier`)
2.  Cluster the full length reads where a cluster corresponds to a gene/gene-family
3.  Make fastq files of each cluster
4.  Correct individual clusters
5.  Optional (join reads from separate clusters back to a single file)

Below shows specific pipeline script to go from raw reads `raw_reads.fq` to corrected full-length reads `all_corrected_reads.fq` (please modify/remove arguments as needed). 

```
#!/bin/bash

# isonano pipeline to get high quality full length reads from transcriots

cdna_classifier.py  raw_reads.fq outfolder/reads_full_length.fq \
                      [-t cores]  [-w outfolder/rescued.fq  -u outfolder/unclassified.fq  -S outfolder/stats.txt] 

isONclust  --ont --fastq outfolder/reads_full_length.fq \
             --outfolder outfolder/clustering  [--t cores] 
isONclust write_fastq --N 1 --clusters outfolder/clustering/final_clusters.csv \
          --fastq reads_full_length.fq --outfolder  outfolder/clustering/fastq_files 
run_isoncorrect --fastq_folder outfolder/clustering/fastq_files  --outfolder /outfolder/correction/ 

# OPTIONAL BELOW 
touch all_corrected_reads.fq
for f in in outfolder/clustering/fastq_files/*/corrected_reads.fastq; 
do 
  cat {f} >> all_corrected_reads.fq
done
```

isONcorrect does not need ONT reads to be full-length (i.e., produced by `pychopper`), but unless you have specific other goals, it is advised to run pychopper for any kind of downstream analysis to guarantee full-length reads. 

### Output

The output of `run_isoncorrect` are one file per cluster with identical headers to the original reads.



CREDITS
----------------

Please cite [1] when using isONcorrect.

1. Kristoffer Sahlin, Botond Sipos, Phillip L James, Daniel Turner, Paul Medvedev (2020)  [Link]().

Bib record: 

LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/isONcorect/blob/master/LICENCE.txt).

