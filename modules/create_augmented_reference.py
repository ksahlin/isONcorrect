from __future__ import print_function
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout




def run_spoa(reads, ref_out_file, spoa_out_file, spoa_path):
    with open(spoa_out_file, "w") as output_file:
        print('Running spoa...', end=' ')
        stdout.flush()
        null = open("/dev/null", "w")
        subprocess.check_call([ spoa_path, reads, "-l", "0", "-r", "2", "--gap", "-2"], stdout=output_file, stderr=null)
        print('Done.')
        stdout.flush()
    output_file.close()
    l = open(spoa_out_file, "r").readlines()
    consensus = l[1].strip()
    msa = [s.strip() for s in l[3:]]
    print(consensus)
    
    r = open(ref_out_file, "w")
    r.write(">{0}\n{1}".format("reference", consensus))
    r.close()

    return consensus, msa

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('reads', type=str, help='Fasta or fastq')
    parser.add_argument('outfile', type=str, help='Fasta or fastq')
    parser.add_argument('--spoa_path', type=str, default='spoa', required=False, help='Path to spoa binary with bwa binary name at the end.')


    args = parser.parse_args()

    run_spoa(args.reads, args.outfile, args.spoa_path)
