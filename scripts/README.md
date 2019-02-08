### Example usage:

For fasta

```
    python simulate_reads.py --coords 0 50 100 150 200 250 --ref /Users/kxs624/tmp/ISONCORRECT/ref_sim_6exons.fa --outfile /Users/kxs624/tmp/ISONCORRECT/SIM_READS/sim.fa --probs 1.0 0.2 1.0 --nr_reads 100
```

for fastq format, simply specify q as last letter in outfile, e.g.,

```
    python simulate_reads.py --coords 0 50 100 150 200 250 --ref /Users/kxs624/tmp/ISONCORRECT/ref_sim_6exons.fa --outfile /Users/kxs624/tmp/ISONCORRECT/SIM_READS/sim.fq --probs 1.0 0.2 1.0 --nr_reads 100
```