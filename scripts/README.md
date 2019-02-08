### Example usage:

For fasta

```
    python simulate_reads.py --coords 0 50 100 150 200 250 --ref test_data/six_exons.fa --outfile <outpath>/sim.fa --probs 1.0 0.2 1.0 --nr_reads 100
```

for fastq format, simply specify q as last letter in outfile, e.g.,

```
    python simulate_reads.py --coords 0 50 100 150 200 250 --ref test_data/six_exons.fa --outfile <outpath>/sim.fq --probs 1.0 0.2 1.0 --nr_reads 100
```

### Test data
The file 100reads.fq are 100 reads simulated from test_data/six_exons.fa. Exon coordinates are (0,49) (50,99), (100,149), (150,200), (200,250), (250,300). Exon probabilities are 1.0, 0.1, 1.0, 1.0, 0.1, 1.0.