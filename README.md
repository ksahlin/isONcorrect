isONcorrect
===========

[![ci](https://github.com/ksahlin/isONcorrect/actions/workflows/ci.yml/badge.svg)](https://github.com/ksahlin/isONcorrect/actions/workflows/ci.yml)

isONcorrect error-corrects Oxford Nanopore cDNA reads. It handles highly variable coverage and exon variation within reads, and leverages regions shared between reads from different isoforms to reach low error rates even for low-abundance transcripts. See the [paper](https://www.nature.com/articles/s41467-020-20340-8).

> ## v0.2.0 — isONcorrect is now written in Rust
>
> Same command line, so existing pipelines do not need editing. About **10x faster**, uses **less
> memory**, and **more accurate** — a bug in region selection meant every previous version corrected
> fewer regions than it should have. **Corrected output therefore differs from earlier releases.**
>
> The Python implementation is still here but is **deprecated**, kept as the reference the Rust
> version is verified against. Details, numbers, and the removed flags: **[CHANGELOG.md](CHANGELOG.md)**.


Install
-------

Needs a Rust toolchain (1.85+) and nothing else — no C/C++ toolchain, no CMake.

```bash
git clone https://github.com/ksahlin/isONcorrect.git
cd isONcorrect
cargo build --release --manifest-path rust/Cargo.toml
```

This builds two binaries, `isONcorrect` and `run_isoncorrect`, in
`rust/target/release/`. Copy them somewhere on your `PATH`.

No Rust? `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

Prefer not to build? Binaries for Linux and macOS (x86_64, arm64) are attached to
[release v0.2.0](https://github.com/ksahlin/isONcorrect/releases/tag/v0.2.0). Linux comes in two
flavours: take `musl` if your glibc is older than 2.34 (CentOS 7, RHEL 8, Ubuntu 20.04) or you would
rather not check, and `gnu` otherwise — see [CHANGELOG.md](CHANGELOG.md).

For the full pipeline you also want pychopper and isONclust:

```bash
conda create -n isoncorrect python=3.9 pip && conda activate isoncorrect
pip install isONclust
conda install -c bioconda "hmmer>=3.0" "pychopper>=2.0"
```

The deprecated Python version, if you need it to reproduce the paper:
`pip install isONcorrect && conda install -c bioconda spoa`.

### Test the install

From the repository root:

```bash
rust/target/release/isONcorrect --fastq test_data/isoncorrect/0.fastq --outfolder /tmp/isoncorrect_test
```

Under a second. Writes `/tmp/isoncorrect_test/corrected_reads.fastq` — 100 reads, same headers as the
input, corrected sequences.


Run
---

One command for the whole pipeline:

```bash
./scripts/correction_pipeline.sh raw_reads.fq outfolder 20   # reads, outdir, cores
```

Or the steps yourself — pychopper for full-length reads, isONclust to group them into genes, then
isONcorrect per cluster:

```bash
cdna_classifier.py raw_reads.fq out/full_length.fq -t 20

isONclust --t 20 --ont --fastq out/full_length.fq --outfolder out/clustering
isONclust write_fastq --N 1 --clusters out/clustering/final_clusters.tsv \
          --fastq out/full_length.fq --outfolder out/clustering/fastq_files

run_isoncorrect --t 20 --fastq_folder out/clustering/fastq_files --outfolder out/correction

cat out/correction/*/corrected_reads.fastq > out/all_corrected_reads.fq
```

Reads need not be full-length, but running pychopper first is advised for downstream analysis.

**Output** is one `corrected_reads.fastq` per cluster with the input headers. Note the quality string
is not real — it is `+` repeated to the length of the sequence, as in the Python version.

**Useful flags** (`--help` for the rest): `--split_wrt_batches` cuts runtime when isONclust produces a
few very large clusters; `--split_mod n --residual i` spreads `run_isoncorrect` across *n* nodes;
`--k 9 --w 10 --max_seqs 1000` reproduces the paper's settings rather than the faster defaults.


Paper data
----------

The result data behind the paper used to be committed here, making the repo ~2.4 GB to clone. It is
now on Zenodo — **<https://zenodo.org/records/21920617>** — and `tools/fetch_data.sh` restores and
checksums it into `data/`. You only need it to regenerate the paper figures under `paper/`.

Because the data was stripped from history, commit SHAs from before the rewrite no longer resolve.


Credits
-------

Please cite:

Sahlin, K., Medvedev, P. Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. *Nat Commun* **12**, 2 (2021). https://doi.org/10.1038/s41467-020-20340-8


Licence
-------

GPL v3.0, see [LICENSE.txt](LICENSE.txt).
