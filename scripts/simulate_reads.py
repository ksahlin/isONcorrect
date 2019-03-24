import os,sys
import random
import itertools
import argparse
import errno
# from collections import deque

'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break

def check_valid_args(args, ref):
    # assert args.start < args.stop and args.start < min(args.coords)
    # assert args.stop > args.start and args.stop > max(args.coords)
    assert args.coords == sorted(args.coords) # has to be in increasing order
    assert len(ref[list(ref.keys())[0]][0]) >= max(args.coords)
    print(args.coords, args.probs)
    assert (len(args.coords)-1) == len(args.probs)



def simulate_reads( args, ref ):

    outfile = open(os.path.join(args.outfolder,"reads.fq"), "w")
    is_fastq = True #if outfile[-1] == "q" else False
    seq, acc = ref[list(ref.keys())[0]]
    seq = seq.upper()

    # exons = [seq[j_start: j_stop] for (j_start, j_stop) in zip(range(0,300, 50), range(50, 301, 50)) ]
    # exons_probs = [1.0, 0.2, 1.0, 1.0, 0.2, 1.0]

    exon_coords = [(start, stop) for start, stop in zip(args.coords[:-1], args.coords[1:]) ]
    exons = [seq[j_start: j_stop] for (j_start, j_stop) in exon_coords ]
    exons_probs = args.probs
    print(exon_coords)
    print(exons)
    print(exons_probs)
    # sys.exit()
    reads = {}
    for i in range(args.nr_reads):
        read = []
        for j, e in enumerate(exons):
            r = random.uniform(0,1)
            if r > exons_probs[j]:
                continue
            for n in e:
                r = random.uniform(0,1)
                if r < 0.84:
                    read.append(n)
                elif 0.84 <= r <= 0.96:
                    pass 
                else:
                    read.append(n)
                    read.append(random.choice("ACGT"))
                    r_ins = random.uniform(0,1)
                    while r_ins >= 0.96:
                        read.append(random.choice("ACGT"))
                        r_ins = random.uniform(0,1)
        if not read:
            continue
        read_seq = "".join([n for n in read])
        reads[str(i)] = read_seq


    if is_fastq:
        for acc, read_seq in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
            outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", "+"*len(read_seq)))
    else:
        for acc, read_seq in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
            outfile.write(">{0}\n{1}\n".format(acc, read_seq))
    
    outfile.close()


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def generate_isoforms(args, ref_path):
    isoforms_out = open(os.path.join(args.outfolder,"isoforms.fa"), "w")
    ref = {acc : seq for acc, (seq, qual) in readfq(open(ref_path, 'r'))}
    # print(ref_path, ref)
    seq = ref[list(ref.keys())[0]]
    exon_coords = [(start, stop) for start, stop in zip(args.coords[:-1], args.coords[1:]) ]
    exons = [seq[j_start: j_stop] for (j_start, j_stop) in exon_coords ]
    for i in range(1, len(exons)+1):
        for j, e in enumerate(itertools.combinations(exons,i)):
            isoform = "".join([ex for ex in e])
            isoforms_out.write(">{0}\n{1}\n".format("isoform_{0}_{1}".format(i,j), isoform))

    isoforms_out.close()

def main(args):
    mkdir_p(args.outfolder)

    if args.sim_genome_len > 0:
        genome_out = open(os.path.join(args.outfolder,"reference.fa"), "w")
        genome_out.write(">{0}\n{1}\n".format("ref", "".join([random.choice("ACGT") for i in range(args.sim_genome_len)]) ))
        genome_out.close()
        generate_isoforms(args,genome_out.name)
        sys.exit()
    else:
        ref = {acc : (seq, qual) for acc, (seq, qual) in readfq(open(args.ref, 'r'))}

    check_valid_args(args, ref)

    simulate_reads(args, ref)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot p-minimizers shared.")
    parser.add_argument('--ref', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('--nr_reads', type=int, default = 200, help='Outfolder path')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    # parser.add_argument('--outfile', type=str, help='Simulated reads file. If ending in "q", fastq format will be output with all quality values being 10, i.e., "+". ')

    parser.add_argument('--coords', nargs='+', type=int, help='Exon coordinates. For example, --coords 0 50 100 150 220 240 gives exons 0-50, 100-150, 220-240. Has to be an even number.')
    parser.add_argument('--probs', nargs='+', type=float, help='Probability for each exon to be sampled. For example, --p 1.0, 0.2, 1.0 includes first and third exon in all reads and second exon is, on average, included in 1/5 of the reads.')
    parser.add_argument('--sim_genome_len', type=int, default = 0, help='Length if simulation of genome')
    
    # parser.add_argument('--start', type=int, help='Start sequencing coordinate, Has to be smaller than smallest jcoord and stop. ')
    # parser.add_argument('--stop', type=int, help='Stop sequencing coordinate, Has to be larger than largest jcoord and start.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)