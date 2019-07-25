import os,sys
import random
import itertools
import argparse
import errno
import math
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



def simulate_reads( args, isoforms ):

    outfile = open(os.path.join(args.outfolder,"reads.fq"), "w")
    is_fastq = True #if outfile[-1] == "q" else False
    middle_exon_frequnecy = args.probs[1]
    no_middel_exon = [(isoforms[0][0] + "_{0}".format(i), isoforms[0][1]) for i in range(int(round(args.nr_reads*(1-args.probs[1]))))]
    has_middel_exon = [(isoforms[1][0] + "_{0}".format(i), isoforms[1][1]) for i in range(int(round(args.nr_reads*args.probs[1])))]
    isoforms_generated = has_middel_exon + no_middel_exon

    isoforms_abundance_out = open(os.path.join(args.outfolder,"isoforms_abundance.fa"), "w")
    for i_acc, isoform in isoforms_generated:
        isoforms_abundance_out.write(">{0}\n{1}\n".format(i_acc, isoform))
    isoforms_abundance_out.close()

    # print(len(has_middel_exon), len(no_middel_exon), args.nr_reads)
    assert len(isoforms_generated) == args.nr_reads
    # seq, acc = isoforms[list(isoforms.items())
    # seq = seq.upper()

    # exons = [seq[j_start: j_stop] for (j_start, j_stop) in zip(range(0,300, 50), range(50, 301, 50)) ]
    # exons_probs = [1.0, 0.2, 1.0, 1.0, 0.2, 1.0]

    # exon_coords = [(start, stop) for start, stop in zip(args.coords[:-1], args.coords[1:]) ]
    # exons = [seq[j_start: j_stop] for (j_start, j_stop) in exon_coords ]
    # exons_probs = args.probs
    # print(exon_coords)
    # print(exons)
    # print(exons_probs)

    # sys.exit()
    reads = {}
    error_lvls = [0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,0.98, 0.99, 0.995]
    for i, (i_acc, isoform) in enumerate(isoforms_generated):
        read = []
        qual = []
        # for j, e in enumerate(exons):
        # r = random.uniform(0,1)
        # if r > middle_exon_frequnecy:
        #     # has_exon = 0
        #     i_acc, isoform = isoforms[1]
        #     # continue
        # else:
        #     i_acc, isoform = isoforms[0]

        #     # has_exon = 1

        was_del = False
        for l, n in enumerate(isoform):
            if l <= 15 or l >= len(isoform) - 15: # no errors first and last 15 bases
                p_correct_reading = 1.0
                p_error = 1.0 - 0.995
            else:
                p_correct_reading = random.choice(error_lvls)
                p_error = 1.0 - p_correct_reading

            # print(round(10**(-math.log(p_correct_reading))), p_error, -math.log(p_error,10)*10)

            r = random.uniform(0,1)
            if r > p_correct_reading:
                error = True
            else:
                error = False

            if error:
                r = random.uniform(0,1)
                if r < 0.6: #deletion
                    was_del = p_error
                    pass 
                elif 0.6 <= r < 0.9:
                    read.append(random.choice("ACGT"))
                    qual.append( round(-math.log(p_error,10)*10) )

                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )

                    r_ins = random.uniform(0,1)
                    while r_ins >= 0.7:
                        read.append(random.choice("ACGT"))
                        r_ins = random.uniform(0,1)
                        qual.append( round(-math.log(0.7,10)*10) )

                # if r < 0.84:
                #     read.append(n)
                # elif 0.84 <= r <= 0.96:
                #     pass 
                # else:
                #     read.append(n)
                #     read.append(random.choice("ACGT"))
                #     r_ins = random.uniform(0,1)
                #     while r_ins >= 0.96:
                #         read.append(random.choice("ACGT"))
                #         r_ins = random.uniform(0,1)
            else:
                if was_del: # adding uncertainty from prevous deleted base
                    read.append(n)
                    qual.append( round(-math.log(was_del,10)*10) )
                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )
                was_del = False



        if not read:
            continue
        read_seq = "".join([n for n in read])
        qual_seq = "".join([chr(q + 33) for q in qual])
        reads[str(i_acc) + "_" + str(i)  ] = (read_seq, qual_seq)

        # print(read_seq)
        # print(qual_seq)


    if is_fastq:
        for acc, (read_seq,qual_seq) in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
            outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))
    else:
        for acc, (read_seq,qual_seq) in sorted(reads.items(), key = lambda x: len(x[1]), reverse = True):
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
    seq = ref[list(ref.keys())[0]]
    exon_coords = [(start, stop) for start, stop in zip(args.coords[:-1], args.coords[1:]) ]
    exons = [seq[j_start: j_stop] for (j_start, j_stop) in exon_coords ]

    # only two
    isoform = "".join([ex for ex in exons])
    isoforms_out.write(">sim|sim|{0}\n{1}\n".format("1", isoform))
    isoform = "".join([exons[0], exons[2]])
    isoforms_out.write(">sim|sim|{0}\n{1}\n".format("2", isoform))
    isoforms_out.close()

    # # all combinations
    # for i in range(1, len(exons)+1):
    #     for j, e in enumerate(itertools.combinations(exons,i)):
    #         isoform = "".join([ex for ex in e])
    #         isoforms_out.write(">{0}\n{1}\n".format("isoform_{0}_{1}".format(i,j), isoform))

    # isoforms_out.close()



def main(args):
    mkdir_p(args.outfolder)

    if args.sim_genome_len > 0:
        genome_out = open(os.path.join(args.outfolder,"reference.fa"), "w")
        genome_out.write(">{0}\n{1}\n".format("ref", "".join([random.choice("ACGT") for i in range(args.sim_genome_len)]) ))
        genome_out.close()
        generate_isoforms(args,genome_out.name)
        sys.exit()
    else:
        isoforms = [(acc,seq) for acc, (seq, qual) in readfq(open(args.isoforms, 'r'))]
    # print(isoforms)
    # check_valid_args(args, ref)

    simulate_reads(args, isoforms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot p-minimizers shared.")
    parser.add_argument('--ref', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('--isoforms', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
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