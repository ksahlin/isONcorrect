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




def simulate_reads( args, isoforms ):

    outfile = open(os.path.join(args.outfolder,"reads.fq"), "w")
    is_fastq = True #if outfile[-1] == "q" else False
    no_mutation = [(isoforms[0][0] + "_{0}".format(i), isoforms[0][1])  for i in range(int(round(args.nr_reads*(1-args.probs))))]
    has_mutation = [(isoforms[1][0] + "_{0}".format(i), isoforms[1][1]) for i in range(int(round(args.nr_reads*args.probs)))]
    isoforms_generated =  no_mutation + has_mutation
    random.shuffle(isoforms_generated)

    isoforms_abundance_out = open(os.path.join(args.outfolder,"isoforms_abundance.fa"), "w")
    for i_acc, isoform in isoforms_generated:
        isoforms_abundance_out.write(">{0}\n{1}\n".format(i_acc, isoform))
    isoforms_abundance_out.close()

    assert len(isoforms_generated) == args.nr_reads
    reads = {}
    error_lvls = [0.85, 0.875, 0.9, 0.92, 0.96, 0.98, 0.99, 0.995] # 7%

    for i, (i_acc, isoform) in enumerate(isoforms_generated):
        read = []
        qual = []
        del_, ins, subs = 0,0,0
        was_del = False
        for l, n in enumerate(isoform):
            p_correct_reading = random.choice(error_lvls)
            p_error = 1.0 - p_correct_reading

            r = random.uniform(0,1)
            if r > p_correct_reading:
                error = True
            else:
                error = False

            if error:
                r = random.uniform(0,1)
                if r < 0.45: #deletion
                    was_del = p_error
                    del_ += 1
                    pass 
                elif 0.45 <= r < 0.8:
                    read.append(random.choice( list(set("ACGT").difference( set(n))) ))
                    qual.append( round(-math.log(p_error,10)*10) )
                    subs += 1
                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )
                    read.append(random.choice("ACGT"))
                    qual.append( round(-math.log(p_error,10)*10) )
                    ins += 1
                    r_ins = random.uniform(0,1)
                    while r_ins >= 0.7:
                        read.append(random.choice("ACGT"))
                        r_ins = random.uniform(0,1)
                        qual.append( round(-math.log(0.7,10)*10) )
                        ins += 1

            else:
                if was_del: # adding uncertainty from prevous deleted base
                    read.append(n)
                    qual.append( round(-math.log(was_del,10)*10) )
                else:
                    read.append(n)
                    qual.append( round(-math.log(p_error,10)*10) )
                was_del = False

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


def generate_transcripts(transcript, ref_path):
    transcripts_out = open(os.path.join(args.outfolder,"isoforms.fa"), "w")
    # only two
    transcripts_out.write(">sim|sim|{0}\n{1}\n".format("1", transcript))
    
    mut_position = random.randrange(75,125)
    nucl = transcript[mut_position]
    new_nucl = set(["A","C","G","T"]) - set([nucl])
    subs = random.choice(list(new_nucl))
    assert nucl != subs
    mutated_transcript = transcript[ : mut_position] + subs + transcript[ mut_position + 1 :]
    
    # nucl2 = transcript[mut_position+1]
    # new_nucl = set(["A","C","G","T"]) - set([nucl2])
    # subs2 = random.choice(list(new_nucl))
    # assert nucl2 != subs2
    # # mutated_transcript = transcript[ : mut_position] + subs + subs2 + transcript[ mut_position + 2 :]
    # mutated_transcript = transcript[ : mut_position] + transcript[ mut_position + 2 :]

    assert mutated_transcript != transcript
    transcripts_out.write(">sim|sim|{0}\n{1}\n".format("2", mutated_transcript))
    transcripts_out.close()


def main(args):
    mkdir_p(args.outfolder)

    if args.sim_mutation:
        generate_transcripts(args.ref, args.outfolder)
        sys.exit()
    else:
        isoforms = [(acc,seq) for acc, (seq, qual) in readfq(open(args.isoforms, 'r'))]

    simulate_reads(args, isoforms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot p-minimizers shared.")
    parser.add_argument('--ref', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('--isoforms', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('--nr_reads', type=int, default = 200, help='Outfolder path')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    parser.add_argument('--probs', type=float, help='Probability frequency for mutation.')
    parser.add_argument('--sim_mutation',  action="store_true", help='Simulate mutation')
        
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)