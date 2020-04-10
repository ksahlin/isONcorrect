import os
import argparse
import numpy as np
import misc_functions
import random
import math
from collections import defaultdict

# log how many reads are simulated for each variant.


def simulate_read(i, transcript_acc, isoform, error_lvls):
    read = []
    qual = []
    del_, ins, subs = 0,0,0
    was_del = False
    for l, n in enumerate(isoform):
        # if l <= 15 or l >= len(isoform) - 15: # no errors first and last 15 bases
        #     p_correct_reading = 1.0
        #     p_error = 1.0 - 0.995
        # else:
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
    assert len(read_seq) == len(qual_seq)
    acc = str(transcript_acc) + "_" +  str(i) #= (read_seq, qual_seq)

    # print(read_seq)
    # print(qual_seq)

    return acc, read_seq, qual_seq, del_, ins, subs 


def main(args):
    sequence_transcripts = {}
    sequence_transcripts = dict(misc_functions.read_fasta(open(args.sequence_material,"r")))

    # just generate all numbers at once and draw from this 5x should be enough
    ont_reads = {}
    reads_generated_log = defaultdict(int)
    errors = []
    tot_del_, tot_ins, tot_subs, tot_len = 0, 0, 0, 0
    read_errors = []
    if args.error_level == 4:
        error_lvls = [0.9, 0.95, 0.96, 0.98, 0.99, 0.995]
    elif  args.error_level == 12:
        error_lvls = [0.75, 0.85, 0.875, 0.91, 0.95, 0.98]        
    else:
        error_lvls = [0.85, 0.875, 0.9, 0.92, 0.96, 0.98, 0.99, 0.995]

    print("mean_error", sum(error_lvls)/len(error_lvls))
    all_transctript_accessions = list(sequence_transcripts.keys())
    if args.subsampling_experiment:
        # 1+2+..+10+20+..+100 = 545
        abundance_vector = []
        for i in list(range(1,11)) +  list(range(20,101,10)):
            for j in range(int(100/i)):
                abundance_vector.append(i)
        print(abundance_vector)
        print(len(abundance_vector))
        for i, acc in enumerate(sequence_transcripts):
            transcript = sequence_transcripts[acc]
            abundance = random.choice(abundance_vector)
            # print(abundance)
            for a in range(abundance):
                read_acc, read, qual, del_, ins, subs  = simulate_read(a, acc, transcript, error_lvls)
                tot_err = (del_ + ins + subs)/float(len(transcript)+ ins)
                ont_reads[read_acc] = (read, qual)
                args.logfile.write("del:{0}, ins:{1}, subs:{2}, tot_err:{3}\n".format(del_, ins, subs, tot_err))
                tot_del_ += del_
                tot_ins += ins
                tot_subs += subs
                tot_len += len(transcript)
                read_errors.append( (del_ + ins +subs)/float(len(transcript) + ins)  )
            if i % 500 == 0:
                print(i, "transcripts simulated from.")                


    elif args.uneven:
        # abundance = [1,2,4,8,16,32]
        abundance_vector = []
        for i in list(range(1,11)) +  list(range(20,101,10)):
            for j in range(int(100/i)):
                abundance_vector.append(i)

        transcript_weights = [ random.choice(abundance_vector) for i in range(len(all_transctript_accessions))]
        accessions = random.choices(all_transctript_accessions, weights=transcript_weights, k = args.read_count)
        for i, acc in enumerate(accessions):
            transcript = sequence_transcripts[acc]
            read_acc, read, qual, del_, ins, subs  = simulate_read(i, acc, transcript, error_lvls)
            tot_err = (del_ + ins + subs)/float(len(transcript))
            ont_reads[read_acc] = (read, qual)
            args.logfile.write("del:{0}, ins:{1}, subs:{2}, tot_err:{3}\n".format(del_, ins, subs, tot_err))
            tot_del_ += del_
            tot_ins += ins
            tot_subs += subs
            tot_len += len(transcript)
            read_errors.append( (del_ + ins +subs)/float(len(transcript) + ins)  )
            if i % 5000 == 0:
                print(i, "reads simulated.")
    else:
        accessions = random.choices(all_transctript_accessions, k = args.read_count)
        for i, acc in enumerate(accessions):
            transcript = sequence_transcripts[acc]
            read_acc, read, qual, del_, ins, subs  = simulate_read(i, acc, transcript, error_lvls)
            read_errors.append( (del_ + ins +subs)/float(len(transcript) + ins)  )
            ont_reads[read_acc] = (read, qual)
            tot_del_ += del_
            tot_ins += ins
            tot_subs += subs
            tot_len += len(transcript)
            if i % 5000 == 0:
                print(i, "reads simulated.")
    print("median error rate (divided by aligmnet length):", sorted(read_errors)[int(len(read_errors)/2)])
    print(tot_del_, tot_ins, tot_subs, (tot_del_ + tot_ins + tot_subs)/float(tot_len))

    # for acc, abundance in misc_functions.iteritems(reads_generated_log):
    #     args.logfile.write("{0}\t{1}\n".format(acc, abundance))

    # n = float(len(errors))
    # mu =  sum(errors) / n
    # sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), errors))) / (n - 1)) ** 0.5
    # min_error = min(errors)
    # max_error = max(errors)
    # errors.sort()
    # if len(errors) %2 == 0:
    #     median_error = (errors[int(len(errors)/2)-1] + errors[int(len(errors)/2)]) / 2.0
    # else:
    #     median_error = errors[int(len(errors)/2)]

    # args.logfile.write("mean error: {0}, sd error:{1}, min_error:{2}, max_error:{3}, median_error:{4}\n".format(mu, sigma, min_error, max_error, median_error))

    outfile = open(args.outfile, "w")
    for acc, (read_seq,qual_seq) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
        outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('sequence_material', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')
    parser.add_argument('read_count', type=int, help='Number of reads to simulate.')
    parser.add_argument('--uneven', action="store_true", help='Transcripts gets relative abundance transcript_weightss from 1,2,..,10,20,...,100.')
    parser.add_argument('--subsampling_experiment', action="store_true", help='Each transcript is sampled with abundance 1,2,..,10,20,...,100')
    parser.add_argument('--error_level', type=int, default = 7, help='Set error level to 4%, 7%, or 12%')


    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)
    misc_functions.mkdir_p(path_)
    args.logfile = open(os.path.join(path_, file_prefix + ".log"), "w")
    main(args)
