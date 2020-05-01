
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

from collections import defaultdict

import parasail
import pysam
import gffutils




def get_error_rates(sam_file, annotated_splice_coordinates_pairs, args): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_accuracies = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            called_intron_lengths = [length for type_, length in read.cigartuples if type_ == 3]


            # Deletion is treated separately because of potential short introns -- if in annotation
            if read.reference_name in annotated_splice_coordinates_pairs:
                del_new = 0
                chr_coordinate_pairs = annotated_splice_coordinates_pairs[read.reference_name]
                read_deletion_coordinate_pairs = get_deletion_sites(read)
                for start_del, stop_del, del_length in read_deletion_coordinate_pairs:
                    if (start_del, stop_del) not in chr_coordinate_pairs:
                        del_new += del_length
                    else:
                        print("True intron masked as deletion, length:", del_length, read.reference_name, start_del, stop_del)

            else:
                del_new = sum([length for type_, length in read.cigartuples if type_ == 2])

            read_accuracies[read.query_name] = round( 100*((ins + del_new + subs)/(ins + del_new + subs + matches)), 4)



    return read_accuracies     


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



def get_deletion_sites(read):
    deletion_sites = []
    ref_pos = read.reference_start
    
    for i, (t,l) in enumerate(read.cigartuples):
        if t == 2:
            deletion_start = ref_pos
            ref_pos += l
            deletion_stop = ref_pos
            deletion_sites.append( (deletion_start, deletion_stop, l) )

        elif t == 7 or t== 0 or t == 8:
            ref_pos += l
        elif t == 3:
            # splice_start = ref_pos
            ref_pos += l
            # splice_stop = ref_pos
            # splice_sites.append( (splice_start, splice_stop) )

        elif t == 1 or t == 4 or t == 5: # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return deletion_sites


def pickle_dump(data, filename):
    with open(os.path.join(args.outfolder,filename), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data


def main(args):

    annotated_ref_isoforms = pickle_load(os.path.join( args.database_path, 'annotated_ref_isoforms.pickle') )
    annotated_splice_coordinates = pickle_load(os.path.join( args.database_path, 'annotated_splice_coordinates.pickle') )
    annotated_splice_coordinates_pairs = pickle_load(os.path.join( args.database_path, 'annotated_splice_coordinates_pairs.pickle') )
    minimum_annotated_intron = pickle_load(os.path.join( args.database_path, 'minimum_annotated_intron.pickle') )

    read_accuracies = get_error_rates(args.samfile, annotated_splice_coordinates_pairs, args)

    # orig, orig_detailed = get_error_rate_stats_per_read(orig_primary_locations, reads, annotated_splice_coordinates_pairs, args)

    outfile = open(args.outfile, "w")
    outfile.write("Read\tAcc\n")
    for i, (read_id, accuracy) in enumerate(read_accuracies.items()):
        outfile.write("{0}\t{1}\n".format(i, accuracy))

    outfile.close()    




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('samfile', type=str, help='Path to the original read file')
    parser.add_argument('database_path', type=str, default='', help='Path to load database from.')
    parser.add_argument('outfile', type=str, help='Output path of results')

    args = parser.parse_args()

    main(args)

