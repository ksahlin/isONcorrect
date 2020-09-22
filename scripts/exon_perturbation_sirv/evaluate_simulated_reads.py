from __future__ import print_function
import os,sys
import argparse
import re
# import ssw
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
from collections import Counter
import parasail
from collections import defaultdict



import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns


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





import edlib

def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    # print(result["cigar"])
    return ed


def get_minimizers_2set_simple(querys, targets):
    best_edit_distances = {}

    for acc1, seq1 in querys.items():
        best_ed = 100000 #len(seq1)
        for acc2, seq2 in targets.items():
            if acc1[:3] == acc2[:3]: #same gene family
                edit_distance = edlib_ed(seq1, seq2, mode="NW", k = 2*len(seq1)) # seq1 = query, seq2 = target
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = best_ed
                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = best_ed
        # print(best_ed)
        # print(seq2)
        # print(seq1)
        # print()


    return best_edit_distances



def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


# def parasail_alignment(read, reference, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -2, opening_penalty = 2, gap_ext = 1, ends_discrepancy_threshold = 0):
#     user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
#     result = parasail.sg_trace_scan_16(read, reference, opening_penalty, gap_ext, user_matrix)
#     if result.saturated:
#         print("SATURATED!")
#         result = parasail.sg_trace_scan_32(read, reference, opening_penalty, gap_ext, user_matrix)
#     if sys.version_info[0] < 3:
#         cigar_string = str(result.cigar.decode).decode('utf-8')
#     else:
#         cigar_string = str(result.cigar.decode, 'utf-8')
    
#     read_alignment, ref_alignment = cigar_to_seq(cigar_string, read, reference)
#     return read_alignment, ref_alignment

# def get_alignments(best_edit_distances, querys, targets):
#     # score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
#     # aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
#     best_edit_distances_parasail = {}
#     for acc1 in best_edit_distances:
#         seq1 = querys[acc1]
#         best_ed = len(seq1)
#         best_edit_distances_parasail[acc1] = {}

#         for acc2 in best_edit_distances[acc1]:
#             seq2 = targets[acc2]            
#             read_alignment, ref_alignment = parasail_alignment(seq1, seq2)
#             insertions = ref_alignment.count("-")
#             deletions = read_alignment.count("-")
#             indels =  insertions + deletions
#             mismatches = len([1 for n1, n2 in zip(read_alignment, ref_alignment) if n1 != n2 and n1 != "-" and n2 != "-"] )

#             insertions_minus_ends = ref_alignment[20:-20].count("-")
#             deletions_minus_ends = read_alignment[20:-20].count("-")
#             indels_minus_ends =  insertions_minus_ends + deletions_minus_ends
#             mismatches_minus_ends = len([1 for n1, n2 in zip(read_alignment[20:-20], ref_alignment[20:-20]) if n1 != n2 and n1 != "-" and n2 != "-"] )


#             sw_ed = mismatches + indels
#             if sw_ed < best_ed:
#                 best_edit_distances_parasail[acc1] = {}
#                 best_edit_distances_parasail[acc1][acc2] = (deletions, insertions, mismatches)
#                 best_ed = sw_ed

#             elif sw_ed == best_ed:
#                 best_edit_distances_parasail[acc1][acc2] = (deletions, insertions, mismatches)

#             # seq1_aln, match_line, seq2_aln = result.alignment

#     return best_edit_distances_parasail


def get_best_match(corrected_reads, reference_transcripts, outfolder, params):
    out_file = open(os.path.join(outfolder, "results.csv"), "w")
    summary_file = open(os.path.join(outfolder, "summary.csv"), "w")

    errors_container = {}

    original_abundances = defaultdict(int)
    for acc in corrected_reads.keys():
        original_abundances[ acc.split("_")[0] ] += 1 # true origin transcript is present in header

    corrected_read_abundances = {}    
    for acc in reference_transcripts:
        corrected_read_abundances[ acc.split("_")[0]]  = 0

    errors_container = {}
    tot_errors = 0
    total_nucleotides = 0
    for acc, read in corrected_reads.items():
        best_ed = 1000
        best_ed_isoform = ''
        for ref_acc, ref in reference_transcripts.items():
            # since reads might be a bit shorter in ends, or contain adapters outside ends of transcript we do "semi-global" alignment with allowing infix.
            edit_distance1 = edlib_ed(read, ref, mode="HW", task="path", k = 2*len(ref)) 
            edit_distance2 = edlib_ed(ref, read, mode="HW", task="path", k = 2*len(ref)) 
            edit_distance = min(edit_distance1, edit_distance2)
            # print(acc,ref_acc, edit_distance)
            ref_len = len(ref)
            if edit_distance < best_ed:
                best_ed = edit_distance
                best_ed_isoform = ref_acc
            elif edit_distance == best_ed:
                print(best_ed_isoform, ref_acc, edit_distance, best_ed,edit_distance1, edit_distance2, "Original error")

        errors_container[acc] = (best_ed,best_ed_isoform)
        total_nucleotides += ref_len
        tot_errors += best_ed



    out_file.write("{0},{1},{2},{3},{4},{5},{6}\n".format("q_acc", "ref_acc", "total_errors", "error_rate", "abundance", "preserved", "minor_isoform"))


    minor_mut_retained = 0
    for q_acc in errors_container:
        q_acc_mod = q_acc.split("_")[0]
        # print(q_acc_mod)
        # print(q_acc_mod.split("|")[2], errors_container[q_acc][1])
        best_ref_ed, best_ref_acc = errors_container[q_acc]
        # switch = 1 if q_acc_mod !=  errors_container[q_acc][1] else 0
        true_abundance = original_abundances[q_acc_mod]
        true_isoform = q_acc_mod.split("|")[2]
        minor_isoform = 1 if 'minor' in q_acc_mod else 0   # true_abundance <= max(list(original_abundances.values())) else 0
        correct_structure = 1 if true_isoform == errors_container[q_acc][1] else 0
        if correct_structure == 1 and minor_isoform == 1:
            # print("here",correct_structure,minor_isoform )
            minor_mut_retained += 1
        out_file.write("{0},{1},{2},{3},{4},{5},{6}\n".format(q_acc, best_ref_acc, best_ref_ed, round(100*best_ref_ed/float(ref_len),4), true_abundance, correct_structure, minor_isoform )) 


    summary_file.write("{0},{1},{2},{3}\n".format(total_nucleotides, tot_errors, round(100*tot_errors/float(total_nucleotides), 3), minor_mut_retained )) 

    # print("TOTAL ERRORS:", sum([ ed for acc, ed in errors_container.items()]))


def main(args):
    corrected_reads = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.corrected_reads, 'r'))}
    reference_transcripts = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.reference_transcripts, 'r'))}
        
    get_best_match(corrected_reads, reference_transcripts, args.outfolder, args)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('corrected_reads', type=str, help='Path to the consensus fasta file')
    parser.add_argument('reference_transcripts', type=str, help='Path to the transcript fasta file')

    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
