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

def plot_heatmap(file_name, data_frame_masked, outfolder, annotation = True):
    out_file1 = os.path.join(outfolder, file_name + ".eps")
    out_file2 = os.path.join(outfolder, file_name + ".png")
    out_file3 = os.path.join(outfolder, file_name + ".pdf")
    ###################
    # import the data directly into a pandas dataframe
    # nba = pd.read_csv("http://datasets.flowingdata.com/ppg2008.csv", index_col='Name  ')
    # remove index title
    # nba.index.name = ""
    # normalize data columns
    # nba_norm = (nba - nba.mean()) / (nba.max() - nba.min())
    # relabel columns
    # labels = ['Games', 'Minutes', 'Points', 'Field goals made', 'Field goal attempts', 'Field goal percentage', 'Free throws made',
    #           'Free throws attempts', 'Free throws percentage','Three-pointers made', 'Three-point attempt', 'Three-point percentage',
    #           'Offensive rebounds', 'Defensive rebounds', 'Total rebounds', 'Assists', 'Steals', 'Blocks', 'Turnover', 'Personal foul']
    # nba_norm.columns = labels
    # set appropriate font and dpi
    sns.set(font_scale=1.2)
    sns.set_style({"savefig.dpi": 100})
    # plot it out
    ax = sns.heatmap(data_frame_masked, cmap=plt.cm.Blues, linewidths=.5, annot=annotation)
    # set the x-axis labels on the top
    ax.xaxis.tick_top()
    # rotate the x-axis labels
    plt.xticks(rotation=90)
    # get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
    fig = ax.get_figure()
    # specify dimensions and save
    fig.set_size_inches(25, 20)
    fig.savefig(out_file1)
    fig.savefig(out_file2)
    fig.savefig(out_file3)
    fig.clf()
    #############################  


def reference_abundances(sampled_transcripts, outfolder, params):
    seqs_seen = set()
    transcript_abundances = defaultdict(int)
    for acc, seq in sorted(sampled_transcripts.items(), key=lambda x: len(x[1])):
        acc_mod = acc.split("_")[0]
        transcript_abundances[acc_mod] += 1

    # sorted_reference_tuples = [ (acc.split("_")[0], transcript_abundances[acc.split("_")[0]]) for acc, seq in sorted(sampled_transcripts.items(), key = lambda x: len(x[1]))]
    sorted_reference_tuples = [ (acc, ab) for acc, ab in transcript_abundances.items()]
    reference_abundances = [0]*len(sorted_reference_tuples)
    for i, (acc1,ab1) in enumerate(sorted_reference_tuples):
        reference_abundances[i] = [0]*len(sorted_reference_tuples)
        for j, (acc2,ab2) in enumerate(sorted_reference_tuples):
            reference_abundances[i][j] = float(ab1)/ab2

    relative_abundance_matrix_data_frame = pd.DataFrame(reference_abundances)
    # msk = relative_abundance_matrix_data_frame > 99
    # relative_abundance_matrix_data_frame_masked = relative_abundance_matrix_data_frame.mask(msk)
    plot_heatmap("relative_abundance", relative_abundance_matrix_data_frame, outfolder)
    # plot_heatmap("relative_abundance", relative_abundance_matrix_data_frame, annotation=True)



def reference_similarity(reference_transcripts, outfolder, params):
    """
        Stats about reference transcripts
    """

    reference_similarities = []
    if not params.no_ref_sim:

        print("calculating reference similarities")
        # do SW
        sorted_reference_tuples = sorted(reference_transcripts.items(), key = lambda x: len(x[1]))
        reference_similarities = [0]*len(sorted_reference_tuples)
        for i, (q_acc, q_seq) in enumerate(sorted_reference_tuples):
            reference_similarities[i] = [0]*len(sorted_reference_tuples)
            # print("ref", i)
            for j, (r_acc, r_seq) in enumerate(sorted_reference_tuples):
                # if len(q_seq) - len(r_seq) > 99 or len(q_seq) - len(r_seq) < 99:
                #     reference_similarities[q_acc][r_acc] = min(len(q_seq), len(r_seq))
                #     continue

                ed = edlib_ed(q_seq, r_seq, mode="NW", task="distance", k=2*max(len(q_seq), len(r_seq) ) )
                reference_similarities[i][j] = ed

                # r_aligned, q_aligned, stats = ssw_alignment( q_acc, r_acc, q_seq, r_seq, i,j, max_discrepancy = 10000 )
                # if stats:
                #     matches, mismatches, indels, deletions, insertions = stats
                #     errors = mismatches + indels
                #     identity = errors/ float(errors + matches)
                #     reference_similarities[i][j] = errors
                #     # print(ed, errors)
                # else:
                #     reference_similarities[i][j] = min(len(q_seq), len(r_seq))



        ref_sim_data_frame = pd.DataFrame(reference_similarities)
        msk = ref_sim_data_frame > 99
        ref_sim_data_frame_masked = ref_sim_data_frame.mask(msk)
        plot_heatmap("similarities", ref_sim_data_frame_masked, outfolder)

    return reference_similarities

import edlib

def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
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


def parasail_alignment(read, reference, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -2, opening_penalty = 2, gap_ext = 1, ends_discrepancy_threshold = 0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(read, reference, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(read, reference, opening_penalty, gap_ext, user_matrix)
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    
    read_alignment, ref_alignment = cigar_to_seq(cigar_string, read, reference)
    return read_alignment, ref_alignment

def get_alignments(best_edit_distances, querys, targets):
    # score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    # aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
    best_edit_distances_parasail = {}
    for acc1 in best_edit_distances:
        seq1 = querys[acc1]
        best_ed = len(seq1)
        best_edit_distances_parasail[acc1] = {}

        for acc2 in best_edit_distances[acc1]:
            seq2 = targets[acc2]            
            read_alignment, ref_alignment = parasail_alignment(seq1, seq2)
            insertions = ref_alignment.count("-")
            deletions = read_alignment.count("-")
            indels =  insertions + deletions
            mismatches = len([1 for n1, n2 in zip(read_alignment, ref_alignment) if n1 != n2 and n1 != "-" and n2 != "-"] )

            insertions_minus_ends = ref_alignment[20:-20].count("-")
            deletions_minus_ends = read_alignment[20:-20].count("-")
            indels_minus_ends =  insertions_minus_ends + deletions_minus_ends
            mismatches_minus_ends = len([1 for n1, n2 in zip(read_alignment[20:-20], ref_alignment[20:-20]) if n1 != n2 and n1 != "-" and n2 != "-"] )


            sw_ed = mismatches + indels
            if sw_ed < best_ed:
                best_edit_distances_parasail[acc1] = {}
                best_edit_distances_parasail[acc1][acc2] = (deletions, insertions, mismatches)
                best_ed = sw_ed

            elif sw_ed == best_ed:
                best_edit_distances_parasail[acc1][acc2] = (deletions, insertions, mismatches)

            # seq1_aln, match_line, seq2_aln = result.alignment

    return best_edit_distances_parasail


def get_best_match(corrected_reads, reference_transcripts, outfolder, params):
    out_file = open(os.path.join(outfolder, "results.csv"), "w")
    summary_file = open(os.path.join(outfolder, "summary.csv"), "w")
    #aligner = ssw.Aligner(gap_open=2, gap_extend=1)
    # do SW
    nr_unique_refs = len(reference_transcripts)
    errors_container = {}
    error_rate_container = {}
    error_types_container = {}
    best_match_container = {}

    original_abundances = defaultdict(int)
    for acc in corrected_reads.keys():
        original_abundances[ acc.split("_")[0] ] += 1 # true origin transcript is present in header

    corrected_read_abundances = {}    
    for acc in reference_transcripts:
        corrected_read_abundances[ acc.split("_")[0]]  = 0

    # print(len(original_abundances))
    # print(len(corrected_read_abundances))
    # print(set(corrected_read_abundances) ^ set(original_abundances) )
    # sys.exit()
    not_FN = set()
    if params.deal_with_ties:
        mutation_present = {}

    sorted_lengths = sorted([(len(q_seq), q_acc) for q_acc, q_seq in corrected_reads.items()])
    # for l in sorted_lengths:
    #     print(l)

    # print("REF LENGHTS")
    sorted_lengths = sorted([len(r_seq) for r_acc, r_seq in reference_transcripts.items()])
    # for l in sorted_lengths:
    #     print(l)
    # pre check exact matches:
    if params.only_exact:
        ref_seq_to_acc = {seq : acc for acc, seq in reference_transcripts.items()}
        ref_seqs = set(reference_transcripts.values())
        exact_matches = set()
        for q_acc, q_seq in corrected_reads.items():
            if q_seq in ref_seqs:
                exact_matches.add(q_acc)
                ref_acc = ref_seq_to_acc[q_seq] #.split("copy")[0]
                errors_container[q_acc] = 0
                best_match_container[q_acc] = ref_acc
                error_rate_container[q_acc] = 1.0
                error_types_container[q_acc] = (0, 0, 0)
                not_FN.add(ref_acc)

        # print(len(ref_seqs))
        # print(len(corrected_reads))
        # print("EXACT MATCHES:", len(exact_matches))

    else:
        # print("Start1")
        best_edit_distances = get_minimizers_2set_simple(corrected_reads, reference_transcripts)
        minimizer_graph_c_to_t = get_alignments(best_edit_distances, corrected_reads, reference_transcripts)
        for i, (q_acc, q_seq) in enumerate(minimizer_graph_c_to_t.items()): 
            best_ed = 200000
            r_acc_max_id = "NONE"
            fewest_errors = len(q_seq)
            best_mismatches, best_insertions, best_deletions = len(q_seq), len(q_seq), len(q_seq)
            q_acc_mod = q_acc.split("_")[0]

            if len(minimizer_graph_c_to_t[q_acc]) > 1:
                if params.deal_with_ties:   
                    mutation_present[q_acc] = 0

                print("TIE!!", q_acc, minimizer_graph_c_to_t[q_acc])
                if q_acc_mod in minimizer_graph_c_to_t[q_acc]:
                    tmp = minimizer_graph_c_to_t[q_acc][q_acc_mod]
                    minimizer_graph_c_to_t[q_acc] = {}
                    minimizer_graph_c_to_t[q_acc] = { q_acc_mod : tmp}
                    print("Done")

            # tie = False
            else:
                if q_acc_mod == list(minimizer_graph_c_to_t[q_acc].keys())[0]:
                    # print("correct")
                    pass
                else:
                    # print("Wrong", q_acc, minimizer_graph_c_to_t[q_acc])
                    pass
                    # continue
                # print(q_acc, minimizer_graph_c_to_t[q_acc])

            for j, (r_acc, r_seq) in enumerate(minimizer_graph_c_to_t[q_acc].items()):
                deletions, insertions, mismatches = minimizer_graph_c_to_t[q_acc][r_acc]
                edit_distance =  deletions + insertions + mismatches
                if edit_distance < best_ed:
                    best_ed = edit_distance
                    r_acc_max_id = r_acc
                    fewest_errors = edit_distance
                    best_mismatches, best_insertions, best_deletions = mismatches, insertions, deletions  


            errors_container[q_acc] = fewest_errors
            best_match_container[q_acc] = r_acc_max_id

            if params.deal_with_ties: 
                if q_acc not in mutation_present:
                    mutation_present[q_acc] = 1 if q_acc_mod == r_acc_max_id else 0

            
            corrected_read_abundances[ r_acc ] += 1

            # print(q_acc, q_acc.split("_")[-1], r_acc_max_id)

            error_rate_container[q_acc] = (best_ed / float(len(reference_transcripts[r_acc_max_id])))
            error_types_container[q_acc] = (best_mismatches, best_insertions, best_deletions)
            not_FN.add(r_acc_max_id)

        print("Stop1!")


    # total discoveries, total perfect matches (1.0 identity), errors for each consensus
    # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(nr_unique_refs, q_acc, best_match_container[q_acc], errors_container[q_acc], error_rate_container[q_acc], *error_types_container[q_acc]))
    print("errors", sorted([ ed for acc, ed in errors_container.items()]))
    print("TOTAL ERRORS:", sum([ ed for acc, ed in errors_container.items()]))

    tot_errors = 0
    all_s = 0
    all_i = 0
    all_d = 0
    total_read_nucleotides = 0
    for q_acc, q_seq in corrected_reads.items():
        q_acc_mod = q_acc.split("_")[0]
        # if q_acc_mod == best_match_container[q_acc]:
        tot_errors += errors_container[q_acc]
        all_s += error_types_container[q_acc][0]
        all_i += error_types_container[q_acc][1]
        all_d += error_types_container[q_acc][2]
        total_read_nucleotides += len(q_seq)


    summary_file.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(total_read_nucleotides, tot_errors, all_s, all_i, all_d, round(100*tot_errors/float(total_read_nucleotides), 3),
                    round(100*all_s/float(total_read_nucleotides), 3), round(100*all_i/float(total_read_nucleotides), 3), round(100*all_d/float(total_read_nucleotides), 3) )) 
    
    sorted_original_abundances = sorted(original_abundances.items(), key = lambda x: x[1], reverse=True)
    print("Corrected abundances:", [corrected_read_abundances[acc] for acc, ab in sorted_original_abundances])
    print("Original abundances:", [ab for acc, ab in sorted_original_abundances])
    # for acc, ab in sorted_original_abundances:
    #     summary_file.write("{0}\t{1}\t{2}\n".format(acc, ab, corrected_read_abundances[acc]))

    if params.deal_with_ties: 
        out_file.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format("q_acc", "ref_acc", "total_errors", "error_rate", "subs", "ins", "del", "switch", "abundance", "mutation_present", "minor"))
    else:
        out_file.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format("q_acc", "ref_acc", "total_errors", "error_rate", "subs", "ins", "del", "switch", "abundance"))

    for q_acc in errors_container:
        q_acc_mod = q_acc.split("_")[0]
        switch = 1 if q_acc_mod !=  best_match_container[q_acc] else 0
        true_abundance = original_abundances[q_acc_mod]

        if params.deal_with_ties: 
            mut_present = mutation_present[q_acc]
            minor = 1 if true_abundance != max(list(original_abundances.values())) else 0
            out_file.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format(q_acc, best_match_container[q_acc], errors_container[q_acc], round(100*error_rate_container[q_acc],4), *error_types_container[q_acc], switch, true_abundance, mut_present, minor )) 
        else:
            out_file.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(q_acc, best_match_container[q_acc], errors_container[q_acc], round(100*error_rate_container[q_acc],4), *error_types_container[q_acc], switch, true_abundance )) 

    print("TOTAL ERRORS:", sum([ ed for acc, ed in errors_container.items()]))


def main(args):
    corrected_reads = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.corrected_reads, 'r'))}
    reference_transcripts = {acc: seq.upper() for acc, (seq, _) in  readfq(open(args.reference_transcripts, 'r'))}
    
    reference_similarities = reference_similarity(reference_transcripts, args.outfolder, args)
    reference_abundances(corrected_reads, args.outfolder, args)
    
    get_best_match(corrected_reads, reference_transcripts, args.outfolder, args)

    out_file_ref_sim = open(os.path.join(args.outfolder, "ref_similaritiy_distr.csv"), "w")

    if not args.no_ref_sim:
        ref_similarities = []
        for i in range(len(reference_similarities)):
            temp_list = reference_similarities[i]
            del temp_list[i]
            min_sim_to_i = min(temp_list) 
            ref_similarities.append(min_sim_to_i)
            # for j in range(len(reference_similarities[i])):
            #     if i <= j:
            #         continue
            #     else:
            #         ref_similarities.append(reference_similarities[i][j])

        n = float( len(ref_similarities) )
        mu =  sum(ref_similarities) / n
        sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), ref_similarities))) / (n - 1)) ** 0.5
        min_error = min(ref_similarities)
        max_error = max(ref_similarities)
        ref_similarities.sort()
        if len(ref_similarities) %2 == 0:
            median_error = (ref_similarities[int(len(ref_similarities)/2)-1] + ref_similarities[int(len(ref_similarities)/2)]) / 2.0
        else:
            median_error = ref_similarities[int(len(ref_similarities)/2)]

        out_file_ref_sim.write("mean_error\tsd_error\tmin_error\tmax_error\tmedian_error\n")
        out_file_ref_sim.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(mu, sigma, min_error, max_error, median_error))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('corrected_reads', type=str, help='Path to the consensus fasta file')
    parser.add_argument('reference_transcripts', type=str, help='Path to the transcript fasta file')
    parser.add_argument('--deal_with_ties',  action='store_true', help='Specific parameter for mutation analysis')

    parser.add_argument('--transcripts_sampled', type=str, help='')
    parser.add_argument('--only_exact',  action='store_true', help='Only output number of exact matches')
    parser.add_argument('--no_ref_sim',  action='store_true', help='Do not compute reference identy levels')

    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
