from __future__ import print_function
import os,sys
import argparse
import re
import ssw
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
# try:
#     import matplotlib
#     matplotlib.use('agg')
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     sns.set_palette("husl", desat=.6)
#     sns.set(font_scale=1.6)
#     plt.rcParams.update({'font.size': 22})
# except:
#     pass


# try:
#     import matplotlib
#     matplotlib.use('Agg')
#     import matplotlib.pyplot as plt
# except ImportError, RuntimeError:
#     pass


def ssw_alignment( x_acc, y_acc, x, y, i,j, max_discrepancy = 200):
    """
        Aligns two sequences with SSW
        x: query
        y: reference

    """
    if i % 10 == 0 and j % 10000 == 0:
        print("processing alignments on all y's where read x_i is participating. i={0}".format(i+1))

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-1)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)

    # for the ends that SSW leaves behind
    bio_matrix = matlist.blosum62
    g_open = -1
    g_extend = -0.5
    ######################################

    result = aligner.align(x, y, revcomp=True)
    y_alignment, match_line, x_alignment = result.alignment

    matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
    deletions = x_alignment.count("-")
    insertions = y_alignment.count("-")
    assert deletions + insertions == indels
    # alignment_length = len(match_line)
    
    start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
    query_end_discrepancy = len(x) - result.query_end - 1
    ref_end_discrepancy = len(y) - result.reference_end - 1
    end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
    # print(query_end_discrepancy, ref_end_discrepancy)
    tot_discrepancy = start_discrepancy + end_discrepancy

    if 0 < start_discrepancy <= max_discrepancy:
        # print("HERE")
        matches_snippet = 0
        mismatches_snippet = 0
        if result.query_begin and result.reference_begin:
            query_start_snippet = x[:result.query_begin]
            ref_start_snippet = y[:result.reference_begin]
            alns = pairwise2.align.globalds(query_start_snippet, ref_start_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            # print(alns)
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            # print(matches_snippet, mismatches_snippet, indels_snippet)
            query_start_alignment_snippet = top_aln[0]
            ref_start_alignment_snippet = top_aln[1]
        elif result.query_begin:
            query_start_alignment_snippet = x[:result.query_begin]
            ref_start_alignment_snippet = "-"*len(query_start_alignment_snippet)
            indels_snippet = len(ref_start_alignment_snippet)
        elif result.reference_begin:
            ref_start_alignment_snippet = y[:result.reference_begin]
            query_start_alignment_snippet = "-"*len(ref_start_alignment_snippet)
            indels_snippet = len(query_start_alignment_snippet)
        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        # print(ref_start_alignment_snippet)
        # print(query_start_alignment_snippet)
        y_alignment = ref_start_alignment_snippet + y_alignment
        x_alignment = query_start_alignment_snippet + x_alignment

    if 0 < end_discrepancy <= max_discrepancy:
        # print("HERE2", query_end_discrepancy, ref_end_discrepancy)
        # print(y_alignment)
        # print(y)
        # print(match_line)
        # print(x_alignment)
        # print(x)
        # print(matches, len(x_alignment))
        matches_snippet = 0
        mismatches_snippet = 0
        if query_end_discrepancy and ref_end_discrepancy:
            query_end_snippet = x[result.query_end+1:]
            ref_end_snippet = y[result.reference_end+1:]
            alns = pairwise2.align.globalds(query_end_snippet, ref_end_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            query_end_alignment_snippet = top_aln[0]
            ref_end_alignment_snippet = top_aln[1]
        elif query_end_discrepancy:
            query_end_alignment_snippet = x[result.query_end+1:]
            ref_end_alignment_snippet = "-"*len(query_end_alignment_snippet)
            indels_snippet = len(ref_end_alignment_snippet)

        elif ref_end_discrepancy:
            ref_end_alignment_snippet = y[result.reference_end+1:]
            query_end_alignment_snippet = "-"*len(ref_end_alignment_snippet)
            indels_snippet = len(query_end_alignment_snippet)

        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        y_alignment = y_alignment + ref_end_alignment_snippet
        x_alignment = x_alignment + query_end_alignment_snippet

    # matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
    deletions = x_alignment.count("-")
    insertions = y_alignment.count("-")
    assert deletions + insertions == indels

    if start_discrepancy > max_discrepancy or end_discrepancy > max_discrepancy:
        # print("REMOVING", start_discrepancy, end_discrepancy)
        return (y_alignment, x_alignment, None)

    else:
        return (y_alignment, x_alignment, (matches, mismatches, indels, deletions, insertions, match_line)) 


def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip()
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip()
        else:
            temp += line.strip().upper()
    if accession:
        yield accession, temp

def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
    if x_cutoff: 
        plt.hist(data, range=[0, x_cutoff], bins = 100)
    else:
        plt.hist(data, bins = 100)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()

def rand_jitter(arr):
    stdev = .003*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def dot_plot(points, args, x_label='x', y_label='y', title='Dotplot', set_marker='o',  name='dot_plot.png'):
    x = map(lambda x: x[0], points)
    y = map(lambda x: x[1], points)
    plt.scatter(rand_jitter(x), rand_jitter(y), marker=set_marker, alpha=0.3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()

def collapse(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def plot_heatmap(file_name, data_frame_masked, annotation = True):
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

def reference_similarity(reference_transcripts, outfolder, params):
    """
        Stats about reference transcripts
    """
    seqs_seen = set()
    transcript_abundances = {}
    transcript_copies = Counter()
    for acc, seq in sorted(reference_transcripts.items(), key=lambda x: len(x[1])):

        try:
            tr_acc, copy_number_str = acc.split("copy")

        except ValueError:
            tr_acc, copy_number_str = acc, "1"  # viral data not simulated

        transcript_copies[tr_acc] += 1
        try :
            copy_number = int(copy_number_str)
        except ValueError:
            copy_number = 1

        if tr_acc not in transcript_abundances:
            transcript_abundances[tr_acc] = copy_number
        elif tr_acc in transcript_abundances and copy_number >  transcript_abundances[tr_acc]:
             transcript_abundances[tr_acc] = copy_number

        if seq in seqs_seen:
            # print("HERE!", len(seq))
            del reference_transcripts[acc]
        else:
            seqs_seen.add(seq)

    print("Number of unique references:", len(reference_transcripts))
    for t_acc, copy_nr in transcript_copies.items():
        print(t_acc, copy_nr)
    # print("abundances:", transcript_abundances)

    if not params.no_ref_sim:
        # relative_abundance_matrix = {}
        # annotation_matrix = {}

        sorted_reference_tuples = sorted(reference_transcripts.items(), key = lambda x: len(x[1]))
        reference_abundances = [0]*len(sorted_reference_tuples)
        for i, (acc1,seq1) in enumerate(sorted_reference_tuples):
            reference_abundances[i] = [0]*len(sorted_reference_tuples)
            # relative_abundance_matrix[acc1] = {}
            # annotation_matrix[acc1] = {}
            for j, (acc2,seq2) in enumerate(sorted_reference_tuples):
                copy_nr_1 = transcript_abundances[acc1.split("copy")[0]]
                copy_nr_2 = transcript_abundances[acc2.split("copy")[0]]
                reference_abundances[i][j] = float(copy_nr_1)/copy_nr_2
                # relative_abundance_matrix[acc1][acc2] = float(copy_nr_1)/copy_nr_2
                # annotation_matrix[acc1][acc2] = str(fractions.Fraction(copy_nr_1, copy_nr_2))

        # print(relative_abundance_matrix)
        # print(annotation_matrix)
        relative_abundance_matrix_data_frame = pd.DataFrame(reference_abundances)
        # msk = relative_abundance_matrix_data_frame > 99
        # relative_abundance_matrix_data_frame_masked = relative_abundance_matrix_data_frame.mask(msk)

        plot_heatmap("relative_abundance", relative_abundance_matrix_data_frame)
        # plot_heatmap("relative_abundance", relative_abundance_matrix_data_frame, annotation=True)


        print("calculating reference similarities")
        #aligner = ssw.Aligner(gap_open=2, gap_extend=1)
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
        plot_heatmap("similarities", ref_sim_data_frame_masked)

    return transcript_abundances, transcript_copies, reference_similarities

import edlib

def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    return ed


def get_minimizers_2set_simple(querys, targets):
    best_edit_distances = {}

    for acc1, seq1 in querys.items():
        best_ed = len(seq1)
        for acc2, seq2 in targets.items():
            edit_distance = edlib_ed(seq1, seq2, mode="NW", k = len(seq1)) # seq1 = query, seq2 = target
            if 0 <= edit_distance < best_ed:
                best_ed = edit_distance
                best_edit_distances[acc1] = {}
                best_edit_distances[acc1][acc2] = best_ed
            elif edit_distance == best_ed:
                best_edit_distances[acc1][acc2] = best_ed
        # print(best_ed)


    return best_edit_distances

def get_ssw_alignments(best_edit_distances, querys, targets):
    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
    best_edit_distances_ssw = {}
    for acc1 in best_edit_distances:
        seq1 = querys[acc1]
        best_ed = len(seq1)
        best_edit_distances_ssw[acc1] = {}

        for acc2 in best_edit_distances[acc1]:
            seq2 = targets[acc2]
            r_aligned, q_aligned, stats = ssw_alignment( acc1, acc2, seq1, seq2, 1111,1111 )
            if stats:
                matches, mismatches, indels, deletions, insertions, match_line = stats
            else:
                continue
            # result = aligner.align(seq1, seq2, revcomp=False)
            # seq2_aln, match_line, seq1_aln = result.alignment
            # matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
            # insertion_count = seq2_aln.count("-")
            # deletion_count = seq1_aln.count("-")

            sw_ed = mismatches + indels
            if sw_ed < best_ed:
                best_edit_distances_ssw[acc1] = {}
                best_edit_distances_ssw[acc1][acc2] = (deletions, insertions, mismatches, match_line)
                best_ed = sw_ed
            elif sw_ed == best_ed:
                best_edit_distances_ssw[acc1][acc2] = (deletions, insertions, mismatches, match_line)

            # seq1_aln, match_line, seq2_aln = result.alignment

    return best_edit_distances_ssw


def get_best_match(consensus_transcripts, reference_transcripts, outfolder, transcript_abundances, transcript_copies, sampled_dict, params):
    out_file = open(os.path.join(outfolder, "results.tsv"), "w")
    #aligner = ssw.Aligner(gap_open=2, gap_extend=1)
    # do SW
    nr_unique_refs = len(reference_transcripts)
    errors_container = {}
    identity_container = {}
    error_types_container = {}
    best_match_container = {}
    not_FN = set()
    # print(consensus_transcripts)
    # if len(consensus_transcripts) == 0:
    #     out_file.write("{0}\t{1}\t{2}\n".format(nr_unique_refs, len(consensus_transcripts), ",".join([ str(a) for a in transcript_abundances.values()])) )
    #     return

    sorted_lengths = sorted([(len(q_seq), q_acc) for q_acc, q_seq in consensus_transcripts.items()])
    # for l in sorted_lengths:
    #     print(l)

    print("REF LENGHTS")
    sorted_lengths = sorted([len(r_seq) for r_acc, r_seq in reference_transcripts.items()])
    # for l in sorted_lengths:
    #     print(l)
    total_read_nucleotides = sum([len(q_seq) for q_acc, q_seq in consensus_transcripts.items()])
    # pre check exact matches:
    if params.only_exact:
        ref_seq_to_acc = {seq : acc for acc, seq in reference_transcripts.items()}
        ref_seqs = set(reference_transcripts.values())
        exact_matches = set()
        for q_acc, q_seq in consensus_transcripts.items():
            if q_seq in ref_seqs:
                exact_matches.add(q_acc)
                ref_acc = ref_seq_to_acc[q_seq] #.split("copy")[0]
                print("Exact", q_acc, "to transcript with copy number:", transcript_copies[ref_acc]) 
                errors_container[q_acc] = 0
                best_match_container[q_acc] = ref_acc
                identity_container[q_acc] = 1.0
                error_types_container[q_acc] = (0, 0, 0)
                not_FN.add(ref_acc)

        print(len(ref_seqs))
        print(len(consensus_transcripts))
        print("EXACT MATCHES:", len(exact_matches))

    else:
        print("Start1")
        best_edit_distances = get_minimizers_2set_simple(consensus_transcripts, reference_transcripts)
        minimizer_graph_c_to_t = get_ssw_alignments(best_edit_distances, consensus_transcripts, reference_transcripts)
        for i, (q_acc, q_seq) in enumerate(minimizer_graph_c_to_t.items()): 
            best_ed = 200000
            r_acc_max_id = "NONE"
            fewest_errors = len(q_seq)
            best_mismatches, best_insertions, best_deletions = len(q_seq), len(q_seq), len(q_seq)
            for j, (r_acc, r_seq) in enumerate(minimizer_graph_c_to_t[q_acc].items()):
                deletions, insertions, mismatches, match_line = minimizer_graph_c_to_t[q_acc][r_acc]
                edit_distance =  deletions + insertions + mismatches
                print(match_line)
                if edit_distance < best_ed:
                    best_ed = edit_distance
                    r_acc_max_id = r_acc
                    fewest_errors = edit_distance
                    best_mismatches, best_insertions, best_deletions = mismatches, insertions, deletions    

            errors_container[q_acc] = fewest_errors
            best_match_container[q_acc] = r_acc_max_id
            identity_container[q_acc] = 1.0 - (best_ed / float(max(len(q_seq), len(reference_transcripts[r_acc_max_id])) ))
            error_types_container[q_acc] = (best_mismatches, best_insertions, best_deletions)
            not_FN.add(r_acc_max_id)

        print("Stop1!")

    if sampled_dict:
        FN = set(sampled_dict.keys()).difference(not_FN)        
    else:
        FN = set(reference_transcripts.keys()).difference(not_FN)

    for ref in FN:
        print("FN:",ref, len(reference_transcripts[ref]) )
    # current logging:
    # first row display info of number of uniqur reference transcripts, and number of inferred transcripts
    # if sampled_dict:
    #     out_file.write("{0}\t{1}\t{2}\n".format(len(sampled_dict), len(consensus_transcripts), ",".join([ str(a) for a in transcript_abundances.values()])) )
    # else:
    #     out_file.write("{0}\t{1}\t{2}\n".format(nr_unique_refs, len(consensus_transcripts), ",".join([ str(a) for a in transcript_abundances.values()])) )

    # total discoveries, total perfect matches (1.0 identity), errors for each consensus
    # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(nr_unique_refs, q_acc, best_match_container[q_acc], errors_container[q_acc], identity_container[q_acc], *error_types_container[q_acc]))

    print("TOTAL ERRORS:", sum([ ed for acc, ed in errors_container.items()]))
    tot_errors = sum([ ed for acc, ed in errors_container.items()])
    all_errors = [error_types_container[acc] for acc in error_types_container]
    all_s = sum([s for s,i,d in all_errors])
    all_i = sum([i for s,i,d in all_errors])
    all_d = sum([d for s,i,d in all_errors])
    out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(total_read_nucleotides, tot_errors, all_s, all_i, all_d, round(100*tot_errors/float(total_read_nucleotides), 3)))
    
    out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format("q_acc", "ref_acc", "total_errors", "identity", "subs", "ins", "del"))

    for q_acc in errors_container:
        # each ro displays values for a consensus transcript
        if  identity_container[q_acc] > params.sim_cutoff:
            ssw_stats = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(q_acc, best_match_container[q_acc], errors_container[q_acc], round(identity_container[q_acc],4), *error_types_container[q_acc])
            # print(ssw_stats, minimizer_graph_c_to_t[q_acc])
            # print()
            out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(q_acc, best_match_container[q_acc], errors_container[q_acc], round(identity_container[q_acc],4), *error_types_container[q_acc]))

    print("TOTAL ERRORS:", sum([ ed for acc, ed in errors_container.items()]))


def main(args):
    consensus_transcripts = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}
    reference_transcripts = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.reference_transcripts, 'r'))}
    
    sampled_dict = {}
    if args.transcripts_sampled:
        for line in open(args.transcripts_sampled, 'r').readlines()[:-1]:
            member, nr_reads = line.strip().split()
            sampled_dict[member] = nr_reads
        print("Number of references sampled in reads:", len(sampled_dict))


    print("Number of references:", len(reference_transcripts))
    print("Number of consensus:", len(consensus_transcripts))

    # if args.only_exact:
    #     consensus_seq_to_acc = {seq.upper(): acc for (acc, seq) in  read_fasta(open(args.consensus_transcripts, 'r'))}
    #     ref_set = set(list(reference_transcripts.values()))
    #     consensus_set = set(list(consensus_transcripts.values()))
    #     # assert len(consensus_set) == len(list(consensus_transcripts.keys()))
    #     TP = consensus_set.intersection(ref_set)

    #     print("Nr unique refs: {0}, Nr unique queries: {1}".format(len(ref_set), len(consensus_set) ))
    #     print("Exact matches: {0}".format(len(TP)))
    #     sys.exit()


    transcript_abundances, transcript_copies, reference_similarities = reference_similarity(reference_transcripts, args.outfolder, args)
    get_best_match(consensus_transcripts, reference_transcripts, args.outfolder, transcript_abundances, transcript_copies, sampled_dict, args)

    out_file_ref_sim = open(os.path.join(args.outfolder, "ref_similaritiy_distr.tsv"), "w")

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
    parser.add_argument('consensus_transcripts', type=str, help='Path to the consensus fasta file')
    parser.add_argument('reference_transcripts', type=str, help='Path to the transcript fasta file')
    parser.add_argument('--transcripts_sampled', type=str, help='Path to sampled transcripts txt file')
    parser.add_argument('--only_exact',  action='store_true', help='Only output number of exact matches')
    parser.add_argument('--no_ref_sim',  action='store_true', help='Do not compute reference identy levels')
    parser.add_argument('--sim_cutoff',  type=float, default= -100.0, help='Output only hits with higher identity than (default: output all)')

    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
