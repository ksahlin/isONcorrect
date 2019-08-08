
#!/usr/bin/env python

from __future__ import print_function
import os,sys
import argparse
import pickle
import subprocess
import re
import pysam 
from collections import defaultdict
import errno

from networkx import nx

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd

import gffutils

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
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



def get_splice_sites(cigar_tuples, first_exon_start):
    splice_sites = []
    ref_pos = first_exon_start
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "=" or t== "D" or  t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos

            splice_sites.append( (splice_start, splice_stop) )

        elif t == "I" or t == "S": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return splice_sites


def is_same_isoform_cigar(q_isoform, ref_isoform):
    # compare cs tag at intron sites
    q_cigar = q_isoform.cigarstring
    q_start = q_isoform.reference_start
    q_end = q_isoform.reference_end
    q_cigar_tuples = []
    result = re.split(r'[=DXSMIN]+', q_cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = q_cigar[i]
        i += 1
        q_cigar_tuples.append((int(length), type_ ))

    ref_cigar = ref_isoform.cigarstring
    ref_start = ref_isoform.reference_start
    ref_end = ref_isoform.reference_end
    ref_cigar_tuples = []
    result = re.split(r'[=DXSMIN]+', ref_cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = ref_cigar[i]
        i += 1
        ref_cigar_tuples.append((int(length), type_ ))
    
    # print(q_cigar_tuples)
    # print(ref_cigar_tuples)
    
    q_splice_sites = get_splice_sites(q_cigar_tuples, q_start)
    all_q_splice_sites = set(q_splice_sites)
    ref_splice_sites = get_splice_sites(ref_cigar_tuples, ref_start) 

    if len(all_q_splice_sites) != len(ref_splice_sites):
        return False
    for r_start, r_stop in ref_splice_sites:
        if (r_start, r_stop) not in all_q_splice_sites:
            return False

    return True

def read_fasta_modify_header(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            accession = accession.replace("(","")
            accession = accession.replace(")","")
            accession = accession.replace(",","")

            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip().split()[0]
            accession = accession.replace("(","")
            accession = accession.replace(")","")
            accession = accession.replace(",","")
        else:
            temp += line.strip()
    yield accession, temp

def cigar_to_quality(q_isoform):
    cigarstring = q_isoform.cigarstring
    cigar_tuples = []
    result = re.split(r'[=DXSMIN]+', cigarstring)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigarstring[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))
    
    unaligned = 0
    difference = 0
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "=" or  t== "M":
            pass
        elif t== "D":
            difference += l
        elif t == "X":
            difference += l
        elif t == "N":
            pass
        elif t == "I":
            unaligned += l
            difference += l        
        elif t == "S": 
            unaligned += l
            difference += l

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()
    
    q_length = float(q_isoform.infer_query_length())

    alignment_id = (q_length - difference)/ q_length
    alignment_coverage = (q_length - unaligned)/q_length

    return alignment_id, alignment_coverage     

def pass_quality_filters(q_isoform):

    ALIGNMENT_START = (90640000,90650000)
    ALIGNMENT_END = (90750000,90760000)
    # check if all our transcripts fulfill the identity and coverage.
    ALIGN_COVERAGE = 0.99
    ALIGN_IDENTITY = 0.95
    if q_isoform.reference_name != "chr4":
        # print("Wrong chromosome", q_isoform.reference_name)
        return False

    alignment_id, alignment_coverage = cigar_to_quality(q_isoform)
    if alignment_id < ALIGN_IDENTITY:
        # print(q_isoform.query_name, "Below identity", alignment_id)  
        return False

    if alignment_coverage < ALIGN_COVERAGE:
        # print(q_isoform.query_name, "Below coverage", alignment_coverage)  
        return False

    if not (ALIGNMENT_START[0] <= q_isoform.reference_start <= ALIGNMENT_START[1]):
        # print(q_isoform.query_name, "Bad start", q_isoform.reference_start)  
        return False

    if not (ALIGNMENT_END[0] <= q_isoform.reference_end <= ALIGNMENT_END[1]):
        # print(q_isoform.query_name, "Bad end", q_isoform.reference_start) 
        return False

    return True



def get_query_splice_sites(q_isoform):
    # compare cs tag at intron sites
    q_cigar = q_isoform.cigarstring
    q_start = q_isoform.reference_start
    q_end = q_isoform.reference_end
    q_cigar_tuples = []
    result = re.split(r'[=DXSMIN]+', q_cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = q_cigar[i]
        i += 1
        q_cigar_tuples.append((int(length), type_ ))    
    q_splice_sites = get_splice_sites(q_cigar_tuples, q_start)
    return q_splice_sites



def detect_isoforms(ref_gff_file, pred_samfile_path, query_fasta, outfolder):
    fn = gffutils.example_filename(ref_gff_file)
    db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    db = gffutils.FeatureDB('test.db', keep_order=True)
    gene = db["PB.1016"]
    # print(transcript)
    ref_isoforms = {}
    for tr in db.children(gene, featuretype='transcript', order_by='start'):
        # print(tr.id, dir(tr)) 
        splice_sites = []
        for j in db.children(tr, featuretype='exon', order_by='start'):
            # print(j, j.start, j.end)
            splice_sites.append(j.start -1)
            splice_sites.append(j.end)
        splice_sites_tmp = splice_sites[1:-1]
        splice_sites = []
        for i in range(0, len(splice_sites_tmp),2):
            splice_sites.append( (splice_sites_tmp[i], splice_sites_tmp[i+1]) )
        # splice_sites = [item for item in zip(splice_sites[:-1], splice_sites[1:])]
        ref_isoforms[tr.id] = splice_sites

    # print(ref_isoforms)

    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)

    query_isoforms = {}
    prev_query = ""
    for q_isoform in pred_samfile.fetch(until_eof=True):
        if q_isoform.query_name != prev_query:
            #apply quality filters such as: is aligned to chrX between start and end, has > 95% id and >=99% aligned bases
            if pass_quality_filters(q_isoform):
                # query_isoforms.append(q_isoform)
                splice_sites = get_query_splice_sites(q_isoform)
                query_isoforms[q_isoform.query_name] = splice_sites

            prev_query = q_isoform.query_name
    print(len(ref_isoforms), len(query_isoforms))
    ref_splices = {  tuple(sp) : acc for acc, sp in ref_isoforms.items()}
    print(len(ref_splices))
    query_splices = {  tuple(sp) : [acc] for acc, sp in query_isoforms.items()}
    for acc, sp in query_isoforms.items():
        if acc not in query_splices[tuple(sp)]:
            query_splices[tuple(sp)].append(acc)

    print(len(query_splices))

    # print(query_isoforms)
    cnt = 0
    cnt_novel = 0
    out = open(os.path.join(outfolder, "transcripts_novel_junctions.fa"), "w")
    for q_sp in query_splices:
        if q_sp in ref_splices:
            print("found")
            cnt += 1
            pass
        else:
            print("Novel", q_sp, query_splices[q_sp])
            cnt_novel +=1
            for acc in query_splices[q_sp]:
                out.write(">{0}\n{1}\n".format(acc, query_fasta[acc]))
    out.close()
    # print(query_splices)
    print("Confirmed", cnt)
    print("Novel:", cnt_novel)
    # print( set(query_splices) & set(ref_splices))
    # print(sorted(ref_splices))
    # print(sorted(query_splices))

    sys.exit()


    counter_old = 0
    counter_new = 0
    ref_to_queries = { ref.query_name : set() for ref in ref_isoforms }
    queries_to_ref = { query.query_name : set() for query in query_isoforms }
    new_isoforms = set()

    for q_isoform in query_isoforms:
        is_new = True
        for ref_isoform in ref_isoforms:

            if is_same_isoform_cigar(q_isoform, ref_isoform) and is_same_isoform_cigar(ref_isoform, q_isoform): # defined as having identical splice sites throughout the alignment
                # print("YO")
                # print(q_isoform.query_name)
                # print(queries_to_ref)
                # print(queries_to_ref[q_isoform.query_name])
                queries_to_ref[q_isoform.query_name].add(ref_isoform.query_name)
                if len(queries_to_ref[q_isoform.query_name]) > 1:
                    print("More than 1 ref")
                    print("Same", q_isoform.query_name, queries_to_ref[q_isoform.query_name] )
                is_new = False
                counter_old += 1
            
        if is_new:
            counter_new += 1
            print("New", q_isoform.query_name, q_isoform.cigarstring)
            new_isoforms.add(q_isoform.query_name)
            queries_to_ref[q_isoform.query_name] = ""

        else:
            assert len(queries_to_ref[q_isoform.query_name]) == 1
            ref = queries_to_ref[q_isoform.query_name].pop()
            queries_to_ref[q_isoform.query_name] = ref

            ref_to_queries[ref].add(q_isoform.query_name)


    print([ len(ref_to_queries[r]) for r in ref_to_queries] )
    for r in sorted(ref_to_queries):
        print(len(ref_to_queries[r]), r)
    total_predictions = len(query_isoforms)
    print(total_predictions, "Total predictions")
    print(counter_old, "predictions had the same isoform structure as ref")
    print(counter_new, "predictions had new isoform structure to ref")
    
    ref_isoform_dict = { r.query_name : r for r in ref_isoforms}

    return queries_to_ref, new_isoforms, query_isoforms, ref_isoform_dict


def group_novel_isoforms(new_isoforms, all_filter_passing_query_isoforms, pred_samfile_path):
    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)
    query_new_isoforms = [q_isoform for q_isoform in all_filter_passing_query_isoforms if q_isoform.query_name in new_isoforms]
    G = nx.Graph()
    for n in new_isoforms:
        G.add_node(n)
    print("nr new:", len(query_new_isoforms))
    for i1 in query_new_isoforms:
        for i2 in query_new_isoforms:
            if i1.query_name == i2.query_name:
                continue
            else:
                if is_same_isoform_cigar(i1, i2) and is_same_isoform_cigar(i2, i1):
                    G.add_edge(i1.query_name, i2.query_name)

    print(len(list(G.nodes())))
    print(len(list(G.edges())))
    maximal_cliques = [cl for cl in nx.find_cliques(G)]


    print([ len(cl) for cl in maximal_cliques] )
    print(sum([ len(cl) for cl in maximal_cliques]) )
    print(len([ len(cl) for cl in maximal_cliques]), "unique splice sites isoforms")

    queries_to_new = { q_acc :  "new_isoform_"+str(i)  for i, cl in enumerate(sorted(maximal_cliques, key=len)) for q_acc in cl }
    return queries_to_new


def main(args):
    query_fasta = {acc : seq for acc, seq, _ in readfq(open(args.query_fasta,"r"))} 
    # ref_fasta = {acc : seq for acc, seq in read_fasta_modify_header(open(args.ref_fasta,"r"))} 

    # outfile = open(os.path.join(args.outfolder, args.prefix + ".fa"), "w")
    # outfile.close()
    queries_to_ref, new_isoforms, all_filter_passing_query_isoforms, ref_isoforms = detect_isoforms(args.ref_gff_file, args.querysamfile, query_fasta, args.outfolder)
    queries_to_new = group_novel_isoforms(new_isoforms, all_filter_passing_query_isoforms, args.querysamfile)
    # new_isoform_tags = get_novelty_feature(new_isoforms, args.querysamfile, args.ref_gff_file, outfile)


    # sam_to_alignment_fasta(queries_to_ref, queries_to_new, all_filter_passing_query_isoforms, ref_isoforms, args.outfolder, query_fasta, ref_fasta)

    # outfile = open( os.path.join(args.outfolder, args.prefix + ".tsv" ), "w" )
    # for q_acc in queries_to_ref:
    #     ref = queries_to_ref[q_acc]
    #     if ref:
    #         outfile.write("{0}\t{1}\t{2}\n".format(q_acc, ref, args.prefix) )
    # outfile.close()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('ref_gff_file', type=str, help='Samfile.')
    parser.add_argument('querysamfile', type=str, help='Samfile.')
    # parser.add_argument('predictions', type=str, help='Fasta file with only filtered isoform hits to FMR region (output of "filter_hits_on_hg19" script).')
    parser.add_argument('--outfolder', type=str, help='outfolder.')  
    # parser.add_argument('prefix', type=str, help='prefix to outfile.') 
    parser.add_argument('--query_fasta', type=str, help='fasta file with queries.')
    # parser.add_argument('--ref_fasta', type=str, help='fasta file with references.')


    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)

    







