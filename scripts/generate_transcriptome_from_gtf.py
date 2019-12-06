import os, sys
import argparse
import gffutils

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
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


def write_isoforms(gff, refs, outfolder):
    db_name = os.path.join(outfolder, 'database.db')
    fn = gffutils.example_filename(gff)
    db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                            sort_attribute_values=True)
    db = gffutils.FeatureDB(db_name, keep_order=True)
    transcripts = {}
    for gene in db.features_of_type('gene'):
        ref_id = gene.seqid
        ref_seq = refs[ref_id]

        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            transcript_seq = []
            for j in db.children(transcript, featuretype='exon', order_by='start'):
                exon_seq = ref_seq[j.start -1: j.end]
                transcript_seq.append(exon_seq)
            transcript_seq = ''.join([e for e in transcript_seq])
            transcripts[transcript.id] = transcript_seq
    return transcripts


def main(args):
    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}
    transcripts = write_isoforms(args.gff, refs, args.outfolder)
    outfile = open(os.path.join(args.outfolder, "sirv_transcriptome.fa"), "w")
    for acc, seq in transcripts.items():
        outfile.write(">{0}\n{1}\n".format(acc, seq))
    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('gff', type=str, help='Samfile.')
    parser.add_argument('refs', type=str, help='Samfile.')
    parser.add_argument('outfolder', type=str, help='outfolder.')  


    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)