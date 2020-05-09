import os
import sys

import csv
import argparse
from collections import defaultdict

def check_same_region(corr_positions,orig_positions, args):
    if args.sirv_iso_cov:
        outfile = open(args.outfile_csv, "w")
        outfile.write("nr_reads,nr_isoforms,is_overcorrected\n")
    cnt_diff, cnt_unaligned = 0, 0
    cnt_sirv_5 = 0
    nr_reads = defaultdict(lambda: defaultdict(int)) #{1:0, 4:0, 8:0}
    nr_isoform_switches = defaultdict(lambda: defaultdict(int)) #{1:0, 4:0, 8:0}

    for read_id in orig_positions:
        (o_id, o_start,o_stop, _) = orig_positions[read_id]
        if read_id in corr_positions:
            (c_id, c_start,c_stop, _) = corr_positions[read_id]
        else:
            cnt_unaligned += 1
            continue
        # print(c_start, c_stop)
        on_same_ref = False
        if args.sirv_iso_cov:
            n_isoforms = int(read_id.split(",")[2])
            read_depth = int(read_id.split(",")[1])
            c_all_id =set( c_id.split(":"))
            o_all_id =set( o_id.split(":"))
            on_same_ref = True if len(c_all_id & o_all_id) > 0 else False
            if not on_same_ref:
                nr_isoform_switches[read_depth][n_isoforms] += 1
                if "SIRV506" in (c_all_id | o_all_id) and "SIRV511" in (c_all_id | o_all_id):
                    cnt_sirv_5 +=1
            nr_reads[read_depth][n_isoforms] += 1

            outfile.write("{0},{1},{2}\n".format(read_depth, n_isoforms, 0 if on_same_ref else 1))

        else:
            on_same_ref = (c_id == o_id)

        if on_same_ref and ( o_start <= c_stop <= o_stop  or o_start <= c_start <= o_stop ): # intercecting
            pass
        else:
            cnt_diff += 1
            # print(c_id,o_id, read_id)
    
    if args.sirv_iso_cov:
        outfile.close()
        print("Depth,1,4,8")
        tot_1, tot_4, tot_8 = 0,0,0
        switch_1, switch_4, switch_8 = 0,0,0
        for cov in nr_reads:
            perc_overcorr_per_depth = []
            for n_iso in [1,4,8]:
                nr_overcorr = nr_isoform_switches[cov][n_iso] if nr_isoform_switches[cov][n_iso] else 0
                n_reads_in_batch = nr_reads[cov][n_iso]
                perc_overcorr_per_depth.append(round(100*nr_overcorr/n_reads_in_batch,1))
                if n_iso == 1:
                    switch_1 += nr_overcorr
                    tot_1 += n_reads_in_batch
                elif n_iso == 4:
                    switch_4 += nr_overcorr
                    tot_4 += n_reads_in_batch
                elif n_iso == 8:
                    switch_8 += nr_overcorr
                    tot_8 += n_reads_in_batch

            print("{0},{1}".format(cov, ",".join([str(p) for p in perc_overcorr_per_depth]) ))

        print("SIRV506 - SIRV511", cnt_sirv_5)
        print("Switching isoforms. Frac overcorrected per n_iso (1,4,8 isoforms): ", round(100*switch_1/float(tot_1), 2), round(100*switch_4/float(tot_4), 2), round(100*switch_8/float(tot_8), 2) )
        print("Number of reads per experiment (1,4,8 isoforms): ", nr_reads)

    return cnt_diff, cnt_unaligned

def check_err_rate(corr_positions,orig_positions):
    cnt_broken = 0
    for read_id in orig_positions:
        (_, _, _, o_err_rate) = orig_positions[read_id]
        if read_id in corr_positions:
            (_, _, _, c_err_rate) = corr_positions[read_id]
            # print(o_err_rate, c_err_rate )
            if c_err_rate > o_err_rate:
                cnt_broken += 1
    return cnt_broken

def main(args):

    corr_positions = {}
    orig_positions = {}
    with open(args.csv, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            if i == 0: # skip header
                continue

            read_id = row[0]

            if args.dros:
                read_type = row[1]
                err_rate = float(row[6])
                chr_id = row[22]
                ref_start = int(row[23])
                ref_stop = int(row[24])
            if args.sirv:
                read_type = row[1]
                err_rate = float(row[6])
                chr_id = row[9]
                ref_start = 0
                ref_stop = 0
            if args.sim:
                read_type = row[8]
                err_rate = float(row[7])
                chr_id = "known"
                ref_start = 0
                ref_stop = 0
            if args.sirv_iso_cov:
                read_id = ",".join(row[:4]) #  row[3] 
                read_type = row[4]
                err_rate = float(row[9])
                chr_id = row[12] # "".join(sorted( [s for s in row[12].split(":")]))
                ref_start = 0
                ref_stop = 0
            # print(chr_id,ref_start, ref_stop)
            if read_type == 'corrected':
                corr_positions[read_id] = (chr_id, ref_start, ref_stop, err_rate)
            else:
                orig_positions[read_id] = (chr_id, ref_start, ref_stop, err_rate)

    cnt_diff, cnt_unaligned = check_same_region(corr_positions,orig_positions, args)
    cnt_broken = check_err_rate(corr_positions,orig_positions)
    print("Total reads change location:", cnt_diff)
    print("Total reads aligned in origina format but unaligned after correction:", cnt_unaligned)
    print("Total reads higher error rate after:", cnt_broken)
    print("Total reads original:", len(orig_positions))
    print("Total reads corrected:", len(corr_positions))


    print("Total aligned orig:", len(orig_positions) )
    print("Total aligned corr:", len(corr_positions) )

    cor_errors = sorted([e for _,_,_,e in corr_positions.values()])
    orig_errors = sorted([e for _,_,_,e in orig_positions.values()])

    print("Median error rate corr:", cor_errors[ int(len(cor_errors)/2)] )
    print("Median error rate orig:", orig_errors[ int(len(orig_errors)/2)] )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot p-minimizers shared.")
    parser.add_argument('csv', type=str, help='Path to CSV input.')
    parser.add_argument('--dros', action="store_true", help='Dros.')
    parser.add_argument('--sirv', action="store_true", help='sirv.')
    parser.add_argument('--sirv_iso_cov', action="store_true", help='sirv.')
    parser.add_argument('--sim', action="store_true", help='sim.')
    parser.add_argument('--outfile_csv', type=str, help='Path to CSV output.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)

