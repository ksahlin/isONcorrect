import os
import sys

import csv
import argparse
# from collections import deque

def check_same_region(corr_positions,orig_positions):
    cnt_diff, cnt_unaligned = 0, 0
    for read_id in orig_positions:
        (o_id, o_start,o_stop, _) = orig_positions[read_id]
        if read_id in corr_positions:
            (c_id, c_start,c_stop, _) = corr_positions[read_id]
        else:
            cnt_unaligned += 1
            continue
        # print(c_start, c_stop)
        if c_id == o_id and ( o_start <= c_stop <= o_stop  or o_start <= c_start <= o_stop ): # intercecting
            pass
        else:
            cnt_diff += 1
            print(c_id,o_id, read_id)
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
                read_id = ",".join(row[:4])
                read_type = row[4]
                err_rate = float(row[9])
                chr_id = row[12]
                ref_start = 0
                ref_stop = 0
            # print(chr_id,ref_start, ref_stop)
            if read_type == 'corrected':
                corr_positions[read_id] = (chr_id, ref_start, ref_stop, err_rate)
            else:
                orig_positions[read_id] = (chr_id, ref_start, ref_stop, err_rate)

    cnt_diff, cnt_unaligned = check_same_region(corr_positions,orig_positions)
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
    parser.add_argument('csv', type=str, help='Path to fasta file with a nucleotide sequence (e.g., gene locus) to simulate isoforms from.')
    parser.add_argument('--dros', action="store_true", help='Dros.')
    parser.add_argument('--sirv', action="store_true", help='sirv.')
    parser.add_argument('--sirv_iso_cov', action="store_true", help='sirv.')
    parser.add_argument('--sim', action="store_true", help='sim.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)

