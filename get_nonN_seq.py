# /home/encode_dream/data/annotations/test_regions.blacklistfiltered.merged.fa
import os
import sys


def parse_seq(seq_path):
    # DNASE.A549.conservative.narrowPeak
    dna_bp = ["a","t","c","g"]
    open_seq = open(seq_path)
    chr = ""
    start = 0
    for k in open_seq:
        k = k.rstrip()
        if k.startswith(">"):
            splitted = k.split(":")
            chr = splitted[0].lstrip(">")
            pos_range = splitted[1]
            start = int(pos_range.split("-")[0])
        else:
            new_start = True
            seq_str = ""
            for bp in k:
                if "n" in bp.lower():
                    new_start = True
                    start += 1
                if bp.lower() in dna_bp:
                    if new_start:
                        seq_str = ""
                        sys.stdout.write("\n>" + chr + ":" + str(start) + "\n" + bp)
                        new_start = False
                    else:
                        sys.stdout.write(bp)
                    start += 1


def main():
    seq_path = "/home/encode_dream/data/annotations/test_regions.blacklistfiltered.merged.fa"
    parse_seq(seq_path)


if __name__ == '__main__':
    main()