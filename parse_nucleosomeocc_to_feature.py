

def read_nucleosome_occ(occ_file, nucleosome_occ_map):
    # gene_expression.A549.biorep1.tsv
    # Sequence        Position        P start P occupied
    # CHR10:60000     0       0.00279 0.00279
    open_occ_file = open(occ_file)
    for i in open_occ_file:
        i = i.rstrip()
        if not i.startswith("CHR"):
            continue
        splitted = i.split("\t")
        base_pos = splitted[0]
        pos_offset = int(splitted[1])
        p_occ = float(splitted[3])
        base_pos_splitted = base_pos.split(":")
        chr = base_pos_splitted[0]
        base_pos_site = int(base_pos_splitted[1])
        actual_pos = base_pos_site + pos_offset
        loci_key = chr.lower() + "_" + str(actual_pos)
        nucleosome_occ_map[loci_key] = p_occ


def calculate_bin_nucleosome_occ(bin_file, nucleosome_occ_map):
    # chr10   600     800
    # chr10   650     850
    open_bin_file = open(bin_file)
    for i in open_bin_file:
        i = i.rstrip()
        if not i.startswith("chr"):
            continue
        splitted = i.split("\t")
        chr = splitted[0]
        pos_start = int(splitted[1])
        pos_end = int(splitted[2])
        p_occ = 0
        for k in range(pos_start, pos_end):
            loci_key = chr.lower() + "_" + str(k)
            if loci_key in nucleosome_occ_map:
                p_occ += nucleosome_occ_map[loci_key]
        p_occ_avg = p_occ/float(200)
        print(chr + "\t" + str(pos_start) + "\t" + str(pos_end) + "\t" + str(p_occ_avg))


def main():
    occ_file = "/home/callsobing/projects/encode_dream/process/annotation/nucleosomeOcc/merge.all.tab"
    bin_file = "/home/encode_dream/data/annotations/train_regions.blacklistfiltered.bed"
    nucleosome_occ_map = {}
    read_nucleosome_occ(occ_file, nucleosome_occ_map)
    calculate_bin_nucleosome_occ(bin_file, nucleosome_occ_map)


if __name__ == '__main__':
    main()