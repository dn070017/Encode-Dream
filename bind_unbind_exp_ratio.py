# gene_expression.A549.biorep1.tsv
# chr	start	stop	A549	H1-hESC	HeLa-S3	HepG2	IMR-90	K562	MCF-7
# chr10	120000	120200	B	B	B	A	A	U	B
# chr10	120050	120250	B	B	B	U	U	U	B
# chr10	120100	120300	B	B	B	U	U	U	B

import os


def merge_gene_expression(tsv_path, dir_path, expression_map, gene_list):
    openi = open(dir_path + tsv_path)
    output_name = ""
    cell_line_expression_maps = {}
    empty_gene_list = True
    if len(gene_list) > 0:
        empty_gene_list = False
    for k in openi:
        if k.startswith("gene_id"):
            continue
        splitted_names = tsv_path.split(".")
        cell_line = splitted_names[1]
        biorep = splitted_names[2]
        output_name = cell_line + "." + biorep
        k = k.rstrip()
        splitted = k.split("\t")
        gene_id = splitted[0]
        tpm = splitted[5]
        cell_line_expression_maps[gene_id] = tpm
        if empty_gene_list is True:
            gene_list.append(gene_id)
    expression_map[output_name] = cell_line_expression_maps


def open_expression_file(promoter_file_path, promoter_expression_map):
    tissue_list = []
    openi = open(promoter_file_path)
    for k in openi:
        k = k.rstrip()
        if "chr\tstart\tstop" in k:
            inner_splitted = k.split("\t")
            for i in range(3, 9):
                tissue_list.append(inner_splitted[i])
        splitted = k.split("\t")
        chr = splitted[0]

        gene = splitted[0]
        second_idx = k.index(splitted[1])
        promoter_expression_map[gene] = k[second_idx:]


def main():
    exp_file = "/home/callsobing/dream10/process/gene_expression/merged.PromoterToGene.expression.txt"
    tf_chip_file = "/home/callsobing/dream10/process/ChIP-seq/label/bind.CTCF.train.labels.tsv"
    promoter_expression_map = {}
    open_expression_file(exp_file, promoter_expression_map)
    for i in os.listdir(dir_path):
        if i.endswith(".tsv"):
            merge_gene_expression(i, dir_path, expression_map, gene_list)
    print "gene_name" + "\t" + "chromosome" + "\t" + "promoter_start" + "\t" + "promoter_end" + "\t" + "strand" + "\t",
    for k in expression_map:
        print k + "\t",
        cell_lines.append(k)
    print "\n",
    for gene_id in gene_list:
        print gene_id + "\t" + promoter_list[gene_id] + "\t",
        for cell_line_id in cell_lines:
            print expression_map[cell_line_id][gene_id] + "\t",
        print "\n",


if __name__ == '__main__':
    main()