openi = open("/home/callsobing/dream10/data/annotations/gencode/genes.list.gtf")

for k in openi:
    ## chr1    HAVANA  gene    11869   14412   .       +       .       gene_id "ENSG00000223972.4";
    if k.startswith("#"):
        continue
    k = k.rstrip()
    splitted = k.split("\t")
    chr = splitted[0]
    pos1 = splitted[3]
    pos2 = splitted[4]
    strand = splitted[6]
    meta = splitted[8]
    # meta => gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene";....
    splitted_meta = meta.split(";")
    for meta_data in splitted_meta:
        if meta_data.startswith("gene_id"):
            gene_id_pair = meta_data.split("\"")
            gene_id = gene_id_pair[1].rstrip("\"")
            if strand is "+":
                print gene_id + "\t" + chr + "\t" + str(int(pos1) - 1000) + "\t" + pos1 + "\t" + "+"
            if strand is "-":
                print gene_id + "\t" + chr + "\t" + pos2 + "\t" + str(int(pos2) + 1000) + "\t" + "-"