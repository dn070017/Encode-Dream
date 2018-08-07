# ARID3A.train.labels.tsv
import os


def parse_binding_site(tsv_path, dir_path, out_path):
    openi = open(dir_path + tsv_path)
    b_file = open(out_path + "bind." + tsv_path, 'w')
    for k in openi:
        k = k.rstrip()
        if k.startswith("chr\tstart"):
            b_file.write(k + "\n")
        if "B" in k:
            b_file.write(k + "\n")
    openi.close()
    b_file.close()

def main():
    dir_path = "/home/callsobing/dream10/data/ChIPseq/label/"
    out_path = "/home/callsobing/dream10/process/ChIP-seq/label/"
    for i in os.listdir(dir_path):
        if i.endswith(".tsv"):
            parse_binding_site(i, dir_path, out_path)


if __name__ == '__main__':
    main()