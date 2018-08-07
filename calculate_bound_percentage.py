import os


def basic_statistics_on_tsv(tsv_path):
    suffix_index = tsv_path.rindex(".train.labels.tsv")
    tf_name = tsv_path[:suffix_index]
    openi = open(tsv_path)
    cell_lines = []
    line_count = 0
    a_count = []
    b_count = []
    for k in openi:
        k = k.rstrip()
        splitted = k.split("\t")
        if k.startswith("chr\tstart\tstop"):
            cell_lines = splitted[3:]
            a_count = [0 for i in range(0, len(splitted) - 3)]
            b_count = [0 for i in range(0, len(splitted) - 3)]
        else:
            line_count += 1
            for l in range(3, len(splitted)):
                status = splitted[l]
                if status is "A":
                    a_count[l - 3] += 1
                if status is "B":
                    b_count[l - 3] += 1
    print "\n" + tf_name + "(" + str(line_count) + ")" + "\t",
    for cell_line_idx in range(0, len(cell_lines) ):
        print cell_lines[cell_line_idx] + ":" + \
              str(round(b_count[cell_line_idx]/float(line_count), 4)) \
              + "/" + str(round(a_count[cell_line_idx]/float(line_count), 4)) + "\t",


for i in os.listdir("/home/callsobing/dream10/data/ChIPseq/label"):
    if i.endswith(".tsv"):
        basic_statistics_on_tsv(i)