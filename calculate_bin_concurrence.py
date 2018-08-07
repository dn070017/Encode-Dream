import os


def calculate_bin_concurrence(tsv_path, dir_path):
    suffix_index = tsv_path.rindex(".train.labels.tsv")
    tf_name = tsv_path[:suffix_index]
    openi = open(dir_path + tsv_path)
    bound_line_count = 0
    bin_counts = []
    for k in openi:
        k = k.rstrip()
        splitted = k.split("\t")
        if k.startswith("chr\tstart\tstop"):
            bin_counts = [0 for i in range(0, len(splitted) - 3)]
        else:
            bound_flag = False
            bound_count = 0
            for l in range(3, len(splitted)):
                status = splitted[l]
                if status is "B":
                    bound_flag = True
                    bound_count += 1
            if bound_flag:
                bound_line_count += 1
                bin_counts[bound_count - 1] += 1

    print "\n" + tf_name + "(" + str(bound_line_count) + ")" + "\t",
    for bin_count_idx in range(0, len(bin_counts)):
        print str(bin_count_idx + 1) + ":" + str(bin_counts[bin_count_idx]/float(bound_line_count)) + "\t",


def main():
    dir_path = "/home/callsobing/dream10/data/ChIPseq/label/"
    for i in os.listdir(dir_path):
        if i.endswith(".tsv"):
            calculate_bin_concurrence(i, dir_path)


if __name__ == '__main__':
    main()