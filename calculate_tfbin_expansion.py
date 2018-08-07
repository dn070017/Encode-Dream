import os


def calculate_tfbin_expansion(tsv_path, dir_path):
    suffix_index = tsv_path.rindex(".train.labels.tsv")
    tf_name = tsv_path[:suffix_index]
    openi = open(dir_path + tsv_path)
    bin_counts = 0
    bin_status = []
    total_bin_length = 0
    bin_length = []
    tmp_bin_length = []
    for k in openi:
        k = k.rstrip()
        splitted = k.split("\t")
        if k.startswith("chr\tstart\tstop"):
            bin_status = [False for i in range(0, len(splitted) - 3)]
            tmp_bin_length = [0 for i in range(0, len(splitted) - 3)]
        else:
            for l in range(3, len(splitted)):
                status = splitted[l]
                if status is not "B" and bin_status[l - 3] is True:
                    bin_status[l - 3] = False
                    bin_length.append(tmp_bin_length[l - 3])
                    tmp_bin_length[l - 3] = 0
                    bin_counts += 1
                if status is not "B":
                    continue
                if bin_status[l - 3] is False:
                    bin_status[l - 3] = True
                total_bin_length += 1
                tmp_bin_length[l - 3] += 1

    print tf_name + "\t" + str(total_bin_length/float(bin_counts)) + "\t" + str(pstdev(bin_length))


def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) # in Python 2 use sum(data)/float(n)


def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss


def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5


def main():
    dir_path = "/home/callsobing/dream10/data/ChIPseq/label/"
    print "TF_name" + "\t" + "Avg_tfbs_bin_len" + "\t" + "Std_tfbs_bin_len"
    for i in os.listdir(dir_path):
        if i.endswith(".tsv"):
            calculate_tfbin_expansion(i, dir_path)


if __name__ == '__main__':
    main()