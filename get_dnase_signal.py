# DNASE.A549.conservative.narrowPeak
import os


def read_bin_sites(train_bin_path, dnase_bin_signal):
    open_bin = open(train_bin_path)
    for k in open_bin:
        k = k.rstrip()
        splitted_line = k.split("\t")
        chr = splitted_line[0]
        start = splitted_line[1]
        end = splitted_line[2]
        key = chr + "_" + start + "_" + end
        dnase_bin_signal[key] = {}
    open_bin.close()


def parse_dnase(dir_path, dnase_bin_signal):
    # DNASE.A549.conservative.narrowPeak
    for i in os.listdir(dir_path):
        if i.endswith(".narrowPeak"):
            splitted_name = i.split(".")
            cell_line_id = splitted_name[1]
            open_peak = open(dir_path + i)
            dnase_base_signal = {}
            for k in open_peak:
                k = k.rstrip()
                splitted_content = k.split("\t")
                chr = splitted_content[0]
                start = splitted_content[1]
                end = splitted_content[2]
                signal = float(splitted_content[6])
                for pos in range(int(start), int(end)):
                    key = chr + "_" + str(pos)
                    dnase_base_signal[key] = signal
            open_peak.close()
            for keys in dnase_bin_signal:
                splitted_key = keys.split("_")
                chr = splitted_key[0]
                start = int(splitted_key[1])
                end = int(splitted_key[2])
                total_signal = 0
                for pos in range(start, end):
                    key = chr + "_" + str(pos)
                    if key in dnase_base_signal:
                        total_signal += dnase_base_signal[key]
                if total_signal > 0:
                    dnase_bin_signal[keys][cell_line_id] = total_signal/float(200)
                else:
                    dnase_bin_signal[keys][cell_line_id] = 0


def main():
    dir_path = "/home/callsobing/dream10/data/DNASE/peak/conservative/"
    train_bin_path = "/home/callsobing/dream10/data/annotations/train_regions.blacklistfiltered.bed"
    dnase_bin_signal = {}
    read_bin_sites(train_bin_path, dnase_bin_signal)
    parse_dnase(dir_path, dnase_bin_signal)
    ## print header
    print "bin_id\t",
    for cell_line in dnase_bin_signal["chr10_600_800"]:
        print cell_line + "\t",
    print "\n",
    ## print data
    for key in dnase_bin_signal:
        print key + "\t",
        for cell_line_signal in dnase_bin_signal[key]:
            print str(dnase_bin_signal[key][cell_line_signal]) + "\t",
        print "\n",


if __name__ == '__main__':
    main()