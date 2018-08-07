# DNASE.A549.conservative.narrowPeak
import os
import numpy as np
from scipy.stats import ttest_ind
from scipy.special import stdtr


def calculate_differential_level(chip_input_label, expression_signal):
    for i in os.listdir(chip_input_label):
        if i.endswith(".tsv"):
            open_chip = open(chip_input_label + i)
            b_count = 0
            b_express_total = 0
            u_count = 0
            u_express_total = 0
            b_expression_array = []
            u_expression_array = []
            cell_lines = []
            for data in open_chip:
                data = data.rstrip()
                splitted = data.split("\t")
                if data.startswith("chr\t"):
                    for count in range(3, len(splitted)):
                        if "IMR-90" in splitted[count]:
                            splitted[count] = "IMR90"
                        cell_lines.append(splitted[count])
                    continue
                chr = splitted[0]
                start = int(splitted[1])
                end = int(splitted[2])
                if ("B" in data) and ("U" in data):
                    for pos in range(start, end):
                        key = chr + "_" + str(pos)
                        if key in expression_signal[cell_lines[0]]:
                            b_inner_express = 0
                            b_inner_count = 0
                            u_inner_express = 0
                            u_inner_count = 0
                            for count in range(0, len(cell_lines)):
                                if "B" in splitted[count + 3]:
                                    b_inner_count += 1
                                    b_inner_express += expression_signal[cell_lines[count]][key]
                                    b_count += 1
                                    b_express_total += expression_signal[cell_lines[count]][key]
                                if "U" in splitted[count + 3]:
                                    u_inner_count += 1
                                    u_inner_express += expression_signal[cell_lines[count]][key]
                                    u_count += 1
                                    u_express_total += expression_signal[cell_lines[count]][key]
                            b_expression_array.append(b_inner_express / float(b_inner_count))
                            u_expression_array.append(u_inner_express / float(u_inner_count))
                            break
            open_chip.close()
            if b_count == 0 or u_count == 0:
                print i + "\t" \
                      + "\tB:" + str(b_count) \
                      + "\tU:" + str(u_count) \
                      + "\tpvalue:" + str("X")
            else:
                tvalue, pvalue = ttest_ind(b_expression_array, u_expression_array, equal_var=False)
                print i + "\t" \
                      + "\tB:" + str(b_express_total / float(b_count)) \
                      + "\tU:" + str(u_express_total / float(u_count)) \
                      + "\tpvalue:" + str(pvalue)


def parse_expresion(expression_file, expression_signal):
    # DNASE.A549.conservative.narrowPeak
    open_expression = open(expression_file)
    first_in = True
    mcf7 = ""
    pan = ""
    gm12878 = ""
    hepg2 = ""
    k562 = ""
    liver = ""
    h1hesc = ""
    imr90 = ""
    ipsc = ""
    hela = ""
    a549 = ""
    pc3 = ""
    hct = ""
    skn = ""
    for k in open_expression:
        k = k.rstrip()
        splitted_content = k.split("\t")
        if first_in:
            mcf7 = splitted_content[5].split(".")[0]
            pan = splitted_content[6].split(".")[0]
            gm12878 = splitted_content[7].split(".")[0]
            hepg2 = splitted_content[8].split(".")[0]
            k562 = splitted_content[9].split(".")[0]
            liver = splitted_content[11].split(".")[0]
            h1hesc = splitted_content[12].split(".")[0]
            imr90 = splitted_content[14].split(".")[0]
            ipsc = splitted_content[20].split(".")[0]
            hela = splitted_content[21].split(".")[0]
            a549 = splitted_content[24].split(".")[0]
            pc3 = splitted_content[26].split(".")[0]
            hct = splitted_content[28].split(".")[0]
            skn = splitted_content[30].split(".")[0]
            expression_signal[mcf7] = {}
            expression_signal[pan] = {}
            expression_signal[gm12878] = {}
            expression_signal[hepg2] = {}
            expression_signal[k562] = {}
            expression_signal[liver] = {}
            expression_signal[h1hesc] = {}
            expression_signal[imr90] = {}
            expression_signal[ipsc] = {}
            expression_signal[hela] = {}
            expression_signal[a549] = {}
            expression_signal[pc3] = {}
            expression_signal[hct] = {}
            expression_signal[skn] = {}
            first_in = False
            continue
        chr = splitted_content[1]
        start = int(splitted_content[2])
        end = int(splitted_content[3])
        for pos in range(start, end):
            key = chr + "_" + str(pos)
            expression_signal[mcf7][key] = (float(splitted_content[5]) + float(splitted_content[19])) / 2
            expression_signal[pan][key] = (float(splitted_content[6]) + float(splitted_content[32])) / 2
            expression_signal[gm12878][key] = (float(splitted_content[7]) + float(splitted_content[18])) / 2
            expression_signal[hepg2][key] = (float(splitted_content[8]) + float(splitted_content[15])) / 2
            expression_signal[k562][key] = (float(splitted_content[9]) + float(splitted_content[10])) / 2
            expression_signal[liver][key] = (float(splitted_content[11]) + float(splitted_content[13])) / 2
            expression_signal[h1hesc][key] = (float(splitted_content[12]) + float(splitted_content[17])) / 2
            expression_signal[imr90][key] = (float(splitted_content[14]) + float(splitted_content[16])) / 2
            expression_signal[ipsc][key] = (float(splitted_content[20]) + float(splitted_content[23])) / 2
            expression_signal[hela][key] = (float(splitted_content[21]) + float(splitted_content[22])) / 2
            expression_signal[a549][key] = (float(splitted_content[24]) + float(splitted_content[25])) / 2
            expression_signal[pc3][key] = (float(splitted_content[26]) + float(splitted_content[27])) / 2
            expression_signal[hct][key] = (float(splitted_content[28]) + float(splitted_content[29])) / 2
            expression_signal[skn][key] = (float(splitted_content[30]) + float(splitted_content[31])) / 2
    open_expression.close()


def main():
    expression_file = "/home/callsobing/dream10/process/gene_expression/merged.PromoterToGene.expression.txt"
    chip_input_label = "/home/callsobing/dream10/process/ChIP-seq/label/"
    expression_signal = {}
    parse_expresion(expression_file, expression_signal)
    calculate_differential_level(chip_input_label, expression_signal)


if __name__ == '__main__':
    main()
