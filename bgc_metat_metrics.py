#!usr/env/bin python3

"""
Author: Catarina Loureiro

Calclute a number of RNA expression values for BGCs
"""

import os
import argparse

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-b', '-bg', dest='bg', help='bedgraph location',\
                 required=True, metavar='<dir>')
        parser.add_argument('-l', '-lst', dest='lst', help='file with bgc list',\
                 required=True, metavar='<file>')
        parser.add_argument('-c', '-cnt', dest='cnt', help='countsfile location',\
                 required=True, metavar='<dir>')
        parser.add_argument('-o', '-out', dest='out', help='output table',\
                required=True, metavar='<file>')
        return parser.parse_args()


def get_bgc_names(bgc_file):
    """
    extract bin names

    bgc_file: str, path
    bgc_list: [str]
    """

    bgc_list = []
    fileobj = open(bgc_file, 'r')
    for line in fileobj:
        line = line.strip()
        bgc = line
        bgc_list.append(bgc)

    return bgc_list

def parse_bedgraph(bgc_name, bedgraph_file, bedgraph_dict):
    """
    parse bedgraph file, process seqs of interest

    bedgraph_file: str, path
    bgc_name: str
    bedgraph_dict: {bgc:[max,min,wavg,percov]}
    """
    fileobj = open(bedgraph_file, 'r')
    maxcov = ''

    for line in fileobj:
        line = line.strip()
        elms = line.split()
        curr_node = elms[0]

        if curr_node != bgc_name:
            continue

        curr_start = int(elms[1])
        curr_end = int(elms[2])
        curr_reads = int(elms[3])

        if not maxcov:
            mincov = curr_reads
            maxcov = curr_reads
            base_count = curr_end - curr_start
            weight_count = base_count * curr_reads
            if curr_reads > 0:
                cov_bases = base_count
                uncov_bases = 0
            else:
                cov_bases = 0
                uncov_bases = base_count
        else:
            if curr_reads < mincov:
                mincov = curr_reads
            elif curr_reads > maxcov:
                maxcov = curr_reads

            curr_base_count = (curr_end - curr_start)
            base_count += curr_base_count

            curr_weight_count = curr_base_count * curr_reads
            weight_count += curr_weight_count

            if curr_reads > 0:
                cov_bases += curr_base_count
            else:
                uncov_bases += curr_base_count

    perc_cov = cov_bases / base_count
    wavg = weight_count / base_count

    bedgraph_dict[bgc_name].extend([maxcov, mincov, wavg, perc_cov])
    return bedgraph_dict


def calc_TPM_RPKM(counts_file, TPM_dict, RPKM_dict, BGC):
    """
    TPM = rate/sum(rate)
    rate = nreads/cluster_length (kb)
    RPKM = read_counts/(cluster_length * sum(read_counts)) * 10^9

    counts_file: str, path
    TPM_dict:{bgc:tpm}
    RPKM_dict:{bgc:rpkm}
    BGC:str
    """

    # TPM: calculate tpm for bgc from native mapping, incl ratesums of other bgcs in that mapping

    rates = {}
    ratesum = 0

    sum_reads = 0
    read_counts = {}
    cluster_lengths = {}

    with open(counts_file, "r") as f:
        for line in f:
            if "*" not in line:
                line = line.strip()
                cluster, length, nreads, nnoreads = line.split("\t")

                read_counts[cluster] = float(nreads)
                cluster_lengths[cluster] = float(length)
                sum_reads += float(nreads)

                try:
                    rate = float(nreads) / float(length)
                    rates[cluster] = rate
                    ratesum += rate
                except(ZeroDivisionError):
                    print('Unable to calculate TPM counts')
                    pass

    TPM_dict[BGC] = rates[BGC] / ratesum
    RPKM_dict[BGC] = read_counts[BGC] / (sum_reads * cluster_lengths[BGC]) * 1000000000
    return None

def write_output(bedgraph_dict, out_file, TPM_dict, RPKM_dict):
    """
    write table with per bin values, each column is a pfam

    bedgraph_dict: {pfam:[max,min,wavg,percov]}
    binname: str
    out_file: str
    all_pfams: list
    TPM_dict:{bgc:tpm}
    RPKM_dict:{bgc:rpkm}
    """
    out_file_headers = (f'#BGC\tmax\tmin\tpercov\twavg_cov\tTPM\tRPKM')
    outobj = open(out_file, 'w')
    outobj.write(f'{out_file_headers}\n')

    for bgc, vals in bedgraph_dict.items():
        vals = [str(val) for val in vals]
        vals.append(str(TPM_dict[bgc]))
        vals.append(str(RPKM_dict[bgc]))
        vals = '\t'.join(vals)
        outobj.write(f'{bgc}\t{vals}\n')

    outobj.close()
    return None


if __name__ == '__main__':

    cmds = get_cmds()

    list_bgcs = get_bgc_names(cmds.lst)

    dict_bedgraph = {bgc: [] for bgc in list_bgcs}
    dict_TPM = {bgc: '' for bgc in list_bgcs}
    dict_RPKM = {bgc: '' for bgc in list_bgcs}

    for bgc in list_bgcs:
        sample = bgc.split('|')[1]
        sample = sample.split('_')[0]

        bg_name = sample + '.bg'
        file_bedgraph = os.path.join(cmds.bg, bg_name)

        cnt_name = sample + '.sorted.count'
        file_counts = os.path.join(cmds.cnt, cnt_name)

        parse_bedgraph(bgc, file_bedgraph, dict_bedgraph)
        calc_TPM_RPKM(file_counts, dict_TPM, dict_RPKM, bgc)

        write_output(dict_bedgraph, cmds.out, dict_TPM, dict_RPKM)