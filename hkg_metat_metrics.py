#!usr/env/bin python3

"""
Author: Catarina Loureiro

Calculate a number of RNA expression metrics for specific genes
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
        parser.add_argument('-l', '-lst', dest='lst', help='file with bin list',\
                 required=True, metavar='<file>')
        parser.add_argument('-p', '-pfam', dest='pfam', help='pfam pairings location',\
                 required=True, metavar='<dir>')
        parser.add_argument('-o', '-out', dest='out', help='output table',\
                required=True, metavar='<file>')

        return parser.parse_args()


def get_bin_names(bin_file):
    """
    extract bin names

    bin_file: str, path
    bin_list: [str]
    """

    bin_list = []
    fileobj = open(bin_file, 'r')
    for line in fileobj:
        line = line.strip()
        bin = line[:-3]
        bin_list.append(bin)

    return bin_list


def get_pfam_node(pairings_file):
    """
    get node:pfam pairs

    pairings_file: str, path
    pairings_dict: {node:pfam}
    """

    pairings_dict = {}
    fileobj = open(pairings_file, 'r')
    for line in fileobj:
        line = line.strip()
        elms = line.split()
        node = elms[0]
        pfam = elms[1]
        pairings_dict[node] = pfam

    return (pairings_dict)

def parse_bedgraph(bedgraph_file, pairings_dict):
    """
    parse bedgraph file, process seqs of interest

    bedgraph_file: str, path
    pairings_dict: {node:pfam}
    bedgraph_dict: {pfam:[max,min,wavg,percov]}
    """
    nodes = pairings_dict.keys()
    bedgraph_dict = {pairings_dict[node]: [] for node in nodes}
    fileobj = open(bedgraph_file, 'r')
    prev_node = ''
    maxcov = ''
    for line in fileobj:
        line = line.strip()
        elms = line.split()
        curr_node = elms[0]

        if curr_node not in nodes:
            continue

        curr_start = int(elms[1])
        curr_end = int(elms[2])
        curr_reads = int(elms[3])
        if not prev_node == curr_node:
            if maxcov:  # if values exist write values to dict & reset
                # print('------> new node in already growing dict')
                perc_cov = cov_bases / base_count
                wavg = weight_count / base_count
                pfam = pairings_dict[prev_node]
                bedgraph_dict[pfam].extend([maxcov, mincov, wavg, perc_cov])

                mincov = curr_reads
                maxcov = curr_reads

                base_count = curr_end - curr_start
                # print(base_count)
                weight_count = base_count * curr_reads

                if curr_reads > 0:
                    cov_bases = base_count
                    uncov_bases = 0
                else:
                    cov_bases = 0
                    uncov_bases = base_count
                prev_node = curr_node
                continue

            else:  # else create growing vals
                # print('-----> dict start')
                mincov = curr_reads
                maxcov = curr_reads
                base_count = curr_end - curr_start
                # print(base_count)
                weight_count = base_count * curr_reads
                if curr_reads > 0:
                    cov_bases = base_count
                    uncov_bases = 0
                else:
                    cov_bases = 0
                    uncov_bases = base_count
                prev_node = curr_node
                continue
                # add to growing vals
        if curr_reads < mincov:
            mincov = curr_reads
        elif curr_reads > maxcov:
            maxcov = curr_reads

        curr_base_count = (curr_end - curr_start)
        # print(curr_base_count)
        base_count += curr_base_count
        # print(base_count)

        curr_weight_count = curr_base_count * curr_reads
        weight_count += curr_weight_count

        if curr_reads > 0:
            cov_bases += curr_base_count
        else:
            uncov_bases += curr_base_count

    perc_cov = cov_bases / base_count
    # print(cov_bases, base_count)
    wavg = weight_count / base_count
    # print(weight_count, base_count)
    pfam = pairings_dict[prev_node]
    bedgraph_dict[pfam].extend([maxcov, mincov, wavg, perc_cov])

    return bedgraph_dict


def calc_TPM(counts_file, pairings_dict):
    """
    TPM = rate/sum(rate)
    rate = nreads/cluster_length (kb)

    counts_file:str, path
    pairings_dict: {node:pfam}
    TPM_dict: {pfam:tpm}
    """
    rates = {}
    ratesum = 0

    with open(counts_file, "r") as f:
        for line in f:
            if "*" not in line:
                line = line.strip()
                node, length, nreads, nnoreads = line.split("\t")
                pfam = pairings_dict[node]
                try:
                    rate = float(nreads) / float(length)
                    rates[pfam] = rate
                    ratesum += rate
                except(ZeroDivisionError):
                    print('Unable to calculate TPM counts')
                    pass

    TPM_dict = {}
    for key in rates:
        try:
            TPM_dict[key] = rates[key] / ratesum
        except(ZeroDivisionError):
            TPM_dict[key] = 0

    return TPM_dict

def calc_RPKM(counts_file, pairings_dict):
        """
        RPKM = read_counts/(cluster_length * sum(read_counts)) * 10^9

        counts_file:str, path
        pairings_dict: {node:pfam}
        RPKM_dict: {pfam:rpkm}
        """
        sum_reads = 0
        read_counts = {}
        node_lengths = {}

        with open(counts_file, "r") as f:
            for line in f:
                if "*" not in line:
                    line = line.strip()
                    node, length, nreads, nnoreads = line.split("\t")
                    pfam = pairings_dict[node]
                    nreads = float(nreads)
                    read_counts[pfam] = nreads
                    node_lengths[pfam] = float(length)
                    sum_reads += nreads
        RPKM_dict = {}
        for key in read_counts:
            try:
                RPKM_dict[key] = read_counts[key] / (sum_reads * node_lengths[key]) * 1000000000
            except(ZeroDivisionError):
                RPKM_dict[key] = 0

        return RPKM_dict


def write_output(bedgraph_dict, binname, out_file, all_pfams, TPM_dict, RPKM_dict):
    """
    write table with per bin values, each column is a pfam

    bedgraph_dict: {pfam:[max,min,wavg,percov]}
    TPM_dict: {pfam:tpm}
    RPKM_dict:{pfam:rpkm}
    binname: str
    out_file: str
    all_pfams: list
    """

    fileobj = open(out_file, 'a')
    new_line = [binname]

    for pfam in all_pfams:
        if pfam in bedgraph_dict.keys():
            new_line.extend([str(val) for val in bedgraph_dict[pfam]])
            new_line.append(str(TPM_dict[pfam]))
            new_line.append(str(RPKM_dict[pfam]))
        else:
            new_line.extend(['-', '-', '-', '-', '-', '-'])
    new_line = '\t'.join(new_line)
    fileobj.write(f'{new_line}\n')
    fileobj.close()

    return None

if __name__ == '__main__':

        cmds = get_cmds()

        list_bins = get_bin_names(cmds.lst)

        pfam_list = ['PF00521', 'PF00204', 'PF00154', 'PF01000', 'PF00562']
        out_file_headers = ['bin']
        for pfam in pfam_list:
                out_file_headers.extend([pfam+'_max',pfam+'_min',pfam+'_percov', pfam+'_wavg_cov',pfam+'_tmp',pfam+'_rpkm'])
        out_file_headers = '\t'.join(out_file_headers)

        outobj = open(cmds.out, 'w')
        outobj.write(f'{out_file_headers}\n')
        outobj.close()

        for bin in list_bins:
                sample = bin.split('_')[0]
                pairing_name = bin+'_prod_trans_hkgcds.pairing.txt'
                file_pairing = os.path.join(cmds.pfam, pairing_name)
                bg_name = bin+'_'+sample+'.sorted.bg'
                file_bedgraph = os.path.join(cmds.bg, bg_name)
                counts_name = bin+'_'+sample+'.hkgcds.sorted.aln.count'
                file_counts = os.path.join(cmds.bg, counts_name)

                print(bin, sample)
                print(pairing_name, bg_name)
                print(file_pairing)
                print(file_bedgraph)
                print(file_counts)

                dict_pairings = get_pfam_node(file_pairing)
                print(dict_pairings)

                dict_bedgraph = parse_bedgraph(file_bedgraph, dict_pairings)
                print(dict_bedgraph)

                dict_TPM = calc_TPM(file_counts,dict_pairings)
                print(dict_TPM)

                dict_RPKM = calc_RPKM(file_counts,dict_pairings)
                print(dict_RPKM)

                write_output(dict_bedgraph, bin, cmds.out, pfam_list,dict_TPM,dict_RPKM)