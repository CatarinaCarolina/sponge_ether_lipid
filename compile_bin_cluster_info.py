#!usr/env/bin python3

"""
Author: Catarina Loureiro

Compile dREP bin cluster elements and rep bin info
"""

import argparse

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-r', '-rep', dest='rep', help='represntative bins',\
                 required=True, metavar='<file>')
        parser.add_argument('-c', '-cluster', dest='cluster', help='bin cluster info',\
                 required=True, metavar='<file>')
        parser.add_argument('-o', '-out', dest='out', help='output file',\
                 required=True, metavar='<file>')
        return parser.parse_args()


def parse_reps(rep_bins):
    """
    A function to parse a text file into list

    rep_bins: str, file path
    rep_list: lst[str]
    """

    file_obj = open(rep_bins, 'r')
    rep_list = []

    for line in file_obj:
        line = line.strip()
        rep_list.append(line)

    return rep_list

def parse_cluster(cluster_file):
    """
    A function to extract bin cluster info and rep

    cluster_file: str, filepath
    dict_cluster_bins: dict{cluster:[bins]}
    """
    dict_cluster_bins = {}
    c_fobj = open(cluster_file, 'r')

    for line in c_fobj:
        line = line.strip()
        if line.startswith('genome'):
            header = line
        else:
            bin, sec_c, thres, meth, algo, pri_c = line.split(',')
            if sec_c not in dict_cluster_bins.keys():
                dict_cluster_bins[sec_c] = [bin]
            else:
                dict_cluster_bins[sec_c].append(bin)

    return dict_cluster_bins


def get_cluster_rep(dict_cluster_bins, rep_list):
    """
    A function to extract the representative of each cluster

    dict_cluster_bins: dict{cluster:[bins]}
    rep_list: list[cluster]
    dict_cluster_rep: dict{cluster:rep}
    """

    dict_cluster_rep = {key: '' for key in dict_cluster_bins.keys()}

    for clu, bins in dict_cluster_bins.items():
        for bin in bins:
            if bin in rep_list:
                rep = bin
        dict_cluster_rep[clu] = rep

    return dict_cluster_rep

def get_cluster_samples(dict_cluster_bins):
    """
    A function to extract the samples the cluster is present in

    dict_cluster_bins: dict{cluster:[bins]}
    dict_cluster_samp: dict{cluster:[samples]}
    """

    dict_cluster_samp = {key: [] for key in dict_cluster_bins.keys()}

    for clu, bins in dict_cluster_bins.items():
        for bin in bins:
            sample = bin.split('_bin')[0]
            if sample not in dict_cluster_samp[clu]:
                dict_cluster_samp[clu].append(sample)

    return dict_cluster_samp


def write_bin_output(file_out, dict_cluster_bins, dict_cluster_rep, dict_cluster_samp, dict_eco_groups):
    """
    A sample to write all elements to an output file

    file_out: str, filepath
    dict_cluster_bins: dict{cluster:[bins]}
    dict_cluster_rep: dict{cluster:rep}
    dict_cluster_samp: dict{cluster:[samples]}
    dict_eco_groups: dict{group:[samples]}
    """
    # bin_cluster bin_rep bins samples

    file_obj = open(file_out, 'w')

    for cluster in dict_cluster_bins.keys():
        rep = dict_cluster_rep[cluster]
        bins = ','.join(dict_cluster_bins[cluster])
        samples = ','.join(dict_cluster_samp[cluster])
        eco_groups = []
        for sample in dict_cluster_samp[cluster]:
            for eco in dict_eco_groups.keys():
                if sample + '_' in dict_eco_groups[eco] and eco not in eco_groups:
                    eco_groups.append(eco)
        eco_groups = ','.join(eco_groups)

        file_obj.write('{}\t{}\t{}\t{}\t{}\n'.format(cluster, rep, samples, bins, eco_groups))

    return None

if __name__ == '__main__':

    cmds = get_cmds()

    bin_reps = cmds.rep
    bin_clu = cmds.cluster
    out_file = cmds.out

    ecogroup2_samples = {'Aplysina':['Aply16_', 'Aply21_', 'Aply22_', 'Aply23_'],\
                        'Geodia_shallow':['gb2_2_', 'gb4_2_', 'gb5_2_', 'gb6_'],\
                        'Geodia_mid':['gb1_','gb7_', 'gb8_2_','gb9_','gb10_'],\
                        'Geodia_deep':['gb126_', 'gb278_', 'gb305_'],\
                        'Petrosia':['Pf4_', 'Pf5_', 'Pf6_','Pf7_', 'Pf8_', 'Pf9_','Pf10_', 'Pf11_', 'Pf12_'],\
                        'Med_SW':['sw_7_', 'sw_8_', 'sw_9_'],\
                        'Atl_SW': ['gb_1_f_', 'gb_2_f_', 'gb_f_3_','gb5_6_f_', 'gb_f_9_', 'gb10_f_'],\
                        'DOM': ['DOM10_', 'DOM14_','DOM33_', 'DOM40_', 'DOM43_']}


    list_reps = parse_reps(bin_reps)
    clu_bin_dict = parse_cluster(bin_clu)
    clu_rep_dict = get_cluster_rep(clu_bin_dict, list_reps)
    clu_samp_dict = get_cluster_samples(clu_bin_dict)

    write_bin_output(out_file,clu_bin_dict, clu_rep_dict,clu_samp_dict, ecogroup2_samples)
