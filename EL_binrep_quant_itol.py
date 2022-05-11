#!usr/env/bin python3

"""
Author: Catarina Loureiro

A script to build itol databases with dREP representative bin quantification
"""

import argparse
import pandas as pd
import numpy as np
import collections

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-q', '-quant', dest='quant', help='all bin quantification',\
                required=True, metavar='<file>')
        parser.add_argument('-d', '-drep', dest='drep', help='drep clusters',\
                 required=True, metavar='<file>')
        parser.add_argument('-i', '-ind', dest='ind', help='output table individual',\
                 required=True, metavar='<file>')
        parser.add_argument('-s', '-scaled', dest='scaled', help='output scaled table individual',\
                required=True, metavar='<file>')

        return parser.parse_args()

def parse_ELbinquant(binquant_file):
        """
        parse binquant file to extract links and individuals

        binquant_file: str, path
        inds_list: list[str], individuals
        binquant_dict: dict{bin:quant}
        """
        inds_list = []
        binquant_dict = {}

        fileobj = open(binquant_file, 'r')
        for line in fileobj:
                line = line.strip()
                elms = line.split('\t')
                bin = elms[0]+'.fa'
                quant = float(elms[1])
                ind = bin.split('bin.')[0][:-1]
                binquant_dict[bin] = quant
                if ind not in inds_list:
                        inds_list.append(ind)


        fileobj.close()
        return (inds_list, binquant_dict)


def scale_ELbinquant(binquant_dict):
        """
        a function to scale the quant values

        binquant_dict: dict{bin:quant}
        scaled_binquant_dict: dict{bin:quant}
        """

        bins = []
        quants = []
        scaled_binquant_dict = {}
        # np.interp(your_arr, (1, max(quants)), (1, 5))

        binquant_dict = collections.OrderedDict(sorted(binquant_dict.items(), key=lambda x: x[1]))

        for bin, quant in binquant_dict.items():
                bins.append(bin)
                quants.append(quant)

        scaled_quants = np.interp(quants, (1, max(quants)), (1, 5))

        for i in range(len(scaled_quants)):
                scaled_binquant_dict[bins[i]] = scaled_quants[i]

        return scaled_binquant_dict


def parse_drep(drep_file):
        """
        parse drep_spsavg file to extract drep cluster structure

        drep_file: str, path
        drep_dict: dict{rep:[bins]}
        """

        drep_dict = {}
        fileobj = open(drep_file, 'r')

        for line in fileobj:
                line = line.strip()
                elms = line.split('\t')
                drep_dict[elms[1]] = elms[2].split(',')

        fileobj.close()
        return drep_dict

def write_output(out_df_file, binquant_dict, inds_list, drep_dict):
        """
        for each rep add abundances to ind-matrix for ELbins

        inds_list: list[str], individuals
        binquant_dict: dict{bin:quant}
        drep_dict: dict{rep:[bins]}
        out_df_file: str, path
        """

        matrix = [[0]*len(inds_list) for i in range(len(drep_dict.keys()))]
        reps = []
        rep_index = 0
        for rep, bins in drep_dict.items():
                reps.append(rep)
                for bin in bins:
                        ind = bin.split('bin.')[0][:-1]
                        try:
                                ind_index = inds_list.index(ind)
                                matrix[rep_index][ind_index] = binquant_dict[bin]
                        except:
                                continue
                rep_index += 1

        df = pd.DataFrame(matrix, columns = inds_list, index = reps)
        df.to_csv(out_df_file, sep = '\t')

        return None

if __name__ == '__main__':

        cmds = get_cmds()
        file_binquant = cmds.quant
        file_drep = cmds.drep
        file_out_ind = cmds.ind
        file_out_scaled = cmds.scaled


        list_inds, dict_quant = parse_ELbinquant(file_binquant)
        dict_scaled_quant = scale_ELbinquant(dict_quant)


        dict_drep = parse_drep(file_drep)
        write_output(file_out_ind,dict_quant,list_inds,dict_drep)

        write_output(file_out_scaled,dict_scaled_quant,list_inds,dict_drep)