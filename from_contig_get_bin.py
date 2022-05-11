#!usr/env/bin python3

"""
Author: Catarina Loureiro

Find contigs in respective bins
"""

import argparse
import subprocess

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-i', '-in_f', dest='in_f', help='tsv with contigs',\
                 required=True, metavar='<file>')
        parser.add_argument('-o', '-out_f', dest='out_f', help='tsv with contigs and bins',\
                required=True, metavar='<file>')
        parser.add_argument('-b', '-bin_loc', dest='bin_loc', help='path to bins folder',\
                required=True, metavar='<str>')

        return parser.parse_args()


def check_output(command, shell=False):
        """
        A function to use subprocess to run cmdline actions

        command: str
        """

        try:
                output = subprocess.check_output(command, shell=shell)

        except subprocess.CalledProcessError as e:
                output = e.output
                return output

        return output.decode('utf-8')


def fetch_contigs(in_file, out_file, bins_folder):
        """
        A function to fetch bin ids from contigs

        in_file: str, filepath
        out_file: str, filepath
        bins_folder: str, folder path
        """

        in_fobj = open(in_file, 'r')
        out_fobj = open(out_file, 'w')

        for line in in_fobj:
                line = line.strip()
                if line.startswith('#'):
                        out_fobj.write('{}\tBin\n'.format(line))
                else:
                        # BGC,contig,gcf_class,GCF,sample = line.split('\t')
                        BCG, rep, contig, sample = line.split('\t')
                        command = 'grep "{}" {}{}*'.format(contig, bins_folder, sample)
                        # print(command)
                        bin_path = check_output(command, shell=True)

                        if str(bin_path).startswith("b''"):
                                bin_id = '-'

                        else:
                                bin_loc = bin_path.split(':>')[0]
                                bin_id = bin_loc.split('/')[-1]

                        # print(BGC, ' ', bin_id)
                        out_fobj.write('{}\t{}\n'.format(line, bin_id))

        return None

if __name__ == '__main__':

        cmds = get_cmds()

        data_tsv = cmds.in_f
        data_bins_tsv = cmds.out_f
        bins_path = cmds.bin_loc

        if not bins_path.endswith('/'):
                bins_path = bins_path+'/'

        fetch_contigs(data_tsv, data_bins_tsv, bins_path)