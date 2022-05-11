#!usr/env/bin python3

"""
Author: Catarina Loureiro

Add rep bin info to bgc table
"""

import argparse

def get_cmds():
        """
        Capture args from the cmdline
        """
        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-b', '-bin', dest='bin', help='bin to rep table',\
                 required=True, metavar='<file>')
        parser.add_argument('-t', '-tsv', dest='tsv', help='bgc info tsv',\
                 required=True, metavar='<file>')
        parser.add_argument('-o', '-out', dest='out', help='updated info tsv',\
                 required=True, metavar='<file>')
        return parser.parse_args()


def build_binrep_dict(binrep_file):
    """
    extract bin:binrep info

    binrep_file: str, filepath
    binrep_dict: dict{str:str}, bin:rep
    """

    fileobj = open(binrep_file, 'r')
    binrep_dict = {}

    for line in fileobj:
        line = line.strip()
        elms = line.split('\t')
        rep = elms[1]
        bins = elms[3]
        if ',' in bins:
            bins = bins.split(',')
            for bin in bins:
                binrep_dict[bin] = rep

        else:
            binrep_dict[bins] = rep

    fileobj.close()
    return binrep_dict

def write_tsv(in_tsv, out_tsv, binrep_dict):
        """
        update tsv with binrep info

        in_tsv, out_tsv: str, filepath
        binrep_dict: dict{str:str}, bin:rep
        """

        inobj = open(in_tsv, 'r')
        outobj = open(out_tsv, 'w')

        for line in inobj:
            line = line.strip()
            if line.startswith('#'):
                outobj.write('{}\tBinRep\n'.format(line))
            else:
                elms = line.split('\t')
                bin = elms[4]
                if bin == '-':
                    outobj.write('{}\t-\n'.format(line))
                else:
                    try:
                        outobj.write('{}\t{}\n'.format(line, binrep_dict[bin]))

                    except:
                        outobj.write('{}\t-\n'.format(line))

        inobj.close()
        outobj.close()
        return None

if __name__ == '__main__':

        cmds = get_cmds()

        dict_binrep = build_binrep_dict(cmds.bin)
        write_tsv(cmds.tsv, cmds.out, dict_binrep)