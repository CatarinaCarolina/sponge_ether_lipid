#!usr/env/bin python3

"""
Author: Catarina Loureiro

Add rep bin classification from gtdbtk to BGC table
"""

import argparse

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-c', '-cla', dest='cla', help='gtdbtk classification',\
                 required=True, metavar='<file>')
        parser.add_argument('-t', '-tsv', dest='tsv', help='bgc info tsv',\
                 required=True, metavar='<file>')
        parser.add_argument('-o', '-out', dest='out', help='updated info tsv',\
                 required=True, metavar='<file>')
        return parser.parse_args()


def build_binrepclass_dict(gtdbtk_file):
    """
    extract binrep:class info

    gtdbtk_file: str, filepath
    gtdbtk_dict: dict{str:str}, binrep:class
    """

    fileobj = open(gtdbtk_file, 'r')
    gtdbtk_dict = {}

    for line in fileobj:
        line = line.strip()
        if line.startswith('user'):
            header = line
            continue
        elms = line.split('\t')
        binrep = elms[0] + '.fa'  ## check
        class_elms = elms[1].split(';')
        class_cut = ';'.join([class_elms[1], class_elms[5], class_elms[6]])

        gtdbtk_dict[binrep] = class_cut

    fileobj.close()
    return gtdbtk_dict

def write_tsv(in_tsv, out_tsv, gtdbtk_dict):
        """
        update tsv with binrep info

        in_tsv, out_tsv: str, filepath
        gtdbtk_dict: dict{str:str}, binrep:class
        """

        inobj = open(in_tsv, 'r')
        outobj = open(out_tsv, 'w')

        for line in inobj:
            line = line.strip()
            if line.startswith('#'):
                outobj.write('{}\tBinRepClass\n'.format(line))
            else:
                elms = line.split('\t')
                binrep = elms[5]
                if binrep == '-':
                    outobj.write('{}\t-\n'.format(line))
                else:
                    try:
                        outobj.write('{}\t{}\n'.format(line, gtdbtk_dict[binrep]))

                    except:
                        outobj.write('{}\t-\n'.format(line))

        inobj.close()
        outobj.close()
        return None

if __name__ == '__main__':

        cmds = get_cmds()
        dict_gtdbtk = build_binrepclass_dict(cmds.cla)
        write_tsv(cmds.tsv, cmds.out, dict_gtdbtk)