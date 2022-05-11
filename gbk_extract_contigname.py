#!usr/env/bin python3

"""
Author: Catarina Loureiro

Extract contig name from gbks
"""

import os
import argparse
from Bio import SeqIO

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-g', '-gbks', dest='gbks', help='folder with input gbks',\
                required=True, metavar='<path>')
        parser.add_argument('-t', '-tsv', dest='tsv', help='tsv to add column to',\
                 required=True, metavar='<file>')
        return parser.parse_args()

def get_gbk_files(indir):
        """
        Get gbk files from indir incl their paths

        indir: str, folder path
        """

        try:
                for dirpath, dirnames, files in os.walk(indir):
                        for f in files:
                                bgc_name = f.split('/')[-1]
                                yield(os.path.join(dirpath, f))
        except:
                pass

def extract_contig(filepath):
        """
        extract contig name from gbk

        gbk_path: str, path
        """

        gbkcontents = SeqIO.parse(filepath, "genbank")
        for record in gbkcontents:
                return record.description

def write_tsv(tsv_file, contigs_dict, sname_dict):
        """
        take contig name information and update tsv

        tsv_file: str, filepath
        contigs_dict, sname_dict: dict{str:str}
        """

        tsv = tsv_file.split('/')[-1]
        outdir = tsv_file[:-(len(tsv))]
        tsv_out = tsv.split('.ts')[0] + '_contigs_sample.tsv'
        tsv_out_path = os.path.join(outdir, tsv_out)

        infile_obj = open(tsv_file, 'r')
        outfile_objs = open(tsv_out_path, 'w')


        for line in infile_obj:
                line = line.strip()
                if line.startswith('#'):
                        header = line
                        outfile_objs.write('{}\t{}\n'.format(header, 'contig'))
                else:
                        elms = line.split('\t')
                        bgc = elms[0]
                        rep = elms[1]
                        contig = contigs_dict[bgc]
                        sample = sname_dict[bgc]
                        outfile_objs.write('{}\t{}\t{}\t{}\n'.format(bgc, rep, contig, sample))

        infile_obj.close()
        outfile_objs.close()

        return None


if __name__ == '__main__':

        cmds = get_cmds()

        cname_dict = {}
        sample_dict = {}

        for filepath in get_gbk_files(cmds.gbks):
                infile = filepath.split('/')[-1]
                header = infile.split('.gb')[0]
                cname_dict[header] = extract_contig(filepath)
                if 'meta' in header:
                        sample = header.split('_meta')[0]
                        if '14B' in sample:
                                sample = sample.replace("14B", "14")
                        sample_dict[header] = sample

                else:
                        sample = header.split('_c0')[0]
                        sample_dict[header] = sample

        write_tsv(cmds.tsv, cname_dict, sample_dict)
