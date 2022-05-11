#!usr/env/bin python3

"""
Author: Catarina Loureiro

Extract aa sequence of A domains in antismash gbks
"""

import os
import argparse
from Bio import SeqIO

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-i', '-index', dest='index', help='text file with list\
                of gbk names', required=True, metavar='<file>')
        parser.add_argument('-g', '-gbks', dest='gbks', help='folder with input gbks',\
                required=True, metavar='<path>')
        parser.add_argument('-o', '-out', dest='out', help='nucleotide fasta',\
                 required=True, metavar='<file>')
        return parser.parse_args()

def parse_index(index_file):
        """
        Parse gbk names

        index_file: str, filepath
        gbk_list: list[str], gbk names
        """

        fileobj = open(index_file, 'r')
        gbk_list = []

        for line in fileobj:
                line = line.strip()
                gbk_filename = line + '.gbk'
                gbk_list.append(gbk_filename)

        fileobj.close()
        return gbk_list



def get_gbk_files(indir, gbk_list):
        """
        Get gbk files from indir incl their paths

        indir: str, folder path
        """

        try:
                for dirpath, dirnames, files in os.walk(indir):
                        for f in files:
                                bgc_name = f.split('/')[-1]
                                if bgc_name in gbk_list:
                                        yield(os.path.join(dirpath, f))
        except:
                pass

def extract_adomain(filepath, fasta_obj):
        """
        trim a gbk file based on orf borders

        gbk_path, outdir: str, path
        borders: list[int]
        """

        infile = filepath.split('/')[-1]
        header = infile.split('.gb')[0]

        gbkcontents = SeqIO.parse(filepath, "genbank")
        for record in gbkcontents:
                for index, feature in enumerate(record.features):
                        if feature.type == 'PFAM_domain' and feature.qualifiers['aSDomain'] == ['AMP-binding']:
                                print(feature)
                                sequence = feature.qualifiers['translation'][0]
                                header = '>{}_{}_AMP-binding_{}_{}'.format(header,feature.qualifiers['locus_tag'][0], feature.qualifiers['protein_start'][0], feature.qualifiers['protein_end'][0])
                                fasta_obj.write('{}\n{}\n'.format(header, sequence))


                        elif feature.type == 'aSDomain' and feature.qualifiers['label'] == ['MXAN_1528_AMP-binding.1']:
                                print(feature)
                                sequence = feature.qualifiers['translation'][0]
                                header = '>{}_{}_AMP-binding_{}_{}'.format(header,feature.qualifiers['locus_tag'][0], feature.qualifiers['protein_start'][0], feature.qualifiers['protein_end'][0])
                                fasta_obj.write('{}\n{}\n'.format(header, sequence))

        return None


if __name__ == '__main__':

        cmds = get_cmds()
        fileobj = open(cmds.out, 'w')

        list_gbks = parse_index(cmds.index)

        for filepath in get_gbk_files(cmds.gbks, list_gbks):
                print(filepath)
                extract_adomain(filepath, fileobj)

        fileobj.close()
