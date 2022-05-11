#!usr/env/bin python3

"""
Author: Catarina Loureiro

Extract coordinates from genes that have hits for specific pfams and make a bedfile
"""

import os
import argparse


def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-d', '-indir', dest='indir', help='dir with hmm outputs',\
                 required=True, metavar='<path>')
        parser.add_argument('-p', '-pfam', dest='pfam', help='txt file with pfams',\
                 required=True, metavar='<file>')
        parser.add_argument('-n', '-nuc', dest='nuc', help='nucleotide fastas dir',\
                required=True, metavar='<path>')
        parser.add_argument('-o', '-out', dest='out', help='out dir for bedfiles',\
                 required=True, metavar='<path>')
        return parser.parse_args()


def get_files(indir):
    """
    Get fa files from indir incl their paths

    indir: str, folder path
    """
    try:
        for dirpath, dirnames, files in os.walk(indir):
            for f in files:
                file_name = f.split('/')[-1]
                if file_name.endswith('Pfam-A.domtable'):
                    yield (os.path.join(dirpath, f))
    except:
        pass


def parse_fasta_header(file_prefix, fasta_dir):
    """
    get query and length

    fasta_dir: str, filepath
    file_prefix: str
    fasta_dict: dict {query:lenght}
    """

    fasta_dict = {}

    nucpath = os.path.join(fasta_dir, file_prefix + '_refined_prodigal_nuc.fa')
    nucobj = open(nucpath, 'r')

    for line in nucobj:
        if not line.startswith('>'):
            continue
        line = line.strip()
        elms = line.split()
        query = elms[0][1:]
        length = int(elms[4]) - int(elms[2]) + 1
        fasta_dict[query] = length
    nucobj.close()
    return fasta_dict

def get_pfams(pfam_file):
    """
    get pfam names

    pfam_file: str, path
    pfam_list:[str]
    """
    pfam_list = []

    fileobj = open(pfam_file, 'r')
    for line in fileobj:
        line = line.strip()
        pfam_list.append(line)

    fileobj.close()
    return pfam_list


def get_pfam_hits(domtable_file, pfam_list):
    """
    get locs and bitscore for pfams in list

    domtable_file: str, filepath
    pfam_list:[str]
    pfam_dict: {pfam:[query,bitscore,start,end]}
    """

    pfam_dict = {}

    fileobj = open(domtable_file, 'r')

    for line in fileobj:
        line = line.strip()
        if line.startswith('#'):
            continue
        elms = line.split()
        pfam = elms[1]
        pfam = pfam.split('.')[0]
        query = elms[3]
        bitscore = elms[13]
        start = 0
        end = elms[5]
        if pfam in pfam_list:
            try:
                old_bitscore = pfam_dict[pfam][1]
            except KeyError:
                pfam_dict[pfam] = [query, bitscore, start, end]
            else:
                if bitscore > old_bitscore:
                    pfam_dict[pfam] = [query, bitscore, start, end]
    fileobj.close()
    return pfam_dict

def write_output(pfam_dict, out_dir, domtable_file, fasta_dict):
        """
        write output bedfile

        pfam_dict: {pfam:[query,bitscore,start,end]}
        fasta_dict: dict {query:length}
        domtable_file: str, filepath
        out_dir: str, dirpath
        """

        domtable_name = domtable_file.split('/')[-1]
        bed_filename = domtable_name.split('Pfam-A')[0] + 'hkgcds.bed'
        bed_filepath = os.path.join(out_dir, bed_filename)

        pair_filename = domtable_name.split('Pfam-A')[0] + 'hkgcds.pairing.txt'
        pair_filepath = os.path.join(out_dir, pair_filename)

        pair_fileobj = open(pair_filepath, 'w')
        fileobj = open(bed_filepath, 'w')
        for pfam, elms in pfam_dict.items():
            query = elms[0]
            end = fasta_dict[query]
            fileobj.write(f'{query}\t0\t{end}\n')
            pair_fileobj.write(f'{query}\t{pfam}\n')
            # fileobj.write(f'{elms[0]}\t{elms[2]}\t{elms[3]}\n')

        pair_fileobj.close()
        fileobj.close()
        return None

if __name__ == '__main__':

        cmds = get_cmds()

        list_pfam = get_pfams(cmds.pfam)
        dir_nuc = cmds.nuc
        print(dir_nuc)

        for file_path in get_files(cmds.indir):
                print(file_path)
                binname = file_path.split('_prod')[0]
                binname = binname.split('/')[-1]
                dict_fasta = parse_fasta_header(binname, dir_nuc)
                dict_pfam = get_pfam_hits(file_path, list_pfam)
                write_output(dict_pfam, cmds.out, file_path, dict_fasta)