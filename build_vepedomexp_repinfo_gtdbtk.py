#!usr/env/bin python3

"""
Author: Catarina Loureiro

Gather bin classificaton and add to BGC attribute table
"""

import argparse
import json

def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-i', '-info', dest='info', help='general info tsv',\
                 required=True, metavar='<file>')
        parser.add_argument('-j', '-json', dest='json', help='vepe reps json',\
                 required=True, metavar='<file>')
        parser.add_argument('-o', '-out', dest='out', help='updated rep tsv',\
                 required=True, metavar='<file>')
        return parser.parse_args()


def parse_allinfo(bgcinfo_file):
    """
    extract bgc, binrep

    bgcinfo_file: str, path
    bcgbin_dict: dict{str:str}
    bingtdbtk_dict: dict{str:str}
    """

    fileobj = open(bgcinfo_file, 'r')
    bgcbin_dict = {}
    bingtdbtk_dict = {}

    for line in fileobj:
        line = line.strip()
        if line.startswith('#'):
            header = line
        else:
            elms = line.split('\t')
            bgc = elms[0]
            binrep = elms[5]
            binclass = elms[6]
            bgcbin_dict[bgc] = binrep
            bingtdbtk_dict[binrep] = binclass

    return (bgcbin_dict, bingtdbtk_dict)

def write_reinfo(reps_dict, bgcbin_dict, bingtdbtk_dict, outfile):
        """
        compile info on reps

        reps_dict: dict{str:[str]}
        bcgbin_dict: dict{str:str}
        bingtdbtk_dict: dict{str:str}
        outfile: str, path
        """

        fileobj = open(outfile, 'w')
        fileobj.write('#repBGC\tSamples\tRepBins\n'.format())

        for rep, bgcs in reps_dict.items():
            samples = []
            binreps = []
            binclasses = []
            if rep.startswith('BGC'):
                fileobj.write('{}\t-\t-\n'.format(rep))
                continue
            for bgc in bgcs:
                # bgc = bgc.replace('14B', '14')
                binrep = bgcbin_dict[bgc]
                binclass = bingtdbtk_dict[binrep]
                if binrep not in binreps:
                    binreps.append(binrep)
                if binclass not in binclasses:
                    binclasses.append(binclass)

                if 'meta' in bgc:
                    sample = bgc.split('_meta')[0]
                    sample = sample.replace('14B', '14')
                    if sample not in samples:
                        samples.append(sample)
                else:
                    sample = bgc.split('_c0')[0]
                    if sample not in samples:
                        samples.append(sample)
            try:
                samples = ','.join(samples)
            except:
                pass

            try:
                binreps = ','.join(binreps)
            except:
                pass

            try:
                binclasses = ','.join(binclasses)
            except:
                pass

            fileobj.write('{}\t{}\t{}\t{}\n'.format(rep, samples, binreps, binclasses))

            return None

if __name__ == '__main__':

        cmds = get_cmds()

        file_json = open(cmds.json, 'r')
        dict_bgcrep = json.load(file_json)
        dict_bgcbin, dict_bingtdbtk = parse_allinfo(cmds.info)

        write_reinfo(dict_bgcrep, dict_bgcbin, dict_bingtdbtk, cmds.out)