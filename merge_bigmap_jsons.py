#!usr/env/bin python3

"""
Author: Catarina Loureiro

A script to take the 2-level jsons outputted by bigmap and merge them
"""

import argparse
import json


def get_cmds():
        """
        Capture args from the cmdline
        """

        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-m', '-mash', dest='mash', help='json from mash redundancy',\
                         required=True, metavar='<file>')
        parser.add_argument('-b', '-bigscape', dest='bigscape', help='json from bigscape\
                         redundancy', required=True, metavar='<file>')
        parser.add_argument('-j', '-json_out', dest='json_out', help='merged json',\
                         required=True, metavar='<file>')
        parser.add_argument('-t', '-tsv_out', dest='tsv_out', help='tsv with bgc gcf info',\
                         required=True, metavar='<file>')
        return parser.parse_args()


def edit_dicts(inner, outer):
    """
    edit the  two dicts

    inner: dict
    outer: dict
    merged: dict
    """
    merged = {}
    total_keys = 0
    total_items = 0

    for key, items in outer.items():
        rep_bgc = key.split('|')[1]
        if rep_bgc not in merged.keys():
            merged[rep_bgc] = []
        for item in items:
            in_bgc = item.split('|')[1]
            merged[rep_bgc].append(in_bgc)

    inner_edit = {}

    for key, items in inner.items():
        rep_bgc = key.split('|')[1]
        if rep_bgc not in inner_edit.keys():
            inner_edit[rep_bgc] = []
        for item in items:
            in_bgc = item.split('|')[1]
            inner_edit[rep_bgc].append(in_bgc)

    return (inner_edit, merged)

def merge_dicts(inner, merged):
    """
    merge two dicts

    inner: dict
    merged: dict
    """
    for i_rep, i_items in inner.items():
        if i_rep in merged.keys():
            for i_item in i_items:
                if i_item not in merged[i_rep]:
                    merged[i_rep].append(i_item)

        elif i_rep not in merged.keys():
            for m_key, m_items in merged.items():
                if i_rep in m_items:
                    for i_item in i_items:
                        if i_item not in merged[m_key]:
                            merged[m_key].append(i_item)
        return merged

def make_tsv(merged, tsv_file):
        """
        make bgc, rep tsv

        merged: dict
        tsv_file: file object
        """

        tsv_file.write('#BGC\tGCF_Rep\n'.format())

        for key, items in merged.items():
            for item in items:
                tsv_file.write('{}\t{}\n'.format(item, key))

        return None

if __name__ == '__main__':

        cmds = get_cmds()

        inner_file = open(cmds.mash, 'r')
        inner_data = json.load(inner_file)
        inner_file.close()

        outer_file = open(cmds.bigscape, 'r')
        outer_data = json.load(outer_file)
        outer_file.close()

        edited_inner, edited_merged = edit_dicts(inner_data, outer_data)
        merge_data = merge_dicts(edited_inner, edited_merged)

        json_file = open(cmds.json_out, 'w')
        json.dump(merge_data, json_file, indent=4)
        json_file.close()

        tsv_fileobj = open(cmds.tsv_out, 'w')
        make_tsv(merge_data, tsv_fileobj)
        tsv_fileobj.close()