# *-* coding:utf-8 *-*
# @Time:2022/1/25 14:30
# @Author:Ruijing Cheng
# @File:get_description2.py
# @Software:PyCharm


import os
import warnings
from progress.bar import Bar
from pathlib import Path
from numpy import *
from Bio import SeqIO
from datetime import datetime
import pandas as pd


def get_description(in_path, out_path):
    """get descriptions of gb files and write a table"""
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # col_name_list = ['accession', 'length', 'ambiguous_bases_percentage', 'date',
    #                  'organelle', 'organism', 'description', 'taxonomy']
    col_name_list = ['accession', 'length', 'gap_percentage', 'date',
                     'organelle', 'organism', 'description', 'taxonomy']
    # todo: add a column of ambiguous bases percentage
    row_name_list = []
    file_list = os.listdir(in_path)
    for my_file in file_list:
        if my_file.split(".")[-1] in ["gb", "gbk", "gbf"]:
            row_name_list.append(my_file)
    descriptions = zeros((len(row_name_list), len(col_name_list)), dtype=int)  # (nrow, ncol)
    df0 = pd.DataFrame(descriptions, index=row_name_list, columns=col_name_list, dtype=str)
    print("Extracting descriptions of %d gb files..." % (len(row_name_list)))
    bar = Bar('Processing', max=len(row_name_list), fill='=', suffix='%(percent)d%%')
    for row_name in row_name_list:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                seq_record = SeqIO.read(Path(in_path) / Path(row_name), "gb")
                # print(Path(row_name).stem)
            # count CDS and rRNA
            df0.loc[row_name]["accession"] = Path(row_name).stem
            df0.loc[row_name]["date"] = seq_record.annotations["date"]
            if "organelle" in seq_record.features[0].qualifiers.keys():
                df0.loc[row_name]["organelle"] = str(seq_record.features[0].qualifiers["organelle"][0])
            else:
                df0.loc[row_name]["organelle"] = "unknown"
            df0.loc[row_name]["organism"] = seq_record.annotations["organism"]
            df0.loc[row_name]["description"] = seq_record.description
            df0.loc[row_name]["taxonomy"] = ""
            for category in seq_record.annotations["taxonomy"]:
                df0.loc[row_name]["taxonomy"] += "%s|" % category
            df0.loc[row_name]["taxonomy"] += seq_record.annotations["organism"]
            seq_str = str(seq_record.seq)
            df0.loc[row_name]["length"] = len(seq_str)
            df0.loc[row_name]["gap_percentage"] = 100*seq_str.count("N")/len(seq_str)
        except Exception as result:
            print("Error: %s || Please check %s" % (result, row_name))
        finally:
            bar.next()
    bar.finish()
    print("Writing results...")
    df0.to_csv(Path(out_path) / Path("description.txt"), index=False, sep="\t")
    print("Summary of description written in description.txt")


if __name__ == "__main__":
    my_in_path = input("Please input the gb files directory path: ")
    my_out_path = input("Please input the output directory path: ")
    time0 = datetime.now()
    get_description(my_in_path, my_out_path)
    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))
