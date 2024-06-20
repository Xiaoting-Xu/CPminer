# -*- coding = utf-8 -*-
# @Time : 2021/9/5 22:35
# @Author : Ruijing Cheng
# @File : select_seq_by_acc.py
# @Software: PyCharm

from Bio import SeqIO
import os
import datetime
import pandas as pd
from pathlib import Path


def get_accession(record):
    parts = record.description.split("|")
    # assert len(parts) == 3
    return parts[0]


def select_seq_by_acc(in_path, out_path, table_path):
    selected_table = pd.read_table(Path(table_path), index_col=0, sep='\t', engine='python')
    file_list = os.listdir(in_path)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    for file in file_list:
        if file.split(".")[-1] in ["fa", "fas", "fasta"]:
            print("Start processing: %s" % file)
            time1 = datetime.datetime.now()
            seq_dict = SeqIO.to_dict(SeqIO.parse(Path(in_path)/Path(file), "fasta"), key_function=get_accession)
            for index in selected_table.index:
                if index in seq_dict.keys():
                    # if not os.path.exists(Path(out_path)/Path(selected_table.loc[index]["Profile"])):
                    #     os.makedirs(Path(out_path)/Path(selected_table.loc[index]["Profile"]))
                    fas_description = ">%s|%s" % (selected_table.loc[index]["species"].replace(" ", "_"), index)
                    # fw = open(Path(out_path)/Path(selected_table.loc[index]["Profile"])/Path(file), "a")
                    fw = open(Path(out_path) / Path(file), "a")
                    fw.write(fas_description)
                    fw.write("\n")
                    fw.write(str(seq_dict[index].seq))
                    fw.write("\n")
                    fw.close()
            time2 = datetime.datetime.now()
            print("Running time: %s seconds" % (time2 - time1))


if __name__ == "__main__":
    # my_in_path = r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta\uni_acc"
    # my_out_path = r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta\uni_sp"
    # my_table_path = r"D:\Ruijing\07_data\04_cp_genome\03_info_table\Viridiplantae_30648\Bryophyta_uni_sp.txt"
    # my_in_path = r"./Asterales_2682_80_CDS_20230531/uni_acc"
    # my_out_path = r"./Asterales_2682_80_CDS_20230531/uni_sp"
    # my_table_path = r"./selected.txt"
    my_in_path = input("input:")
    my_out_path = input("output: ")
    my_table_path = input("table of selected sequences: ")
    select_seq_by_acc(my_in_path, my_out_path, my_table_path)
