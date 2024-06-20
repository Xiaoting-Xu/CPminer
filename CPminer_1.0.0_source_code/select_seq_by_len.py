# -*- codeing = utf-8 -*-
# @Time : 2021/8/30 15:43
# @Author : Ruijing Cheng
# @File : select_seq_by_len.py
# @Software : PyCharm

import os
import datetime
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import numpy as np


def create_seq_dict(fasta_path):
    def get_accession(record):
        parts = record.description
        return parts
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"), key_function=get_accession)
    return seq_dict


def create_length_df(seq_dict):
    data = np.zeros((len(seq_dict), 3), dtype=str)
    df = pd.DataFrame(data, columns=["accession", "description", "length"], dtype=str)
    i = 0
    for key in seq_dict.keys():
        df.loc[i]["accession"] = key.split("|")[0]  # accession
        df.loc[i]["description"] = key  # description
        df.loc[i]["length"] = len(seq_dict[key].seq)
        i += 1
    df["length"] = pd.to_numeric(df["length"])
    return df


def select_seq_by_len(in_path, out_path):
    """select sequence you want"""
    file_list = os.listdir(in_path)
    if "filtered.csv" in file_list:
        df0 = pd.read_table(Path(in_path) / Path("filtered.csv"), index_col=0, sep=',', engine='python')
    else:
        print("Could not find filtered.csv")
        return
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    for file in file_list:
        if Path(file).suffix in [".fa", ".fas", ".fasta"]:
            print("Start processing: %s" % str(Path(file)))
            try:
                time1 = datetime.datetime.now()
                seq_dict = create_seq_dict(Path(in_path)/Path(file))
                seq_df = create_length_df(seq_dict)
                out_file = open(Path(out_path)/Path(file), "w")
                for name, group in seq_df.groupby(seq_df["accession"]):
                    id_max = group["length"].idxmax()
                    out_file.write(">"+group.loc[id_max, "description"]+"\n")
                    out_file.write(str(seq_dict[group.loc[id_max, "description"]].seq)+"\n")
                    df0.loc[name][Path(file).stem] = group.loc[id_max, "length"]
                out_file.close()
                df0.to_csv(Path(out_path) / Path("length.csv"), index=True, sep=",")
                time2 = datetime.datetime.now()
                print("Running time: %s seconds" % (time2-time1))
            except Exception as result:
                print(result)


if __name__ == "__main__":
    time0 = datetime.datetime.now()
    # my_in_path = r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta\filtered"
    # my_out_path = r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta\uni_acc"
    my_in_path = input("The input directory: ")
    my_out_path = input("The output directory: ")
    select_seq_by_len(my_in_path, my_out_path)
    time3 = datetime.datetime.now()
    print("Total running time: %s seconds" % (time3-time0))
