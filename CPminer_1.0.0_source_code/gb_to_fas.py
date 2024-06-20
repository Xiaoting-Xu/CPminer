# -*- coding = utf-8 -*-
# @Time : 2022/1/22 19:57
# @Author : Ruijing Cheng
# @File : gb_to_fas.py
# @Software: PyCharm

import os
import warnings
from pathlib import Path
from Bio import SeqIO


def gb_to_fas(in_path, in_file):
    """convert gb files to fasta files"""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            record = SeqIO.read(Path(in_path) / Path(in_file), "gb")
        sub_folder = record.annotations["taxonomy"][6]
        out_path = Path(in_path) / Path(sub_folder)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_file_path = Path(out_path)/Path(record.id+".fasta")
        fw = open(out_file_path, "w")
        fw.write(">")
        fw.write(record.annotations["taxonomy"][6])
        fw.write("|")
        fw.write(record.annotations["organism"])
        fw.write("|")
        fw.write(record.id)
        fw.write("\n")
        fw.write(str(record.seq))
        fw.write("\n")
        fw.close()
    except Exception as result:
        print("Error: %s, please check %s" % (result, in_file))


if __name__ == "__main__":
    my_in_path = input("input path: ")
    file_list = os.listdir(my_in_path)
    for my_in_file in file_list:
        if my_in_file.split(".")[-1] == "gb":
            gb_to_fas(my_in_path, my_in_file)
