# *-* coding:utf-8 *-*
# @Time:2022/6/16 16:21
# @Author:Ruijing Cheng
# @File:trim_start_end3.py
# @Software:PyCharm


import os
from Bio import AlignIO
from pathlib import Path


def trim_start_end(in_path, out_path, file):
    """从碱基gap比例小于2.5%的位点开始保留"""
    alignments = AlignIO.read(Path(in_path)/Path(file), "fasta")
    n_row = len(alignments)
    n_col = len(alignments[0])

    i = 0
    while True:
        n_gap = alignments[:, i].count("-")
        if n_gap/n_row < 0.025:
            break
        i += 1  # 18

    j = -1
    while True:
        n_gap = alignments[:, j].count("-")
        if n_gap/n_row < 0.025:
            break
        j -= 1  # -19

    fw = open(Path(out_path) / Path(file), "w")
    for n in range(n_row):
        fw.write(">"+alignments[n].description+"\n")
        fw.write(str(alignments[n].seq[i:j]+alignments[n].seq[j])+"\n")
    fw.close()


if __name__ == "__main__":
    # my_in_path = r"D:\Ruijing\07_data\04_cp_genome\Embryophyta_12248_trimmed"
    # my_out_path = r"D:\Ruijing\07_data\04_cp_genome\Embryophyta_12247_trimmed_start_end3"

    my_in_path = input("input: ")
    my_out_path = input("output: ")
    if not os.path.exists(my_out_path):
        os.makedirs(my_out_path)
    file_list = os.listdir(my_in_path)
    for my_file in file_list:
        if Path(my_file).suffix in [".fasta", ".fas", ".fa"]:
            try:
                trim_start_end(my_in_path, my_out_path, my_file)
            except Exception as result:
                print(result)