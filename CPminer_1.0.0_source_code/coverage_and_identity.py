# *-* coding:utf-8 *-*
# @Time:2022/4/13 10:28
# @Author:Ruijing Cheng
# @File:coverage_and_identity.py
# @Software:PyCharm


from Bio import SeqIO
from Bio import AlignIO
import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd


def coverage_and_identity(in_path, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    file_list = os.listdir(in_path)
    for file in file_list:
        if Path(file).suffix in [".fa", ".fas", ".fasta"]:
            try:
                print("Start processing: "+file)
                alignments = AlignIO.read(in_path/Path(file), "fasta")
                n_row = len(alignments)
                n_col = len(alignments[0])
                ref_len = n_col - alignments[0].seq.count("-")
                sum_mat = np.zeros((n_row, 3), dtype=str)
                sum_table = pd.DataFrame(sum_mat, columns=["description", "coverage", "identity"], dtype=str)
                col_list = []
                for i in range(n_col):
                    if alignments[0][i] != "-":
                        col_list.append(i)
                for i in range(n_row):
                    match = 0
                    mismatch = 0
                    for j in col_list:
                        if alignments[i][j] != "-":
                            if alignments[0][j] == alignments[i][j]:
                                match += 1
                            else:
                                mismatch += 1
                    coverage = (match + mismatch) / ref_len
                    identity = match / (match + mismatch)
                    sum_table.loc[i][0] = alignments[i].description
                    sum_table.loc[i][1] = coverage
                    sum_table.loc[i][2] = identity
                sum_table.to_csv(out_path/Path(Path(file).stem+"_coverage_and_identity.txt"), index=False, sep="\t")
            except Exception as result:
                print(result)


if __name__ == "__main__":
    my_in_path = Path(input("Please type the input directory: "))
    my_out_path = Path(input("Please type the output directory: "))
    coverage_and_identity(my_in_path, my_out_path)


