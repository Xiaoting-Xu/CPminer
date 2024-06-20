# *-* coding:utf-8 *-*
# @Time:2021/12/2 22:22
# @Author:Ruijing Cheng
# @File:seq_align.py
# @Software:PyCharm

import os
from datetime import datetime
from pathlib import Path


def seq_align(in_path, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    cmd_str = "mafft --thread -1 --auto --reorder  %s/%s > %s/%s"
    # cmd_str = "mafft --thread -1 --auto %s/%s > %s/%s"
    file_name_list = os.listdir(in_path)
    for file_name in file_name_list:
        t0 = datetime.now()
        try:
            if file_name.split(".")[-1] in ["fasta", "fas", "fa"]:
                print("Processing: "+file_name)
                out_file_name = "msa_%s.fasta" % Path(file_name).stem
                tmp_cmd = cmd_str % (in_path, file_name, out_path, out_file_name)
                print("Command: %s" % tmp_cmd)
                os.system(tmp_cmd)
                t1 = datetime.now()
                print("Running time: %s Seconds" % (t1 - t0))
        except Exception as result:
            print("Error: %s || Please check %s" % (result, file_name))


if __name__ == "__main__":

    # my_in_path = "D:\\Ruijing\\07_data\\04_cp_genome\\04_coding_sequence\\CDS_Angiosperm\\renamed"
    # my_out_path = "D:\\Ruijing\\07_data\\04_cp_genome\\05_multi_seq_align\\CDS_Angiosperm"
    my_in_path = input("input directory: ")
    my_out_path = input("output directory: ")
    # my_in_path = r"./Asterales_2682_80_CDS_20230531/uni_sp"
    # my_out_path = r"./Asterales_2682_80_CDS_20230531/uni_sp_msa"
    time0 = datetime.now()
    print("Start processing...")
    seq_align(my_in_path, my_out_path)
    time1 = datetime.now()
    print("Total running time: %s Seconds" % (time1 - time0))
