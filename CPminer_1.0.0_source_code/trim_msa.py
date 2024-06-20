# *-* coding:utf-8 *-*
# @Time:2022/4/26 21:05
# @Author:Ruijing Cheng
# @File:trim_msa.py
# @Software:PyCharm


import os
from pathlib import Path
from datetime import datetime


def trim_msa(in_path, out_path, cmd):
    cmd_str = "%s -in %s -out %s -fasta -automated1"
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    file_list = os.listdir(in_path)
    for file in file_list:
        t0 = datetime.now()
        try:
            if Path(file).suffix in [".fa", ".fas", ".fasta"]:
                in_file = in_path/Path(file)
                out_file = out_path/Path(file)
                tmp_cmd = cmd_str % (cmd, str(in_file), str(out_file))
                print("Command: %s" % tmp_cmd)
                os.system(tmp_cmd)
                t1 = datetime.now()
                print("Running time: %s Seconds" % (t1 - t0))
        except Exception as result:
            print("Error: %s || Please check %s" % (result, file))


if __name__ == "__main__":
    my_in_path = Path(r"gymnosperm_uni_sp_with_reference")
    my_out_path = Path(r"gymnosperm_uni_sp_with_reference_trimmed")
    my_cmd = r"D:\Public\Chengruijing\02_software\trimal.v1.2rev59\trimAl\bin\trimal"
    trim_msa(my_in_path, my_out_path, my_cmd)
