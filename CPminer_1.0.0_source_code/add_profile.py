# *-* coding:utf-8 *-*
# @Time:2022/4/2 17:19
# @Author:Ruijing Cheng
# @File:add_profile.py
# @Software:PyCharm


import os
from datetime import datetime
from pathlib import Path


def add_profile(profile1, profile2, out_path):

    if not os.path.exists(out_path):
        os.makedirs(out_path)
    #
    # mafft --addprofile msa1 msa2 > profile
    # msa1 must form a monophyletic cluster
    cmd_str = r"D:\Public\Chengruijing\02_software\mafft-7.487-win64-signed\mafft-win\mafft --thread -1 --reorder --add %s %s > %s"
    # cmd_str = "mafft --thread -1 --auto %s/%s > %s/%s"
    msa_list = os.listdir(profile1)
    for msa in msa_list:
        t0 = datetime.now()
        try:
            if msa.split(".")[-1] in ["fasta", "fas", "fa"]:
                print("Processing: "+msa)
                # out_file = "profile_%s" % msa
                tmp_cmd = cmd_str % (Path(profile1)/Path(msa), Path(profile2)/Path("msa_"+msa), Path(out_path)/Path("msa_"+msa))
                print("Command: %s" % tmp_cmd)
                os.system(tmp_cmd)
                t1 = datetime.now()
                print("Running time: %s Seconds" % (t1 - t0))
        except Exception as result:
            print("Error: %s || Please check %s" % (result, msa))


if __name__ == "__main__":
    # my_profile1 = "./Saxifraga_YLX"
    # my_profile2 = "./angiosperm_all_with_reference0.75"
    # my_out_path = "./angiosperm_all_with_reference0.75_add"
    # my_in_path = input("input directory: ")
    # my_out_path = input("output directory: ")

    my_profile1 = input("New sequences to be added: ")
    my_profile2 = input("Aligned seuquences: ")
    my_out_path = input("Output directory: ")
    time0 = datetime.now()
    print("Start processing...")
    add_profile(my_profile1, my_profile2, my_out_path)
    time1 = datetime.now()
    print("Total running time: %s Seconds" % (time1 - time0))