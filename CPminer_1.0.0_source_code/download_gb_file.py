# *-* coding:utf-8 *-*
# @Time:2022/1/5 10:02
# @Author:Ruijing Cheng
# @File:download_gb_file.py
# @Software:PyCharm


from Bio import Entrez
from Bio import SeqIO
import os
from pathlib import Path
from datetime import datetime
import func_timeout.exceptions
from func_timeout import func_set_timeout


@func_set_timeout(600)
def my_efetch(accession):
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=accession)
    return handle


def download_gb_file(in_path, out_path):
    """download chloroplast genome from Genbank according to given accessions"""
    Entrez.email = "ruijing.cheng@qq.com"
    fr = open(in_path, "r")
    accession_list = fr.read().splitlines()
    fr.close()
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    else: 
        file_list = os.listdir(out_path)
        file_list = [Path(file).stem for file in file_list if Path(file).suffix == ".gb"]
        accession_list = list(set(accession_list) - set(file_list))
    accession_list.sort()
    while len(accession_list) != 0:
        accession = accession_list[0]
        try:
            t0 = datetime.now()
            # hd1 = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")
            hd1 = my_efetch(accession=accession)
            seq = SeqIO.read(hd1, "gb")
            out_file_path = Path(out_path)/Path(accession+".gb")
            out_file = open(out_file_path, "w")
            SeqIO.write(seq, out_file, "gb")
            out_file.close()
            t1 = datetime.now()
            print("%s downloaded in %s seconds" % (accession, t1 - t0))
            accession_list.pop(0)
        except func_timeout.exceptions.FunctionTimedOut:
            print("Time out, try downloading %s again..." % accession)
        except Exception as result:
            print("Error: %s, try downloading %s again..." % (result, accession))


if __name__ == "__main__":
    in_str = input("Please input the accession list file path: ")
    out_str = input("Please input the output directory path: ")
    print("Start downloading...")
    time0 = datetime.now()
    download_gb_file(in_str, out_str)
    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))
