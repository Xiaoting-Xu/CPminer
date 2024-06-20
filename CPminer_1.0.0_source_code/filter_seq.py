# *-* coding:utf-8 *-*
# @Time:2022/4/2 15:55
# @Author:Ruijing Cheng
# @File:filter_seq.py.py
# @Software:PyCharm

from Bio import SeqIO
from pathlib import Path
import os
import pandas as pd


def filter_seq(in_path, out_path, ref_len_dict):
    """remove duplicate keys and sequences too long or too short"""
    file_list = os.listdir(in_path)
    if "cds_num.csv" in file_list:
        df0 = pd.read_table(Path(in_path) / Path("cds_num.csv"), index_col=0, sep=',', engine='python')
    else:
        print("Could not find cds_num.csv")
        return
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    for file in file_list:
        if Path(file).suffix in [".fa", ".fas", ".fasta"]:
            try:
                print("Start processing: %s" % str(Path(in_path)/Path(file)))
                key_list = []
                for record in SeqIO.parse(Path(in_path)/Path(file), "fasta"):
                    if record.description in key_list:
                        df0.loc[record.description.split("|")[0]][Path(file).stem] -= 1
                        if df0.loc[record.description.split("|")[0]][Path(file).stem] == 0:
                            df0.loc[record.description.split("|")[0]]["cds_num"] -= 1
                        continue
                    else:
                        key_list.append(record.description)
                    if 0.5 * ref_len_dict[Path(file).stem] <= len(record.seq) <= 2 * ref_len_dict[Path(file).stem]:
                        fw0 = open(Path(out_path) / Path(file), "a")
                        fw0.write(">%s" % record.description)
                        fw0.write("\n")
                        fw0.write(str(record.seq))
                        fw0.write("\n")
                        fw0.close()
                    else:
                        df0.loc[record.description.split("|")[0]][Path(file).stem] -= 1
                        if df0.loc[record.description.split("|")[0]][Path(file).stem] == 0:
                            df0.loc[record.description.split("|")[0]]["cds_num"] -= 1
            except Exception as result:
                print(result)
    df0.to_csv(Path(out_path) / Path("filtered.csv"), index=True, sep=",")


if __name__ == "__main__":
    angiosperm_ref_len_dict = {'accD': 1599, 'atpA': 1524, 'atpB': 1503, 'atpE': 405, 'atpF': 555, 'atpH': 246,
                               'atpI': 747, 'ccsA': 942, 'cemA': 690, 'clpP': 609, 'infA': 234, 'matK': 1506,
                               'ndhA': 1092, 'ndhB': 1479, 'ndhC': 363, 'ndhD': 1503, 'ndhE': 303, 'ndhF': 2241,
                               'ndhG': 534, 'ndhH': 1182, 'ndhI': 543, 'ndhJ': 477, 'ndhK': 762, 'petA': 963,
                               'petB': 648, 'petD': 483, 'petG': 114, 'petL': 96, 'petN': 90, 'psaA': 2253,
                               'psaB': 2205, 'psaC': 246, 'psaI': 111, 'psaJ': 129, 'psbA': 1053, 'psbB': 1527,
                               'psbC': 1422, 'psbD': 1062, 'psbE': 252, 'psbF': 120, 'psbH': 222, 'psbI': 111,
                               'psbJ': 123, 'psbK': 186, 'psbL': 117, 'psbM': 108, 'psbN': 132, 'psbT': 108,
                               'psbZ': 189, 'rbcL': 1428, 'rpl14': 369, 'rpl16': 408, 'rpl2': 822, 'rpl20': 360,
                               'rpl22': 375, 'rpl23': 288, 'rpl32': 174, 'rpl33': 207, 'rpl36': 114, 'rpoA': 1005,
                               'rpoB': 3219, 'rpoC1': 2043, 'rpoC2': 4110, 'rps11': 417, 'rps12': 258, 'rps14': 303,
                               'rps15': 264, 'rps16': 237, 'rps18': 306, 'rps19': 279, 'rps2': 711, 'rps3': 657,
                               'rps4': 606, 'rps7': 468, 'rps8': 399, 'rrn16': 1490, 'rrn23': 2814, 'rrn4.5': 103,
                               'rrn5': 121, 'ycf1': 5385, 'ycf15': 147, 'ycf2': 6915, 'ycf3': 507, 'ycf4': 708,
                               'ycf68': 234}
    gymnosperm_ref_len_dict = {'accD': 1041, 'atpA': 1524, 'atpB': 1479, 'atpE': 417, 'atpF': 555, 'atpH': 246,
                               'atpI': 747, 'ccsA': 966, 'cemA': 786, 'clpP': 609, 'infA': 249, 'matK': 1500,
                               'ndhA': 1107, 'ndhB': 1485, 'ndhC': 363, 'ndhD': 1503, 'ndhE': 303, 'ndhF': 2208,
                               'ndhG': 543, 'ndhH': 1182, 'ndhI': 558, 'ndhJ': 522, 'ndhK': 813, 'petA': 963,
                               'petB': 648, 'petD': 507, 'petG': 114, 'petL': 129, 'petN': 90, 'psaA': 2253,
                               'psaB': 2205, 'psaC': 246, 'psaI': 111, 'psaJ': 135, 'psbA': 1062, 'psbB': 1527,
                               'psbC': 1422, 'psbD': 1062, 'psbE': 252, 'psbF': 120, 'psbH': 228, 'psbI': 111,
                               'psbJ': 123, 'psbK': 177, 'psbL': 117, 'psbM': 105, 'psbN': 132, 'psbT': 108,
                               'psbZ': 264, 'rbcL': 1428, 'rpl14': 369, 'rpl16': 417, 'rpl2': 831, 'rpl20': 348,
                               'rpl22': 420, 'rpl23': 276, 'rpl32': 210, 'rpl33': 201, 'rpl36': 114, 'rpoA': 1026,
                               'rpoB': 3219, 'rpoC1': 2046, 'rpoC2': 4101, 'rps11': 393, 'rps12': 258, 'rps14': 303,
                               'rps15': 270, 'rps16': 255, 'rps18': 228, 'rps19': 279, 'rps2': 708, 'rps3': 657,
                               'rps4': 606, 'rps7': 471, 'rps8': 399, 'rrn16': 1475, 'rrn23': 2813, 'rrn4.5': 103,
                               'rrn5': 121, 'ycf1': 5082, 'ycf15': 363, 'ycf2': 7308, 'ycf3': 510, 'ycf4': 555,
                               'ycf68': 240}
    # my_in_path = Path(r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta")
    # my_out_path = Path(r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta\filtered")
    my_in_path = input("The input directory: ")
    my_out_path = input("The output direcotory: ")
    filter_seq(my_in_path, my_out_path, angiosperm_ref_len_dict)
