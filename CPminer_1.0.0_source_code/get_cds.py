# -*- coding = utf-8 -*-
# @Time : 2022/1/23 15:02
# @Author : Ruijing Cheng
# @File : get_cds.py
# @Software: PyCharm

import os
import warnings
# from progress.bar import Bar
from pathlib import Path
import numpy as np
from Bio import SeqIO
from datetime import datetime
import pandas as pd


def get_cds(in_path, out_path):
    """Extract CDS and rRNA sequences and save as fasta files, the number of the same gene found are written"""
    if not os.path.exists(Path(out_path)):
        os.makedirs(Path(out_path))
    cds_list = ['accD', 'atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI', 'ccsA', 'cemA',
                'clpP', 'infA', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI',
                'ndhJ', 'ndhK', 'petA', 'petB', 'petD', 'petG', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI',
                'psaJ', 'psbA', 'psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbI', 'psbJ', 'psbK', 'psbL',
                'psbM', 'psbN', 'psbT', 'psbZ', 'rbcL', 'rpl14', 'rpl16', 'rpl2', 'rpl20', 'rpl22', 'rpl23',
                'rpl32', 'rpl33', 'rpl36', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps14', 'rps15',
                'rps16', 'rps18', 'rps19', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4.5',
                'rrn5', 'ycf1', 'ycf15', 'ycf2', 'ycf3', 'ycf4', 'ycf68']
    alternative_name = {'16S ribosomal RNA': 'rrn16', '23S ribosomal RNA': 'rrn23',
                        '4.5S ribosomal RNA': 'rrn4.5', '5S ribosomal RNA': 'rrn5',
                        'rrn16S': 'rrn16', 'rrn23S': 'rrn23',
                        'rrn4.5S': 'rrn4.5', 'rrn5S': 'rrn5'
                        }
    gb_file_list = []
    file_list = os.listdir(Path(in_path))
    for my_file in file_list:
        if Path(my_file).suffix in [".gb", ".gbf"]:
            gb_file_list.append(Path(my_file).stem)
    locations = np.zeros((len(gb_file_list), len(cds_list)+1), dtype=int)  # (nrow, ncol)
    df0 = pd.DataFrame(locations, index=gb_file_list, columns=cds_list+['cds_num'], dtype=int)
    df1 = pd.DataFrame(locations, index=gb_file_list, columns=cds_list+['cds_num'], dtype=str)
    for cds in cds_list:
        fw = open((Path(out_path) / Path("%s.fasta" % cds)), "w")
        fw.close()
    print("Start extracting %d CDS and rRNA in %d gb files..." % (len(cds_list), (len(gb_file_list))))
    bar = Bar('Processing', max=len(gb_file_list), fill='=', suffix='%(percent)d%%')
    for gb_file_stem in gb_file_list:
        gene_name = "unknown"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                seq_record = SeqIO.read(Path(in_path)/Path(gb_file_stem + ".gb"), "gb")
            for feature in seq_record.features:
                try:
                    if feature.type in ["CDS", "rRNA"]:
                        gene_name = "unknown"
                        if "gene" in feature.qualifiers.keys():
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name in alternative_name.keys():
                                gene_name = alternative_name[gene_name]
                        elif "product" in feature.qualifiers.keys():
                            product_name = feature.qualifiers["product"][0]
                            if product_name in alternative_name.keys():
                                gene_name = alternative_name[product_name]
                        if gene_name in cds_list:
                            fas_seq = ""
                            for part in feature.location.parts:  # 'NoneType' object has no attribute 'parts'
                                start = part.nofuzzy_start
                                end = part.nofuzzy_end
                                tmp_seq = seq_record.seq[start:end]
                                if part.strand == -1:
                                    tmp_seq = tmp_seq.reverse_complement()
                                fas_seq += tmp_seq
                            if df0.loc[gb_file_stem][gene_name] == 0:
                                df0.loc[gb_file_stem]["cds_num"] += 1
                                df1.loc[gb_file_stem][gene_name] = str(feature.location)
                            else:
                                df1.loc[gb_file_stem][gene_name] += "|%s" % str(feature.location)
                            df0.loc[gb_file_stem][gene_name] += 1
                            fw = open((Path(out_path)/Path("%s.fasta" % gene_name)), "a")
                            # fas_description = ">%s" % seq_record.id
                            fas_description = ">%s|%s" % (gb_file_stem, str(feature.location))
                            fw.write(fas_description)
                            fw.write("\n")
                            fw.write(str(fas_seq))
                            fw.write("\n")
                            fw.close()
                except Exception as result:
                    print("Error: %s, please check %s %s" % (result, gb_file_stem, gene_name))
            df1.loc[gb_file_stem]["cds_num"] = df0.loc[gb_file_stem]["cds_num"]
        except Exception as result:
            print("Error: %s, please check %s %s" % (result, gb_file_stem, gene_name))
        finally:
            bar.next()
    bar.finish()
    print("Writing summary...")
    df0.to_csv(Path(out_path) / Path("cds_num.csv"), index=True, sep=",")
    df1.to_csv(Path(out_path) / Path("cds_loc.txt"), index=True, sep="\t")
    print("Summary of cds number written in cds_num.csv")
    print("Summary of cds location written in cds_loc.txt")


if __name__ == "__main__":
    # my_in_path = r"D:\Ruijing\07_data\04_cp_genome\02_gb_file\Bryophyta"
    # my_out_path = r"D:\Ruijing\07_data\04_cp_genome\04_coding_sequence\Bryophyta"
    my_in_path = input("Input path: ")
    my_out_path = input("Output path: ")
    time0 = datetime.now()
    get_cds(my_in_path, my_out_path)
    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))
