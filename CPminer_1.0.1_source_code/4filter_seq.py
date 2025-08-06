
# -*- codeing = utf-8 -*-
# @Time : 2021/8/30 15:43
# @Author : Ruijing Cheng
# @File : select_seq_by_len.py
# @Software : PyCharm


from Bio import SeqIO
from pathlib import Path
import os, sys
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import numpy as np
from datetime import datetime


def create_seq_dict(fasta_path):
    def get_accession(record):
        parts = record.description
        return parts
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"), key_function=get_accession)
    return seq_dict


def create_length_df(seq_dict):
    data = np.zeros((len(seq_dict), 3), dtype=str)
    df = pd.DataFrame(data, columns=["accession", "description","length"], dtype=str)
    i = 0
    for key in seq_dict.keys():
        df.loc[i, "accession"] = key.split("|")[0]  # accession
        df.loc[i, "description"] = key  # accession
        df.loc[i, "length"] = len(seq_dict[key].seq)
        i += 1
    df["length"] = pd.to_numeric(df["length"])
    return df


def select_seq_by_len(in_path, out_path, ref_len_dict, lower_bound, upper_bound, threads=3):
    file_list = [f for f in os.listdir(in_path) 
                if Path(f).suffix.lower() in [".fa", ".fas", ".fasta"]]
    if not file_list:
        print("No fasta files found")
        return

    df0 = pd.read_csv(Path(in_path) / "cds_num.csv", index_col='filename', sep=',')
    os.makedirs(out_path, exist_ok=True)
    
    # 存储各线程返回的结果: {文件名: (更新数据, 异常信息)}
    results = {}
    
    def process_file(file):
        print(f"Filter record of : {file} by length")
        gene = Path(file).stem
        try:
            # 检查基因是否在参考字典中
            if gene not in ref_len_dict:
                raise KeyError(f"基因名 {gene} 不在参考字典中")
            
            seq_dict = create_seq_dict(Path(in_path)/file)
            seq_df = create_length_df(seq_dict)
            ref_len = ref_len_dict[gene]
            updates = {}  # 存储更新数据: {accession: length}
            
            with open(Path(out_path)/file, "w") as out_file:
                for name, group in seq_df.groupby("accession"):
                    idx_max = group["length"].idxmax()
                    seq_rec = seq_dict[group.loc[idx_max, "description"]]
                    seq_len = len(seq_rec.seq)
                    
                    if lower_bound * ref_len <= seq_len <= upper_bound * ref_len:
                        out_file.write(f">{seq_rec.description}\n{seq_rec.seq}\n")
                        updates[name] = seq_len
            
            return updates, None
        except Exception as e:
            return None, str(e)

    # 使用线程池并行处理
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_file = {
            executor.submit(process_file, file): file 
            for file in file_list
        }
        
        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                updates, error = future.result()
                results[file] = (updates, error)
            except Exception as e:
                results[file] = (None, str(e))

    # 主线程更新DataFrame
    for file, (updates, error) in results.items():
        gene = Path(file).stem
        if error:
            print(f"处理文件 {file} 出错: {error}")
            continue
            
        if gene not in df0.columns:
            df0[gene] = np.nan  # 添加新列
            
        for acc, length in updates.items():
            if acc in df0.index:
                df0.loc[acc, gene] = length

    df0.to_csv(Path(out_path) / "length.csv", index=True, sep=",")

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Filter sequences with CDS lengths within a defined percentile range(default: 50%-200\%) of the orthologous CDS length distribution from reference genomes.
        Author: Ruijing Cheng
        Organization: 
        GitHub Site:
        Usage Example:
          添加示例
        """
        return f"{description}\n\n{help_text}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    # 输入输出参数
    parser.add_argument("-i", "--input_dir", help="Input directory containing CDS files (.fasta)", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for filtered CDS sequences (.fasta)", required=True)
    
    # 参考基因组选择
    parser.add_argument("-r", "--ref_orthologous", choices=['Ang', 'Gym'],
                        help="Reference orthologous CDS lengths:\n"
                             "Ang: Angiosperm (flowering plants)\n"
                             "Gym: Gymnosperm (conifers, cycads, etc.)",required=True)
    
    # 分位数范围参数
    parser.add_argument("-lb", "--lower_bound", type=float, default=0.5,
                        help="Lower percentile bound for orthologous CDS lengths (default: %(default)s)")
    parser.add_argument("-ub", "--upper_bound", type=float, default=2.0,
                        help="Upper percentile bound for orthologous CDS lengths (default: %(default)s)")

    # 线程参数
    parser.add_argument("-t", "--threads", type=int, default=3,
                        help="Number of processing threads to use (default: %(default)s)")
    
    args = parser.parse_args()
    
    # ===== 参数验证与处理 =====
    # 1. 检查分位数范围有效性
    if not (0 <= args.lower_bound <= 1) or not (1 <= args.upper_bound <= 3):
        print("Error: Percentile bounds must be between 0 and 3")
        sys.exit(1)
    if args.lower_bound >= args.upper_bound:
        print("Error: Lower bound must be less than upper bound")
        sys.exit(1)

    # 2. 验证线程数
    if args.threads < 1:
        print("Warning: Thread count reset to minimum value 1")
        args.threads = 1
    
    # 3. 处理目录路径
    # 确保输入目录存在
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    # 创建输出目录（如果不存在）
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 4. 打印配置摘要
    print("\n" + "="*50)
    print(f"{'Configuration Summary':^50}")
    print("="*50)
    print(f"Input directory:     {args.input_dir}")
    print(f"Output directory:    {args.output_dir}")
    print(f"Reference type:      {'Angiosperm' if args.ref_orthologous == 'Ang' else 'Gymnosperm'}")
    print(f"Percentile range:    {args.lower_bound} - {args.upper_bound}")
    print(f"Processing threads:  {args.threads}")
    print("="*50 + "\n")

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
    
    time0 = datetime.now()

    if args.ref_orthologous == 'Ang':
        select_seq_by_len(args.input_dir, args.output_dir, angiosperm_ref_len_dict, args.lower_bound, args.upper_bound, args.threads)
    if args.ref_orthologous == 'Gym':
        select_seq_by_len(args.input_dir, args.output_dir, gymnosperm_ref_len_dict, args.lower_bound, args.upper_bound, args.threads)

    time1 = datetime.now()
    print(f"Total running time: {time1 - time0}")


"""
==================================================
              Configuration Summary
==================================================
Input directory:     H:\cpminer2\CPminer命令行代码\cds
Output directory:    H:\cpminer2\CPminer命令行代码\cds_filter
Reference type:      Angiosperm
Percentile range:    0.5 - 2.0
Processing threads:  10
==================================================
......
Total running time: 0:00:17.725347
"""



