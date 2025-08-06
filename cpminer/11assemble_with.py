
from Bio import SeqIO
from pathlib import Path
import os, sys
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import numpy as np
from datetime import datetime
import shutil




def create_seq_dict(fasta_path):
    def get_accession(record):
        parts = record.description
        return parts
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"), key_function=get_accession)
    return seq_dict


def create_length_df(seq_dict):
    """ 提取fasta信息转换为df """
    data = np.zeros((len(seq_dict), 5), dtype=str)
    df = pd.DataFrame(data, columns=["accession", "organism", "description","length","seq"], dtype=str)
    i = 0
    for key in seq_dict.keys():
        df.loc[i, "accession"] = key.split(":")[0]  # accession
        df.loc[i, "organism"] = key.split("|")[1]  # accession
        df.loc[i, "description"] = key  # accession
        df.loc[i, "length"] = len(seq_dict[key].seq)
        df.loc[i, "seq"] = seq_dict[key].seq
        i += 1
    df["length"] = pd.to_numeric(df["length"])
    return df

# 检测df有没有重复值
# 去重保留最长的
# 根据要求进行替换
import pandas as pd

def filter(in_file_path):
    df = create_length_df(create_seq_dict(in_file_path))
    # 检查organism列是否有重复值
    if df['organism'].duplicated().any():
        print("存在重复的organism值，进行去重处理")
        
        # 按organism分组，并在每个分组内：
        # 1. 保留length最大的一行
        # 2. 如果length相同，则保留第一行
        result = df.sort_values('length', ascending=False)  \
            .drop_duplicates('organism')    \
            .sort_index()
    else:
        print("没有重复的organism值")
        result = df

    return result

def replace(in_file_path1, in_file_path2, out_file_dir, mode):
    """ 合并df """
    df1 = create_length_df(create_seq_dict(in_file_path1))
    df2 = create_length_df(create_seq_dict(in_file_path2))

    mode = int(mode)
    if mode == 1:
        # 将两个DataFrame垂直拼接
        combined = pd.concat([df1, df2])
        # 按organism去重，保留第一个出现的（即df1中的行）
        result = combined.drop_duplicates(subset='organism', keep='first')
        print("1")

    if mode == 2:
        # 将两个DataFrame垂直拼接
        combined = pd.concat([df2, df1])
        # 按organism去重，保留第一个出现的（即df1中的行）
        result = combined.drop_duplicates(subset='organism', keep='first')

    if mode == 3:
        # 将两个DataFrame垂直拼接
        combined = pd.concat([df2, df1])
        # 按organism去重，保留第一个出现的（即df1中的行）
        result = combined.loc[combined.groupby("organism")['length'].idxmax()]

    with open(Path(out_file_dir) / Path(in_file_path1).name, "w") as out_file:
        for index, row in result.iterrows():
            out_file.write(f">{row['description']}\n{row['seq']}\n")

def merge(in_dir1, in_dir2, out_dir = None, mode = None):
    """判断是否有重复基因，重复合并，不是重复copy"""
    fasta_files1 = [file for file in os.listdir(in_dir1) if file.endswith(".fasta")]
    fasta_files2 = [file for file in os.listdir(in_dir2) if file.endswith(".fasta")]

    print(fasta_files1)
    print(fasta_files2)

    for file in set(fasta_files1) - set(fasta_files2):
        shutil.copy2(Path(in_dir1) / file, Path(out_dir) / file)  # 保留元数据

    for file in set(fasta_files2) - set(fasta_files1):
        shutil.copy2(Path(in_dir2) / file, Path(out_dir) / file)  # 保留元数据

    common_files = set(fasta_files1) & set(fasta_files2)
    for file in common_files:
        replace(Path(in_dir1) / file, Path(in_dir2) / file, out_dir, mode)
        
        







class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: assemble other gene to the composite matrix
        Author: Ruijing Cheng
        Organization: 
        GitHub Site:
        Usage Example:
          添加示例
        """
        return f"{description}\n\n{help_text}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sequence assembly tool')

    # 使用单破折号短选项和双破折号长选项
    parser.add_argument("-c", "--input_CDS_dir", 
                        help="Input directory path of CDS sequence(.fasta)", 
                        required=True,
                        metavar="DIR")

    parser.add_argument("-g", "--input_other_gene_dir", 
                        help="Input directory for other gene sequence(.fasta)", 
                        required=True,
                        metavar="DIR")

    parser.add_argument("-o", "--output_dir", 
                        help="Output directory for assembled sequences",
                        required=True,
                        metavar="OUTPUT_DIR")

    parser.add_argument("-m", "--mode", 
                        help="Sequence replacement method",
                        required=True)

    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'Configuration Summary':^50}")
    print("="*50)
    print(f"CDS input directory:    {args.input_CDS_dir}")
    print(f"Other gene directory:   {args.input_other_gene_dir}")
    print(f"Output directory:       {args.output_dir if args.output_dir else 'Not specified'}")
    print(f"Mode:   {args.mode}")
    print("="*50 + "\n")

    cds_dir = os.path.abspath(args.input_CDS_dir)
    other_dir = os.path.abspath(args.input_other_gene_dir)
    output_dir = os.path.abspath(args.output_dir)
    mode = args.mode

    # 验证输入目录是否存在
    if not os.path.isdir(cds_dir):
        raise ValueError(f"CDS directory not found: {cds_dir}")
    if not os.path.isdir(other_dir):
        raise ValueError(f"Other genes directory not found: {other_dir}")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    merge(cds_dir, other_dir, output_dir, mode) 

    

"""
流程：
获取输入参数：叶绿体基因文件夹，小片段文件夹，balstn地址，替换模式，
验证文件夹是否有效：
    叶绿体文件夹：验证是否有.fasta文件，（验证基因是否都为叶绿体基因）

将替换文件夹的预处理提取出来， 

"""