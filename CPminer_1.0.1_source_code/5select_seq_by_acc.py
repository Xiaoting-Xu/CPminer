# -*- coding = utf-8 -*-
# @Time : 2021/9/5 22:35
# @Author : Ruijing Cheng
# @File : select_seq_by_acc.py
# @Software: PyCharm

import concurrent.futures
from Bio import SeqIO
import os
import datetime
import pandas as pd
from pathlib import Path
import argparse
from datetime import datetime



def process_file(file, in_path, out_path, selected_table):
    """处理单个文件的辅助函数"""
    if file.split(".")[-1].lower() not in ["fa", "fas", "fasta"]:
        return  # 跳过非fasta文件
    
    # print(f"Start processing: {file}")
    time1 = datetime.now()
    output_file = Path(out_path) / file
    seq_dict = SeqIO.to_dict(
        SeqIO.parse(Path(in_path) / file, "fasta"),
        key_function=lambda record: record.description.split("|")[0]
    )
    
    with open(output_file, "w") as fw:
        for index in selected_table.index:
            if index in seq_dict:
                species_name = selected_table.loc[index]["organism"].replace(" ", "_")
                fas_description = f">{index}|{species_name}"
                fw.write(fas_description + "\n")
                fw.write(str(seq_dict[index].seq) + "\n")
    
    time2 = datetime.now()
    print(f"Finished processing {file} in {(time2 - time1).total_seconds():.2f} seconds")

def select_seq_by_acc(in_path: str, out_path: str, selected_table: pd.DataFrame, num_threads: int = 3):
    """使用多线程处理文件"""
    selected_table = selected_table.set_index('filename')
    file_list = os.listdir(in_path)
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    # 使用线程池处理文件
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for file in file_list:
            futures.append(
                executor.submit(
                    process_file,
                    file, in_path, out_path, selected_table
                )
            )
        
        # 等待所有任务完成
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error processing file: {e}")

def make_tab(in_path, table_path):
    df =  pd.read_csv(Path(in_path) / "length.csv", sep=',')
    df1 = df[["filename", "organism", "cds_num"]].copy()
    df1['length'] = df.loc[:, 'accD':'ycf68'].sum(axis=1)

    if args.file_organism_name:
        df_rename = pd.read_csv(Path(table_path))
        df1 = pd.merge(df1, df_rename, on="organism", how="right")
        print("Organism name standardization completed")
    
    duplicate_counts = df1['organism'].value_counts().loc[lambda x: x > 1]
    if len(duplicate_counts):
            while True:
                reannotation_choice = input(f"{len(duplicate_counts)} organism(s) have duplicate values. Continue to retain the longest? (y/n): ")
                # 标准化输入并检查
                normalized_input = reannotation_choice.strip().lower()
                if normalized_input in ('y', 'yes'):
                    print("Proceeding to filter...")
                    return df1.loc[df1.groupby('organism')['length'].idxmax()]
                    # 执行重新注释的逻辑
                elif normalized_input in ('n', 'no'):
                    print("Program execution terminated.")
                    return None
                else:
                    print("Invalid input.")
    else:
        print("No duplicate samples exist.")
        return df1.loc[df1.groupby('organism')['length'].idxmax()]


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Correct the organism's scientific names and select the representative chloroplast genome for each species based on the maximum cumulative CDS length.
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
    parser.add_argument("-i", "--input_dir", help="Input directory of filtered CDS files (.fasta)", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for CDS sequences (.fasta) of unique organisms", required=True)
    parser.add_argument("-f", "--file_organism_name", help="Input path of standardized organism name file(.csv)", required=False)
    parser.add_argument("-t", "--threads", type=int, default=3, help="Number of threads to use %(default)s", required=False)

    args = parser.parse_args()
    
    # 验证输入目录存在
    if not os.path.isdir(args.input_dir):
        raise NotADirectoryError(f"Input directory does not exist: {args.input_dir}")

    print("\n" + "="*50)
    print(f"{' CONFIGURATION SUMMARY ':=^50}")
    print(f"{'Input directory:':<20} {os.path.abspath(args.input_dir)}")
    print(f"{'Output directory:':<20} {os.path.abspath(args.output_dir)}")
    print(f"{'Organism names:':<20} {args.file_organism_name or 'Not specified'}")
    print(f"{'Threads:':<20} {args.threads}")
    print("="*50 + "\n")

    if args.file_organism_name:
        if not os.path.isfile(args.file_organism_name):
            raise FileNotFoundError(f"Organism name file not found: {args.file_organism_name}")

    
    time0 = datetime.now()
    # 矫正名称
    df_organism = make_tab(args.input_dir, args.file_organism_name)
    if not df_organism is None:
        # 将线程数传递给函数
        select_seq_by_acc(
            args.input_dir, 
            args.output_dir, 
            df_organism,
            num_threads=args.threads  # 添加线程数参数
        )

    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))


"""

==================================================
============= CONFIGURATION SUMMARY ==============
Input directory:     H:\cpminer2\CPminer命令行代码\cds_filter
Output directory:    H:\cpminer2\CPminer命令行代码\cds_uni
Organism names:      Not specified
Threads:             3
==================================================

No duplicate samples exist.
Finished processing atpA.fasta in 0.13 seconds
.....
Total running time: 0:00:03.856980 seconds

"""




