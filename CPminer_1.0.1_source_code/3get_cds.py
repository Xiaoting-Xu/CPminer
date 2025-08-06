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
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed




def process_gb_file(gb_file_name, in_path, out_path, cds_list, alternative_name, gene_locks):
    """处理单个GB文件的CDS提取（线程安全）"""
    # print(f"Extract CDS of : {gb_file_name}")
    # 初始化结果存储
    df0_row = {gene: 0 for gene in cds_list}
    df0_row['cds_num'] = 0
    df1_row = {gene: "" for gene in cds_list}
    df1_row['cds_num'] = ""
    
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            seq_record = SeqIO.read(Path(in_path) / Path(gb_file_name), "gb")
            
        for feature in seq_record.features:
            try:
                if feature.type in ["CDS", "rRNA"]:
                    gene_name = "unknown"
                    if "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name in alternative_name:
                            gene_name = alternative_name[gene_name]
                    elif "product" in feature.qualifiers:
                        product_name = feature.qualifiers["product"][0]
                        if product_name in alternative_name:
                            gene_name = alternative_name[product_name]
                    
                    # 获取基因位置并提取序列
                    if gene_name in cds_list:
                        fas_seq = ""
                        for part in feature.location.parts:
                            start = part.start
                            end = part.end
                            tmp_seq = seq_record.seq[start:end]
                            if part.strand == -1:
                                tmp_seq = tmp_seq.reverse_complement()
                            fas_seq += tmp_seq
                        
                        # 更新统计信息
                        if df0_row[gene_name] == 0:
                            df0_row["cds_num"] += 1
                            df1_row[gene_name] = str(feature.location)
                        else:
                            df1_row[gene_name] += "|%s" % str(feature.location)
                        df0_row[gene_name] += 1
                        
                        # 使用锁安全地写入文件
                        with gene_locks[gene_name]:
                            fas_description = f">{Path(gb_file_name).stem}|{feature.location}"
                            with open(Path(out_path) / f"{gene_name}.fasta", "a") as fw:
                                fw.write(f"{fas_description}\n{fas_seq}\n")
            except Exception as e:
                print(f"Error in feature: {e} - {gb_file_name} {gene_name}")
    except Exception as e:
        print(f"Error processing file: {e} - {gb_file_name}")
    
    # 设置cds_num为总数
    df1_row["cds_num"] = str(df0_row["cds_num"])
    return gb_file_name, df0_row, df1_row

def get_cds(in_path, out_path, threads=3):
    """多线程提取CDS和rRNA序列"""
    if not os.path.exists(Path(out_path)):
        os.makedirs(Path(out_path))
    
    # CDS基因列表
    cds_list = ['accD', 'atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI', 'ccsA', 'cemA',
                'clpP', 'infA', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI',
                'ndhJ', 'ndhK', 'petA', 'petB', 'petD', 'petG', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI',
                'psaJ', 'psbA', 'psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbI', 'psbJ', 'psbK', 'psbL',
                'psbM', 'psbN', 'psbT', 'psbZ', 'rbcL', 'rpl14', 'rpl16', 'rpl2', 'rpl20', 'rpl22', 'rpl23',
                'rpl32', 'rpl33', 'rpl36', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps14', 'rps15',
                'rps16', 'rps18', 'rps19', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4.5',
                'rrn5', 'ycf1', 'ycf15', 'ycf2', 'ycf3', 'ycf4', 'ycf68']

    # 基因名映射
    alternative_name = {'16S ribosomal RNA': 'rrn16', '23S ribosomal RNA': 'rrn23',
                        '4.5S ribosomal RNA': 'rrn4.5', '5S ribosomal RNA': 'rrn5',
                        'rrn16S': 'rrn16', 'rrn23S': 'rrn23',
                        'rrn4.5S': 'rrn4.5', 'rrn5S': 'rrn5'}
    
    # 创建初始空文件（确保文件存在）
    for cds in cds_list:
        with open(Path(out_path) / f"{cds}.fasta", "w"):
            pass
    
    # 获取GB文件列表
    gb_file_list = []
    for my_file in os.listdir(Path(in_path)):
        if Path(my_file).suffix.lower() in [".gb", ".gbf", ".gbk"]:
            gb_file_list.append(Path(my_file).name)
    
    print(f"Start extracting {len(cds_list)} CDS and rRNA in {len(gb_file_list)} gb files...")
    
    # 创建基因文件锁（每个基因一个锁）
    gene_locks = {gene: threading.Lock() for gene in cds_list}
    
    # 初始化结果DataFrame
    print(len(gb_file_list), len(cds_list)+1)
    locations = np.zeros((len(gb_file_list), len(cds_list)+1), dtype=int)
    df0 = pd.DataFrame(locations, index=[Path(x).stem for x in gb_file_list], columns=cds_list+['cds_num'], dtype=int)
    df1 = pd.DataFrame(locations, index=[Path(x).stem for x in gb_file_list], columns=cds_list+['cds_num'], dtype=str)
    
    # 设置线程数（默认为CPU核心数）
    if threads is None:
        threads = min(os.cpu_count(), len(gb_file_list))
    threads = max(1, min(threads, len(gb_file_list)))
    print(f"Using {threads} threads for processing")
    
    # 使用线程池处理文件
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {}
        for gb_name in gb_file_list:
            future = executor.submit(
                process_gb_file,
                gb_name,
                in_path,
                out_path,
                cds_list,
                alternative_name,
                gene_locks
            )
            futures[future] = os.path.splitext(gb_name)[0]
        
        # 收集处理结果
        for future in as_completed(futures):
            gb_name = futures[future]
            print(gb_name)
            try:
                _, df0_row, df1_row = future.result()
                # 更新结果DataFrame
                for gene in cds_list + ['cds_num']:
                    df0.loc[gb_name, gene] = df0_row[gene]
                    df1.loc[gb_name, gene] = df1_row[gene]

            except Exception as e:
                print(f"Error collecting results for {gb_name}: {e}")
    
    # 保存结果
    print("Writing summary...")
    df_info = pd.read_csv(Path(in_path) / "gb_info.csv", usecols=['filename', 'organism'])
    df0 = pd.merge(df_info, df0, left_on='filename', right_index=True, how='right')
    df0.to_csv(Path(out_path) / "cds_num.csv", index=False, sep=",")
    df1.to_csv(Path(out_path) / "cds_loc.txt", index=True, sep="\t")
    print("Summary of cds number written in cds_num.csv")
    print("Summary of cds location written in cds_loc.txt")
    
    # 清理空文件
    deleted_files = []
    for root, dirs, files in os.walk(out_path):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.getsize(file_path) == 0:
                try:
                    os.remove(file_path)
                    deleted_files.append(file_path)
                    print(f"Deleted empty file: {file_path}")
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")
    
    print(f"\nOperation completed, deleted {len(deleted_files)} empty files")

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Coding sequences (CDSs) are extracted from all GenBank files (.gb extension) within the input gb_folder_path.Each CDS gene is subsequently written to an individual FASTA file, where the filename corresponds to the standardized gene
        nomenclature. These segmented sequence files are systematically stored within the designated cds_folder_path.
        Author: Ruijing Cheng
        Organization: 
        GitHub Site:
        Usage Example:
          添加示例
        """
        return f"{description}\n\n{help_text}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--input_dir", help="Input directory path of gb file.", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for CDS sequence(.fasta)", required=True)
    parser.add_argument("-t", "--threads", type=int, default=3, help="Number of threads to use %(default)s", required=False)
    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'Configuration Summary':^50}")
    print("="*50)
    print(f"Input directory:     {args.input_dir}")
    print(f"Output directory:    {args.output_dir}")
    print(f"Processing threads:  {args.threads}")
    print("="*50 + "\n")



    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    threads = args.threads


    time0 = datetime.now()
    get_cds(input_dir, output_dir, threads)
    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))

"""
==================================================
              Configuration Summary
==================================================
Input directory:     H:\cpminer2\CPminer命令行代码\testFile
Output directory:    H:\cpminer2\CPminer命令行代码\cds
Processing threads:  3
==================================================

Start extracting 85 CDS and rRNA in 728 gb files...
Using 3 threads for processing
Writing summary...
Summary of cds number written in cds_num.csv
Summary of cds location written in cds_loc.txt

Operation completed, deleted 0 empty files
Total running time: 0:00:27.125098 seconds
"""