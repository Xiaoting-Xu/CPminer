# *-* coding:utf-8 *-*
# @Time:2022/4/13 10:28
# @Author:Ruijing Cheng
# @File:coverage_and_identity.py
# @Software:PyCharm


from Bio import SeqIO
from Bio import AlignIO
import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from functools import reduce
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from datetime import datetime
import os
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import os
from pathlib import Path
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import concurrent.futures
import numpy as np
import multiprocessing
import argparse

def process_single_file(file, in_path, out_path, coverage_threshold, identity_threshold):
    """处理单个文件的函数，用于多线程执行"""
    try:
        print(f"Processing: {file}")
        alignments = AlignIO.read(in_path / file, "fasta")
        
        n_row = len(alignments)
        n_col = len(alignments[0])
        ref_len = n_col - str(alignments[0].seq).count("-")
        
        # 创建统计表格
        sum_table = pd.DataFrame(index=range(n_row), 
                                 columns=["description", "coverage", "identity", "filter_res"])
        col_list = [j for j in range(n_col) if alignments[0][j] != "-"]
        filtered_records = []
        
        for i in range(n_row):
            match = mismatch = 0
            for j in col_list:
                base = alignments[i][j]
                if base != "-":
                    if alignments[0][j] == base:
                        match += 1
                    else:
                        mismatch += 1
            
            total = match + mismatch
            coverage = total / ref_len if ref_len > 0 else 0
            identity = match / total if total > 0 else 0
            
            # 标记过滤结果
            filter_pass = (coverage >= coverage_threshold) and (identity >= identity_threshold)
            filter_res = "PASS" if filter_pass else "FAIL"
            
            sum_table.loc[i] = [
                alignments[i].description,
                f"{coverage:.4f}",
                f"{identity:.4f}",
                filter_res
            ]
            
            # 保留高质量序列
            if filter_pass:
                record = SeqRecord(
                    Seq(str(alignments[i].seq)),
                    id=alignments[i].id,
                    description=alignments[i].description
                )
                filtered_records.append(record)
        
        # 写入过滤后的FASTA
        if filtered_records:
            SeqIO.write(filtered_records, out_path / file, "fasta")
        
        # 保存单个文件的统计CSV
        stats_path = out_path / f"{Path(file).stem}_stats.csv"
        sum_table.to_csv(stats_path, sep=",", index=False)
        
        return stats_path  # 返回统计文件路径用于后续合并
    
    except Exception as e:
        print(f"Error processing {file}: {str(e)}")
        return None

def coverage_and_identity(in_path, out_path, coverage_threshold=0.75, identity_threshold=0.75, max_workers=None):
    """
    使用多线程处理FASTA文件，计算覆盖度和一致性并过滤序列
    
    参数:
    in_path: 输入目录路径
    out_path: 输出目录路径
    coverage_threshold: 覆盖度阈值 (默认0.75)
    identity_threshold: 一致性阈值 (默认0.75)
    max_workers: 最大线程数 (默认None，使用系统核心数)
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)
    
    # 获取所有FASTA文件
    fasta_files = [f for f in os.listdir(in_path) if f.endswith((".fa", ".fas", ".fasta"))]
    if not fasta_files:
        print("未找到FASTA文件")
        return
    
    # 使用线程池并行处理文件
    stats_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        futures = {executor.submit(
            process_single_file, 
            file, in_path, out_path, 
            coverage_threshold, identity_threshold
        ): file for file in fasta_files}
        
        # 收集结果
        for future in concurrent.futures.as_completed(futures):
            file = futures[future]
            try:
                result = future.result()
                if result:
                    stats_files.append(result)
                    print(f"完成处理: {file}")
            except Exception as e:
                print(f"处理 {file} 时出错: {str(e)}")
    



class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Retain sequences with both coverage and similarity to the reference sequence > threshold (default: 75%).
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
    parser.add_argument("-i", "--input_dir", help="Input directory of aligned CDS files (.fasta)", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for fileter CDS sequences (.fasta) and statistics summary of coverage&identity (.csv)", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=0.75, help="threshold of coverage&identity. default: %(default)s", required=False)
    parser.add_argument("-p", "--para_file", type=int, default=3, help="Number of parallel files to process, default: %(default)s", required=False)

    args = parser.parse_args()

        # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'CDS Filter Configuration':^50}")
    print("="*50)
    print(f"{'Input directory:':<30}{args.input_dir}")
    print(f"{'Output directory:':<30}{args.output_dir}")
    print(f"{'Coverage&Identity threshold:':<30}{args.threshold}")
    print(f"{'Parallel files:':<30}{args.para_file}")
    print("="*50 + "\n")

    # 确保输出目录存在
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if os.path.exists(args.input_dir):
        file_name_list = os.listdir(args.input_dir)
        if not file_name_list:  # 检查目录是否为空
            raise FileNotFoundError(f"input_dir is empty: {args.input_dir}")
    else:
        raise FileNotFoundError(f"input_dir not exists: {args.input_dir}")
    
    time0 = datetime.now()
    coverage_and_identity(args.input_dir, args.output_dir, args.threshold, args.threshold, args.para_file)
    time1 = datetime.now()
    print("Total running time: %s seconds" % (time1 - time0))

# Total running time: 0:02:04.451261 seconds