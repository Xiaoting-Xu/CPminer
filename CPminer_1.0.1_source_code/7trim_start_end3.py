# *-* coding:utf-8 *-*
# @Time:2022/6/16 16:21
# @Author:Ruijing Cheng
# @File:trim_start_end3.py
# @Software:PyCharm


import os
from Bio import AlignIO
from pathlib import Path
import os
import sys
import subprocess
from datetime import datetime
from pathlib import Path
import shutil
import multiprocessing
import argparse



def is_trimAL_installed(trimal_path: str) -> bool:
    """检查是否安装并可用"""
    if not os.path.exists(trimal_path):
        return False
    try:
        result = subprocess.run([trimal_path, "--version"], 
                              capture_output=True, 
                              text=True)
        return result.returncode == 0
    except:
        return False


def trim_start_end(boundaries_threshold, in_file_path, out_file_path, ):
    """从碱基gap比例小于2.5%的位点开始保留"""
    
    alignments = AlignIO.read(Path(in_file_path), "fasta")
    n_row = len(alignments)
    n_col = len(alignments[0])

    i = 0
    while True:
        n_gap = alignments[:, i].count("-")
        if n_gap/n_row < boundaries_threshold:
            break
        i += 1  # 18

    j = -1
    while True:
        n_gap = alignments[:, j].count("-")
        if n_gap/n_row < boundaries_threshold:
            break
        j -= 1  # -19

    fw = open(Path(out_file_path), "w")
    for n in range(n_row):
        fw.write(">"+alignments[n].description+"\n")
        fw.write(str(alignments[n].seq[i:j]+alignments[n].seq[j])+"\n")
    fw.close()



def trim_msa(trimal_path, in_file_path, out_file_path):
    print("trimal:"+in_file_path)

    t0 = datetime.now()
    cmd_str = "%s -in %s -out %s -fasta -automated1"
    tmp_cmd = cmd_str % (trimal_path, str(in_file_path), str(out_file_path))
    print("Command: %s" % tmp_cmd)
    os.system(tmp_cmd)
    t1 = datetime.now()
    print("Running time: %s Seconds" % (t1 - t0))


def process_bound_trimal(args):
    """包装函数用于多进程调用"""
    return trim_start_end(*args)
def process_interior_trimal(args):
    """包装函数用于多进程调用"""
    return trim_msa(*args)


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Removes unaligned terminal regions where gap proportions exceed the predefined threshold, And remove internally unaligned fragments where gap proportions exceed the threshold in the sequence matrix.
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
    parser.add_argument("-i", "--input_dir", help="Input directory of CDS files (.fasta)", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory of trimAL CDS sequences (.fasta)", required=True)
    # 阈值
    parser.add_argument("-bt", "--boundaries_threshold", type=float, default=0.025,
                        help="Lower percentile bound for orthologous CDS lengths (default: %(default)s)")
    parser.add_argument("-tp", "--trimal_path", type=str, help="trimAL path", required=False)
    parser.add_argument("-p", "--para_file", type=int, default=3, help="Number of parallel files to process, default: %(default)s", required=False)

    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'TrimAL Configuration':^50}")
    print("="*50)
    print(f"{'Input directory:':<25}{args.input_dir}")
    print(f"{'Output directory:':<25}{args.output_dir}")
    print(f"{'Boundaries threshold:':<25}{args.boundaries_threshold}")
    print(f"{'TrimAL path:':<25}{args.trimal_path}")
    print(f"{'Parallel files:':<25}{args.para_file}")
    
    print("="*50 + "\n")

    # 确保输出目录存在
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    
    if not os.path.exists(args.input_dir):
        raise FileNotFoundError(f"input_dir not exists: {args.input_dir}")
    
    file_name_list = os.listdir(args.input_dir)
    if not file_name_list:  # 检查目录是否为空
        raise FileNotFoundError(f"input_dir is empty: {args.input_dir}")

    # 准备任务列表
    tasks = []
    for in_file_name in file_name_list:
        if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
            in_file_path = os.path.join(args.input_dir, in_file_name)
            out_file_path = os.path.join(args.output_dir, Path(in_file_name).name)
            tasks.append((args.boundaries_threshold, in_file_path, out_file_path))

    t0 = datetime.now()
    # 并行处理
    if args.para_file > 1:
        print(f"Starting parallel processing with {args.para_file} workers...")
        with multiprocessing.Pool(processes=args.para_file) as pool:
            results = pool.map(process_bound_trimal, tasks)
        
        # 检查结果
        success_count = sum(1 for r in results if r)
        print(f"\nAlignment boundaries completed: {success_count}/{len(tasks)} files succeeded")
        if success_count < len(tasks):
            print("Warning: Some files failed. Check error messages above.")
    else:
        # 单线程处理
        success_count = 0
        for task in tasks:
            if process_bound_trimal(task):
                success_count += 1
        print(f"\nAlignment completed: {success_count}/{len(tasks)} files succeeded")
    t1 = datetime.now()
    print("Total running time: %s seconds" % (t1 - t0))

    if args.trimal_path:
        if not is_trimAL_installed(args.trimal_path):
            raise FileNotFoundError(f"trimAL not found or not working at {args.trimal_path}")

        file_name_list = os.listdir(args.output_dir)
        # 准备任务列表
        tasks = []
        for in_file_name in file_name_list:
            if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
                in_file_path = os.path.join(args.output_dir, in_file_name)
                out_file_path = os.path.join(args.output_dir, Path(in_file_name).name)
                tasks.append((args.trimal_path, in_file_path, out_file_path))

        t0 = datetime.now()
        # 并行处理
        if args.para_file > 1:
            print(f"Starting parallel processing with {args.para_file} workers...")
            with multiprocessing.Pool(processes=args.para_file) as pool:
                results = pool.map(process_interior_trimal, tasks)
            
            # 检查结果
            success_count = sum(1 for r in results if r)
            print(f"\nTrimAL process completed: {success_count}/{len(tasks)} files succeeded")
            if success_count < len(tasks):
                print("Warning: Some files failed trimAL. Check error messages above.")
        else:
            # 单线程处理
            success_count = 0
            for task in tasks:
                if process_bound_trimal(task):
                    success_count += 1
            print(f"\nTrimAL process completed: {success_count}/{len(tasks)} files succeeded")
        t1 = datetime.now()
        print("Total running time: %s seconds" % (t1 - t0))







    

# "H:\cpminer\mafft\mafft\mafft-win\mafft.bat"
# trim Total running time: 0:01:30.484482 seconds
# sr 3s
