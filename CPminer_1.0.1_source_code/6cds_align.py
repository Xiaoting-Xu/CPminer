import os
import sys
import subprocess
from datetime import datetime
from pathlib import Path
import shutil
import multiprocessing
import argparse

def is_mafft_installed(mafft_path: str) -> bool:
    """检查MAFFT是否安装并可用"""
    if not os.path.exists(mafft_path):
        return False
    try:
        result = subprocess.run([mafft_path, "--version"], 
                              capture_output=True, 
                              text=True)
        return result.returncode == 0
    except:
        return False

def run_mafft_safely(mafft_path, mafft_option, input_file, output_file,):
    """安全运行MAFFT多序列比对"""
    try:
        cmd_str = str(Path(mafft_path)) + " %s %s > %s"
        tmp_cmd = cmd_str % (mafft_option, input_file, output_file)
        print(f"Aligning: {input_file}")
        os.system(tmp_cmd)
        return True
    except Exception as result:
        print("Error: %s || Please check %s" % (result, input_file))
        return False

def add_profile(mafft_path, profile1, profile2, out_path):
    """将新序列添加到已有的多序列比对结果中"""
    cmd_str = str(Path(mafft_path)) + " --thread 1 --reorder --add %s %s > %s"
    t0 = datetime.now()
    print(f"Merging: {profile1} and {profile2}")
    tmp_cmd = cmd_str % (Path(profile1), Path(profile2), Path(out_path))
    os.system(tmp_cmd)
    t1 = datetime.now()
    print(f"Merging completed in {(t1 - t0).total_seconds():.2f} seconds")

def process_file_wrapper(args):
    """包装函数用于多进程调用"""
    return run_mafft_safely(*args)

# 自定义帮助信息格式化类
class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        # 创建装饰性标题
        help_header = "\n" + "=" * 50 + "\n"
        help_header += f"{'MAFFT Alignment Help':^50}\n"
        help_header += "=" * 50 + "\n"
        
        # 获取标准帮助信息
        standard_help = """options:
  -h, --help            Show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Input directory of CDS files (.fasta)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory of aligned CDS sequences (.fasta)
  -p PARA_FILE, --para_file PARA_FILE
                        Number of parallel files to process, default: 3
  -mp MAFFT_PATH, --mafft_path MAFFT_PATH
                        Mafft path
  -mo MAFFT_OPTION, --mafft_option MAFFT_OPTION
                        Mafft process option, default: --auto --reorder
"""
        mafft_help = """
mafft options:
    --op # :         Gap opening penalty, default: 1.53
    --ep # :         Offset (works like gap extension penalty), default: 0.0
    --maxiterate # : Maximum number of iterative refinement, default: 0
    --clustalout :   Output: clustal format, default: fasta
    --reorder :      Outorder: aligned, default: input order
    --quiet :        Do not report progress
    --thread # :     Number of threads (if unsure, --thread -1)
    --dash :         Add structural information (Rozewicki et al, submitted)
"""
# 完善这里的内容，描述和示例
        
        # 组合自定义标题和标准帮助内容
        return help_header + standard_help + "\n" + mafft_help + "=" * 50 + "\n"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomHelpFormatter,  # 应用自定义格式化器
        add_help=False  # 禁用默认的help，以便手动添加
    )
    
    # 手动添加帮助参数（使用自定义格式化）
    parser.add_argument(
        "-h", "--help", 
        action="help",  # 特殊action触发帮助
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    
    # 输入输出参数
    parser.add_argument("-i", "--input_dir", help="Input directory of CDS files (.fasta)", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory of aligned CDS sequences (.fasta)", required=True)
    parser.add_argument("-p", "--para_file", type=int, default=3, help="Number of parallel files to process, default: %(default)s", required=False)
    parser.add_argument("-mp", "--mafft_path", type=str, help="Mafft path", required=True)
    parser.add_argument("-mo", "--mafft_option", type=str, default="--auto --reorder", help="Mafft process option, default: %(default)s", required=False)

    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'MAFFT Alignment Configuration':^50}")
    print("="*50)
    print(f"{'Input directory:':<25}{args.input_dir}")
    print(f"{'Output directory:':<25}{args.output_dir}")
    print(f"{'Parallel files:':<25}{args.para_file}")
    print(f"{'mafft path:':<25}{args.mafft_path}")
    print(f"{'mafft option:':<25}{args.mafft_option}")
    print("="*50 + "\n")

    if not is_mafft_installed(args.mafft_path):
        raise FileNotFoundError(f"MAFFT not found or not working at {args.mafft_path}")

    # 确保输出目录存在
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    
    if os.path.exists(args.input_dir):
        file_name_list = os.listdir(args.input_dir)
        if not file_name_list:  # 检查目录是否为空
            raise FileNotFoundError(f"input_dir is empty: {args.input_dir}")

        # 准备任务列表
        tasks = []
        for in_file_name in file_name_list:
            if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
                in_file_path = os.path.join(args.input_dir, in_file_name)
                out_file_path = os.path.join(args.output_dir, Path(in_file_name).name)
                tasks.append((args.mafft_path, args.mafft_option, in_file_path, out_file_path))

        t0 = datetime.now()
        # 并行处理
        if args.para_file > 1:
            print(f"Starting parallel processing with {args.para_file} workers...")
            with multiprocessing.Pool(processes=args.para_file) as pool:
                results = pool.map(process_file_wrapper, tasks)
            
            # 检查结果
            success_count = sum(1 for r in results if r)
            print(f"\nAlignment completed: {success_count}/{len(tasks)} files succeeded")
            if success_count < len(tasks):
                print("Warning: Some files failed alignment. Check error messages above.")
        else:
            # 单线程处理
            success_count = 0
            for task in tasks:
                if process_file_wrapper(task):
                    success_count += 1
            print(f"\nAlignment completed: {success_count}/{len(tasks)} files succeeded")
        t1 = datetime.now()
        print("Total running time: %s seconds" % (t1 - t0))

    else:
        raise FileNotFoundError(f"input_dir not exists: {args.input_dir}")
    

# "H:\cpminer\mafft\mafft\mafft-win\mafft.bat"
"""
Alignment completed: 89/89 files succeeded
Total running time: 0:05:05.063862 seconds
"""
