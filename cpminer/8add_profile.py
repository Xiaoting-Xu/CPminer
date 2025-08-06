# *-* coding:utf-8 *-*
# @Time:2022/4/2 17:19
# @Author:Ruijing Cheng
# @File:add_profile.py
# @Software:PyCharm


import os
from datetime import datetime
from pathlib import Path
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
import tempfile


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



def check_file(file_path):
    """检查文件是否存在且非空"""
    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在: {file_path}")
        return False
    
    # 检查是否为文件（非目录）
    if not os.path.isfile(file_path):
        print(f"❌ 路径不是文件: {file_path}")
        return False
    
    # 检查文件大小
    file_size = os.path.getsize(file_path)
    if file_size == 0:
        print(f"⚠️ 文件存在但为空 (0字节): {file_path}")
        return False
    
    print(f"✅ 文件存在且非空: {file_path} (大小: {file_size} 字节)")
    return True



def add_profile(mafft_path, mafft_thread, profile1, profile2):
    """将新序列添加到已有的多序列比对结果中"""

    if not check_file(profile1) and not check_file(profile2):
        return

    with tempfile.NamedTemporaryFile(
        mode='w', 
        suffix=".fasta", 
        delete=False,
        dir=os.path.dirname(profile1) or None
    ) as tmp_file:
        temp_output = tmp_file.name

    cmd_str = str(Path(mafft_path)) + " --thread %s --reorder --add %s %s > %s"
    t0 = datetime.now()
    print(f"Merging: {profile1} and {profile2}")
    tmp_cmd = cmd_str % ( mafft_thread, profile1, profile2, temp_output)
    print(tmp_cmd)
    
    
    if not os.system(tmp_cmd) and check_file(temp_output):
        t1 = datetime.now()
        print(f"Merging completed in {(t1 - t0).total_seconds():.2f} seconds")

        try:
            shutil.move(temp_output, profile1)
            return True
        except Exception as e:
            return False




    




def process_file_wrapper(args):
    """包装函数用于多进程调用"""
    return add_profile(*args)


def find_executable_path(executable_name):
    """查找可执行文件的完整路径"""
    
    # 方法2：手动检查 PATH 环境变量（兼容性更强）
    path_values = os.environ.get('PATH', '').split(os.pathsep)
    #print(path_values)
    for dir_path in path_values:
        if not dir_path:
            continue
        full_path = os.path.join(dir_path, executable_name)
        if os.path.isfile(full_path):
            return os.path.abspath(full_path)
    return None

def is_path_installed(path: str) -> bool:
    """检查perl是否安装并可用"""
    if not os.path.exists(path):
        return False
    try:
        result = subprocess.run([path, "-version"], capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False






class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Align matching FASTA files with MAFFT and merge into consolidated files preserving original headers.
        Author: Ruijing Cheng
        Organization: 
        GitHub Site:
        Usage Example:
        8add_profile.py: error: the following arguments are required: -i/--input_dir, -a/--add_profile, -o/--output_dir, -mp/--mafft_path
        """
        return f"{description}\n\n{help_text}"




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    
    # 输入输出参数
    parser.add_argument("-i", "--input_dir", help="Input directory of align CDS files (.fasta)", required=True)
    parser.add_argument("-a", "--add_profile", help="Input directory of alingn CDS files (.fasta) need to merge", required=True)
    parser.add_argument("-p", "--para_file", type=int, default=3, help="Number of parallel files to process, default: %(default)s", required=False)
    parser.add_argument("-t", "--thread", type=int, default=1, help="Mafft thread process option, default: %(default)s", required=False)

    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'MAFFT Alignment Configuration':^50}")
    print("="*50)
    print(f"{'Input directory:':<25}{args.input_dir}")
    print(f"{'Add Profile:':<25}{args.add_profile}")
    print(f"{'Parallel files:':<25}{args.para_file}")
    print(f"{'mafft thread:':<25}{args.thread}")
    print("="*50 + "\n")

    # 检测mafft路径

    mafft_path = find_executable_path('mafft') or find_executable_path('mafft.bat')
    # 检测perl路径是否有效
    if mafft_path and is_path_installed(mafft_path):
        print(f"Mafft path:{mafft_path}")
    else:
        mafft_path = input("No valid Mafft runtime path was found in the system variables. Please enter the Mafft installation directory(./mafft/mafft.bat):")
        while True:
            if is_path_installed(mafft_path):
                break
            else:
                mafft_path = input("The entered Mafft path is invalid. Please re-enter it(./mafft/mafft.bat):")


    if os.path.exists(args.add_profile):
        file_name_list = os.listdir(args.add_profile)
        if not file_name_list:  # 检查目录是否为空
            raise FileNotFoundError(f"add_profile is empty: {args.add_profile}")
    else:
        raise FileNotFoundError(f"add_profile not exists: {args.add_profile}")


    if os.path.exists(args.input_dir):
        file_name_list = os.listdir(args.input_dir)
        if not file_name_list:  # 检查目录是否为空
            raise FileNotFoundError(f"input_dir is empty: {args.input_dir}")

        # 准备任务列表
        tasks = []
        for in_file_name in file_name_list:
            if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
                in_file_path1 = os.path.join(args.input_dir, in_file_name)
                in_file_path2 = os.path.join(args.add_profile, in_file_name)
                if  os.path.exists(in_file_path2):
                    print(in_file_path2)
                    tasks.append((mafft_path, args.thread, in_file_path1, in_file_path2))
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