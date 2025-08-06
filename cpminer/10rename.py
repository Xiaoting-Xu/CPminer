"""
    小片段文件夹：默认其中文件名称为：基因名称.fasta
        检验文件夹是否存在，检验文件夹是否为空(必须包含.fasta文件)，
        fasta文件必须对其，record的长度必须相同，不相同报错
        文件命名会给出说明文档，进行一个自动化矫正
        叶绿体基因只认可说明文档中的内容，对叶绿体基因进行比对，比对不符合报错
        在说明文档中说明序列中给出record.id的格式要求，
        是否为替换文档添加名称矫正功能，经行矫正时如何确定我提取到了合适位置

"""

import re
import sys
from Bio import SeqIO
import argparse
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
import multiprocessing

def check_sequence_lengths(fasta_file):
    """
    检查FASTA文件中所有序列是否长度一致
    """
    lengths = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.add(len(record.seq))
        if len(lengths) > 1:
            return (False, None)
    
    return True if lengths else False


def parse_header_format(format_string):
    """
    解析用户输入的格式字符串，提取字段顺序和分隔符
    :param format_string: 用户输入的格式字符串（如 "accession:location|organism|description"）
    :return: 元组 (字段顺序列表, 分隔符列表)
    """
    # 预定义所有可能的字段名
    field_names = [
        'accession', 
        'location', 
        'organism', 
        'description', 
        'other information'
    ]
    
    # 按长度降序排序以确保优先匹配长字段名（如"other information"）
    sorted_fields = sorted(field_names, key=len, reverse=True)
    pattern = '|'.join(map(re.escape, sorted_fields))
    
    # 分割格式字符串
    tokens = re.split(f'({pattern})', format_string)
    
    # 移除首尾空字符串
    if tokens and tokens[0] == '':
        tokens = tokens[1:]
    if tokens and tokens[-1] == '':
        tokens = tokens[:-1]
    
    # 提取字段顺序和分隔符
    fields = tokens[0::2]   # 偶数索引位置是字段名
    separators = tokens[1::2]  # 奇数索引位置是分隔符
    
    # 验证字段有效性
    for field in fields:
        if field not in field_names:
            raise ValueError(f"未知字段名: '{field}'。可用字段: {field_names}")
    
    # 验证分隔符数量
    if len(separators) != len(fields) - 1:
        raise ValueError("字段与分隔符数量不匹配")
    
    return fields, separators

def parse_record_id(record_id, fields, separators):
    """
    根据解析的格式解析单个记录ID
    :param record_id: FASTA记录ID字符串
    :param fields: 字段顺序列表
    :param separators: 分隔符列表
    :return: 字典 {字段名: 对应值}
    """
    parts = []
    remaining = record_id
    
    # 使用分隔符逐步分割字符串
    for sep in separators:
        if not remaining:
            raise ValueError("记录ID过早结束")
            
        # 查找分隔符位置
        index = remaining.find(sep)
        if index == -1:
            raise ValueError(f"未找到分隔符: '{sep}'")
        
        # 提取当前字段并更新剩余字符串
        parts.append(remaining[:index])
        remaining = remaining[index + len(sep):]
    
    parts.append(remaining)  # 添加最后一个字段
    
    # 组合字段名和值
    return dict(zip(fields, parts))


def export_fasta_info(filename, header_format):
    """
    处理FASTA文件并提取所需字段
    :param filename: FASTA文件路径
    :param header_format: 用户提供的格式字符串
    :return: 包含解析结果的DataFrame
    """
    try:
        fields, separators = parse_header_format(header_format)
        print(f"解析格式成功: 字段={fields}, 分隔符={separators}")
    except ValueError as e:
        print(f"格式解析错误: {e}")
        return pd.DataFrame()  # 返回空DataFrame
    
    # 确保所需字段存在
    required_fields = {'accession', 'location', 'organism'}
    missing_fields = required_fields - set(fields)
    if missing_fields:
        print(f"错误: 格式中缺少必要字段: {missing_fields}")
        return pd.DataFrame()
    
    # 存储解析结果的列表
    results = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line[1:].strip()
                try:
                    parsed = parse_record_id(header, fields, separators)
                    # 只提取需要的三个字段
                    result = {
                        'accession': parsed['accession'],
                        'location': parsed['location'],
                        'organism': parsed['organism']
                    }
                    results.append(result)
                except Exception as e:
                    print(f"解析失败: '{header}'\n原因: {str(e)}")
    
    # 创建DataFrame
    df = pd.DataFrame(results)
    return df


def correct_fasta_info(in_file_path, header_format, out_file_path, df: pd.DataFrame):
    """
    处理FASTA文件并提取所需字段
    :param filename: FASTA文件路径
    :param header_format: 用户提供的格式字符串
    :return: 包含解析结果的DataFrame
    """
    try:
        fields, separators = parse_header_format(header_format)
        print(f"解析格式成功: 字段={fields}, 分隔符={separators}")
    except ValueError as e:
        print(f"格式解析错误: {e}")
    
    # 确保所需字段存在
    required_fields = {'accession', 'location', 'organism'}
    missing_fields = required_fields - set(fields)
    if missing_fields:
        print(f"错误: 格式中缺少必要字段: {missing_fields}")

    seq_dict = SeqIO.to_dict(
        SeqIO.parse(in_file_path, "fasta"),
        key_function=lambda record: record.description
    )

    with open(out_file_path, 'w') as file:
        for key, record in seq_dict.items():  # 正确遍历字典项
            parsed = parse_record_id(record.description, fields, separators)
            accession = parsed['accession']
            acc_rows = df[df['accession'] == accession]
            if not acc_rows.empty:
                file.write(f">{acc_rows['accession'].values[0]}:{acc_rows['location'].values[0]}|{acc_rows['organism'].values[0]}\n")
                file.write(record.seq)

def process_file_export(args):
    """包装函数用于多进程调用"""
    return export_fasta_info(*args)

def process_file_correct(args):
    """包装函数用于多进程调用"""
    return correct_fasta_info(*args)



def get_valid_input():
    while True:
        user_input = input("please input your choice(0/1/2): ")
        if user_input in ['0', '1', '2']:
            return int(user_input)
        print("Input invalid! Please try again.")

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Preprocessing of the FASTA files to be assemble
        Author: Li Dan
        Organization: 
        GitHub Site:
        Usage Example:
          添加示例
        """
        return f"{description}\n\n{help_text}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--input", help="Input path of file or directory(.fasta).", required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for CDS sequence(.fasta)", required=True)
    parser.add_argument("-p", "--para_file", type=int, default=3, help="When handling multiple files, choose the thread count for parallel execution. %(default)s", required=False)
    args = parser.parse_args()



    if os.path.isdir(args.input):

        file_name_list = os.listdir(args.input)
        if not file_name_list:  # 检查目录是否为空
            raise FileNotFoundError(f"input_dir is empty: {args.input}")

        print("Select operation:  \n1 - Export records information \n2 - Fix organism annotations \n0 - stop execution ")
        choice = get_valid_input()

        if choice == 0:
            print("Execution stopped.")
            sys.exit(1)
        
        if choice == 1:
            header_format = input("Please enter the format of record.id, for example:\n\t \"accession:location|organism|description\"\n your input:")

            # 准备任务列表
            tasks = []
            for in_file_name in file_name_list:
                if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
                    in_file_path = os.path.join(args.input, in_file_name)
                    tasks.append((in_file_path, header_format))

            t0 = datetime.now()
            # 并行处理

            print(f"Starting parallel processing with {args.para_file} workers...")
            with multiprocessing.Pool(processes=args.para_file) as pool:
                results = pool.map(process_file_export, tasks)
            
            # 检查结果
            success_count = sum(1 for r in results if not r.empty)
            print(f"\nAlignment completed: {success_count}/{len(tasks)} files succeeded")
            # 从每个结果中提取所需列并合并
            combined_df = pd.concat(
                [df[['accession', 'organism']] for df in results],
                ignore_index=True
            )

            # 可选：删除重复行（如果需要）
            combined_df = combined_df.drop_duplicates()

            # 检查name列是否有重复值
            has_duplicates = combined_df['accession'].duplicated().any()
            if has_duplicates:
                print("accession have duplicates")
                # 查看重复值的具体情况
                duplicates = combined_df[combined_df['accession'].duplicated(keep=False)]
                raise ValueError("repeat：\n", duplicates)
            
            # 检查name列是否有重复值
            has_duplicates = combined_df['organism'].duplicated().any()
            if has_duplicates:
                print("organism have duplicates")
                # 查看重复值的具体情况
                duplicates = combined_df[combined_df['organism'].duplicated(keep=False)]
                raise ValueError("repeat：\n", duplicates)

            if success_count < len(tasks):
                print("Warning: Some files failed alignment. Check error messages above.")

            output_csv = Path(args.output_dir)/"seq_info.csv"
            combined_df.to_csv(output_csv, index=False)
            print(f"sequence information have save in {output_csv}")


        if choice == 2:
            header_format = input("Please enter the format of record.id, for example:\n\t \"accession:location|organism|description\"\n your input:")
            organism_info = input("Please enter organism information file(.csv): ")
            organism_df = pd.read_csv(organism_info)

            # 准备任务列表
            tasks = []
            for in_file_name in file_name_list:
                if in_file_name.split(".")[-1].lower() in ["fasta", "fas", "fa"]:
                    in_file_path = os.path.join(args.input, in_file_name)
                    out_file_path = os.path.join(args.output_dir, Path(in_file_name).name)
                    tasks.append((in_file_path, header_format, out_file_path, organism_df))


            # 并行处理
            print(f"Starting parallel processing with {args.para_file} workers...")
            with multiprocessing.Pool(processes=args.para_file) as pool:
                results = pool.map(process_file_correct, tasks)
            # 检查结果
            success_count = sum(1 for r in results if r)
            print(f"\nAlignment completed: {success_count}/{len(tasks)} files succeeded")
            if success_count < len(tasks):
                print("Warning: Some files failed alignment. Check error messages above.")


    elif os.path.isfile(args.input):
        if not check_sequence_lengths(args.input):
            raise FileExistsError(f"The sequences in the FASTA file are not aligned: {args.input}")
        print("Select FASTA operation:  \n1 - Display sequence records  \n2 - Fix organism annotations \n0 - stop execution ")
        choice = get_valid_input()

        if choice == 0:
            print("Execution stopped.")
            sys.exit(1)
        
        if choice == 1:
            header_format = input("Please enter the format of record.id, for example: \"accession:location|organism|description\"")
            result_df = export_fasta_info(args.input, header_format)
            result_df.to_csv(Path(args.output_dir)/"seq_info.csv", index=False)
            print(f"sequence information have save in {args.output_dir}seq_info.csv")

        if choice == 2:
            header_format = input("Please enter the format of record.id, for example: \"accession:location|organism|description\"")
            organism_info = input("Please enter organism information file(.csv): ")
            organism_df = pd.read_csv(organism_info)
            process_file_correct(args.input, header_format, args.output_dir, organism_df )

    else:
        print(f"'{args.input}' does not exist or is neither a file nor directory.")



