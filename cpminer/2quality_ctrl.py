# -*- coding = utf-8 -*-

import argparse
from pathlib import Path
from Bio import SeqIO
import os
from datetime import datetime
import pandas as pd
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

def process_single_file(file_path, out_folder_path, cds_threshold, ambig_threshold):
    """处理单个GB文件的辅助函数"""
    print(f"Collect information of : {Path(file_path).name}")
    results = {
        'filename': Path(file_path).stem,
        'organism': None,
        'unclear_bases': None,
        'unclear_ratio': None,
        'CDS': 0,
        'tRNA': 0,
        'rRNA': 0,
        'gene_count': 0,
        'sequence_length': None,
        'error': None
    }
    is_problematic = False  # 标记文件是否不符合阈值

    try:
        # 解析文件并获取第一个记录
        record = next(SeqIO.parse(file_path, 'genbank'))
        
        # 序列统计
        gb_seq = str(record.seq).upper()
        results['sequence_length'] = len(gb_seq)
        un_atcg = sum(1 for c in gb_seq if c not in {'A', 'T', 'C', 'G'})
        results['organism'] = record.annotations["organism"]
        results['unclear_bases'] = un_atcg
        results['unclear_ratio'] = un_atcg / len(gb_seq) if len(gb_seq) > 0 else 0
        
        # 特征统计
        feature_counts = {'CDS': 0, 'tRNA': 0, 'rRNA': 0}
        for feature in record.features:
            if feature.type == 'gene':
                results['gene_count'] += 1
            elif feature.type in feature_counts:
                feature_counts[feature.type] += 1
        
        # 更新特征计数
        results.update(feature_counts)
        
        # 检查阈值
        if results['CDS'] < cds_threshold or results['unclear_ratio'] > ambig_threshold:
            is_problematic = True
            try:
                os.remove(file_path)  # 删除原GB文件
                # 保存为FASTA
                out_file = Path(out_folder_path) / (Path(file_path).stem + ".fasta")
                with open(out_file, "w") as fw:
                    fw.write(f">{record.id}|{record.annotations.get('organism', '')}|{record.description}\n")
                    fw.write(str(record.seq) + "\n")
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                results['error'] = str(e)
    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        results['error'] = str(e)
    
    return results, is_problematic

def generate_genome_report(in_folder_path, out_folder_path, cds_threshold=80, ambig_threshold=0.2, threads=3):
    """多线程生成基因组统计报告"""
    gb_files = [
        Path(in_folder_path) / f for f in os.listdir(in_folder_path) 
        if f.lower().endswith((".gb", ".gbf", ".gbk", ".genbank"))
    ]
    print(f"Found {len(gb_files)} GenBank files in {in_folder_path}")
    
    # 确保输出目录存在
    os.makedirs(out_folder_path, exist_ok=True)
    
    all_results = []
    error_count = 0
    
    # 设置线程数
    if threads is None:
        # 默认使用CPU核心数，但不超过文件数量
        threads = min(multiprocessing.cpu_count(), len(gb_files))
    else:
        # 确保线程数不超过文件数量
        threads = min(threads, len(gb_files))
    
    print(f"Using {threads} threads for processing")
    
    # 使用线程池处理文件
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                process_single_file, 
                file_path, 
                out_folder_path,
                cds_threshold,
                ambig_threshold
            ): file_path for file_path in gb_files
        }
        
        for future in as_completed(futures):
            file_path = futures[future]
            try:
                results, is_problematic = future.result()
                all_results.append(results)
                if is_problematic:
                    error_count += 1
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                all_results.append({
                    'filename': Path(file_path).stem,
                    'error': str(e)
                })
    
    # 创建DataFrame并保存
    df = pd.DataFrame(all_results)
    columns = [
        'filename', 'organism', 'sequence_length', 'gene_count', 
        'CDS', 'tRNA', 'rRNA', 
        'unclear_bases', 'unclear_ratio', 'error'
    ]
    df = df.reindex(columns=columns)  # 确保列顺序正确
    df.to_csv(Path(in_folder_path) / "gb_info.csv", index=False)
    
    return error_count

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: For each gb file within the gb_folder_path, CDS counts and ambiguous base proportions are quantified.\nThese metrics subsequently inform filtration criteria, where gb records are evaluated against predefined thresholds \nto determine retention eligibility.
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
    parser.add_argument("-o", "--output_dir", help="Output directory for fasta files (only created for problematic GenBank annotations)", required=True)
    parser.add_argument("-c", "--cds_threshold", type=int, default=80,
                        help="Minimum CDS count threshold (default: %(default)s)", required=False)
    parser.add_argument("-a", "--ambig_threshold", type=float, default=0.2,
                        help="Maximum ambiguous base proportion (default: %(default)s)", required=False)
    parser.add_argument("-p", "--perl_path", help="Execution path for PGA.pl (Perl Genome Annotator)", required=False)
    # 添加线程数参数
    parser.add_argument("-t", "--threads", type=int, default=3,
                        help="Number of threads to use %(default)s", required=False)
    args = parser.parse_args()

    # 打印配置摘要
    print("\n" + "="*50)
    print(f"{'Configuration Summary':^50}")
    print("="*50)
    print(f"Input directory:      {args.input_dir}")
    print(f"Output directory:     {args.output_dir}")
    print(f"CDS threshold:        {args.cds_threshold}")
    print(f"Ambiguity threshold:  {args.ambig_threshold:.1%}")
    print(f"Perl path:            {args.perl_path or 'Not specified'}")
    print(f"Processing threads:   {args.threads}")
    print("="*50 + "\n")

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    cds_threshold = args.cds_threshold
    ambig_threshold = args.ambig_threshold
    perl_path = args.perl_path
    threads = args.threads

    time0 = datetime.now()
    count_error_gb = generate_genome_report(
        input_dir, 
        output_dir, 
        cds_threshold, 
        ambig_threshold, 
        threads  # 传递线程数参数
    )
    if count_error_gb :
        while True:
            reannotation_choice = input(f"Annotations are problematic in {count_error_gb} GenBank (.gb) files. Need re-annotation? (y/n): ")
            # 标准化输入并检查
            normalized_input = reannotation_choice.strip().lower()
            if normalized_input in ('y', 'yes'):
                print("Proceeding with re-annotation...")
                # 执行重新注释的逻辑

                break
            elif normalized_input in ('n', 'no'):
                print("Skipping re-annotation.")
                # 跳过重新注释的逻辑
                break
            else:
                print("Invalid input. Skipping re-annotation.")

    else:
        print("All gb files meet the criteria.")

    time1 = datetime.now()
    print(f"Total running time: {time1 - time0}")

"""
==================================================
              Configuration Summary
==================================================
Input directory:      H:\cpminer2\CPminer命令行代码\testFile
Output directory:     H:\cpminer2\CPminer命令行代码\fasta
CDS threshold:        80
Ambiguity threshold:  20.0%
Perl path:            Not specified
Processing threads:   3
==================================================

Found 939 GenBank files in H:\cpminer2\CPminer命令行代码\testFile
Using 3 threads for processing
Found 939 GenBank files in H:\cpminer2\CPminer命令行代码\testFile
Using 3 threads for processing
Annotations are problematic in 211 GenBank (.gb) files. Need re-annotation? (y/n): n
Skipping re-annotation.
Total running time: 0:00:13.973249
"""