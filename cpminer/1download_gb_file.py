# *-* coding:utf-8 *-*
# @Time:2022/1/5 10:02
# @Author:Ruijing Cheng
# @File:download_gb_file.py
# @Software:PyCharm

import argparse
from pathlib import Path
import os
import threading
import time
from datetime import datetime
from Bio import Entrez, SeqIO
import func_timeout
from func_timeout import func_set_timeout
import queue

@func_set_timeout(600)
def my_efetch(accessions):
    """带超时设置的API请求函数，支持批量下载"""
    if isinstance(accessions, list):
        accessions = ",".join(accessions)
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=accessions)
    return handle

def download_batch(batch_accessions, out_path, request_lock, last_request_time):
    """下载一批基因组文件"""
    max_retries = 3
    retry_delay = 5
    
    for attempt in range(max_retries):
        try:
            # 速率控制
            with request_lock:
                current_time = time.time()
                elapsed = current_time - last_request_time[0]
                if elapsed < 0.34:  # 限制每秒3次请求
                    time.sleep(0.34 - elapsed)
                last_request_time[0] = time.time()
            
            # 执行批量下载
            t0 = datetime.now()
            with my_efetch(batch_accessions) as handle:
                records = list(SeqIO.parse(handle, "gb"))
            
            # 保存所有记录
            for record in records:
                accession = record.id.split('.')[0]  # 使用基础accession作为文件名
                out_file = Path(out_path) / f"{accession}.gb"
                with open(out_file, "w") as f:
                    SeqIO.write(record, f, "gb")
            
            t1 = datetime.now()
            print(f"Batch [{', '.join(batch_accessions)}] downloaded in {t1 - t0}")
            return True
        
        except func_timeout.exceptions.FunctionTimedOut:
            print(f"Batch timeout ({attempt+1}/{max_retries}) for {batch_accessions}")
        except Exception as e:
            print(f"Batch error ({attempt+1}/{max_retries}) for {batch_accessions}: {str(e)}")
        
        time.sleep(retry_delay * (attempt + 1))
    
    print(f"Failed to download batch after {max_retries} attempts: {batch_accessions}")
    return False

def download_worker(out_path, request_lock, last_request_time, task_queue, batch_size=3):
    """工作线程函数，处理下载任务"""
    while True:
        try:
            accessions = task_queue.get(timeout=10)  # 10秒超时
            if not accessions:
                break
                
            download_batch(accessions, out_path, request_lock, last_request_time)
        except queue.Empty:
            break
        finally:
            task_queue.task_done()

def download_gb_file(email, in_path, out_path, max_threads=10, batch_size=3):
    """多线程批量下载基因组文件"""
    Entrez.email = email
    
    # 读取accession列表
    with open(in_path, "r") as fr:
        accession_list = fr.read().splitlines()
    
    # 创建输出目录
    os.makedirs(out_path, exist_ok=True)
    
    # 过滤已存在文件
    existing_files = {p.stem for p in Path(out_path).glob("*.gb")}
    accession_list = list(set(accession_list) - existing_files)
    
    print(f"Starting download of {len(accession_list)} files using {max_threads} threads (batch size: {batch_size})")
    
    # 创建任务队列
    task_queue = queue.Queue()
    
    # 将accession列表分成批次
    for i in range(0, len(accession_list), batch_size):
        batch = accession_list[i:i+batch_size]
        task_queue.put(batch)
    
    # 共享变量用于请求速率控制
    last_request_time = [time.time()]
    request_lock = threading.Lock()
    
    # 创建工作线程
    threads = []
    for _ in range(max_threads):
        t = threading.Thread(
            target=download_worker,
            args=(out_path, request_lock, last_request_time, task_queue, batch_size)
        )
        t.daemon = True
        t.start()
        threads.append(t)
    
    # 等待所有任务完成
    task_queue.join()
    
    # 通知线程退出
    for _ in range(max_threads):
        task_queue.put(None)
    
    # 等待所有线程结束
    for t in threads:
        t.join()
    
    print("All downloads completed")

class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Function: Download chloroplast genomes from GenBank using batch requests
        另外，我们设置最大重试次数为5，如果5次都失败，则记录失败并放弃。
        Author: Ruijing Cheng
        Organization: 
        GitHub Site:
        Usage Example:
          添加示例
        """
        return f"{description}\n\n{help_text}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-e", "--your_email", help="Your email address", required=False)
    parser.add_argument("-i", "--input_file", 
                        help="Input path of accession list file. All accessions must be from the nucleotide database", 
                        required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory for GenBank files", required=True)
    parser.add_argument("-t", "--threads", type=int, default=10, 
                        help="Number of download threads (default: 10)")
    parser.add_argument("-b", "--batch_size", type=int, default=3, 
                        help="Number of accessions per batch request (default: 3)")
    
    args = parser.parse_args()

    email = args.your_email or "your_email@example.com"
    in_path = Path(args.input_file)
    out_path = Path(args.output_dir)
    max_threads = args.threads
    batch_size = args.batch_size

    print(f"Input file path: {in_path}")
    print(f"Output directory: {out_path}")
    print(f"Using {max_threads} threads with batch size {batch_size}")

    print("Start downloading...")
    time0 = datetime.now()
    download_gb_file(email, in_path, out_path, max_threads, batch_size)
    time1 = datetime.now()
    print(f"Total running time: {time1 - time0}")