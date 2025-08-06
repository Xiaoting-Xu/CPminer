import argparse
import os
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import glob
from pathlib import Path

def concatenate_alignments(fasta_files, output_phy, partition_file):
    """
    将多个比对文件合并为超矩阵并生成分区文件
    """
    # 收集所有物种和基因数据
    all_species = set()
    gene_lengths = {}
    valid_files = []
    
    print("Processing FASTA files:")
    # 第一遍：获取所有物种和每个基因的长度
    for file in fasta_files:
        try:
            # 确保文件存在且可读
            if not os.path.isfile(file):
                print(f"  Warning: '{file}' is not a valid file, skipping")
                continue
                
            print(f"  Reading: {os.path.basename(file)}")
            alignment = AlignIO.read(file, "fasta")
            gene_name = Path(file).stem
            gene_lengths[gene_name] = alignment.get_alignment_length()
            
            for record in alignment:
                all_species.add(str(record.id).split('|')[-1])
            
            valid_files.append(file)
        except Exception as e:
            print(f"  Error processing {file}: {str(e)}")
            continue
    
    if not all_species:
        print("Error: No valid species found in any alignment files")
        sys.exit(1)
    
    print("\nGene lengths:")
    for gene, length in gene_lengths.items():
        print(f"  {gene}: {length} bp")
    
    total_length = sum(gene_lengths.values())
    print(f"\nTotal concatenated length: {total_length} bp")
    print(f"Number of species: {len(all_species)}")
    
    # 创建物种列表并按字母排序
    species_list = sorted(list(all_species))
    
    # 第二遍：为每个物种创建拼接序列
    concatenated_seqs = {species: ["-" * total_length] for species in species_list}
    current_position = 0
    partitions = []
    
    print("\nConcatenating alignments:")
    for file in valid_files:
        gene_name = Path(file).stem
        try:
            print(f"  Processing: {gene_name}")
            alignment = AlignIO.read(file, "fasta")
            aln_dict = {record.id: str(record.seq) for record in alignment}
            
            gene_length = gene_lengths[gene_name]
            start = current_position + 1
            end = current_position + gene_length
            partitions.append(f"DNA, {gene_name} = {start}-{end}")
            
            for species in species_list:
                seq = concatenated_seqs[species][0]
                if species in aln_dict:
                    gene_seq = aln_dict[species]
                else:
                    gene_seq = "-" * gene_length
                    
                # 构建新的拼接序列
                new_seq = seq[:current_position] + gene_seq + seq[current_position+gene_length:]
                concatenated_seqs[species] = [new_seq]
            
            current_position += gene_length
        except Exception as e:
            print(f"  Error processing {file}: {str(e)}")
            continue
    
    # 创建多序列比对对象
    records = []
    for species in species_list:
        records.append(SeqRecord(seq=concatenated_seqs[species][0], id=species, description=""))
    
    super_alignment = MultipleSeqAlignment(records)
    
    # 写入PHYLIP文件
    try:
        with open(output_phy, "w") as phy_handle:
            phy_handle.write(f" {len(species_list)} {current_position}\n")
            for record in super_alignment:
                # PHYLIP格式：最多10字符的ID，然后序列
                phy_id = record.id[:10].ljust(10)
                phy_handle.write(f"{phy_id}{record.seq}\n")
        print(f"\nSupermatrix saved to: {output_phy}")
    except Exception as e:
        print(f"\nError writing PHYLIP file: {str(e)}")
        sys.exit(1)
    
    # 写入分区文件
    try:
        with open(partition_file, "w") as part_handle:
            for line in partitions:
                part_handle.write(line + "\n")
        print(f"Partition file saved to: {partition_file}")
    except Exception as e:
        print(f"Error writing partition file: {str(e)}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Concatenate multiple sequence alignments into a supermatrix',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-w', '--work_dir', required=True, 
                        help='Working directory path')
    parser.add_argument('-o', '--output_phy', default="supermatrix.phy", 
                        help='Output PHYLIP file name')
    parser.add_argument('-pf', '--partition_file', default="partitionFile.txt", 
                        help='Output partition file name')
    parser.add_argument('-p', '--patterns', nargs='+', default=["*.fasta"], 
                        help='File patterns to match (e.g., "*.fasta")')
    
    args = parser.parse_args()
    
    # 设置工作目录
    try:
        os.chdir(args.work_dir)
        print(f"Working directory set to: {os.getcwd()}")
    except Exception as e:
        print(f"Error setting working directory: {str(e)}")
        sys.exit(1)
    
    # 获取所有匹配的文件
    fasta_files = []
    for pattern in args.patterns:
        try:
            matched_files = glob.glob(pattern)
            if not matched_files:
                print(f"Warning: No files found for pattern '{pattern}'")
            
            for file in matched_files:
                if os.path.isfile(file):
                    fasta_files.append(file)
                else:
                    print(f"Warning: '{file}' is not a file, skipping")
        except Exception as e:
            print(f"Error with pattern '{pattern}': {str(e)}")
    
    if not fasta_files:
        print("Error: No valid FASTA files found")
        sys.exit(1)
    
    print(f"\nFound {len(fasta_files)} alignment files to concatenate")
    
    # 执行拼接
    concatenate_alignments(
        fasta_files=fasta_files,
        output_phy=args.output_phy,
        partition_file=args.partition_file
    )

if __name__ == "__main__":
    main()