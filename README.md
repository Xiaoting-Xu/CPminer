# CPminer
CPminer is a Python-based software designed for efficient downloading of chloroplast genomes. Its workflow includes automated and batch processing steps such as chloroplast genome downloading, CDS extraction and filtering, sequence alignment, and trimming Each step allows for user intervention, enabling efficient construction of chloroplast genome phylogenetic datasets. CPminer, developed using Python (v3.9.5), primarily leverages Biopython (version 1.79) for handling nucleotide sequence objects. The Entrez.search function from the Biopython library is encapsulated for automated batch sequence downloading and supports resumption from interruptions. The SeqIO function is employed for rapid parsing of chloroplast genome GenBank files. CPminer also integrates widely used phylogenetic analysis software, including the sequence alignment tool MAFFT and the alignment trimming tool trimAl. To address the challenges of aligning large sets of sequences simultaneously, MAFFT’s --add parameter is used, allowing profile alignment of sub-datasets classified by different taxa to enhance the quality of the sequence alignment matrix. The alignment trimming module also includes parameters tailored for chloroplast gene alignment matrices.
