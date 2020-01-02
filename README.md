ML-Barcoding
Depository contains scripts used for sequence identification using Machine learnign approach as implemented in WEKA. In brief, these
python scripts generated the input file for WEKA machine learning tool. 

1). Mismatch kernal
  usage:python3 kmerMismatchV4.py <input_file> <value of L> <value of K> <value of M>
  This will generate two files: (i) Similarity matrix text file (ii) .csv file with leaf kmers. The .csv file can be used directly in WEKA machine learning software to test data against various available classifiers.

2). Gappy kernel
  usage:python3 kmerGappyV4.py <input_file> <value of L> <value of K> <value of M>
  This will generate .csv file of sparse matrix, which contain kmer frequencies. This .csv file can also be used with WEKA GUI interface. 
  
  
