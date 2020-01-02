import sys
from collections import OrderedDict
from strkernel.mismatch_kernel import MismatchKernel
from strkernel.mismatch_kernel import preprocess
from Bio import SeqIO
from Bio.Seq import Seq

def GenerateMatrices(file1,L,K,M):
    dest=open(file1+".mismatch_similarity_matrix.txt",'w')
    dest1=open(file1+".mismatch_kernel.csv",'w')
    data=list(SeqIO.parse(file1,"fasta"))
    data_processed=preprocess(data)
    mismatch_kernell=MismatchKernel(l=int(L),k=int(K),m=int(M)).get_kernel(data_processed)
    similarity_mat1=mismatch_kernell.kernel
    NameDict=OrderedDict()
    NameDict=ExtractInfoFasta(file1)
    for i in range(len(similarity_mat1)):
        for j in range(len(similarity_mat1[i])):
            dest.write(str(similarity_mat1[i][j]))
            dest.write("\t")
        dest.write("\n")
    ke=list(mismatch_kernell.leaf_kmers.keys())
    dest1.write("id,"+",".join(ke)+",species")
    dest1.write("\n")
    k_ke=[]
    for i in range(len(ke)):
        tmp=list(mismatch_kernell.leaf_kmers[ke[i]].keys())
        for i in tmp:
            if i not in k_ke:
                k_ke.append(i)
    k_ke.sort()
    for k_k in k_ke:
        li=[NameDict[k_k][0]]
        for i in ke:
            if mismatch_kernell.leaf_kmers[i].get(k_k,"-9")=="-9":
                li.append("0")
            else:li.append(str(mismatch_kernell.leaf_kmers[i][k_k]))
        li.append(NameDict[k_k][1])
        dest1.write(",".join(li))
        dest1.write("\n")
    dest.close()
    dest1.close()

def ExtractInfoFasta(file1):
    TmpDict=OrderedDict()
    cnt=-1
    with open(file1) as fast:
        for lin in fast:
            if lin.startswith(">"):
                cnt+=1
                a=lin.lstrip(">").split("|")
                TmpDict[cnt]=[a[0],a[1]]
    return TmpDict


inp=sys.argv[1]
inp1=sys.argv[2]
inp2=sys.argv[3]
inp3=sys.argv[4]
GenerateMatrices(inp,inp1,inp2,inp3)
    
