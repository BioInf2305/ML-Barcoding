import sys
import tempfile
import scipy.io
from collections import OrderedDict
from strkernel.gappy_kernel import gappypair_kernel as gk
from strkernel.gappy_trie import gappypair_kernel as gt
from Bio import SeqIO
from Bio.Seq import Seq

def MakeKmerGappy(file1,K,T,G):
    ti= tempfile.NamedTemporaryFile(delete=False,mode='w+b')
    dest=open(file1+"_gappy_kmers.csv","w")
    sequences = list(SeqIO.parse(file1, "fasta"))
    X1=gk(sequences,k=int(K),t=int(T),g=int(G))
    scipy.io.mmwrite(ti,X1)
    ti.seek(0)
    cnt=0
    KmerDict=OrderedDict()
    FirstRow=[]
    NameDict=OrderedDict()
    NameDict=ExtractInfoFasta(file1)
    for line in ti:
        line=line.decode()
        line=line.rstrip()
        if cnt<3:
            cnt+=1
        else:
            a=line.split()
            if a[1] not in FirstRow:
                FirstRow.append(a[1])
            if a[0] not in KmerDict:
                KmerDict[a[0]]={}
            KmerDict[a[0]][a[1]]=str(float(a[2]))
    dest.write("id,"+",".join(FirstRow)+",species")
    dest.write("\n")
    Ke=list(KmerDict.keys())
    for i in Ke:
        Li=[]
        Li.append(NameDict[int(i)][0])
        for j in FirstRow:
            if KmerDict[i].get(j,"-9")=="-9":
                Li.append("0")
            else:
                Li.append(str(int(float(KmerDict[i][j]))))
        Li.append(NameDict[int(i)][1])
        dest.write(",".join(Li))
        dest.write("\n")
    dest.close()
    ti.close()

def ExtractInfoFasta(file1):
    TmpDict=OrderedDict()
    cnt=0
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
MakeKmerGappy(inp,inp1,inp2,inp3)
