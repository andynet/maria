from Bio import SeqIO
import numpy as np
from pydivsufsort import divsufsort

def parese_msa(filename):
    pos_to_col = {}
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    for seq_id, fasta in enumerate(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        col = 0
        for i in range(len(sequence)):
            if sequence[i] != '-':
                pos_to_col[ (seq_id, col) ] = i
                col += 1
    return pos_to_col

def parese_msa_bis(msa):
    pos_to_col={}
    for i in range(len(msa)):
        pos=0
        msa[i]+='$'
        for j in range(len(msa[i])):
            if msa[i][j] != '-':
                pos_to_col[(i,pos)]=j
                pos+=1
    return pos_to_col

def findP(EP,p):
    for i in range(len(EP)):
        if EP[i] > p:
            return p - EP[i-1]

def finfSeqn(EP,p):
    for u in range(len(EP)):
        if EP[u]>p:
            return u-1

def constructo_col(filename):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    msa=[]
    for i in fasta_sequences:
        msa.append(str(i.seq))

    T='$'.join(msa).replace('-','')+'$#'
    #sa=genereta_sa(T)
    sa=list(divsufsort(T))
    EP=[0]
    for i in range(len(T)):
        if T[i] == '$':
            EP.append(i+1)


    P=[]
    seqn=[]
    bwt=''
    for i in range(len(T)):
        p=sa[i]-1
        bwt+=T[p]
        P.append(findP(EP,p))
        seqn.append(finfSeqn(EP,p))

    C=[]
    psot_clo=parese_msa_bis(msa)
    for i in range(len(T)):
        if seqn[i]==-1:
            C.append(len(msa)+1)
            continue
        C.append(psot_clo[(seqn[i],P[i])])

    #print(C)
    return (C, seqn, bwt, sa)

def rl_encode(C: [int], R: [int]) -> ([int], [int]):
    # TODO: Alessia
    pass

def rindex_query(T: str, P: str) -> int:
    n = len(T)
    m = len(P)
    for i in range(n-m):
        j = 0
        while j < m and T[i+j] == P[j]:
            j += 1
        if j == m:
            return i
    return -1

def lce(T: str, i: int, j: int) -> (int, bool):
    n = len(T)
    k = 0
    while i+k < n and j+k < n and T[i+k] == T[j+k]:
        k += 1
    if i+k >= n:        return (k, True)    # suffix Ti ends sooner
    if j+k >= n:        return (k, False)   # suffix Tj ends sooner
    if T[i+k] < T[j+k]: return (k, True)    # suffix Ti < suffix Tj
    if T[j+k] < T[i+k]: return (k, False)   # suffix Tj < suffix Ti

    print("I should never end here.")

def binsearch(occ1: int, rle_C: ([int], [int]), T: str, upper: bool):
    # TODO: Goobi
    pass

def doc_listing(C: [int], i: int, j: int) -> [int]:
    res = [C[i]]
    for k in range(i+1,j):
        if C[k-1] != C[k]: res.append(C[k])
    return res

def get_cols(occ1: int, rle_C: ([int], [int]), T: str) -> [int]:
    upper = binsearch(occ1, rle_C, T, True)     # upper=True
    lower = binsearch(occ1, rle_C, T, False)    # upper=False
    C, R = rle_C
    return doc_listing(C, upper, lower)

def search(T: str, P: str) -> list[int]:
    C, R, bwt, sa = constructo_col('mul_aln_meningitidis.5.fa') # TODO: return sa
    rle_C = rl_encode(C, R)

    occ1 = rindex_query(T, P)
    cols = get_cols(occ1, rle_C, T)

if __name__ == '__main__':
    T = "ABBAABBABABABABBABA"
    P = "BABA"
    pos = rindex_query(T, P)
    print(pos)
    i, j = 3, 11
    print(T[i:])
    print(T[j:])
    (lenght, sign) = lce(T, i, j)
    print(lenght, sign)

    i, j = 3, 15
    l = [4,4,5,5,5,5,7,8,8,4,4,8,8,8,8,8,8,1,1]
    cols = doc_listing(l, i, j)
    print(l[i:j])
    print(cols)
