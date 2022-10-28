import sys
from Bio import SeqIO
import numpy as np
from pydivsufsort import divsufsort

def parse_msa(msa):
    pos_to_col={}
    for i in range(len(msa)):
        pos=0
        msa[i]+='$'
        for j in range(len(msa[i])):
            if msa[i][j] != '-':
                pos_to_col[(i,pos)]=j
                pos+=1
    return pos_to_col

def map_msa_to_T(msa):
    msa_to_T = {}
    pos = 0
    for i in range(len(msa)):
        col = 0
        msa[i] += '$'
        for j in range(len(msa[i])):
            if msa[i][j] != '-':
                msa_to_T[(i,col)] = pos
                col += 1
                pos += 1
    return msa_to_T

def find_P(EP, p):
    for i in range(len(EP)):
        if EP[i] > p:
            return p - EP[i-1]

def find_row(EP, p):
    for u in range(len(EP)):
        if EP[u]>p:
            return u-1

def load_msa(filename):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    msa=[]
    for i in fasta_sequences:
        msa.append(str(i.seq))
    return msa

def construct_col(msa):
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
        P.append(find_P(EP,p))
        seqn.append(find_row(EP,p))

    C=[]
    pos_to_col=parse_msa(msa)
    for i in range(len(T)):
        if seqn[i]==-1:
            C.append(len(msa)+1)
            continue
        C.append(pos_to_col[(seqn[i],P[i])])

    return (C, seqn, bwt, sa)

def rl_encode(C: [int], R: [int]) -> ([int], [int]):
    new_C = [C[0]]
    new_R = [R[0]]

    for i in range(1, len(C)):
        if C[i] != new_C[-1]:
            new_C.append(C[i-1])
            new_R.append(R[i-1])
            new_C.append(C[i])
            new_R.append(R[i])

    new_C.append(C[-1])
    new_R.append(R[-1])

    return (new_C, new_R)

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


def check(occ1: int, len1: int, rle_C: ([int], [int]), T: str, msa2t: dict[int, int], upper: bool, middle: int, offset: int) -> bool:
    pos_act, pos_off = msa2t[ (rle_C[1][middle], rle_C[0][middle]) ], msa2t[ (rle_C[1][middle+offset], rle_C[0][middle+offset]) ]
    offset_lce = lce(T, pos_off, occ1)[0]
    middle_lce = lce(T, pos_act, occ1)[0]
    offset_sign = lce(T, pos_off, occ1)[1]
    middle_sign = lce(T, pos_act, occ1)[1]
    if middle_lce >= len1 and offset_lce < len1:
        return "done"
    if upper: # might wanna convert this to a switch case
        if middle_lce == offset_lce and middle_sign:
            return "down"
        elif (middle_lce >= len1) or (not middle_sign):
            return "up"
        else:
            print("I should not end up here", file=sys.stderr)
            return "up"
    else:
        if middle_lce == offset_lce and (not middle_sign):
            return "up"
        elif (middle_lce >= len1) or middle_sign:
            return "down"
        else:
            print("I should not end up here", file=sys.stderr)
            return "up"
    print("I should not end up here", file=sys.stderr)
    return "up"

def binsearch(occ1: int, len1: int, rle_C: ([int], [int]), T: str, msa2t: dict[int, int], upper: bool) -> int:
    # we represent the half-open interval <start, end)
    N = len(rle_C[0])
    start, end, offset  = 1, N, -1
    if not upper:
        start, end, offset = 0, N-1, 1
    middle = None
    while start < end:
        middle = (end - start)//2
        res = check(occ1, len1, rle_C, T, msa2t, upper, middle, offset)
        if res == "done":
            break
        elif res == "up":
            end = middle
        else:
            start = middle + 1
    return middle

def doc_listing(C: [int], i: int, j: int) -> [int]:
    res = [C[i]]
    for k in range(i+1,j):
        if C[k-1] != C[k]: res.append(C[k])
    return res

def get_cols(occ1: int, len1: int, rle_C: ([int], [int]), T: str, msa: [[chr]]) -> [int]:
    msa2t = map_msa_to_T(msa)
    upper = binsearch(occ1, len1, rle_C, T, msa2t, True)     # upper=True
    lower = binsearch(occ1, len1, rle_C, T, msa2t, False)    # upper=False
    C, R = rle_C
    return doc_listing(C, upper, lower)

def search(T: str, P: str) -> list[int]:
    pass

if __name__ == '__main__':
    filename = 'data/test.fa'
    msa = load_msa(filename)

    T='$'.join(msa).replace('-','')+'$#'
    P='AT'

    C, R, bwt, sa = construct_col(msa)
    rle_C = rl_encode(C, R)

    occ1 = rindex_query(T, P)
    cols = get_cols(occ1, len(P), rle_C, T, msa)

if __name__ == '__main__':
    pass
