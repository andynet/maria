# from unittest import TestCase

from create_col import *

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
