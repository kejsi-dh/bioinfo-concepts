import numpy as np

def submatrix(a, b, v, u, match, mismatch):
    # sequence lengths
    submatrix.I = len(a)
    submatrix.J = len(b)

    # size of sub matrix
    r = submatrix.J + 1
    c = submatrix.I + 1

    S = np.zeros((r,c), dtype=int)
    submatrix.E = np.zeros((r,c), dtype=int)
    submatrix.F = np.zeros((r,c), dtype=int)

    # initialization
    S[0,0] = 0
    for i in range (0,r):
      submatrix.E[i,0] = v+u * i
      submatrix.F[i,0] = -10**6

    for j in range (0,c):
      submatrix.F[0,j] = v+u * j
      submatrix.E[0,j] = -10**6

    for i in range(1,r):
      for j in range(1,c):
        submatrix.E[i,j] = max(submatrix.E[i-1,j] + u, S[i-1,j] + v + u)
        submatrix.F[i,j] = max(submatrix.F[i,j-1] + u, S[i,j-1] + v + u)
        if(a[j-1] == b[i-1]):
          S[i,j] = max(S[i-1, j-1] + match, submatrix.E[i,j], submatrix.F[i,j])
        else: S[i,j] = max(S[i-1, j-1] + mismatch, submatrix.E[i,j], submatrix.F[i,j])

    return S


def traceback(a, b, S):
    # alignments
    a_ = []
    b_ = []

    i = submatrix.J
    j = submatrix.I
    count = 0

    while (i > 0 or j > 0) and (count <= 10**6):
        count = count + 1
        if (i > 0 and j > 0 and S[i,j] == S[i-1, j-1] + match and b[i-1] == a[j-1]) or (
                i > 0 and j > 0 and S[i,j] == S[i-1, j-1] + mismatch and b[i-1] != a[j-1]):
          a_.append(a[j-1])
          b_.append(b[i-1])
          i = i-1
          j = j-1

        elif i > 0 and S[i,j] == submatrix.E[i,j]:
          a_.append('-')
          b_.append(b[i-1])
          i = i-1

        elif j > 0 and S[i,j] == submatrix.F[i,j]:
          a_.append(a[j-1])
          b_.append('-')
          j = j-1

    a_.reverse()
    b_.reverse()
    return a_, b_


v = -0.5 # open gap penalty
u = -0.5 # extended gap penalty
match = 1.0
mismatch = -1.0

a = 'AGTAGTGTATAGTAC'
b = 'ATGTAGTACATG'

S = submatrix(a, b, v, u, match, mismatch)
print(S)

a_, b_ = traceback(a, b, S)
print('First sequence after traceback: ')
print(a_)
print('Second sequence after traceback: ')
print(b_)
