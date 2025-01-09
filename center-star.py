import pandas as pd
import numpy as np

def distancematrix(seq_list, penalty):
    seqnr = len(seq_list)
    distance_m = np.zeros((seqnr, seqnr))
    for i in range(seqnr):
        for j in range(seqnr):
            lenx = len(seq_list[j])
            leny = len(seq_list[i])
            A = np.zeros((lenx + 1, leny + 1))
            A[0, :] = np.arange(0, penalty * (leny + 1), step=penalty)
            A[:, 0] = np.arange(0, penalty * (lenx + 1), step=penalty)
            for J in range(1, leny + 1):
                for I in range(1, lenx + 1):
                    seq1 = seq_list[i]
                    seq2 = seq_list[j]
                    if seq2[I - 1] == seq1[J - 1]: sub_penalty = 0
                    else: sub_penalty = 1
                    subPenalty = A[I-1, J-1] + sub_penalty
                    A[I,J] = min(subPenalty, A[I-1, J] + penalty, A[I, J-1] + penalty)

            distance_m[i,j] = A[-1, -1]

    dm_pandas = pd.DataFrame(distance_m, index=seq_list, columns=seq_list)
    return distance_m, dm_pandas

def center_string(distance_m):
    colnr = distance_m.shape[1]
    rowsum = [sum(distance_m[:,i]) for i in range(colnr)]
    minsum = min(rowsum)
    pos_minsum = rowsum.index(minsum)
    return pos_minsum

def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1):
    n = len(a) + 1
    m = len(b) + 1
    score_matrix = np.zeros((n,m), dtype=int)
    traceback_matrix = np.zeros((n,m), dtype=int)

    for i in range(n):
        score_matrix[i][0] = gap * i
    for j in range(m):
        score_matrix[0][j] = gap * j

    for i in range(1,n):
        for j in range(1,m):
            match_score = score_matrix[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch)
            del_score = score_matrix[i-1][j] + gap
            in_score = score_matrix[i][j-1] + gap

            score_matrix[i][j] = max(match_score, del_score, in_score)

            if score_matrix[i][j] == match_score:
                traceback_matrix[i][j] = 1  # diagonal: match/mismatch
            elif score_matrix[i][j] == del_score:
                traceback_matrix[i][j] = 2  # up: deletion
            else:
                traceback_matrix[i][j] = 3  # left: insertion

    # traceback
    aligned_a = ""
    aligned_b = ""
    i, j = len(a), len(b)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 1:
            aligned_a = a[i-1] + aligned_a
            aligned_b = b[j-1] + aligned_b
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or traceback_matrix[i][j] == 2):
            aligned_a = a[i-1] + aligned_a
            aligned_b = "-" + aligned_b
            i -= 1
        elif j > 0 and (i == 0 or traceback_matrix[i][j] == 3):
            aligned_a = "-" + aligned_a
            aligned_b = b[j-1] + aligned_b
            j -= 1

    return aligned_a, aligned_b

def msa(center_str, seq_list):
    seq_list.remove(center_str)
    first_al = needleman_wunsch(center_str, seq_list[0])
    msa = [first_al[0], first_al[1]]
    seq_list.pop(0)

    for seq in seq_list:
        aligned_pair_center, seq_aligned = needleman_wunsch(center_str, seq)
        aligned_msa_center = msa[0]
        pairwise_blanks = 0
        msa_blanks = 0
        for char in aligned_pair_center[::-1]:
            if char == '-':
                pairwise_blanks = pairwise_blanks + 1
            else:
                break
        for char in aligned_msa_center[::-1]:
            if char == '-':
                msa_blanks = msa_blanks + 1
            else:
                break

        msa_ind = 0 # msa indices
        pair_ind = 0 # pairwise indices
        blankcol_pair = []
        blankcol_msa = []
        while msa_ind < len(aligned_msa_center) and pair_ind < len(aligned_pair_center):
            if aligned_pair_center[pair_ind] == '-' and aligned_msa_center[msa_ind] == '-':
                msa_ind += 1
                pair_ind += 1
            elif aligned_pair_center[pair_ind] == aligned_msa_center[msa_ind]:
                msa_ind += 1
                pair_ind += 1
            elif aligned_pair_center[pair_ind] == '-' and aligned_msa_center[msa_ind] != '-':
                pair_ind += 1
                blankcol_msa.append(msa_ind)
            elif aligned_pair_center[pair_ind] != '-' and aligned_msa_center[msa_ind] == '-':
                msa_ind += 1
                blankcol_pair.append(pair_ind)

        for i, el in enumerate(msa):
            updated_row = list(el)
            offset = 0
            for index in blankcol_msa:
                updated_row.insert(index + offset, '-')
                offset += 1
            if pairwise_blanks > msa_blanks:
                for i in range(pairwise_blanks - msa_blanks):
                    updated_row.append('-')
            msa[i] = ''.join(updated_row)

        updated_seq = list(seq_aligned)
        offset = 0
        for index in blankcol_pair:
            updated_seq.insert(index + offset, '-')
            offset += 1
        if pairwise_blanks < msa_blanks:
            for i in range(msa_blanks - pairwise_blanks):
                updated_seq.append('-')
        msa.append(''.join(updated_seq))

    return '\n'.join(msa)

seq_list = ["CGAGGT", "GTAG", "TTGGCCA", "GTCGGA", "CGAAGTT"]
penalty = 1

# SEQUENCE DISTANCE MATRIX
dm_pandas, distancem = distancematrix(seq_list, penalty)
print("1. Sequence distance matrix:")
print(distancem)

# CENTER STRING
centerstr = center_string(dm_pandas)
print("\n")
str2 = "nd sequence: "
str3 = seq_list[centerstr]
print("2. Center string -- %d"%(centerstr + 1) + str2 + str3)

# MSA WITH CENTER STRING
print("\n")
center_str = seq_list[centerstr]
msa = msa(center_str, seq_list)
print("3. (Near-optimal) MSA:")
print(msa)
