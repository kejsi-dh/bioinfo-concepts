def min(a, b):
    return a if a < b else b

def overlapping_pair(a, b):
    max_overlap = 0
    merged_str = ""
    len1, len2 = len(a), len(b)
    for i in range(1, min(len1, len2) + 1):
        if a[len1 - i:] == b[:i]:
            if i > max_overlap:
                max_overlap = i
                merged_str = a + b[i:]

    for i in range(1, min(len1, len2) + 1):
        if a[:i] == b[len2 - i:]:
            if i > max_overlap:
                max_overlap = i
                merged_str = b + a[i:]

    return max_overlap, merged_str

def shortest_superstr(strs):
    while len(strs) > 1:
        max_overlap = 0
        l, r = 0, 0
        merged_str = ""

        for i in range(len(strs)):
            for j in range(i + 1, len(strs)):
                overlap, merged = overlapping_pair(strs[i], strs[j])

                if overlap > max_overlap:
                    max_overlap = overlap
                    merged_str = merged
                    l, r = i, j

        if max_overlap == 0:
            strs[0] += strs.pop()
        else:
            strs[l] = merged_str
            strs.pop(r)

    return strs[0]

strings = ["TACGTA", "ACGTAC", "TACGAT", "CTGACG", "GTACGT"]
print("The Shortest Superstring is:", shortest_superstr(strings))
