import numpy as np

def levenshtein_distance(a: str, b: str) -> int:
    if not a:
        return len(b)
    if not b:
        return len(a)
    if a[0] == b[0]:
        return levenshtein_distance(a[1:], b[1:])

    insert = levenshtein_distance(a, b[1:])
    delete = levenshtein_distance(a[1:], b)
    replace = levenshtein_distance(a[1:], b[1:])

    return 1 + min(insert, delete, replace)

def dist_matrix(str_list: list[str]) -> np.ndarray:
    size = len(str_list)
    dist_matrix = np.zeros((size, size), dtype=int)

    for i, str1 in enumerate(str_list):
        for j, str2 in enumerate(str_list):
            if i <= j:
                dist = levenshtein_distance(str1, str2)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

    return dist_matrix

def compute_median(dist_matrix: np.ndarray) -> int:
    return int(np.median(dist_matrix, axis=0).argmin())

strings = ["ACCG", "TCCG", "CGGTA", "TAACG"]
print(f"Levenshtein Distance of Strings '{strings[0]}' and '{strings[1]}' is "
      f"{levenshtein_distance(strings[0], strings[1])}.")
dist_matrix = dist_matrix(strings)
print("Distance Matrix:")
print(dist_matrix)

median_index = compute_median(dist_matrix)
median_str = strings[median_index]
print("Median String:")
print(median_str)
