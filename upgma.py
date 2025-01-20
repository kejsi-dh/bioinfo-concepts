import numpy as np

def upgma(dist_matrix, labels):
    dist_matrix = np.array(dist_matrix)
    n = len(labels)

    clusters = {i: labels[i] for i in range(n)}
    tree = {labels[i]: {} for i in range(n)}

    while len(clusters) > 1:
        min_dist = float("inf")
        x, y = -1, -1
        cluster_indices = list(clusters.keys())
        for i in range(len(cluster_indices)):
            for j in range(i + 1, len(cluster_indices)):
                a, b = cluster_indices[i], cluster_indices[j]
                if dist_matrix[i, j] < min_dist:
                    min_dist = dist_matrix[i, j]
                    x, y = a, b

        new_cluster_label = f"({clusters[x]},{clusters[y]})"
        tree[new_cluster_label] = {clusters[x]: float(min_dist / 2), clusters[y]: float(min_dist / 2)}

        print(f"Tree after merging {clusters[x]} and {clusters[y]}:")
        for k, v in tree.items():
            print(f"  {k}: {v}")
        print()

        new_dists = []
        for k in cluster_indices:
            if k != x and k != y:
                new_dist = (dist_matrix[cluster_indices.index(x), cluster_indices.index(k)] +
                            dist_matrix[cluster_indices.index(y), cluster_indices.index(k)]) / 2
                new_dists.append(new_dist)

        new_dist_matrix = []
        for i in range(len(cluster_indices)):
            if cluster_indices[i] != x and cluster_indices[i] != y:
                new_row = [dist_matrix[i, j] for j in range(len(cluster_indices))
                           if cluster_indices[j] != x and cluster_indices[j] != y]
                new_dist_matrix.append(new_row)

        new_dist_matrix = np.array(new_dist_matrix)
        if len(new_dist_matrix) > 0:
            new_dist_matrix = np.vstack([new_dist_matrix, new_dists])
            new_dists.append(0)
            new_dist_matrix = np.column_stack([new_dist_matrix, new_dists])

        dist_matrix = new_dist_matrix
        clusters[x] = new_cluster_label
        del clusters[y]

    return tree

dist_matrix = [
    [0.0, 5.0, 9.0, 9.0, 8.0],
    [5.0, 0.0, 10.0, 10.0, 9.0],
    [9.0, 10.0, 0.0, 8.0, 7.0],
    [9.0, 10.0, 8.0, 0.0, 3.0],
    [8.0, 9.0, 7.0, 3.0, 0.0]
]
labels = ["A", "B", "C", "D", "E"]

phylogenetic_tree = upgma(dist_matrix, labels)
print("Phylogenetic Tree:")
for k, v in phylogenetic_tree.items():
    print(f"{k}: {v}")
