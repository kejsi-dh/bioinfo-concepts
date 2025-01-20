from collections import defaultdict

def build_db_graph(reads, k):
    graph = defaultdict(list)

    for read in reads:
        for i in range(len(read) - k + 1):
            kmer1 = read[i:i + k - 1]
            kmer2 = read[i + 1:i + k]
            graph[kmer1].append(kmer2)

    return graph

def print_db_graph(graph):
    print("De Bruijn Graph:")
    for node, edges in graph.items():
        print(f"{node} -> {', '.join(edges)}")

reads = ["AAG", "AGA", "GAT", "TAA"]
k = 3
de_bruijn_graph = build_db_graph(reads, k)
print_db_graph(de_bruijn_graph)
