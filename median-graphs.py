import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations

def hamming_distance(a, b):
    if len(a) != len(b):
        raise ValueError("Sequences must be of equal length.")
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def construct_median_graph(seqs):
    graph = nx.Graph()

    for seq in seqs:
        graph.add_node(seq)

    for seq1, seq2 in combinations(seqs, 2):
        if hamming_distance(seq1, seq2) == 1:
            graph.add_edge(seq1, seq2)
    return graph

def find_median_seq(graph):
    min_distance = float('inf')
    median_seq = None

    for node in graph.nodes:
        total_distance = sum(nx.shortest_path_length(graph, source=node, target=other)
                             for other in graph.nodes if other != node)
        if total_distance < min_distance:
            min_distance = total_distance
            median_seq = node

    return median_seq

seqs = ["ACGT", "ACGA", "CCGT", "ACGG", "TCGT"]
if len(set(len(seq) for seq in seqs)) != 1:
    raise ValueError("All sequences must be of equal length")
median_graph = construct_median_graph(seqs)
median_seq = find_median_seq(median_graph)
print("Median Graph:")
print("Nodes:", median_graph.nodes)
print("Edges:", median_graph.edges)
print("Median Sequence:", median_seq)

plt.figure(figsize=(5, 5))
pos = nx.spring_layout(median_graph)
nx.draw(
    median_graph,
    pos,
    with_labels=True,
    node_color="yellow",
    node_size=2000,
    font_size=10,
    font_weight="bold",
    edge_color="gray",
)
plt.title("Median Graph")
plt.show()
