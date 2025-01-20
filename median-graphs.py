import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def construct_median_graph(sequences):
    graph = nx.Graph()

    for seq in sequences:
        graph.add_node(seq)

    for seq1, seq2 in combinations(sequences, 2):
        if hamming_distance(seq1, seq2) == 1:
            graph.add_edge(seq1, seq2)

    return graph

def find_median_sequence(graph):
    min_distance = float('inf')
    median_sequence = None

    for node in graph.nodes:
        total_distance = sum(nx.shortest_path_length(graph, source=node, target=other)
                             for other in graph.nodes if other != node)
        if total_distance < min_distance:
            min_distance = total_distance
            median_sequence = node

    return median_sequence

sequences = [ "ACGT", "ACGA", "CCGT", "ACGG", "TCGT"]
if len(set(len(seq) for seq in sequences)) != 1:
    raise ValueError("All sequences must be of equal length")
median_graph = construct_median_graph(sequences)
median_sequence = find_median_sequence(median_graph)
print("Median Graph:")
print("Nodes:", median_graph.nodes)
print("Edges:", median_graph.edges)
print("Median Sequence:", median_sequence)

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
