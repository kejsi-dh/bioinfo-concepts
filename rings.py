import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

# compute ring layers of a graph around a node up to depth
def compute_ring(graph: nx.Graph, node: int, max_depth: int):
    layers = 0
    node_set, in_edge, out_edge = set(), set(), set()
    ring_layers = []
    dist = [np.inf] * graph.number_of_nodes()
    dist[node] = 0
    discovered = set()
    queue = [node]

    while queue:
        next_node = queue.pop(0)
        if dist[next_node] > layers:
            ring_layers.append((node_set, in_edge, out_edge))
            layers += 1
            node_set, in_edge, out_edge = set(), set(), set()

        node_set.add(next_node)
        for nbor in graph.neighbors(next_node):
            edge = (next_node, nbor)
            if edge in discovered or (nbor, next_node) in discovered:
                continue
            discovered.add(edge)

            if dist[nbor] == np.inf:
                dist[nbor] = layers + 1
                if dist[nbor] < max_depth:
                    queue.append(nbor)

            if dist[nbor] == layers:
                in_edge.add(edge)
            else:
                out_edge.add(edge)

    ring_layers.append((node_set, in_edge, out_edge))
    return ring_layers

def test_graph():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (0, 2), (0, 4), (1, 2), (1, 5), (2, 3), (3, 4)])

    pos = {0: (0, 0), 1: (-1, 0.3), 2: (2, 0.17), 3: (4, 0.255), 4: (5, 0.03), 5: (3, 0.5)}

    options = {
        "font_size": 12,
        "node_size": 1000,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 2,
        "width": 2,
    }
    nx.draw_networkx(G, pos, **options)
    plt.gca().margins(0.20)
    plt.axis("off")
    plt.show()

    return G

test_graph = test_graph()
start_node = 0
depth = 3

print(f"Computing layers of ring around node n={start_node} of depth {depth}...")
ring_layers = compute_ring(test_graph, start_node, depth)

print("Computed layers of ring:")
for index, layer in enumerate(ring_layers):
    node_set, in_edge, out_edge = layer
    print(f"\tLayer {index}: N = {node_set}, OE = {out_edge}, IE = {in_edge}")
