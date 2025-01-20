import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import networkx as nx
import matplotlib.pyplot as plt

def synthetic_ge_data(num_genes: int, num_samples: int, random_state: int = 42):
    np.random.seed(random_state)
    data = np.random.rand(num_samples, num_genes)
    gene_names = [f"Gene_{i+1}" for i in range(num_genes)]
    return pd.DataFrame(data, columns=gene_names)

def infer_grn(expr_data: pd.DataFrame, n_estimators: int = 100):
    genes = expr_data.columns
    num_genes = len(genes)
    grn = nx.DiGraph()
    grn.add_nodes_from(genes)

    for target_gene in genes:
        X = expr_data.drop(columns=[target_gene])
        y = expr_data[target_gene]

        model = RandomForestRegressor(n_estimators=n_estimators, random_state=42)
        model.fit(X, y)

        importances = model.feature_importances_

        for regulator, importance in zip(X.columns, importances):
            if importance > 0:
                grn.add_edge(regulator, target_gene, weight=importance)

    return grn

def visualize_grn(grn: nx.DiGraph):
    pos = nx.spring_layout(grn)
    edge_weights = nx.get_edge_attributes(grn, 'weight')

    plt.figure(figsize=(12, 8))
    nx.draw(
        grn, pos, with_labels=True, node_color='yellow', edge_color='gray',
        node_size=3000, font_size=10, font_weight='bold'
    )
    nx.draw_networkx_edge_labels(grn, pos, edge_labels={k: f"{v:.2f}" for k, v in edge_weights.items()})
    plt.title("Inferred GRN")
    plt.show()

num_genes = 5
num_samples = 50
expression_data = synthetic_ge_data(num_genes, num_samples)

print("Synthetic Gene Expression Data:")
print(expression_data.head())
print("\nInferring Gene Regulatory Network...")
grn = infer_grn(expression_data)
print("\nVisualizing Gene Regulatory Network...")
visualize_grn(grn)
