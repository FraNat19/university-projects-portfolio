"""Minimum-cut hierarchies for weighted graphs.

Gomory-Hu cut trees built with Gusfield's procedure, and the fragility
analysis they support on dependency networks.

Francesco Natali
"""

import itertools
import time
from dataclasses import dataclass

import networkx as nx
import pandas as pd

CAPACITY = "capacity"


@dataclass
class Boundary:
    value: float
    left: frozenset
    right: frozenset

    @property
    def size(self):
        return min(len(self.left), len(self.right))


def gomory_hu_tree(G, capacity=CAPACITY):
    nodes = list(G.nodes)
    if len(nodes) < 2:
        return nx.Graph((n, n) for n in nodes) if nodes else nx.Graph()

    root = nodes[0]
    parent = {n: root for n in nodes}
    weight = {n: 0.0 for n in nodes}
    calls = 0

    for s in nodes[1:]:
        t = parent[s]
        value, (side_s, _) = nx.minimum_cut(G, s, t, capacity=capacity)
        calls += 1
        weight[s] = value

        for n in nodes:
            if n != s and n in side_s and parent[n] == t:
                parent[n] = s

        if parent[t] in side_s:
            parent[s] = parent[t]
            parent[t] = s
            weight[s] = weight[t]
            weight[t] = value

    tree = nx.Graph()
    tree.add_nodes_from(nodes)
    for n in nodes:
        if n != root:
            tree.add_edge(n, parent[n], weight=weight[n])
    tree.graph["mincut_calls"] = calls
    return tree


def min_cut_value(tree, u, v):
    if u == v:
        return 0.0
    path = nx.shortest_path(tree, u, v)
    return min(tree[a][b]["weight"] for a, b in zip(path, path[1:]))


def all_pairs_min_cut(tree):
    return {
        (u, v): min_cut_value(tree, u, v)
        for u, v in itertools.combinations(sorted(tree.nodes, key=str), 2)
    }


def boundaries(G, tree=None, capacity=CAPACITY):
    tree = gomory_hu_tree(G, capacity) if tree is None else tree
    found = []
    for u, v, data in tree.edges(data=True):
        cut = tree.copy()
        cut.remove_edge(u, v)
        left = frozenset(nx.node_connected_component(cut, u))
        right = frozenset(tree.nodes) - left
        found.append(Boundary(data["weight"], left, right))
    return sorted(found, key=lambda b: b.value)


def weakest_boundaries(G, k=5, min_size=1, capacity=CAPACITY):
    return [b for b in boundaries(G, capacity=capacity) if b.size >= min_size][:k]


def cut_hierarchy(G, capacity=CAPACITY):
    tree = gomory_hu_tree(G, capacity)
    edges = sorted(tree.edges(data=True), key=lambda e: -e[2]["weight"])
    levels = [{"value": None, "clusters": [frozenset(tree.nodes)]}]
    working = tree.copy()
    for u, v, data in edges[::-1]:
        working.remove_edge(u, v)
        clusters = [frozenset(c) for c in nx.connected_components(working)]
        levels.append({"value": data["weight"], "clusters": clusters})
    return levels


def pairwise_baseline(G, capacity=CAPACITY):
    calls = 0
    values = {}
    for u, v in itertools.combinations(sorted(G.nodes, key=str), 2):
        value, _ = nx.minimum_cut(G, u, v, capacity=capacity)
        calls += 1
        values[(u, v)] = value
    return values, calls


def verify(G, capacity=CAPACITY, tol=1e-9):
    tree = gomory_hu_tree(G, capacity)
    reference, _ = pairwise_baseline(G, capacity)
    mismatches = []
    for (u, v), expected in reference.items():
        got = min_cut_value(tree, u, v)
        if abs(got - expected) > tol * max(1.0, expected):
            mismatches.append((u, v, got, expected))
    return mismatches, tree.graph["mincut_calls"], len(reference)


def compare_cost(sizes=(8, 12, 16, 20, 24, 28), seed=11, capacity=CAPACITY):
    import random

    rng = random.Random(seed)
    rows = []
    for n in sizes:
        G = nx.gnp_random_graph(n, 0.45, seed=seed)
        if not nx.is_connected(G):
            G = nx.connected_watts_strogatz_graph(n, 4, 0.3, seed=seed)
        for u, v in G.edges:
            G[u][v][capacity] = rng.randint(1, 20)

        t0 = time.perf_counter()
        tree = gomory_hu_tree(G, capacity)
        t_tree = (time.perf_counter() - t0) * 1000

        t0 = time.perf_counter()
        _, calls_pairwise = pairwise_baseline(G, capacity)
        t_pairwise = (time.perf_counter() - t0) * 1000

        rows.append(
            {
                "nodes": n,
                "edges": G.number_of_edges(),
                "gusfield_calls": tree.graph["mincut_calls"],
                "pairwise_calls": calls_pairwise,
                "gusfield_ms": round(t_tree, 2),
                "pairwise_ms": round(t_pairwise, 2),
                "speedup": round(t_pairwise / t_tree, 2),
            }
        )
    return pd.DataFrame(rows)


def service_graph():
    edges = [
        ("gateway", "auth", 900), ("gateway", "catalog", 1400),
        ("gateway", "checkout", 1100), ("auth", "user-db", 850),
        ("catalog", "search", 700), ("catalog", "product-db", 1300),
        ("search", "product-db", 260), ("checkout", "payments", 640),
        ("checkout", "inventory", 780), ("payments", "ledger", 590),
        ("payments", "fraud", 180), ("fraud", "ledger", 120),
        ("inventory", "product-db", 810), ("inventory", "warehouse", 240),
        ("warehouse", "shipping", 520), ("shipping", "carrier-api", 95),
        ("ledger", "reporting", 150), ("reporting", "analytics-db", 430),
        ("analytics-db", "search", 110),
    ]
    G = nx.Graph()
    for a, b, c in edges:
        G.add_edge(a, b, **{CAPACITY: float(c)})
    return G


def fragility_report(G, k=6, capacity=CAPACITY):
    rows = []
    for b in weakest_boundaries(G, k=k, capacity=capacity):
        smaller = b.left if len(b.left) <= len(b.right) else b.right
        rows.append(
            {
                "cut_capacity": b.value,
                "isolated_services": len(smaller),
                "detaches": ", ".join(sorted(smaller)),
            }
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("=== correctness against exhaustive pairwise minimum cuts ===")
    total_bad = 0
    for n, p, seed in ((10, 0.4, 1), (14, 0.35, 2), (18, 0.3, 3), (22, 0.28, 4)):
        G = nx.gnp_random_graph(n, p, seed=seed)
        if not nx.is_connected(G):
            G = nx.connected_watts_strogatz_graph(n, 4, 0.3, seed=seed)
        import random

        rng = random.Random(seed)
        for u, v in G.edges:
            G[u][v][CAPACITY] = rng.randint(1, 15)
        bad, calls, pairs = verify(G)
        total_bad += len(bad)
        print(f"  n={n:3}  pairs={pairs:4}  gusfield calls={calls:3}  mismatches={len(bad)}")
    print(f"  total mismatches: {total_bad}")

    print("\n=== cost: Gusfield versus exhaustive pairwise ===")
    df = compare_cost()
    print(df.to_string(index=False))
    df.to_csv("cut_hierarchy_cost.csv", index=False)

    print("\n=== fragility of a service dependency graph ===")
    S = service_graph()
    print(f"  {S.number_of_nodes()} services, {S.number_of_edges()} links")
    print(fragility_report(S).to_string(index=False))

    print("\n=== hierarchy levels ===")
    for level in cut_hierarchy(S)[:5]:
        label = "root" if level["value"] is None else f"cut {level['value']:.0f}"
        sizes = sorted((len(c) for c in level["clusters"]), reverse=True)
        print(f"  {label:>10}: {len(level['clusters'])} clusters, sizes {sizes}")
