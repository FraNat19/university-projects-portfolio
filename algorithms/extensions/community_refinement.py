import time
from collections import defaultdict

import networkx as nx
import pandas as pd


def as_partition(communities):
    return {n: i for i, c in enumerate(communities) for n in c}


def as_communities(partition):
    groups = defaultdict(set)
    for node, cid in partition.items():
        groups[cid].add(node)
    return [groups[k] for k in sorted(groups)]


def disconnected_communities(G, communities):
    bad = []
    for i, c in enumerate(communities):
        if len(c) < 2:
            continue
        sub = G.subgraph(c)
        parts = list(nx.connected_components(sub))
        if len(parts) > 1:
            bad.append({"community": i, "size": len(c), "fragments": len(parts),
                        "largest_fragment": max(len(p) for p in parts)})
    return bad


def split_disconnected(G, communities):
    refined = []
    for c in communities:
        sub = G.subgraph(c)
        for part in nx.connected_components(sub):
            refined.append(set(part))
    return refined


def local_move(G, communities, weight="weight", max_passes=10):
    partition = as_partition(communities)
    m2 = 2.0 * G.size(weight=weight)
    if m2 == 0:
        return communities

    deg = dict(G.degree(weight=weight))
    tot = defaultdict(float)
    for n, cid in partition.items():
        tot[cid] += deg[n]

    improved = True
    passes = 0
    while improved and passes < max_passes:
        improved = False
        passes += 1
        for n in G.nodes():
            cur = partition[n]
            links = defaultdict(float)
            for v, data in G[n].items():
                if v != n:
                    links[partition[v]] += float(data.get(weight, 1.0))

            tot[cur] -= deg[n]
            base = links.get(cur, 0.0) - tot[cur] * deg[n] / m2
            best_c, best_gain = cur, base

            for cid, w in links.items():
                if cid == cur:
                    continue
                gain = w - tot[cid] * deg[n] / m2
                if gain > best_gain:
                    best_c, best_gain = cid, gain

            tot[best_c] += deg[n]
            if best_c != cur:
                partition[n] = best_c
                improved = True

    return [c for c in as_communities(partition) if c]


def refined_louvain(G, weight="weight", seed=42, max_rounds=5):
    communities = nx.community.louvain_communities(G, weight=weight, seed=seed)
    history = [{"round": 0, "stage": "louvain",
                "communities": len(communities),
                "modularity": nx.community.modularity(G, communities, weight=weight),
                "disconnected": len(disconnected_communities(G, communities))}]

    for r in range(1, max_rounds + 1):
        bad = disconnected_communities(G, communities)
        if not bad:
            break
        communities = split_disconnected(G, communities)
        history.append({"round": r, "stage": "split",
                        "communities": len(communities),
                        "modularity": nx.community.modularity(G, communities, weight=weight),
                        "disconnected": len(disconnected_communities(G, communities))})
        communities = local_move(G, communities, weight=weight)
        history.append({"round": r, "stage": "local_move",
                        "communities": len(communities),
                        "modularity": nx.community.modularity(G, communities, weight=weight),
                        "disconnected": len(disconnected_communities(G, communities))})

    return communities, pd.DataFrame(history)


def compare(G, weight="weight", seed=42):
    rows = []

    t0 = time.perf_counter()
    louvain = nx.community.louvain_communities(G, weight=weight, seed=seed)
    t_louvain = (time.perf_counter() - t0) * 1000
    rows.append({"method": "louvain", "communities": len(louvain),
                 "modularity": nx.community.modularity(G, louvain, weight=weight),
                 "disconnected": len(disconnected_communities(G, louvain)),
                 "ms": t_louvain})

    t0 = time.perf_counter()
    refined, _ = refined_louvain(G, weight=weight, seed=seed)
    t_refined = (time.perf_counter() - t0) * 1000
    rows.append({"method": "refined", "communities": len(refined),
                 "modularity": nx.community.modularity(G, refined, weight=weight),
                 "disconnected": len(disconnected_communities(G, refined)),
                 "ms": t_refined})

    t0 = time.perf_counter()
    greedy = list(nx.community.greedy_modularity_communities(G, weight=weight))
    t_greedy = (time.perf_counter() - t0) * 1000
    rows.append({"method": "greedy_modularity", "communities": len(greedy),
                 "modularity": nx.community.modularity(G, greedy, weight=weight),
                 "disconnected": len(disconnected_communities(G, greedy)),
                 "ms": t_greedy})

    return pd.DataFrame(rows)


def benchmark(seeds=range(20), n=500, m=3):
    rows = []
    for s in seeds:
        G = nx.barabasi_albert_graph(n, m, seed=s)
        df = compare(G, seed=s)
        df["graph_seed"] = s
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def ring_of_cliques(cliques=12, size=6):
    G = nx.Graph()
    for c in range(cliques):
        base = c * size
        for i in range(size):
            for j in range(i + 1, size):
                G.add_edge(base + i, base + j)
    for c in range(cliques):
        G.add_edge(c * size, ((c + 1) % cliques) * size)
    return G


def graph_families():
    fams = {
        "barabasi_albert_2000_2": [nx.barabasi_albert_graph(2000, 2, seed=s) for s in range(15)],
        "barabasi_albert_5000_1": [nx.barabasi_albert_graph(5000, 1, seed=s) for s in range(10)],
        "powerlaw_cluster_2000": [nx.powerlaw_cluster_graph(2000, 2, 0.3, seed=s) for s in range(10)],
        "stochastic_block_10x100": [
            nx.stochastic_block_model(
                [100] * 10,
                [[0.08 if i == j else 0.006 for j in range(10)] for i in range(10)],
                seed=s,
            )
            for s in range(10)
        ],
        "random_geometric_1500": [nx.random_geometric_graph(1500, 0.045, seed=s) for s in range(10)],
        "erdos_renyi_1500": [nx.gnp_random_graph(1500, 0.003, seed=s) for s in range(10)],
    }
    lfr = []
    for s in range(8):
        try:
            lfr.append(
                nx.LFR_benchmark_graph(1000, 2.8, 1.4, 0.45, average_degree=12,
                                       min_community=30, seed=s)
            )
        except Exception:
            continue
    if lfr:
        fams["lfr_1000_mu045"] = lfr
    return fams


def pathology_scan(families=None):
    families = families or graph_families()
    rows = []
    for name, graphs in families.items():
        affected = 0
        fragments = 0
        for i, G in enumerate(graphs):
            G = nx.Graph(G)
            G.remove_edges_from(nx.selfloop_edges(G))
            communities = nx.community.louvain_communities(G, seed=i)
            bad = disconnected_communities(G, communities)
            if bad:
                affected += 1
                fragments += sum(b["fragments"] - 1 for b in bad)
        rows.append({"family": name, "graphs": len(graphs),
                     "graphs_affected": affected, "extra_fragments": fragments})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("=== disconnected-community scan across graph families ===")
    scan = pathology_scan()
    print(scan.to_string(index=False))
    print(f"total graphs: {scan['graphs'].sum()}   affected: {scan['graphs_affected'].sum()}")
    scan.to_csv("extensions/pathology_scan.csv", index=False)

    print("\n=== ring of cliques (12 x 6) ===")
    R = ring_of_cliques()
    print(compare(R, seed=1).round(4).to_string(index=False))

    print("\n=== Barabasi-Albert n=500 m=3, 20 seeds ===")
    b = benchmark()
    summary = b.groupby("method").agg(
        communities=("communities", "mean"),
        modularity=("modularity", "mean"),
        disconnected=("disconnected", "mean"),
        seeds_with_disconnected=("disconnected", lambda s: int((s > 0).sum())),
        ms=("ms", "mean"),
    )
    print(summary.round(4).to_string())
    b.to_csv("extensions/community_benchmark.csv", index=False)
