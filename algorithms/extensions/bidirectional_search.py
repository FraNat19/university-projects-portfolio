import heapq
import math
import time
from dataclasses import dataclass, field

import networkx as nx
import pandas as pd


@dataclass
class SearchResult:
    path: list
    cost: float
    expanded: int
    elapsed_ms: float
    meta: dict = field(default_factory=dict)


def _reconstruct(parent, target):
    path = []
    node = target
    while node is not None:
        path.append(node)
        node = parent.get(node)
    return path[::-1]


def dijkstra(G, source, target, weight="weight"):
    t0 = time.perf_counter()
    dist = {source: 0.0}
    parent = {source: None}
    settled = set()
    heap = [(0.0, source)]
    expanded = 0
    while heap:
        d, u = heapq.heappop(heap)
        if u in settled:
            continue
        settled.add(u)
        expanded += 1
        if u == target:
            break
        for v, data in G[u].items():
            if v in settled:
                continue
            nd = d + float(data.get(weight, 1.0))
            if nd < dist.get(v, math.inf):
                dist[v] = nd
                parent[v] = u
                heapq.heappush(heap, (nd, v))
    elapsed = (time.perf_counter() - t0) * 1000
    if target not in dist:
        return SearchResult([], math.inf, expanded, elapsed)
    return SearchResult(_reconstruct(parent, target), dist[target], expanded, elapsed)


def bidirectional_dijkstra(G, source, target, weight="weight"):
    t0 = time.perf_counter()
    if source == target:
        return SearchResult([source], 0.0, 0, 0.0)

    Gr = G.reverse(copy=False) if G.is_directed() else G
    dist = [{source: 0.0}, {target: 0.0}]
    parent = [{source: None}, {target: None}]
    settled = [set(), set()]
    heaps = [[(0.0, source)], [(0.0, target)]]
    adj = [G, Gr]

    best = math.inf
    meeting = None
    expanded = 0

    while heaps[0] and heaps[1]:
        side = 0 if heaps[0][0][0] <= heaps[1][0][0] else 1
        d, u = heapq.heappop(heaps[side])
        if u in settled[side]:
            continue
        settled[side].add(u)
        expanded += 1

        if u in settled[1 - side]:
            total = dist[0].get(u, math.inf) + dist[1].get(u, math.inf)
            if total < best:
                best = total
                meeting = u

        if heaps[0] and heaps[1] and heaps[0][0][0] + heaps[1][0][0] >= best:
            break

        for v, data in adj[side][u].items():
            if v in settled[side]:
                continue
            nd = d + float(data.get(weight, 1.0))
            if nd < dist[side].get(v, math.inf):
                dist[side][v] = nd
                parent[side][v] = u
                heapq.heappush(heaps[side], (nd, v))
                other = dist[1 - side].get(v)
                if other is not None and nd + other < best:
                    best = nd + other
                    meeting = v

    elapsed = (time.perf_counter() - t0) * 1000
    if meeting is None:
        return SearchResult([], math.inf, expanded, elapsed)

    forward = _reconstruct(parent[0], meeting)
    backward = _reconstruct(parent[1], meeting)[::-1]
    return SearchResult(forward + backward[1:], best, expanded, elapsed)


def haversine(a, b):
    lat1, lon1 = a
    lat2, lon2 = b
    r = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(h))


def _coord(G, node):
    d = G.nodes[node]
    return (d["y"], d["x"])


def astar(G, source, target, weight="weight", scale=1.0):
    t0 = time.perf_counter()
    goal = _coord(G, target)
    h = lambda n: haversine(_coord(G, n), goal) * scale
    g = {source: 0.0}
    parent = {source: None}
    settled = set()
    heap = [(h(source), source)]
    expanded = 0
    while heap:
        _, u = heapq.heappop(heap)
        if u in settled:
            continue
        settled.add(u)
        expanded += 1
        if u == target:
            break
        gu = g[u]
        for v, data in G[u].items():
            if v in settled:
                continue
            ng = gu + float(data.get(weight, 1.0))
            if ng < g.get(v, math.inf):
                g[v] = ng
                parent[v] = u
                heapq.heappush(heap, (ng + h(v), v))
    elapsed = (time.perf_counter() - t0) * 1000
    if target not in g:
        return SearchResult([], math.inf, expanded, elapsed)
    return SearchResult(_reconstruct(parent, target), g[target], expanded, elapsed)


def bidirectional_astar(G, source, target, weight="weight", scale=1.0):
    t0 = time.perf_counter()
    if source == target:
        return SearchResult([source], 0.0, 0, 0.0)

    src, dst = _coord(G, source), _coord(G, target)
    hf = lambda n: haversine(_coord(G, n), dst) * scale
    hb = lambda n: haversine(_coord(G, n), src) * scale
    pf = lambda n: (hf(n) - hb(n)) / 2.0
    pot = [pf, lambda n: -pf(n)]

    Gr = G.reverse(copy=False) if G.is_directed() else G
    adj = [G, Gr]
    g = [{source: 0.0}, {target: 0.0}]
    parent = [{source: None}, {target: None}]
    settled = [set(), set()]
    heaps = [[(pot[0](source), source)], [(pot[1](target), target)]]

    best = math.inf
    meeting = None
    expanded = 0

    while heaps[0] and heaps[1]:
        side = 0 if heaps[0][0][0] <= heaps[1][0][0] else 1
        _, u = heapq.heappop(heaps[side])
        if u in settled[side]:
            continue
        settled[side].add(u)
        expanded += 1

        if u in settled[1 - side]:
            total = g[0].get(u, math.inf) + g[1].get(u, math.inf)
            if total < best:
                best = total
                meeting = u

        if heaps[0] and heaps[1] and heaps[0][0][0] + heaps[1][0][0] >= best:
            break

        gu = g[side][u]
        for v, data in adj[side][u].items():
            if v in settled[side]:
                continue
            ng = gu + float(data.get(weight, 1.0))
            if ng < g[side].get(v, math.inf):
                g[side][v] = ng
                parent[side][v] = u
                heapq.heappush(heaps[side], (ng + pot[side](v), v))
                other = g[1 - side].get(v)
                if other is not None and ng + other < best:
                    best = ng + other
                    meeting = v

    elapsed = (time.perf_counter() - t0) * 1000
    if meeting is None:
        return SearchResult([], math.inf, expanded, elapsed)

    forward = _reconstruct(parent[0], meeting)
    backward = _reconstruct(parent[1], meeting)[::-1]
    return SearchResult(forward + backward[1:], best, expanded, elapsed)


def grid_graph(side, spacing_deg=0.001):
    G = nx.Graph()
    for i in range(side):
        for j in range(side):
            n = i * side + j
            G.add_node(n, y=45.0 + i * spacing_deg, x=9.0 + j * spacing_deg)
    for i in range(side):
        for j in range(side):
            n = i * side + j
            for di, dj in ((0, 1), (1, 0)):
                a, b = i + di, j + dj
                if a < side and b < side:
                    m = a * side + b
                    w = haversine(_coord(G, n), _coord(G, m))
                    G.add_edge(n, m, weight=w)
    return G


ALGORITHMS = {
    "dijkstra": dijkstra,
    "bidirectional_dijkstra": bidirectional_dijkstra,
    "astar": astar,
    "bidirectional_astar": bidirectional_astar,
}


def scaling_experiment(sides=(10, 20, 30, 40, 50, 60), repeats=15, seed=42):
    import random

    rng = random.Random(seed)
    rows = []
    for side in sides:
        G = grid_graph(side)
        nodes = list(G.nodes)
        pairs = [(rng.choice(nodes), rng.choice(nodes)) for _ in range(repeats)]
        for name, fn in ALGORITHMS.items():
            costs, times, exps = [], [], []
            for s, t in pairs:
                if s == t:
                    continue
                r = fn(G, s, t)
                costs.append(r.cost)
                times.append(r.elapsed_ms)
                exps.append(r.expanded)
            rows.append(
                {
                    "nodes": G.number_of_nodes(),
                    "edges": G.number_of_edges(),
                    "algorithm": name,
                    "mean_ms": sum(times) / len(times),
                    "mean_expanded": sum(exps) / len(exps),
                    "mean_cost": sum(costs) / len(costs),
                }
            )
    return pd.DataFrame(rows)


def verify_optimality(sides=(10, 20, 30), repeats=20, seed=7, tol=1e-6):
    import random

    rng = random.Random(seed)
    mismatches = []
    for side in sides:
        G = grid_graph(side)
        nodes = list(G.nodes)
        for _ in range(repeats):
            s, t = rng.choice(nodes), rng.choice(nodes)
            if s == t:
                continue
            reference = nx.shortest_path_length(G, s, t, weight="weight")
            for name, fn in ALGORITHMS.items():
                got = fn(G, s, t).cost
                if abs(got - reference) > tol * max(1.0, reference):
                    mismatches.append((side, s, t, name, got, reference))
    return mismatches


if __name__ == "__main__":
    bad = verify_optimality()
    print(f"optimality mismatches: {len(bad)}")
    for row in bad[:10]:
        print(row)

    df = scaling_experiment()
    pivot_time = df.pivot(index="nodes", columns="algorithm", values="mean_ms")
    pivot_exp = df.pivot(index="nodes", columns="algorithm", values="mean_expanded")
    print("\nmean runtime (ms)")
    print(pivot_time.round(3).to_string())
    print("\nmean nodes expanded")
    print(pivot_exp.round(1).to_string())
    print("\nspeedup vs dijkstra (time)")
    print(pivot_time.rdiv(pivot_time["dijkstra"], axis=0).round(2).to_string())
    print("\nexpansion reduction vs dijkstra")
    print(pivot_exp.rdiv(pivot_exp["dijkstra"], axis=0).round(2).to_string())
    df.to_csv("extensions/scaling_benchmark.csv", index=False)
