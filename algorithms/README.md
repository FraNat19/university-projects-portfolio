<div align="center">

# 🧭 Graph Algorithms

**Routing on real street networks · Social network analysis · Empirical complexity benchmarking**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/NetworkX-2C3E50?style=for-the-badge" alt="NetworkX">
<img src="https://img.shields.io/badge/OSMnx-7EBC6F?style=for-the-badge" alt="OSMnx">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Personal project extending the Algorithms &amp; Data Structures course · Sapienza University of Rome</sub>

</div>

---

## 📋 Overview

Three projects that take textbook algorithms and measure what they actually do, plus a set of extensions written afterwards to chase down a result the original benchmark could not explain.

The thread running through all of it: **asymptotic complexity tells you the shape of the curve, not where you are on it.** Every project here found a case where the theory-optimal choice lost on the machine.

| | Notebook | What it measures |
|---|---|---|
| 1 | [`Algorithm_comparison.ipynb`](Algorithm_comparison.ipynb) | Sorting and traversal, timed and curve-fitted |
| 2 | [`Algorithm_maps.ipynb`](Algorithm_maps.ipynb) | Shortest paths on OpenStreetMap street networks |
| 3 | [`Network_Analysis.ipynb`](Network_Analysis.ipynb) | Centrality and community detection on scale-free graphs |
| + | [`extensions/`](extensions/) | Bidirectional search and partition validation, written after the fact |

📊 **[`benchmark_dashboard.html`](benchmark_dashboard.html)** — every measurement below, plotted and interactive. Download and open it; it is self-contained, no dependencies.

🗺️ **[`route_map.html`](route_map.html)** — a computed route rendered on a Folium map.

---

## ⏱️ Project 1 — Empirical complexity benchmarking

Merge sort, quick sort and heap sort implemented from scratch, timed against CPython's built-in Timsort across array sizes from 100 to 20,000, then fitted to their theoretical curve.

```python
def fit_complexity_curve(sizes, times):
    model = lambda n, a, b: a * n * np.log2(n) + b
    params, _ = curve_fit(model, sizes, times)
    return params, r2_score(times, model(sizes, *params))
```

| n | Merge | Quick | Heap | Timsort |
|---:|---:|---:|---:|---:|
| 100 | 0.296 | 0.132 | 0.272 | **0.009** |
| 1,000 | 3.640 | 1.615 | 5.509 | **0.114** |
| 5,000 | 26.998 | 11.778 | 56.893 | **1.075** |
| 20,000 | 48.383 | 108.558 | 58.427 | **3.514** |

<sub>Mean milliseconds over random integer arrays.</sub>

**Fitted models:** merge sort `T(n) = 0.000159·n·log₂n + 3.65` at R² = 0.876, quick sort `T(n) = 0.000359·n·log₂n − 4.70` at R² = 0.933.

Two things worth saying plainly. The fits are good but not clean, and the reason is visible in the raw data: merge sort is timed *faster* at n = 10,000 than at n = 5,000, which is noise on a shared Colab runtime rather than anything algorithmic. And Timsort beats all three hand-written sorts by one to two orders of magnitude at every size, because it runs as compiled C while the others run in the interpreter. Same complexity class, hundredfold difference in constant.

### Traversal on Erdős–Rényi graphs

| Nodes | Edges | BFS | DFS | Dijkstra | MST |
|---:|---:|---:|---:|---:|---:|
| 50 | 94.8 | 0.300 | 0.164 | 0.681 | 1.686 |
| 500 | 1,546.0 | 1.238 | 0.985 | 3.573 | 10.154 |
| 2,000 | 7,607.0 | 2.704 | 2.488 | 10.200 | 27.078 |

BFS and DFS track each other within noise, as two O(V+E) traversals should. Dijkstra pays for the priority queue, and minimum spanning tree costs roughly three times Dijkstra throughout.

---

## 🗺️ Project 2 — Multi-criteria route optimization

Real street networks pulled from OpenStreetMap with OSMnx, then routed with algorithms implemented from scratch rather than called from a library.

| Algorithm | Complexity | Notes |
|---|---|---|
| Dijkstra | O((V+E) log V) | Uniform expansion, guaranteed optimal |
| A* | O((V+E) log V) | Haversine heuristic, admissible on great-circle distance |
| Multi-objective | O((V+E) log V) | Weighted blend of distance, travel time and toll cost |

The multi-objective variant is the practical one: it takes a weight vector over distance, time and cost, so the same graph answers "fastest", "shortest" and "cheapest" without rebuilding. Traffic is layered on by scaling edge travel times for rush hour and night conditions.

### The result that needed explaining

Benchmarked on 20 random queries over Piedmont, California (352 nodes, 944 edges):

```
Average Dijkstra time: 1.06 ms
Average A* time:       2.93 ms
Average speedup:       0.41x
```

A* came out **slower than Dijkstra**. That is not what the textbook promises, and the notebook offered a plausible explanation — heuristic overhead exceeding the search saved — without evidence. Establishing whether that was actually true is what the extensions below are for.

---

## 🕸️ Project 3 — Social network analysis

Scale-free networks from the Barabási–Albert model, analysed for structure and community.

| Algorithm | Complexity | Purpose |
|---|---|---|
| Degree centrality | O(V) | Raw connection count |
| Betweenness centrality | O(VE) | Bridges between communities |
| PageRank | O(V+E) | Rank weighted by neighbour quality |
| Louvain | O(E log V) | Modularity optimisation |
| Girvan–Newman | O(E²V) | Hierarchical splitting by edge betweenness |

On n = 500, m = 3:

```
Nodes:             500          Avg clustering:    0.056
Edges:           1,491          Communities:          13
Density:         0.012          Modularity:         0.38
```

The degree distribution follows the expected power law, with hub nodes at degree 66 alongside a long tail at degree 3. Modularity of 0.38 indicates real community structure rather than an artefact of the partitioning. Smaller communities come out denser internally (0.09 against 0.06), which is what you would expect from niche clusters.

---

## 🔬 Extensions

Written after the course, to settle the A* question and to check the community partitions properly. Both modules run standalone and print their own results.

### `extensions/bidirectional_search.py`

Bidirectional Dijkstra and bidirectional A*, the latter with the consistent potential `p(v) = (h_t(v) − h_s(v)) / 2` so the two frontiers can meet without losing optimality. Every implementation is checked against `networkx.shortest_path_length` before it is benchmarked:

```
optimality mismatches: 0
```

Benchmarked on square grid graphs, 15 random source-target pairs per size:

| Nodes | A* expansion reduction | A* wall-clock speedup |
|---:|---:|---:|
| 100 | 2.74× | 1.07× |
| 900 | 3.39× | 1.69× |
| 2,500 | 3.89× | 1.92× |
| 3,600 | 3.92× | 2.04× |

**This settles it.** A* expands three to four times fewer nodes than Dijkstra *at every size, including the smallest*. The algorithmic advantage was there all along. What changes with scale is whether that advantage survives contact with the clock: each expansion avoided costs one Haversine evaluation, and at 100 nodes that trade is nearly break-even at 1.07×. By 3,600 nodes it is 2.04× and still climbing.

So the original finding was real but the stated reason was only half right. A* was not doing more work; it was doing less work more expensively. On the 352-node network in Project 2, the trade had not yet paid off.

Bidirectional A* is the slowest of the four in wall clock until the largest size tested, for the same reason amplified: its potential function evaluates the Haversine twice per node.

### `extensions/community_refinement.py`

Louvain can return communities that are internally disconnected (Traag, Waltman & van Eck, 2019). This module implements the diagnostic, a splitting refinement and a local-move pass, then tests whether the problem actually occurs.

```
family                    graphs  affected  extra_fragments
barabasi_albert_2000_2        15         0                0
barabasi_albert_5000_1        10         0                0
powerlaw_cluster_2000         10         0                0
stochastic_block_10x100       10         0                0
random_geometric_1500         10         0                0
erdos_renyi_1500              10         0                0
lfr_1000_mu045                 7         0                0
total graphs: 72   affected: 0
```

**The refinement is a no-op here, and that is the finding.** Across 72 graphs from seven families, the NetworkX Louvain implementation never produced a disconnected community, so there was nothing to repair. The pathology is documented and real; it did not reproduce on any topology tested. Reporting that is more useful than shipping a fix for a problem that is not occurring.

The comparison the module did produce is worth having. Over 20 Barabási–Albert graphs:

| Method | Communities | Modularity | Time |
|---|---:|---:|---:|
| Louvain | 12.9 | **0.3830** | **261 ms** |
| Greedy modularity | 11.5 | 0.3818 | 1,037 ms |

Louvain wins on both axes: marginally better modularity, four times faster.

---

## ⚙️ Running it

```bash
pip install networkx pandas numpy matplotlib seaborn scipy osmnx folium python-louvain

python extensions/bidirectional_search.py
python extensions/community_refinement.py
```

Both write their results to CSV next to the source. The notebooks were written in Colab; `Algorithm_maps.ipynb` downloads its street network at runtime and needs a working connection.

---

## 📖 References

- Dijkstra, E.W. (1959). A note on two problems in connexion with graphs. *Numerische Mathematik*, 1, 269–271.
- Hart, P.E., Nilsson, N.J. & Raphael, B. (1968). A formal basis for the heuristic determination of minimum cost paths. *IEEE Transactions on Systems Science and Cybernetics*, 4(2), 100–107.
- Goldberg, A.V. & Harrelson, C. (2005). Computing the shortest path: A* search meets graph theory. *SODA*, 156–165.
- Blondel, V.D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *J. Stat. Mech.*, P10008.
- Traag, V.A., Waltman, L. & van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
- Barabási, A.-L. & Albert, R. (1999). Emergence of scaling in random networks. *Science*, 286(5439), 509–512.
