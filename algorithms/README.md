# Graph Algorithms: Routing Optimization & Social Network Analysis

> Personal project extending the Algorithms & Data Structures course | Sapienza University of Rome

This repository contains two independent projects implementing classic graph algorithms on real-world open-source datasets. The goal was to go beyond theoretical concepts and apply them to practical problems: GPS routing and social network analysis.

---

## Project 1: Multi-Criteria Route Optimization

Implementation of shortest-path algorithms on real street networks from OpenStreetMap.

### What it does

- Downloads actual city street networks using OSMnx
- Implements Dijkstra's algorithm from scratch (guaranteed optimal)
- Implements A* algorithm with Haversine heuristic (goal-directed speedup)
- Supports multi-objective optimization (distance, time, toll cost)
- Simulates real-time traffic conditions (rush hour, night time)
- Generates interactive maps with Folium

### Algorithms Implemented

| Algorithm | Complexity | Description |
|-----------|------------|-------------|
| **Dijkstra** | O((V+E) log V) | Classic shortest path, explores all directions equally |
| **A*** | O((V+E) log V) | Heuristic-guided search, typically 2-10× faster |
| **Multi-objective** | O((V+E) log V) | Weighted combination of distance, time, and cost |

### Key Findings

The benchmark on 20 random queries revealed an interesting result: **A* was slower than Dijkstra on small networks** (~350 nodes). This happens because the Haversine heuristic calculation overhead exceeds the time saved by exploring fewer nodes. On larger networks (1000+ nodes), A* consistently outperforms Dijkstra.

```
Average Dijkstra time: 1.06 ms
Average A* time:       2.93 ms
Average speedup:       0.41× (A* slower on small graphs!)
```

### Data Source

- **OpenStreetMap** via OSMnx library
- Test network: Piedmont, California (352 nodes, 944 edges)
- Scalable to any city worldwide

### Output

- Interactive HTML map with route visualization
- Performance benchmark CSV

---

## Project 2: Social Network Analysis

Implementation of community detection and centrality algorithms on synthetic social networks.

### What it does

- Generates scale-free networks using Barabási-Albert model
- Computes multiple centrality metrics to identify influencers
- Detects communities using Louvain and Girvan-Newman algorithms
- Analyzes degree distribution to confirm power-law behavior
- Creates publication-quality network visualizations

### Algorithms Implemented

| Algorithm | Complexity | Purpose |
|-----------|------------|---------|
| **Degree Centrality** | O(V) | Identify nodes with most connections |
| **Betweenness Centrality** | O(VE) | Find bridge nodes connecting communities |
| **PageRank** | O(V+E) | Rank nodes by connection quality |
| **Louvain** | O(E log V) | Fast community detection via modularity optimization |
| **Girvan-Newman** | O(E²V) | Hierarchical community detection via edge removal |

### Key Findings

The Barabási-Albert model (n=500, m=3) produced a network with clear scale-free properties:

```
Nodes:              500
Edges:              1,491
Density:            0.012 (sparse, typical for social networks)
Avg Clustering:     0.056
Communities found:  13 (Louvain)
Modularity:         0.38 (strong community structure)
```

The degree distribution follows a power-law P(k) ∝ k⁻³, confirming the "rich get richer" preferential attachment mechanism. A few hub nodes (degree 66) coexist with many peripheral nodes (degree 3).

### Community Analysis

| Community | Size | Density | Top Influencer |
|-----------|------|---------|----------------|
| 0 | 60 | 0.057 | Node 7 |
| 6 | 57 | 0.058 | Node 6 |
| 7 | 49 | 0.058 | Node 10 |
| ... | ... | ... | ... |

Smaller communities show higher internal density (0.09 vs 0.06), reflecting tighter cohesion in niche groups.

### Data Source

- **Synthetic network**: Barabási-Albert model (configurable)
- Can load real data from Kaggle Facebook dataset

### Output

- Network visualization with community coloring
- Centrality metrics CSV
- Degree distribution plots

---

## 🛠️ Tech Stack

| Component | Technology |
|-----------|------------|
| Language | Python 3.10+ |
| Graph Library | NetworkX |
| Street Networks | OSMnx |
| Visualization | Matplotlib, Folium, Seaborn |
| Community Detection | python-louvain |
| Data Processing | Pandas, NumPy |


## 📚 Theoretical Background

### Dijkstra vs A*

Both algorithms find optimal shortest paths, but differ in exploration strategy:

- **Dijkstra**: Explores nodes in order of distance from source (uniform expansion)
- **A***: Explores nodes in order of f(n) = g(n) + h(n), where h(n) is a heuristic estimate to target

A* is faster when the heuristic is admissible (never overestimates) and the graph is large enough to benefit from guided search.

### Community Detection

Communities are groups of nodes with dense internal connections and sparse external connections. The **modularity** metric Q ∈ [-0.5, 1] measures partition quality:

- Q ≈ 0: Random partition
- Q > 0.3: Significant community structure
- Q > 0.7: Very strong communities

### Scale-Free Networks

Real social networks exhibit power-law degree distributions P(k) ∝ k⁻ᵞ with γ ≈ 2-3. This means:

- Few "hub" nodes with many connections
- Many peripheral nodes with few connections
- Network is robust to random failures but vulnerable to targeted attacks on hubs


## 📖 References

- Dijkstra, E.W. (1959). A note on two problems in connexion with graphs
- Hart, P.E., Nilsson, N.J., Raphael, B. (1968). A Formal Basis for the Heuristic Determination of Minimum Cost Paths
- Blondel, V.D. et al. (2008). Fast unfolding of communities in large networks (Louvain)
- Girvan, M., Newman, M.E.J. (2002). Community structure in social and biological networks
- Barabási, A.L., Albert, R. (1999). Emergence of scaling in random networks

---

## 👤 Author

**Francesco Natali**  
Master's Degree in Statistical Methods and Applications  
Sapienza University of Rome

*Course: Algorithms & Data Structures (A.Y. 2024/2025)*
