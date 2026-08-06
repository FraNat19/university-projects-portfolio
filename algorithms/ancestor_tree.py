# Core Ancestor Tree implementation (from Leonardo Agate's MSc thesis, ATGrafovero notebook)
import networkx as nx
import random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CUT_CALLS = {"count": 0}

def total_capacity(G):
    return sum(float(d.get("capacity", 0.0)) for _, _, d in G.edges(data=True))

def mincut_fst_inplace(Gu, s, t, x, y, M):
    if x == y:
        return (0.0, ({x}, set(Gu.nodes()) - {x}))
    CUT_CALLS["count"] += 2
    Gtemp = Gu.copy()
    for node in [s, t, x, y]:
        if node not in Gtemp: Gtemp.add_node(node)
    Gtemp.add_edge(x, s, capacity=M); Gtemp.add_edge(y, t, capacity=M)
    try: v1, part1 = nx.minimum_cut(Gtemp, s, t, capacity="capacity")
    except nx.NetworkXError: v1, part1 = float("inf"), (set(), set())
    Gtemp = Gu.copy()
    for node in [s, t, x, y]:
        if node not in Gtemp: Gtemp.add_node(node)
    Gtemp.add_edge(x, t, capacity=M); Gtemp.add_edge(y, s, capacity=M)
    try: v2, part2 = nx.minimum_cut(Gtemp, s, t, capacity="capacity")
    except nx.NetworkXError: v2, part2 = float("inf"), (set(), set())
    return (v1, part1) if v1 <= v2 else (v2, part2)

class ATNode:
    def __init__(self, members):
        self.members = set(members); self.value = None
        self.left = None; self.right = None; self.parent = None

def build_ancestor_tree_fst(G, s, t, seed=42, tol=1e-9):
    random.seed(seed)
    all_nodes = list(G.nodes)
    root = ATNode(all_nodes); leaves = [root]
    M = total_capacity(G) + 1.0; cuts_done = 0
    def find_insertion_point(node, cut_value):
        x = node
        while (x.parent is not None and x.parent.value is not None and cut_value < x.parent.value - tol):
            x = x.parent
        return x
    while cuts_done < len(all_nodes) - 1:
        leaf = next((L for L in leaves if len(L.members) > 1), None)
        if leaf is None: break
        members = list(leaf.members)
        best_val = float("inf"); best_part = None
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                x, y = members[i], members[j]
                val, part = mincut_fst_inplace(G, s, t, x, y, M)
                if val < best_val - tol:
                    best_val, best_part = val, part
        if best_part is None:
            leaves.remove(leaf); continue
        S, T = best_part
        left_members = {z for z in leaf.members if z in S}
        right_members = leaf.members - left_members
        if not left_members or not right_members:
            leaves.remove(leaf); continue
        insertion = find_insertion_point(leaf, best_val)
        if insertion.value is not None and best_val < insertion.value - tol:
            new_parent = ATNode(insertion.members); new_parent.value = best_val
            new_parent.left = ATNode(left_members); new_parent.right = ATNode(right_members)
            new_parent.left.parent = new_parent; new_parent.right.parent = new_parent
            parent = insertion.parent
            if parent:
                if parent.left is insertion: parent.left = new_parent
                else: parent.right = new_parent
            new_parent.parent = parent; insertion.parent = new_parent
            if root is insertion: root = new_parent
            cuts_done += 1; continue
        leaf.value = best_val
        Lnode = ATNode(left_members); Rnode = ATNode(right_members)
        Lnode.parent = leaf; Rnode.parent = leaf
        leaf.left = Lnode; leaf.right = Rnode
        if leaf in leaves: leaves.remove(leaf)
        leaves.extend([Lnode, Rnode]); cuts_done += 1
    return root

def lca_value(root, u, v):
    def path_to_root(x):
        path = []
        def dfs(node):
            if node is None: return
            if x in node.members:
                path.append(node)
                if node.left and x in node.left.members: dfs(node.left)
                elif node.right and x in node.right.members: dfs(node.right)
        dfs(root); return path
    pu, pv = path_to_root(u), path_to_root(v)
    if not pu or not pv: return None
    lca = None
    for na, nb in zip(pu, pv):
        if na is nb: lca = na
        else: break
    while lca and lca.value is None: lca = lca.parent
    return None if lca is None else lca.value

# ---- Vitality ----
def maxflow_value(G, s, t):
    return nx.maximum_flow_value(G, s, t, capacity="capacity")

def vitality_bruteforce(G, s, t):
    base = maxflow_value(G, s, t); vit = {}
    for (u, v) in list(G.edges()):
        H = G.copy(); H.remove_edge(u, v)
        try: f = maxflow_value(H, s, t)
        except nx.NetworkXError: f = 0.0
        vit[(u, v)] = round(base - f, 6)
    return base, vit

def vitality_ancestor_tree(G, s, t, seed=42):
    root = build_ancestor_tree_fst(G, s, t, seed=seed)
    base = maxflow_value(G, s, t); vit = {}
    for (u, v, d) in G.edges(data=True):
        ce = float(d.get("capacity", 0.0))
        mc = lca_value(root, u, v)
        vit[(u, v)] = 0.0 if mc is None else round(max(0.0, base - (mc - ce)), 6)
    return base, vit, root
