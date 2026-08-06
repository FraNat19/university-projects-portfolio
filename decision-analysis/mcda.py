"""Multi-criteria decision analysis methods, reimplemented from the course models.

PROMETHEE II, AHP with consistency checking, and the Fishbein multi-attribute
attitude model. Each function reproduces one of the spreadsheet models in this
folder, so the Excel results can be checked rather than trusted.

Francesco Natali
"""

import itertools
from dataclasses import dataclass

import numpy as np
import pandas as pd

RANDOM_INDEX = {1: 0.0, 2: 0.0, 3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24,
                7: 1.32, 8: 1.41, 9: 1.45, 10: 1.49}


# ------------------------------------------------------------------ PROMETHEE

@dataclass
class Criterion:
    name: str
    weight: float
    maximise: bool
    q: float = 0.0
    p: float = 0.0


def linear_preference(d, q, p):
    if d <= q:
        return 0.0
    if d >= p:
        return 1.0
    return (d - q) / (p - q)


def unicriterion_flows(values, criterion):
    n = len(values)
    sign = 1.0 if criterion.maximise else -1.0
    pref = np.zeros((n, n))
    for i, j in itertools.permutations(range(n), 2):
        d = sign * (values[i] - values[j])
        pref[i, j] = linear_preference(d, criterion.q, criterion.p)
    positive = pref.sum(axis=1) / (n - 1)
    negative = pref.sum(axis=0) / (n - 1)
    return positive - negative, pref


def promethee2(table, criteria):
    alternatives = list(table.index)
    n = len(alternatives)
    net = {}
    weighted = np.zeros(n)
    aggregated = np.zeros((n, n))

    for c in criteria:
        flows, pref = unicriterion_flows(table[c.name].to_numpy(float), c)
        net[c.name] = flows
        weighted += c.weight * flows
        aggregated += c.weight * pref

    phi_plus = aggregated.sum(axis=1) / (n - 1)
    phi_minus = aggregated.sum(axis=0) / (n - 1)

    result = pd.DataFrame(net, index=alternatives)
    result["phi+"] = phi_plus
    result["phi-"] = phi_minus
    result["net_from_flows"] = phi_plus - phi_minus
    result["net_from_weights"] = weighted
    result["rank"] = result["net_from_flows"].rank(ascending=False).astype(int)
    return result.sort_values("rank")


# ------------------------------------------------------------------------ AHP

def priority_vector(matrix):
    m = np.asarray(matrix, dtype=float)
    geometric = np.prod(m, axis=1) ** (1.0 / m.shape[0])
    return geometric / geometric.sum()


def consistency(matrix):
    m = np.asarray(matrix, dtype=float)
    n = m.shape[0]
    w = priority_vector(m)
    lam = float((m @ w / w).mean())
    ci = (lam - n) / (n - 1) if n > 1 else 0.0
    ri = RANDOM_INDEX.get(n, 1.49)
    cr = ci / ri if ri else 0.0
    return {"lambda_max": lam, "CI": ci, "RI": ri, "CR": cr, "acceptable": cr < 0.10}


def ahp(criteria_matrix, criteria_names, alternative_matrices, alternative_names):
    weights = priority_vector(criteria_matrix)
    local = pd.DataFrame(
        {name: priority_vector(alternative_matrices[name]) for name in criteria_names},
        index=alternative_names,
    )
    global_scores = local.mul(weights, axis=1).sum(axis=1)
    checks = {name: consistency(alternative_matrices[name]) for name in criteria_names}
    checks["_criteria"] = consistency(criteria_matrix)
    return pd.Series(weights, index=criteria_names), local, global_scores, checks


def ahp_hierarchy(goal_matrix, goal_names, subtrees, alternative_matrices, alternative_names):
    top = pd.Series(priority_vector(goal_matrix), index=goal_names)
    leaf_weights = {}
    checks = {"GOAL": consistency(goal_matrix)}

    for node, weight in top.items():
        if node in subtrees:
            names, matrix = subtrees[node]
            child = pd.Series(priority_vector(matrix), index=names)
            checks[node] = consistency(matrix)
            for name, w in child.items():
                leaf_weights[name] = weight * w
        else:
            leaf_weights[node] = weight

    leaf_weights = pd.Series(leaf_weights)
    local = pd.DataFrame(
        {name: priority_vector(alternative_matrices[name]) for name in leaf_weights.index},
        index=alternative_names,
    )
    for name in leaf_weights.index:
        checks[name] = consistency(alternative_matrices[name])

    scores = local.mul(leaf_weights, axis=1).sum(axis=1)
    return leaf_weights, local, scores, checks


# ------------------------------------------------------------------------ MAA

def fishbein(importance, ratings):
    scores = ratings.mul(importance, axis=0)
    totals = scores.sum()
    out = pd.DataFrame({"score": totals, "rank": totals.rank(ascending=False).astype(int)})
    return out.sort_values("rank"), scores


# ------------------------------------------------------------------- datasets

def phones_promethee():
    table = pd.DataFrame(
        {
            "price": [1489, 239, 100, 1399, 147],
            "weight": [227, 200, 200, 226, 191],
            "videocamera": [10, 7, 3, 8, 5],
            "screen": [6.9, 6.7, 6.88, 6.8, 6.8],
        },
        index=["iPhone 16 Pro Max", "Galaxy A16 5G", "Xiaomi Redmi A3 Pro",
               "Huawei Pura 70 Ultra", "Honor 200 Smart"],
    )
    criteria = [
        Criterion("price", 0.35, False, 100, 250),
        Criterion("weight", 0.15, False, 4, 9),
        Criterion("videocamera", 0.25, True, 1, 2),
        Criterion("screen", 0.25, True, 0.0, 0.04),
    ]
    return table, criteria


def phones_maa():
    importance = pd.Series(
        [5, 8, 6, 7, 6, 3, 5, 3, 9],
        index=["Network", "Display", "Platform", "Memory", "Main Camera",
               "Selfie Camera", "Sounds", "Communications", "Price"],
    )
    ratings = pd.DataFrame(
        {
            "iPhone 16 Pro Max":     [8, 8, 9, 8, 8, 7, 8, 7, 4],
            "Galaxy A16 5G":         [8, 7, 7, 6, 9, 7, 7, 7, 7],
            "Xiaomi Redmi A3 Pro":   [5, 5, 7, 6, 9, 5, 6, 7, 9],
            "Huawei Pura 70 Ultra":  [8, 8, 7, 9, 9, 7, 5, 7, 2],
            "Honor 200 Smart":       [7, 7, 7, 7, 9, 3, 5, 7, 7],
            "Asus Zenfone 11 Ultra": [7, 6, 7, 7, 9, 8, 4, 7, 4],
            "Oppo K12 Plus":         [7, 6, 7, 8, 9, 7, 3, 7, 7],
            "LG W31":                [3, 3, 4, 3, 4, 4, 1, 6, 8],
        },
        index=importance.index,
    )
    return importance, ratings


def campus_parking_ahp():
    goal_matrix = [[1, 2, 5],
                   [0.5, 1, 3],
                   [0.2, 1 / 3, 1]]
    goal_names = ["CUSTOMER", "CONTROL", "INVESTMENT"]

    subtrees = {
        "CUSTOMER": (["SERVICE", "COST", "CONVENIENCE"],
                     [[1, 2, 5], [0.5, 1, 3], [0.2, 1 / 3, 1]])
    }

    alternatives = ["Large area outside", "Garage", "Small area inside"]
    matrices = {
        "SERVICE":     [[1, 0.4, 2 / 3], [2.5, 1, 5 / 3], [1.5, 0.6, 1]],
        "COST":        [[1, 2 / 3, 2], [1.5, 1, 3], [0.5, 1 / 3, 1]],
        "CONVENIENCE": [[1, 2 / 3, 2], [1.5, 1, 3], [0.5, 1 / 3, 1]],
        "CONTROL":     [[1, 1 / 7, 1 / 3], [7, 1, 5], [3, 0.2, 1]],
        "INVESTMENT":  [[1, 1 / 3, 4], [3, 1, 12], [0.25, 1 / 12, 1]],
    }
    return goal_matrix, goal_names, subtrees, matrices, alternatives


if __name__ == "__main__":
    pd.set_option("display.width", 120)

    print("=" * 74)
    print("PROMETHEE II - five phones, four criteria")
    table, criteria = phones_promethee()
    res = promethee2(table, criteria)
    print(res.round(4).to_string())

    gap = (res["net_from_flows"] - res["net_from_weights"]).abs().max()
    print(f"\nmax |phi+ - phi-  minus  weighted unicriterion sum| = {gap:.2e}")
    print("the two routes to the net flow agree" if gap < 1e-9 else "ROUTES DISAGREE")

    excel_from_flows = {"iPhone 16 Pro Max": 0.0938, "Galaxy A16 5G": 0.0022,
                        "Xiaomi Redmi A3 Pro": 0.1415, "Huawei Pura 70 Ultra": -0.375,
                        "Honor 200 Smart": 0.1375}
    excel_weighted = {"iPhone 16 Pro Max": 0.1313, "Galaxy A16 5G": -0.0227,
                      "Xiaomi Redmi A3 Pro": 0.0228, "Huawei Pura 70 Ultra": -0.525,
                      "Honor 200 Smart": 0.2625}
    comparison = pd.DataFrame({
        "python": res["net_from_flows"],
        "excel_phi_difference": pd.Series(excel_from_flows),
        "excel_weighted_sum": pd.Series(excel_weighted),
    })
    comparison["diff_vs_phi"] = (comparison["python"] - comparison["excel_phi_difference"]).abs()
    comparison["diff_vs_weighted"] = (comparison["python"] - comparison["excel_weighted_sum"]).abs()
    print("\ncomparison with the spreadsheet")
    print(comparison.round(4).to_string())

    print("\n" + "=" * 74)
    print("AHP - campus parking, two-level hierarchy")
    gm, gn, subs, mats, alts = campus_parking_ahp()
    leaf_weights, local, scores, checks = ahp_hierarchy(gm, gn, subs, mats, alts)
    print("\nglobal criterion weights")
    print(leaf_weights.round(4).to_string())
    print("\nlocal priorities per criterion")
    print(local.round(4).to_string())
    print("\noverall scores")
    print(scores.sort_values(ascending=False).round(4).to_string())

    excel_utilities = {"Large area outside": 0.1988, "Garage": 0.5898,
                       "Small area inside": 0.2055}
    delta = (scores - pd.Series(excel_utilities)).abs().max()
    print(f"\nlargest gap against the spreadsheet utilities: {delta:.4f}")

    print("\nconsistency")
    print(pd.DataFrame(checks).T.round(4).to_string())

    print("\n" + "=" * 74)
    print("Fishbein multi-attribute attitude model - eight phones")
    imp, rat = phones_maa()
    ranking, _ = fishbein(imp, rat)
    print(ranking.to_string())
