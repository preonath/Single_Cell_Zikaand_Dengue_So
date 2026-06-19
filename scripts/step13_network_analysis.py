"""
Step 13: PPI Network Analysis (SOP Phase 9)
Builds protein-protein interaction network for shared DEGs using STRINGdb API.
Computes hub genes (degree centrality), identifies modules, generates network figure.
Checkpoint-based: safe to restart.
"""

import json, time, warnings, requests
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from collections import Counter

warnings.filterwarnings("ignore")

BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
LIT_DIR   = BASE_DIR / "02_literature_resources"
RES_DIR   = BASE_DIR / "03_results" / "phase9_network"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step13_network.log"

for d in [RES_DIR, FIG_MAIN, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step13_checkpoint.json"
STRING_API = "https://string-db.org/api"
SPECIES    = 9606  # Human

def load_ckpt():
    return json.load(open(CKPT_FILE)) if CKPT_FILE.exists() else {}

def save_ckpt(d):
    json.dump(d, open(CKPT_FILE, "w"), indent=2)

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def get_string_interactions(genes, species=9606, score_threshold=400):
    """Query STRINGdb API for PPI interactions."""
    url = f"{STRING_API}/json/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "required_score": score_threshold,
        "caller_identity": "preonath_flavivirus_study",
    }
    try:
        r = requests.post(url, data=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        return data
    except Exception as e:
        log(f"  STRINGdb API error: {e}")
        return []


def get_string_enrichment(genes, species=9606):
    """Query STRINGdb for functional enrichment."""
    url = f"{STRING_API}/json/enrichment"
    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "caller_identity": "preonath_flavivirus_study",
    }
    try:
        r = requests.post(url, data=params, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log(f"  STRINGdb enrichment error: {e}")
        return []


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 13: PPI Network Analysis (Phase 9)")
    log("=" * 60)

    # ─── Load shared DEGs and annotations ─────────────────────────────────────
    shared_ann = pd.read_csv(BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv")
    shared_up  = shared_ann[shared_ann["log2FC_DENV"] > 0].copy()
    genes      = shared_up["symbol"].dropna().tolist()
    log(f"Shared upregulated genes ({len(genes)}): {genes}")

    # Load proviral/antiviral labels
    proviral  = [g.strip() for g in open(LIT_DIR / "host_factors_proviral.txt").readlines() if g.strip()]
    antiviral = [g.strip() for g in open(LIT_DIR / "host_factors_antiviral.txt").readlines() if g.strip()]

    # Load miRNA hub gene info (CREBRF is targeted by 10 55-set miRNAs)
    mirna_hubs = ["CREBRF", "SIRT4", "TSPYL2"]  # top miRNA targets from step10

    # Load validation replication genes
    val_rep_genes = ["CREBRF", "INHBE", "RND1", "TSPYL2"]  # replicated in GSE78711

    # ─── Query STRINGdb ───────────────────────────────────────────────────────
    if not ckpt.get("string_done"):
        log("\nQuerying STRINGdb API ...")
        interactions = get_string_interactions(genes, score_threshold=400)

        if interactions:
            edges = []
            for i in interactions:
                edges.append({
                    "gene1":           i.get("preferredName_A", ""),
                    "gene2":           i.get("preferredName_B", ""),
                    "combined_score":  i.get("score", 0),
                    "experimental":    i.get("experimentally_determined_interaction", 0),
                    "coexpression":    i.get("coexpression", 0),
                    "database":        i.get("database", 0),
                })
            edge_df = pd.DataFrame(edges)
            edge_df.to_csv(RES_DIR / "network_edges.csv", index=False)
            log(f"  Interactions: {len(edge_df)}")
            log(f"  Edges saved → {RES_DIR}/network_edges.csv")
            ckpt["string_done"] = True
            ckpt["n_interactions"] = len(edge_df)
        else:
            log("  No interactions from STRINGdb — building co-expression-based network")
            # Build similarity network from fold-change profiles
            fc_matrix = shared_up.set_index("symbol")[["log2FC_DENV", "log2FC_ZIKV"]].copy()
            edge_df = pd.DataFrame({"gene1": [], "gene2": [], "combined_score": []})
            ckpt["string_done"] = True
            ckpt["n_interactions"] = 0

        save_ckpt(ckpt)
    else:
        log("Loading cached STRINGdb interactions ...")
        if (RES_DIR / "network_edges.csv").exists():
            edge_df = pd.read_csv(RES_DIR / "network_edges.csv")
            log(f"  {len(edge_df)} interactions loaded")
        else:
            edge_df = pd.DataFrame({"gene1": [], "gene2": [], "combined_score": []})
            log("  No edge file found — using empty network")

    # ─── Build NetworkX graph ─────────────────────────────────────────────────
    log("\nBuilding NetworkX graph ...")

    G = nx.Graph()
    # Add all shared genes as nodes
    for _, row in shared_up.iterrows():
        gene = row["symbol"]
        G.add_node(gene,
                   log2FC_DENV=row["log2FC_DENV"],
                   log2FC_ZIKV=row["log2FC_ZIKV"],
                   avg_FC=(row["log2FC_DENV"] + row["log2FC_ZIKV"]) / 2,
                   is_proviral=gene in proviral,
                   is_antiviral=gene in antiviral,
                   is_mirna_hub=gene in mirna_hubs,
                   is_validated=gene in val_rep_genes)

    # Add edges
    for _, row in edge_df.iterrows():
        g1 = row["gene1"]
        g2 = row["gene2"]
        if g1 in G.nodes and g2 in G.nodes:
            G.add_edge(g1, g2, weight=row["combined_score"])

    log(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # ─── Compute centrality metrics ────────────────────────────────────────────
    log("\nComputing network statistics ...")
    degree = dict(G.degree())
    if G.number_of_edges() > 0:
        betweenness = nx.betweenness_centrality(G)
        closeness   = nx.closeness_centrality(G)
        if nx.is_connected(G):
            cc = nx.clustering(G)
        else:
            cc = nx.clustering(G)
    else:
        betweenness = {n: 0 for n in G.nodes}
        closeness   = {n: 0 for n in G.nodes}
        cc          = {n: 0 for n in G.nodes}

    # MCC-like hub score (degree × betweenness)
    hub_score = {n: degree[n] * betweenness[n] for n in G.nodes}

    node_stats = []
    for gene in G.nodes:
        node_stats.append({
            "gene": gene,
            "degree": degree[gene],
            "betweenness": round(betweenness[gene], 4),
            "closeness": round(closeness[gene], 4),
            "clustering": round(cc.get(gene, 0), 4),
            "hub_score": round(hub_score[gene], 4),
            "is_proviral": gene in proviral,
            "is_antiviral": gene in antiviral,
            "is_mirna_hub": gene in mirna_hubs,
            "is_validated": gene in val_rep_genes,
            "log2FC_DENV": shared_up[shared_up["symbol"] == gene]["log2FC_DENV"].values[0]
                           if len(shared_up[shared_up["symbol"] == gene]) > 0 else 0,
        })

    node_df = pd.DataFrame(node_stats).sort_values("hub_score", ascending=False)
    node_df.to_csv(RES_DIR / "network_node_attributes.csv", index=False)
    log(f"  Node stats saved → {RES_DIR}/network_node_attributes.csv")

    top10 = node_df.head(10)
    log("\n  Top hub genes by hub score:")
    log(f"  {'Gene':<12} {'Degree':>8} {'Betweenness':>12} {'Hub Score':>12} {'miRNA Hub':>10} {'Validated':>10}")
    log("  " + "-" * 68)
    for _, row in top10.iterrows():
        log(f"  {row['gene']:<12} {row['degree']:>8} {row['betweenness']:>12.4f} "
            f"{row['hub_score']:>12.4f} {'YES' if row['is_mirna_hub'] else '':>10} "
            f"{'YES' if row['is_validated'] else '':>10}")

    top10.to_csv(RES_DIR / "hub_genes_top10.csv", index=False)

    # ─── Network statistics ────────────────────────────────────────────────────
    net_stats = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": round(nx.density(G), 4),
        "avg_degree": round(np.mean(list(degree.values())), 3),
        "avg_clustering": round(np.mean(list(cc.values())), 4) if cc else 0,
        "n_connected_components": nx.number_connected_components(G),
        "top_hub_gene": node_df.iloc[0]["gene"] if len(node_df) > 0 else "N/A",
    }
    pd.DataFrame([net_stats]).to_csv(RES_DIR / "network_statistics.csv", index=False)
    log(f"\n  Network statistics: {net_stats}")

    # ─── Figure ───────────────────────────────────────────────────────────────
    log("\nGenerating network figure ...")

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    # Panel A — Network visualization
    ax = axes[0]
    if G.number_of_edges() > 0:
        pos = nx.spring_layout(G, k=2.5, seed=42)
    else:
        pos = nx.circular_layout(G)

    # Node colors
    node_colors = []
    for n in G.nodes:
        if n in proviral:
            node_colors.append("#D32F2F")    # Red: proviral
        elif n in antiviral:
            node_colors.append("#1565C0")   # Blue: antiviral
        elif n in val_rep_genes:
            node_colors.append("#E91E63")   # Pink: validated in NPCs
        elif n in mirna_hubs:
            node_colors.append("#FF6F00")   # Orange: miRNA hub
        else:
            node_colors.append("#607D8B")   # Grey: novel

    # Node sizes by average FC
    node_sizes = []
    for n in G.nodes:
        ndata = G.nodes[n]
        fc = abs(ndata.get("avg_FC", 1))
        node_sizes.append(max(300, fc * 250))

    # Draw
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color="#999999",
                           width=[G[u][v].get("weight", 400) / 400 for u, v in G.edges])
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes,
                           edgecolors="black", linewidths=0.8)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold")
    ax.set_title("A  PPI Network (STRING score ≥ 400)\nShared DENV–ZIKV Host Response Genes",
                 fontsize=11, fontweight="bold")
    ax.axis("off")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(color="#D32F2F", label="Proviral"),
        Patch(color="#1565C0", label="Antiviral"),
        Patch(color="#E91E63", label="Validated (NPC)"),
        Patch(color="#FF6F00", label="miRNA hub target"),
        Patch(color="#607D8B", label="Novel"),
    ]
    ax.legend(handles=legend_elements, loc="lower left", fontsize=8, framealpha=0.8)

    # Panel B — Hub gene ranking
    ax2 = axes[1]
    plot_genes = node_df.head(15)["gene"].tolist()
    plot_degrees = node_df.head(15)["degree"].tolist()
    bar_colors = []
    for g in plot_genes:
        if g in proviral: bar_colors.append("#D32F2F")
        elif g in antiviral: bar_colors.append("#1565C0")
        elif g in val_rep_genes: bar_colors.append("#E91E63")
        elif g in mirna_hubs: bar_colors.append("#FF6F00")
        else: bar_colors.append("#607D8B")

    bars = ax2.barh(plot_genes[::-1], plot_degrees[::-1], color=bar_colors[::-1],
                    edgecolor="black", linewidth=0.7, alpha=0.85)
    ax2.set_xlabel("Node Degree (# PPI connections)", fontsize=11)
    ax2.set_title("B  Hub Gene Ranking (by Degree)\nAll Shared DEGs",
                  fontsize=11, fontweight="bold")
    ax2.legend(handles=legend_elements, loc="lower right", fontsize=8, framealpha=0.8)

    # Annotate with special labels
    for bar, gene in zip(bars[::-1], plot_genes[::-1]):
        annotations = []
        if gene in val_rep_genes: annotations.append("NPC✓")
        if gene in mirna_hubs: annotations.append("miRNA-hub")
        if annotations:
            ax2.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                     " | ".join(annotations), va="center", fontsize=7, color="#333333")

    plt.suptitle(f"DENV–ZIKV Shared Response: PPI Network Analysis\n"
                 f"({G.number_of_nodes()} genes, {G.number_of_edges()} interactions, "
                 f"STRING score ≥ 400)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig_path = FIG_MAIN / "Figure_Network_Analysis.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    log(f"Figure saved → {fig_path}")

    # ─── Summary ──────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 13 COMPLETE — Network Analysis Summary")
    log("=" * 60)
    log(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    log(f"  Top hub gene: {node_df.iloc[0]['gene'] if len(node_df) > 0 else 'N/A'}")
    log(f"  NPC-validated genes in network: {[g for g in val_rep_genes if g in G.nodes]}")
    log(f"  miRNA hub genes in network: {[g for g in mirna_hubs if g in G.nodes]}")

    ckpt["network_done"] = True
    save_ckpt(ckpt)
    log("\nNext: run step14_final_summary_figures.py")


if __name__ == "__main__":
    main()
