#!/usr/bin/env python3
"""
AMR gene profiling across Salmonella enterica genomes.

Reads per-genome NCBI AMRFinderPlus TSV outputs, combines them, and produces
summary tables plus publication-ready figures:
  1. AMR gene count per genome (barplot, MDR highlighted)
  2. AMR gene count per drug class (barplot)
  3. Presence/absence heatmap of genomes x drug classes

Usage: python3 amr_analysis.py   (run from the repo root)
Inputs:  results/amrfinder_raw/*.tsv, data/genome_metadata.csv
Outputs: results/amr_genes_combined.csv, results/per_genome_summary.csv,
         results/per_class_summary.csv, figures/*.png
"""

import glob
import os
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

REPO = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(REPO, "results", "amrfinder_raw")
FIG_DIR = os.path.join(REPO, "figures")
RES_DIR = os.path.join(REPO, "results")
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams.update({"figure.dpi": 150, "font.size": 10})

# ---------------------------------------------------------------- load data
frames = []
for path in sorted(glob.glob(os.path.join(RAW_DIR, "*.tsv"))):
    label = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, sep="\t")
    df["genome"] = label
    frames.append(df)

if not frames:
    raise SystemExit("No AMRFinderPlus TSV files found in results/amrfinder_raw/")

all_hits = pd.concat(frames, ignore_index=True)

# Keep resistance genes only (drop virulence/stress/biocide-only markers
# that AMRFinderPlus also reports, keeping the focus on AMR)
amr = all_hits[all_hits["Type"] == "AMR"].copy()
amr["Class"] = amr["Class"].fillna("Other")

meta = pd.read_csv(os.path.join(REPO, "data", "genome_metadata.csv"))
meta_map = dict(zip(meta["label"], meta["accession"]))
serovar_map = dict(zip(meta["label"], meta["serovar"]))
amr["accession"] = amr["genome"].map(meta_map)
amr["serovar"] = amr["genome"].map(serovar_map)

# ------------------------------------------------------- combined table
combined = amr[[
    "genome", "accession", "serovar",
    "Element symbol", "Element name", "Class", "Subclass",
    "Type", "Scope", "% Coverage of reference", "% Identity to reference",
]].rename(columns={
    "Element symbol": "gene_symbol",
    "Element name": "gene_name",
    "Type": "element_type",
    "Scope": "scope",
    "% Coverage of reference": "coverage_pct",
    "% Identity to reference": "identity_pct",
})
combined.to_csv(os.path.join(RES_DIR, "amr_genes_combined.csv"), index=False)

# ------------------------------------------------------- per-genome summary
# include every genome from the metadata, even those with zero AMR hits
all_labels = meta["label"].tolist()
per_genome = (
    amr.groupby("genome")
    .agg(n_genes=("Element symbol", "size"),
         n_classes=("Class", "nunique"))
    .reindex(all_labels)
    .fillna(0)
    .reset_index()
    .rename(columns={"index": "genome"})
)
per_genome["n_genes"] = per_genome["n_genes"].astype(int)
per_genome["n_classes"] = per_genome["n_classes"].astype(int)
per_genome["accession"] = per_genome["genome"].map(meta_map)
per_genome["serovar"] = per_genome["genome"].map(serovar_map)
per_genome["MDR"] = per_genome["n_classes"] >= 3  # resistant to 3+ drug classes
per_genome = per_genome.sort_values("n_genes", ascending=False)
per_genome.to_csv(os.path.join(RES_DIR, "per_genome_summary.csv"), index=False)

# ------------------------------------------------------- per-class summary
per_class = (
    amr.groupby("Class")
    .agg(n_genes=("Element symbol", "size"),
         n_genomes=("genome", "nunique"))
    .reset_index()
    .sort_values("n_genes", ascending=False)
)
per_class.to_csv(os.path.join(RES_DIR, "per_class_summary.csv"), index=False)

# ------------------------------------------------------- figure 1: genes per genome
fig, ax = plt.subplots(figsize=(10, 5))
order = per_genome["genome"]
colors = ["#c0392b" if m else "#3498db" for m in per_genome["MDR"]]
bars = ax.bar(order, per_genome["n_genes"], color=colors)
ax.set_ylabel("AMR genes detected")
ax.set_xlabel("Genome")
ax.set_title("AMR gene count per Salmonella genome (red = MDR, 3+ drug classes)")
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=30, ha="right")
for bar, val in zip(bars, per_genome["n_genes"]):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
            str(val), ha="center", va="bottom", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "01_amr_genes_per_genome.png"))
plt.close(fig)

# ------------------------------------------------------- figure 2: genes per drug class
fig, ax = plt.subplots(figsize=(10, 5))
bars = ax.bar(per_class["Class"], per_class["n_genes"], color="#2ecc71")
ax.set_ylabel("AMR gene hits (across all genomes)")
ax.set_xlabel("Drug class")
ax.set_title("AMR gene hits per drug class (10 Salmonella genomes)")
ax.set_xticks(range(len(per_class)))
ax.set_xticklabels(per_class["Class"], rotation=30, ha="right")
for bar, val in zip(bars, per_class["n_genes"]):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
            str(val), ha="center", va="bottom", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "02_genes_per_drug_class.png"))
plt.close(fig)

# ------------------------------------------------------- figure 3: presence/absence heatmap
pivot = (
    amr.drop_duplicates(["genome", "Class"])
    .assign(present=1)
    .pivot_table(index="genome", columns="Class", values="present", fill_value=0)
)
# order genomes by total genes, classes by total hits; keep zero-hit genomes
pivot = pivot.reindex(index=per_genome["genome"],
                      columns=per_class["Class"]).fillna(0)
fig, ax = plt.subplots(figsize=(12, 5))
sns.heatmap(pivot, cmap=["#ecf0f1", "#2980b9"], cbar=False,
            linewidths=0.5, linecolor="white", ax=ax)
ax.set_title("AMR drug-class presence/absence across Salmonella genomes")
ax.set_xlabel("Drug class")
ax.set_ylabel("Genome")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "03_presence_absence_heatmap.png"))
plt.close(fig)

# ------------------------------------------------------- console summary (for README)
n_genomes = per_genome.shape[0]
n_hits = len(amr)
n_mdr = int(per_genome["MDR"].sum())
top = per_genome.iloc[0]
print(f"GENOMES: {n_genomes}")
print(f"TOTAL_AMR_HITS: {n_hits}")
print(f"MDR_GENOMES: {n_mdr}")
print(f"MOST_RESISTANT: {top['genome']} ({top['accession']}, {top['n_genes']} genes, {top['n_classes']} classes)")
print("TOP_CLASSES:")
for _, r in per_class.head(3).iterrows():
    print(f"  {r['Class']}: {r['n_genes']} hits in {r['n_genomes']} genomes")
print("PER_GENOME:")
for _, r in per_genome.iterrows():
    print(f"  {r['genome']}: {r['n_genes']} genes, {r['n_classes']} classes, MDR={r['MDR']}")
print("Wrote figures to figures/ and tables to results/")
