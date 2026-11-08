from pathlib import Path
import importlib.metadata as md
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ---------------------------------------------------------------- 0. Parameters
COUNTS_PATH = Path("data/raw/counts.csv.gz")
METADATA_PATH = Path("data/raw/metadata.csv.gz")
RESULTS = Path("results")
LOGS = Path("logs")
CONDITION_COLUMN = "site"
PLATE_COLUMN = "plate"
CONDITION_TEST = "Primary_Solid_Tumor"
CONDITION_REF = "Solid_Tissue_Normal"
MIN_COUNT = 10
MIN_SAMPLES = 20
ALPHA = 0.05
LFC_THRESHOLD = 1

RED, BLUE, GREY = "red", "blue", "grey"

RESULTS.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)



# ------------------------------------------------------------- [1/8] Read Data
print("[1/8] Reading input data (TCGA-BRCA)...")
counts = pd.read_csv(COUNTS_PATH, index_col=0)
metadata = pd.read_csv(METADATA_PATH, index_col=0)
n_genes_input = counts.shape[0]
print("  counts:", counts.shape)

# ------------------------------------------------------------ [2/8] Align Samples
print("[2/8] Aligning samples between counts and metadata...")
metadata = metadata.loc[counts.columns].copy()
assert (counts.columns == metadata.index).all()

metadata[PLATE_COLUMN] = metadata[PLATE_COLUMN].astype(str).astype("category")
metadata[CONDITION_COLUMN] = pd.Categorical(
    metadata[CONDITION_COLUMN], categories=[CONDITION_REF, CONDITION_TEST]
)

# ------------------------------------------------- [3/8] DESeq2 Object and Filtering
print("[3/8] Creating DESeq2 object and filtering low count genes...")
keep = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
counts = counts.loc[keep]
print("  genes after filter:", counts.shape[0])

# NOTE (R vs Python): unlike R/DESeq2 where genes are rows and samples are columns, 
# PyDESeq2 expects samples on rows and genes on columns, so I have to transpose counts here.
counts_t = counts.T.astype(int)

dds = DeseqDataSet(
    counts=counts_t,
    metadata=metadata,
    design=f"~ {CONDITION_COLUMN} + {PLATE_COLUMN}",
)


# ------------------------------------------------------------ [4/8] Run DESeq
print("[4/8] Running DESeq...")
dds.deseq2()

try:
    dds.plot_dispersions(save_path=str(RESULTS / "plot_dispersion.png"))
except Exception as e:
    print("  Could not plot dispersion:", e)
plt.close("all")

# ------------------------------------------------- [5/8] Results, Tables, and Plots
print("[5/8] Extracting results and saving tables...")
stat = DeseqStats(dds, contrast=[CONDITION_COLUMN, CONDITION_TEST, CONDITION_REF], alpha=ALPHA)
stat.summary()
res_df = stat.results_df.copy()
res_df["gene_id"] = res_df.index
res_df = res_df.sort_values("padj", na_position="last").reset_index(drop=True)

sig_df = res_df[res_df["padj"].notna() & (res_df["padj"] < ALPHA)]
de_df = res_df[
    res_df["padj"].notna()
    & (res_df["padj"] < ALPHA)
    & (res_df["log2FoldChange"].abs() > LFC_THRESHOLD)
].copy()
up_df = de_df[de_df["log2FoldChange"] > 0]
down_df = de_df[de_df["log2FoldChange"] < 0]

top10_up = up_df.head(10)
top10_down = down_df.head(10)

de_df["top10_label"] = "other"
de_df.loc[de_df["gene_id"].isin(top10_up["gene_id"]), "top10_label"] = "top10_up"
de_df.loc[de_df["gene_id"].isin(top10_down["gene_id"]), "top10_label"] = "top10_down"

res_df.to_csv(RESULTS / "deseq2_all_results.tsv", sep="\t", index=False, na_rep="NA")
sig_df.to_csv(RESULTS / "deseq2_sig_padj_0.05.tsv", sep="\t", index=False, na_rep="NA")
de_df.to_csv(RESULTS / "deseq2_de_padj_0.05_lfc_1.tsv", sep="\t", index=False, na_rep="NA")
top10_up.to_csv(RESULTS / "top10_up_genes.tsv", sep="\t", index=False, na_rep="NA")
top10_down.to_csv(RESULTS / "top10_down_genes.tsv", sep="\t", index=False, na_rep="NA")

print(f"  significant: {len(sig_df)} | DE: {len(de_df)} | up: {len(up_df)} | down: {len(down_df)}")

top10 = pd.concat([top10_up, top10_down]).sort_values("log2FoldChange")
fig, ax = plt.subplots(figsize=(7, 6))
ax.barh(top10["gene_id"], top10["log2FoldChange"],
        color=[RED if v > 0 else BLUE for v in top10["log2FoldChange"]])
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("log2FoldChange")
fig.tight_layout()
fig.savefig(RESULTS / "plot_top10_genes.png", dpi=150)
plt.close(fig)

volcano = res_df[res_df["padj"].notna()].copy()
volcano["mlog10_padj"] = -np.log10(volcano["padj"].clip(lower=1e-300))
volcano["status"] = "Not significant"
volcano.loc[volcano["gene_id"].isin(up_df["gene_id"]), "status"] = "Up"
volcano.loc[volcano["gene_id"].isin(down_df["gene_id"]), "status"] = "Down"

fig, ax = plt.subplots(figsize=(7, 6))
for status, color in [("Not significant", GREY), ("Up", RED), ("Down", BLUE)]:
    sub = volcano[volcano["status"] == status]
    ax.scatter(sub["log2FoldChange"], sub["mlog10_padj"], s=5, c=color, label=status)
for x in (-LFC_THRESHOLD, LFC_THRESHOLD):
    ax.axvline(x, ls="--", c="black", lw=0.8)
ax.axhline(-np.log10(ALPHA), ls="--", c="black", lw=0.8)
for _, row in volcano.head(10).iterrows():
    ax.text(row["log2FoldChange"], row["mlog10_padj"], row["gene_id"], fontsize=7, ha="center", va="bottom")
ax.set_xlabel("log2FoldChange")
ax.set_ylabel("-log10(padj)")
ax.legend(title="status")
fig.tight_layout()
fig.savefig(RESULTS / "plot_volcano.png", dpi=150)
plt.close(fig)

# NOTE (R vs Python): unlike R's plotMA function which handles limits automatically, 
# here in matplotlib I had to manually clip values and use triangle markers for out-of-bounds points.
ma = res_df[res_df["baseMean"] > 0].copy()
ma["signif"] = ma["padj"].notna() & (ma["padj"] < ALPHA)
ma["y"] = ma["log2FoldChange"].clip(-4, 4)
ma["clipped"] = ma["log2FoldChange"].abs() > 4

fig, ax = plt.subplots(figsize=(8, 6))
for signif, color in [(False, GREY), (True, BLUE)]:
    for clipped, marker in [(False, "o"), (True, "^")]:
        sub = ma[(ma["signif"] == signif) & (ma["clipped"] == clipped)]
        ax.scatter(sub["baseMean"], sub["y"], s=4, c=color, marker=marker)
ax.axhline(0, c="black", lw=1)
ax.set_xscale("log")
ax.set_ylim(-4.3, 4.3)
ax.set_xlabel("mean of normalized counts")
ax.set_ylabel("log fold change")
ax.set_title("DESeq2 MA plot")
fig.tight_layout()
fig.savefig(RESULTS / "plot_ma.png", dpi=150)
plt.close(fig)

# ------------------------------------------------------- [6/8] Top Gene Context
print("[6/8] Biological contextualization: top gene in contrast...")
top_gene = de_df.iloc[0]
(RESULTS / "top_gene_context.txt").write_text(
    f"gene_id: {top_gene['gene_id']}\n"
    f"log2FoldChange: {round(top_gene['log2FoldChange'], 4)}\n"
    f"padj: {top_gene['padj']:.4g}\n"
)

# ------------------------------------------------------- [7/8] Enrichment Analysis
print("[7/8] GO and Reactome enrichment...")
# NOTE (R vs Python): unlike clusterProfiler in R which often requires conversion to ENTREZID, 
# gseapy can take gene symbols directly as input, which saved me some conversion headaches.
universe = de_df["gene_id"].unique().tolist()
GO_LIB = "GO_Biological_Process_2023"
REACTOME_LIB = "Reactome_Pathways_2024"


def run_enrichment(genes, library, out_name):
    genes = list(pd.unique(genes))
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=library,
        organism="human",
        background=universe,
        outdir=None,
    )
    table = enr.results
    table = table[table["Adjusted P-value"] < 0.05]
    table.to_csv(RESULTS / f"{out_name}.tsv", sep="\t", index=False, na_rep="NA")
    try:
        gp.dotplot(table, column="Adjusted P-value", title=out_name,
                   cutoff=0.05, top_term=10, ofname=str(RESULTS / f"plot_{out_name}.png"))
    except Exception as e:
        print(f"  No dotplot for {out_name}:", e)
    plt.close("all")


run_enrichment(up_df["gene_id"], GO_LIB, "go_bp_up")
run_enrichment(down_df["gene_id"], GO_LIB, "go_bp_down")
run_enrichment(up_df["gene_id"], REACTOME_LIB, "reactome_up")
run_enrichment(down_df["gene_id"], REACTOME_LIB, "reactome_down")

# ----------------------------------------------------------- [8/8] Session Info and Logs
print("[8/8] Saving session traceability...")
packages = ["pydeseq2", "gseapy", "pandas", "numpy", "matplotlib"]
session_lines = [f"Python {sys.version}"] + [f"{p}=={md.version(p)}" for p in packages]
(LOGS / "sessionInfo.txt").write_text("\n".join(session_lines) + "\n")

summary_lines = [
    f"n_samples: {counts.shape[1]}",
    f"n_genes_input: {n_genes_input}",
    f"n_genes_after_filter: {counts.shape[0]}",
    f"n_significant_padj_lt_{ALPHA}:{len(sig_df)}",
    f"n_de_absLog2FC_gt_{LFC_THRESHOLD}_padj_lt_{ALPHA}:{len(de_df)}",
    f"n_up: {len(up_df)}",
    f"n_down: {len(down_df)}",
    "top10_up_file: results/top10_up_genes.tsv",
    "top10_down_file: results/top10_down_genes.tsv",
    "top10_plot_file: results/plot_top10_genes.png",
    "volcano_plot_file: results/plot_volcano.png",
]
(LOGS / "run_summary.txt").write_text("\n".join(summary_lines) + "\n")

print("Analysis completed. Check results/ and logs/ folders.")