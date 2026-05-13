#!/usr/bin/env bash
set -u

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT" || exit 1

mkdir -p scripts logs

LOG="logs/105_organize_final_package_and_functional_metabolic_figures.log"

cat > scripts/105_organize_final_package_and_functional_metabolic_figures.py <<'PY'
import os
import re
import json
import shutil
import math
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/mnt/f/Marcos/rossana")
os.chdir(ROOT)

FINAL = ROOT / "04_results" / "FINAL_SEPSIS_MANUSCRIPT_PACKAGE"
TABLES = FINAL / "tables"
FIGURES = FINAL / "figures"
REPORTS = FINAL / "reports"
SCRIPTS = FINAL / "scripts_used"
LOGS = FINAL / "logs_used"
SUPP = FINAL / "supplementary"
MANIFEST_DIR = FINAL / "manifest"

for d in [FINAL, TABLES, FIGURES, REPORTS, SCRIPTS, LOGS, SUPP, MANIFEST_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def msg(x):
    print(x, flush=True)


def exists_nonempty(path):
    p = Path(path)
    return p.exists() and p.is_file() and p.stat().st_size > 0


def read_csv_safe(path):
    path = Path(path)
    if not exists_nonempty(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        try:
            return pd.read_csv(path, sep=";")
        except Exception:
            return pd.DataFrame()


def clean_gene(x):
    if pd.isna(x):
        return ""
    x = str(x).strip().upper()
    x = re.sub(r"[^A-Z0-9]", "", x)
    return x


def find_col(df, patterns):
    if df is None or df.empty:
        return None
    cols = list(df.columns)
    for pat in patterns:
        rgx = re.compile(pat, re.I)
        for c in cols:
            if rgx.search(str(c)):
                return c
    return None


def detect_gene_col(df):
    return find_col(df, [
        r"^gene_symbol$",
        r"gene.*symbol",
        r"symbol",
        r"gene$",
        r"external_gene_name",
        r"preferredName",
        r"queryItem"
    ])


def detect_term_col(df):
    return find_col(df, [
        r"description",
        r"term",
        r"name",
        r"pathway",
        r"category",
        r"process"
    ])


def copy_nonempty(src, dst_dir, new_name=None):
    src = Path(src)
    if exists_nonempty(src):
        dst = Path(dst_dir) / (new_name if new_name else src.name)
        shutil.copy2(src, dst)
        return str(dst)
    return None


def save_table(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def finite_numeric(s):
    return pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan)


def wrap_labels(labels, width=45):
    out = []
    for lab in labels:
        lab = str(lab)
        words = lab.split()
        lines = []
        cur = ""
        for w in words:
            if len(cur) + len(w) + 1 <= width:
                cur = (cur + " " + w).strip()
            else:
                if cur:
                    lines.append(cur)
                cur = w
        if cur:
            lines.append(cur)
        out.append("\n".join(lines))
    return out


msg("============================================================")
msg("PIPELINE 105 - FINAL PACKAGE + FUNCTIONAL/METABOLIC FIGURES")
msg("============================================================")

paths = {
    "final_prioritization": "04_results/repaired_gse13205_gpl570/Table_103G_integrated_prioritization_GPL570_DE_AUC.csv",
    "de_auc_summary_candidate": "04_results/final_manuscript/tables/Table_Final_03_GSE13205_DE_AUC_summary.csv",
    "de": "04_results/repaired_gse13205_gpl570/Table_103E_GSE13205_differential_expression_GPL570.csv",
    "auc": "04_results/repaired_gse13205_gpl570/Table_103F_GSE13205_univariate_AUC_GPL570.csv",
    "probe_mapping": "04_results/repaired_gse13205_gpl570/Table_103C_GSE13205_GPL570_feature_matching_audit.csv",
    "go": "04_results/enrichment/Table_03_GO_terms_selected_genes.csv",
    "go_summary": "04_results/enrichment/Table_04_GO_ontology_summary.csv",
    "kegg": "04_results/enrichment/Table_05_KEGG_pathways_selected_genes.csv",
    "functional_proximity": "04_results/networks/Table_07_gene_gene_functional_proximity_GO_KEGG.csv",
    "gene_pathway_edges": "04_results/networks/Table_06_gene_pathway_network_edges.csv",
    "string_mapping": "04_results/repaired/string/Table_R01_STRING_mapping_repaired_api.csv",
    "string_interactions": "04_results/repaired/string/Table_R02_STRING_interactions_repaired_api.csv",
    "string_centrality": "04_results/repaired/string/Table_R03_STRING_centrality_repaired_api.csv",
    "expression_summary": "04_results/tables/Table_17_selected_gene_expression_summary.csv",
    "final_clean_xlsx": "04_results/final_summary/clean_final_gene_candidates_for_manuscript.xlsx",
    "final_summary_xlsx": "04_results/final_summary/final_sepsis_gene_prioritization_summary.xlsx",
}

data = {k: read_csv_safe(v) for k, v in paths.items()}

# ============================================================
# 1. COPIAR TABELAS ÚTEIS E NÃO VAZIAS
# ============================================================

copied_tables = []

table_candidates = {
    "Table_01_final_integrated_prioritization_GPL570_DE_AUC.csv": paths["final_prioritization"],
    "Table_02_GSE13205_differential_expression.csv": paths["de"],
    "Table_03_GSE13205_univariate_AUC.csv": paths["auc"],
    "Table_04_GSE13205_probe_mapping_GPL570.csv": paths["probe_mapping"],
    "Table_05_GO_terms_selected_genes.csv": paths["go"],
    "Table_06_GO_ontology_summary.csv": paths["go_summary"],
    "Table_07_KEGG_pathways_selected_genes.csv": paths["kegg"],
    "Table_08_gene_gene_functional_proximity_GO_KEGG.csv": paths["functional_proximity"],
    "Table_09_gene_pathway_network_edges.csv": paths["gene_pathway_edges"],
    "Table_10_STRING_mapping_repaired.csv": paths["string_mapping"],
    "Table_11_STRING_interactions_repaired.csv": paths["string_interactions"],
    "Table_12_STRING_centrality_repaired.csv": paths["string_centrality"],
    "Table_13_expression_summary_selected_genes.csv": paths["expression_summary"],
}

for new_name, src in table_candidates.items():
    df = read_csv_safe(src)
    if not df.empty:
        out = copy_nonempty(src, TABLES, new_name)
        if out:
            copied_tables.append(out)

for k in ["final_clean_xlsx", "final_summary_xlsx"]:
    out = copy_nonempty(paths[k], TABLES)
    if out:
        copied_tables.append(out)

# ============================================================
# 2. COPIAR FIGURAS FINAIS EXISTENTES, SE ESTIVEREM BOAS
# ============================================================

copied_figures = []

figure_sources = [
    "04_results/final_manuscript/figures/Figure_Final_01_integrated_gene_prioritization.png",
    "04_results/final_manuscript/figures/Figure_Final_02_GSE13205_AUC_FDR_bubble.png",
    "04_results/final_manuscript/figures/Figure_Final_03_GSE13205_logFC_direction.png",
    "04_results/final_manuscript/figures/Figure_Final_04_top_gene_expression_boxplots.png",
    "04_results/final_manuscript/figures/Figure_Final_05_STRING_network_centrality.png",
    "04_results/final_manuscript/figures/Figure_Final_06_evidence_matrix_top_candidates.png",
    "04_results/figures/Figure_06_top_functional_proximity_pairs_GO_KEGG.png",
    "04_results/figures/Figure_05_integrated_gene_candidate_prioritization.png",
]

for src in figure_sources:
    out = copy_nonempty(src, FIGURES)
    if out:
        copied_figures.append(out)

# ============================================================
# 3. PADRONIZAR TABELA MASTER
# ============================================================

prior = data["final_prioritization"].copy()

if prior.empty:
    msg("[ERRO] Tabela final de priorização não encontrada.")
    msg("Esperado: 04_results/repaired_gse13205_gpl570/Table_103G_integrated_prioritization_GPL570_DE_AUC.csv")
    raise SystemExit(1)

gene_col = detect_gene_col(prior)
if gene_col is None:
    msg("[ERRO] Não encontrei coluna de gene na priorização.")
    raise SystemExit(1)

prior["gene_clean"] = prior[gene_col].map(clean_gene)
prior = prior[prior["gene_clean"] != ""].copy()
prior = prior.drop_duplicates(subset=["gene_clean"], keep="first").copy()

score_col = find_col(prior, [r"integrated_score_final", r"integrated.*score", r"score"])
auc_col = find_col(prior, [r"auc.*directional", r"best_auc", r"auc"])
padj_col = find_col(prior, [r"adj.*p", r"fdr", r"BH"])
direction_col = find_col(prior, [r"direction"])
logfc_col = find_col(prior, [r"logfc", r"logFC", r"mean.*diff"])

if score_col:
    prior["score_numeric"] = finite_numeric(prior[score_col])
else:
    prior["score_numeric"] = 0

if auc_col:
    prior["auc_numeric"] = finite_numeric(prior[auc_col])
else:
    prior["auc_numeric"] = np.nan

if padj_col:
    prior["padj_numeric"] = finite_numeric(prior[padj_col])
else:
    prior["padj_numeric"] = np.nan

if logfc_col:
    prior["logfc_numeric"] = finite_numeric(prior[logfc_col])
else:
    prior["logfc_numeric"] = np.nan

sort_cols = ["score_numeric"]
ascending = [False]
if "auc_numeric" in prior.columns:
    sort_cols.append("auc_numeric")
    ascending.append(False)
if "padj_numeric" in prior.columns:
    sort_cols.append("padj_numeric")
    ascending.append(True)

prior = prior.sort_values(sort_cols, ascending=ascending, na_position="last")

master_cols = []
for c in [
    gene_col,
    "evidence_group",
    "analysis_priority",
    score_col,
    auc_col,
    padj_col,
    direction_col,
    logfc_col,
    "has_GO_annotation",
    "has_KEGG_pathway",
    "has_STRING_interaction_repaired",
    "has_GSE13205_GPL570_probe_mapping",
    "has_GSE13205_GPL570_DE",
    "has_GSE13205_GPL570_AUC"
]:
    if c and c in prior.columns and c not in master_cols:
        master_cols.append(c)

master = prior[master_cols].copy()
master.insert(0, "rank_final", range(1, len(master) + 1))
master.rename(columns={gene_col: "gene_symbol"}, inplace=True)

save_table(master, TABLES / "Table_14_FINAL_clean_ranked_gene_candidates_for_manuscript.csv")

# ============================================================
# 4. CRIAR TABELAS FUNCIONAIS/METABÓLICAS
# ============================================================

go = data["go"].copy()
kegg = data["kegg"].copy()
prox = data["functional_proximity"].copy()
string_cent = data["string_centrality"].copy()
string_edges = data["string_interactions"].copy()
de = data["de"].copy()
auc = data["auc"].copy()

# ---- KEGG gene-pathway edge table
kegg_edges = pd.DataFrame()
if not kegg.empty:
    kgene = detect_gene_col(kegg)
    kterm = detect_term_col(kegg)
    if kgene and kterm:
        kegg_edges = kegg[[kgene, kterm]].copy()
        kegg_edges.columns = ["gene_symbol", "kegg_pathway"]
        kegg_edges["gene_symbol"] = kegg_edges["gene_symbol"].map(clean_gene)
        kegg_edges["kegg_pathway"] = kegg_edges["kegg_pathway"].astype(str)
        kegg_edges = kegg_edges[(kegg_edges["gene_symbol"] != "") & (kegg_edges["kegg_pathway"] != "")]
        kegg_edges = kegg_edges.drop_duplicates()

save_table(kegg_edges, TABLES / "Table_15_functional_KEGG_gene_pathway_edges.csv")

# ---- GO gene-term edge table
go_edges = pd.DataFrame()
if not go.empty:
    ggene = detect_gene_col(go)
    gterm = detect_term_col(go)
    ont_col = find_col(go, [r"ontology", r"ont", r"category"])
    if ggene and gterm:
        keep = [ggene, gterm]
        if ont_col:
            keep.append(ont_col)
        go_edges = go[keep].copy()
        if ont_col:
            go_edges.columns = ["gene_symbol", "go_term", "ontology"]
        else:
            go_edges.columns = ["gene_symbol", "go_term"]
            go_edges["ontology"] = ""
        go_edges["gene_symbol"] = go_edges["gene_symbol"].map(clean_gene)
        go_edges["go_term"] = go_edges["go_term"].astype(str)
        go_edges = go_edges[(go_edges["gene_symbol"] != "") & (go_edges["go_term"] != "")]
        go_edges = go_edges.drop_duplicates()

save_table(go_edges, TABLES / "Table_16_functional_GO_gene_term_edges.csv")

# ---- Functional axes by keyword
axes = {
    "Immune/inflammatory response": [
        "immune", "inflamm", "cytokine", "interleukin", "defense", "response to bacter",
        "response to lipopolysaccharide", "toll", "nf-kappa", "leukocyte"
    ],
    "Leukocyte/neutrophil activation and migration": [
        "neutrophil", "granulocyte", "leukocyte migration", "chemotaxis", "myeloid",
        "phagocyt", "degranulation"
    ],
    "Acute-phase and antimicrobial response": [
        "acute-phase", "acute phase", "antimicrobial", "defensin", "bactericidal",
        "complement", "serum amyloid"
    ],
    "Oxidative/metabolic stress": [
        "oxidative", "glutathione", "reactive oxygen", "redox", "metabolic",
        "metabolism", "oxidation", "mitochond", "response to stress"
    ],
    "Tissue remodeling / extracellular matrix": [
        "extracellular matrix", "collagen", "metalloproteinase", "remodel",
        "wound", "cell adhesion", "matrix organization"
    ],
    "Vascular / endothelial / permeability axis": [
        "vascular", "endothelial", "angiogenesis", "permeability", "blood vessel",
        "vasculature"
    ],
    "Signal transduction / receptor regulation": [
        "receptor", "signaling", "signal transduction", "kinase", "phosphorylation",
        "jak", "stat", "mapk"
    ],
}

axis_records = []

def add_axis_records(edge_df, gene_col_name, term_col_name, source):
    if edge_df.empty:
        return
    for _, row in edge_df.iterrows():
        gene = clean_gene(row.get(gene_col_name, ""))
        term = str(row.get(term_col_name, ""))
        low = term.lower()
        for axis, kws in axes.items():
            hit = any(kw.lower() in low for kw in kws)
            if hit:
                axis_records.append({
                    "gene_symbol": gene,
                    "functional_axis": axis,
                    "source": source,
                    "matched_term": term
                })

add_axis_records(go_edges, "gene_symbol", "go_term", "GO")
add_axis_records(kegg_edges, "gene_symbol", "kegg_pathway", "KEGG")

axis_df = pd.DataFrame(axis_records)
if not axis_df.empty:
    axis_df = axis_df[axis_df["gene_symbol"] != ""].drop_duplicates()
    axis_summary = (
        axis_df.groupby(["gene_symbol", "functional_axis"])
        .agg(n_supporting_terms=("matched_term", "nunique"),
             sources=("source", lambda x: "; ".join(sorted(set(x)))))
        .reset_index()
        .sort_values(["functional_axis", "n_supporting_terms"], ascending=[True, False])
    )
else:
    axis_summary = pd.DataFrame(columns=["gene_symbol", "functional_axis", "n_supporting_terms", "sources"])

save_table(axis_df, TABLES / "Table_17_functional_axis_gene_term_long.csv")
save_table(axis_summary, TABLES / "Table_18_functional_axis_gene_summary.csv")

# ---- Metabolic bridge
metabolic_keywords = [
    "metabolic", "metabolism", "glutathione", "oxidative", "redox", "lipid",
    "arachidonic", "tryptophan", "amino acid", "fatty acid", "peroxisome",
    "carbon", "glycolysis", "mitochond", "oxidation", "biosynthesis"
]

metabolic_records = []
for source_name, edge_df, term_col in [
    ("GO", go_edges, "go_term"),
    ("KEGG", kegg_edges, "kegg_pathway")
]:
    if edge_df.empty:
        continue
    for _, row in edge_df.iterrows():
        term = str(row.get(term_col, ""))
        low = term.lower()
        if any(k in low for k in metabolic_keywords):
            metabolic_records.append({
                "gene_symbol": clean_gene(row.get("gene_symbol", "")),
                "source": source_name,
                "metabolic_or_stress_term": term
            })

metabolic_df = pd.DataFrame(metabolic_records)
if not metabolic_df.empty:
    metabolic_df = metabolic_df[metabolic_df["gene_symbol"] != ""].drop_duplicates()
    metabolic_summary = (
        metabolic_df.groupby("gene_symbol")
        .agg(
            n_metabolic_stress_terms=("metabolic_or_stress_term", "nunique"),
            sources=("source", lambda x: "; ".join(sorted(set(x)))),
            example_terms=("metabolic_or_stress_term", lambda x: " | ".join(list(dict.fromkeys(map(str, x)))[:5]))
        )
        .reset_index()
        .sort_values("n_metabolic_stress_terms", ascending=False)
    )
else:
    metabolic_summary = pd.DataFrame(columns=["gene_symbol", "n_metabolic_stress_terms", "sources", "example_terms"])

save_table(metabolic_df, TABLES / "Table_19_metabolic_stress_gene_term_long.csv")
save_table(metabolic_summary, TABLES / "Table_20_metabolic_stress_gene_summary.csv")

# ============================================================
# 5. FIGURAS NOVAS FUNCIONAIS/METABÓLICAS
# ============================================================

def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=400, bbox_inches="tight")
    plt.close()
    copied_figures.append(str(path))

# Figure 07: KEGG pathway landscape
if not kegg_edges.empty:
    counts = (
        kegg_edges.groupby("kegg_pathway")
        .agg(n_genes=("gene_symbol", "nunique"))
        .reset_index()
        .sort_values("n_genes", ascending=False)
        .head(15)
    )
    if not counts.empty:
        plt.figure(figsize=(10, max(4, 0.45 * len(counts))))
        y = np.arange(len(counts))
        plt.barh(y, counts["n_genes"])
        plt.yticks(y, wrap_labels(counts["kegg_pathway"], 48))
        plt.xlabel("Number of selected genes")
        plt.ylabel("KEGG pathway")
        plt.title("Functional pathway landscape of selected sepsis-related genes")
        plt.gca().invert_yaxis()
        save_fig(FIGURES / "Figure_07_KEGG_pathway_landscape_top_pathways.png")

# Figure 08: Gene x KEGG pathway matrix
if not kegg_edges.empty:
    top_pathways = (
        kegg_edges.groupby("kegg_pathway")["gene_symbol"]
        .nunique()
        .sort_values(ascending=False)
        .head(12)
        .index
        .tolist()
    )
    top_genes = master["gene_symbol"].map(clean_gene).head(20).tolist()
    mat_df = kegg_edges[
        kegg_edges["kegg_pathway"].isin(top_pathways)
        & kegg_edges["gene_symbol"].isin(top_genes)
    ].copy()
    if not mat_df.empty:
        mat = pd.crosstab(mat_df["gene_symbol"], mat_df["kegg_pathway"])
        mat = mat.reindex(index=[g for g in top_genes if g in mat.index])
        plt.figure(figsize=(12, max(4, 0.35 * mat.shape[0])))
        plt.imshow(mat.values, aspect="auto")
        plt.xticks(np.arange(mat.shape[1]), wrap_labels(mat.columns, 25), rotation=70, ha="right")
        plt.yticks(np.arange(mat.shape[0]), mat.index)
        plt.xlabel("KEGG pathway")
        plt.ylabel("Gene")
        plt.title("Gene–KEGG pathway evidence matrix")
        cbar = plt.colorbar()
        cbar.set_label("Evidence present")
        save_fig(FIGURES / "Figure_08_gene_KEGG_pathway_matrix.png")

# Figure 09: Functional axis matrix
if not axis_summary.empty:
    top_genes_axis = master["gene_symbol"].map(clean_gene).head(25).tolist()
    axis_plot = axis_summary[axis_summary["gene_symbol"].isin(top_genes_axis)].copy()
    if not axis_plot.empty:
        mat = axis_plot.pivot_table(
            index="gene_symbol",
            columns="functional_axis",
            values="n_supporting_terms",
            aggfunc="sum",
            fill_value=0
        )
        mat = mat.reindex(index=[g for g in top_genes_axis if g in mat.index])
        plt.figure(figsize=(13, max(4, 0.35 * mat.shape[0])))
        plt.imshow(mat.values, aspect="auto")
        plt.xticks(np.arange(mat.shape[1]), wrap_labels(mat.columns, 26), rotation=60, ha="right")
        plt.yticks(np.arange(mat.shape[0]), mat.index)
        plt.xlabel("Functional biological axis")
        plt.ylabel("Gene")
        plt.title("Functional axis map integrating GO and KEGG evidence")
        cbar = plt.colorbar()
        cbar.set_label("Number of supporting GO/KEGG terms")
        save_fig(FIGURES / "Figure_09_functional_axis_gene_matrix_GO_KEGG.png")

# Figure 10: Metabolic/oxidative stress bridge
if not metabolic_summary.empty:
    met = metabolic_summary.head(15).copy()
    plt.figure(figsize=(9, max(4, 0.45 * len(met))))
    y = np.arange(len(met))
    plt.barh(y, met["n_metabolic_stress_terms"])
    plt.yticks(y, met["gene_symbol"])
    plt.xlabel("Number of metabolism/oxidative-stress related terms")
    plt.ylabel("Gene")
    plt.title("Metabolic and oxidative-stress bridge among selected candidates")
    plt.gca().invert_yaxis()
    save_fig(FIGURES / "Figure_10_metabolic_oxidative_stress_gene_bridge.png")

# Figure 11: Integrated score vs AUC
if "score_numeric" in prior.columns and "auc_numeric" in prior.columns:
    scat = prior.dropna(subset=["score_numeric", "auc_numeric"]).copy()
    scat = scat[scat["gene_clean"] != ""]
    if not scat.empty:
        plt.figure(figsize=(8.5, 6))
        plt.scatter(scat["score_numeric"], scat["auc_numeric"], s=60)
        for _, row in scat.sort_values(["score_numeric", "auc_numeric"], ascending=False).head(12).iterrows():
            plt.text(row["score_numeric"], row["auc_numeric"], row["gene_clean"], fontsize=8)
        plt.xlabel("Final integrated evidence score")
        plt.ylabel("Directional AUC in GSE13205")
        plt.title("Functional prioritization versus transcriptomic discrimination")
        save_fig(FIGURES / "Figure_11_integrated_score_vs_GSE13205_AUC.png")

# Figure 12: Directional DE among strongest candidates
if not prior.empty and "logfc_numeric" in prior.columns:
    deplot = prior.dropna(subset=["logfc_numeric"]).copy()
    deplot = deplot[deplot["gene_clean"] != ""].copy()
    deplot["abs_logfc"] = deplot["logfc_numeric"].abs()
    deplot = deplot.sort_values("abs_logfc", ascending=False).head(20)
    if not deplot.empty:
        plt.figure(figsize=(9, max(4, 0.42 * len(deplot))))
        y = np.arange(len(deplot))
        plt.barh(y, deplot["logfc_numeric"])
        plt.axvline(0, linewidth=0.8)
        plt.yticks(y, deplot["gene_clean"])
        plt.xlabel("Mean expression difference: sepsis minus control")
        plt.ylabel("Gene")
        plt.title("Direction and magnitude of expression differences in GSE13205")
        plt.gca().invert_yaxis()
        save_fig(FIGURES / "Figure_12_GSE13205_expression_direction_top_candidates.png")

# ============================================================
# 6. EXCEL CONSOLIDADO
# ============================================================

xlsx_path = TABLES / "FINAL_SEPSIS_MANUSCRIPT_TABLES_COMPLETE.xlsx"

sheets = {
    "final_ranked_candidates": master,
    "functional_axis_summary": axis_summary,
    "metabolic_stress_summary": metabolic_summary,
    "KEGG_gene_pathway_edges": kegg_edges,
    "GO_gene_term_edges": go_edges,
    "STRING_centrality": string_cent,
    "STRING_interactions": string_edges,
    "GSE13205_DE": de,
    "GSE13205_AUC": auc,
    "functional_proximity": prox,
}

with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    for sheet, df in sheets.items():
        if df is not None and not df.empty:
            df.to_excel(writer, sheet_name=sheet[:31], index=False)

copied_tables.append(str(xlsx_path))

# ============================================================
# 7. COPIAR SCRIPTS E LOGS ÚTEIS
# ============================================================

script_patterns = [
    "scripts/99_run_post_gene_sepsis_pipeline.sh",
    "scripts/100_audit_post_gene_sepsis_outputs.sh",
    "scripts/101_repair_string_auc_differential_sepsis.sh",
    "scripts/102_fix_GSE13205_phenotype_and_rerun_DE_AUC.sh",
    "scripts/103*",
    "scripts/104_make_final_figures_tables_and_results_text.sh",
    "scripts/105_organize_final_package_and_functional_metabolic_figures.sh",
    "scripts/105_organize_final_package_and_functional_metabolic_figures.py",
]

for pat in script_patterns:
    for src in ROOT.glob(pat):
        if exists_nonempty(src):
            copy_nonempty(src, SCRIPTS)

log_patterns = [
    "logs/02_post_gene_sepsis_functional_clinical_pipeline.log",
    "logs/101_repair_string_auc_differential_sepsis.log",
    "logs/102*",
    "logs/103*",
    "logs/104_make_final_figures_tables_and_results_text.log",
    "logs/105_organize_final_package_and_functional_metabolic_figures.log",
]

for pat in log_patterns:
    for src in ROOT.glob(pat):
        if exists_nonempty(src):
            copy_nonempty(src, LOGS)

# ============================================================
# 8. RELATÓRIOS E LEGENDAS
# ============================================================

n_total = len(master)
n_go = len(go_edges)
n_kegg = len(kegg_edges)
n_axis = len(axis_summary)
n_met = len(metabolic_summary)
n_figs = len(list(FIGURES.glob("*.png")))
n_tables = len(list(TABLES.glob("*")))

top_integrated = master["gene_symbol"].head(10).astype(str).tolist()
top_metabolic = metabolic_summary["gene_symbol"].head(10).astype(str).tolist() if not metabolic_summary.empty else []

report = f"""============================================================
FINAL SEPSIS MANUSCRIPT PACKAGE - EXECUTIVE REPORT
============================================================

Generated package:
{FINAL}

This package contains only non-empty and interpretable outputs from the post-selection sepsis gene prioritization workflow.

Main retained evidence layers:
1. Final integrated gene prioritization after GPL570/GSE13205 repair.
2. GO functional annotation.
3. KEGG pathway annotation.
4. Gene-gene functional proximity based on GO/KEGG.
5. STRING/PPI mapping, interactions and centrality.
6. Corrected GSE13205 transcriptomic validation.
7. Differential expression and univariate directional AUC.
8. Functional/metabolic axis analysis generated in this pipeline.

Package summary:
- Clean ranked genes: {n_total}
- GO gene-term edges: {n_go}
- KEGG gene-pathway edges: {n_kegg}
- Functional axis gene summaries: {n_axis}
- Metabolic/oxidative-stress linked genes: {n_met}
- Tables/files in final table folder: {n_tables}
- PNG figures in final figure folder: {n_figs}

Top integrated candidates:
{", ".join(top_integrated)}

Top metabolism/oxidative-stress linked candidates:
{", ".join(top_metabolic) if top_metabolic else "No metabolism/oxidative-stress keyword-linked candidates detected from available GO/KEGG terms."}

Scientific interpretation:
The final analytical package supports a biologically coherent prioritization of sepsis-related genes. The strongest evidence comes from the convergence of GO/KEGG functional annotation, STRING/PPI recovery, corrected GPL570 probe mapping, GSE13205 differential expression and directional AUC analysis. The new functional/metabolic figures should be used as mechanistic support, not as direct metabolomics evidence, because they are inferred from gene-pathway and gene-term annotations rather than measured metabolite concentrations.

Important caution:
GSE13205 is based on vastus lateralis muscle biopsies. Therefore, transcriptomic discrimination should be described as tissue-specific exploratory validation, not as direct blood-based diagnostic validation.

Recommended figures for manuscript:
Figure 1. Integrated prioritization of sepsis-related gene candidates.
Figure 2. Transcriptomic discrimination of selected genes in GSE13205.
Figure 3. Direction of expression differences between sepsis and control samples.
Figure 4. Evidence matrix for top prioritized candidates.
Figure 5. STRING-based network centrality of selected candidates.
Figure 6. Functional axis map integrating GO and KEGG evidence.
Figure 7. Metabolic and oxidative-stress bridge among selected candidates.
Figure 8. Gene–KEGG pathway evidence matrix.

Main folders:
Tables:
{TABLES}

Figures:
{FIGURES}

Reports:
{REPORTS}

Scripts used:
{SCRIPTS}
"""

(REPORTS / "FINAL_PACKAGE_EXECUTIVE_REPORT.txt").write_text(report, encoding="utf-8")

legends = """Figure legends for the final manuscript package

Figure 1. Integrated prioritization of sepsis-related gene candidates.
Bar plot showing top-ranked genes according to the final integrated evidence score. The score summarizes support from functional annotation, KEGG pathway membership, STRING/PPI recovery, GPL570 probe mapping, differential expression and univariate AUC analysis.

Figure 2. Transcriptomic discrimination of selected genes in GSE13205.
Bubble plot showing directional AUC values for selected genes after corrected phenotype mapping in GSE13205. Bubble size represents FDR-adjusted statistical support when available.

Figure 3. Direction of expression differences between sepsis and control samples.
Bar plot showing mean expression differences calculated as septic minus control expression. Positive values indicate higher expression in sepsis, whereas negative values indicate lower expression in sepsis.

Figure 4. Evidence matrix for top prioritized candidates.
Tile plot summarizing the presence or absence of each evidence layer for top-ranked genes.

Figure 5. STRING-based network centrality of selected candidates.
Bar plot showing degree centrality values derived from recovered STRING protein-protein interaction data.

Figure 6. Functional axis map integrating GO and KEGG evidence.
Heatmap showing the number of supporting GO/KEGG terms connecting each gene to major functional axes, including immune/inflammatory response, leukocyte activation, acute-phase response, oxidative/metabolic stress, extracellular matrix remodeling and signaling regulation.

Figure 7. Metabolic and oxidative-stress bridge among selected candidates.
Bar plot showing genes linked to metabolism, oxidative stress, redox biology, glutathione-related pathways or other metabolic/stress-response annotations. This figure should be interpreted as pathway-level functional inference, not direct metabolomics.

Figure 8. Gene–KEGG pathway evidence matrix.
Heatmap showing the presence of selected genes across the most represented KEGG pathways, highlighting pathway convergence among candidates.
"""

(REPORTS / "Figure_legends_final_package.txt").write_text(legends, encoding="utf-8")

results_pt = f"""Após a organização final dos resultados, o pacote analítico consolidou {n_total} candidatos gênicos ranqueados com evidências integradas. A estrutura final combinou priorização funcional, anotação GO, mapeamento em vias KEGG, proximidade funcional gene-gene, recuperação de interações proteína-proteína pelo STRING, remapeamento de sondas GPL570 e validação transcriptômica exploratória no GSE13205 por expressão diferencial e AUC direcional.

Além das figuras finais previamente geradas, este pipeline acrescentou uma camada funcional mais robusta, com mapas gene-via KEGG, matriz gene-eixo funcional e uma análise de ponte metabólica/estresse oxidativo baseada em termos GO/KEGG. Essa abordagem permite discutir os genes não apenas como marcadores isolados, mas como componentes de eixos biológicos complementares da sepse, incluindo resposta imune/inflamatória, ativação leucocitária, fase aguda, remodelamento tecidual, sinalização celular e estresse oxidativo/metabólico.

Os principais candidatos integrados foram: {", ".join(top_integrated)}. A interpretação dos achados transcriptômicos deve considerar que o GSE13205 deriva de biópsias de músculo vasto lateral, de modo que os resultados representam evidência tecido-específica e exploratória, e não validação diagnóstica direta em sangue.
"""

(REPORTS / "Results_text_PT_final_package.txt").write_text(results_pt, encoding="utf-8")

results_en = f"""After final organization of the analytical outputs, the package retained {n_total} ranked gene candidates with integrated evidence. The final framework combined functional prioritization, GO annotation, KEGG pathway mapping, gene-gene functional proximity, STRING protein-protein interaction recovery, GPL570 probe remapping, and exploratory transcriptomic validation in GSE13205 through differential expression and directional AUC analysis.

In addition to the previously generated final figures, this pipeline added a more robust functional layer, including KEGG gene-pathway maps, gene-functional axis matrices, and a metabolism/oxidative-stress bridge based on GO/KEGG terms. This approach allows the candidates to be interpreted not only as isolated markers but also as components of complementary biological axes of sepsis, including immune/inflammatory response, leukocyte activation, acute-phase response, tissue remodeling, cellular signaling, and oxidative/metabolic stress.

The top integrated candidates were: {", ".join(top_integrated)}. Transcriptomic findings should be interpreted considering that GSE13205 is based on vastus lateralis muscle biopsies; therefore, these results represent tissue-specific exploratory evidence rather than direct blood-based diagnostic validation.
"""

(REPORTS / "Results_text_EN_final_package.txt").write_text(results_en, encoding="utf-8")

# Manifest
manifest = {
    "final_package": str(FINAL),
    "tables_copied_or_generated": sorted([str(p) for p in TABLES.glob("*")]),
    "figures_copied_or_generated": sorted([str(p) for p in FIGURES.glob("*")]),
    "reports_generated": sorted([str(p) for p in REPORTS.glob("*")]),
    "scripts_archived": sorted([str(p) for p in SCRIPTS.glob("*")]),
    "logs_archived": sorted([str(p) for p in LOGS.glob("*")]),
}

(MANIFEST_DIR / "final_package_manifest.json").write_text(
    json.dumps(manifest, indent=2, ensure_ascii=False),
    encoding="utf-8"
)

manifest_table = []
for folder_name, folder in [
    ("tables", TABLES),
    ("figures", FIGURES),
    ("reports", REPORTS),
    ("scripts_used", SCRIPTS),
    ("logs_used", LOGS)
]:
    for p in sorted(folder.glob("*")):
        if p.is_file():
            manifest_table.append({
                "section": folder_name,
                "file": p.name,
                "path": str(p),
                "size_kb": round(p.stat().st_size / 1024, 2)
            })

manifest_df = pd.DataFrame(manifest_table)
save_table(manifest_df, MANIFEST_DIR / "final_package_manifest.csv")

msg("============================================================")
msg("PIPELINE 105 FINALIZADO")
msg("============================================================")
msg(f"Pasta final: {FINAL}")
msg(f"Relatório: {REPORTS / 'FINAL_PACKAGE_EXECUTIVE_REPORT.txt'}")
msg(f"Excel consolidado: {xlsx_path}")
msg(f"Manifesto: {MANIFEST_DIR / 'final_package_manifest.csv'}")
msg("============================================================")
PY

python3 scripts/105_organize_final_package_and_functional_metabolic_figures.py 2>&1 | tee "$LOG"

echo ""
echo "============================================================"
echo "PIPELINE 105 FINALIZADO"
echo "Pasta final:"
echo "04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE"
echo ""
echo "Relatório principal:"
echo "04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE/reports/FINAL_PACKAGE_EXECUTIVE_REPORT.txt"
echo ""
echo "Excel consolidado:"
echo "04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE/tables/FINAL_SEPSIS_MANUSCRIPT_TABLES_COMPLETE.xlsx"
echo ""
echo "Figuras:"
echo "04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE/figures"
echo ""
echo "Manifesto:"
echo "04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE/manifest/final_package_manifest.csv"
echo ""
echo "Log:"
echo "$LOG"
echo "============================================================"
