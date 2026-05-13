#!/usr/bin/env bash
set -u

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT" || exit 1

mkdir -p scripts logs 04_results/final_manuscript/tables 04_results/final_manuscript/figures 04_results/final_manuscript/reports

LOG="logs/104_make_final_figures_tables_and_results_text.log"

cat > scripts/104_make_final_figures_tables_and_results_text.R <<'RSCRIPT'
options(stringsAsFactors = FALSE)

root <- "/mnt/f/Marcos/rossana"
setwd(root)

outdir <- "04_results/final_manuscript"
tabdir <- file.path(outdir, "tables")
figdir <- file.path(outdir, "figures")
repdir <- file.path(outdir, "reports")

dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(repdir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n")
}

safe_read <- function(path) {
  if (!file.exists(path)) {
    log_msg("Arquivo não encontrado:", path)
    return(data.frame())
  }
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      log_msg("Erro lendo:", path, conditionMessage(e))
      data.frame()
    }
  )
}

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (p in c("ggplot2", "dplyr", "readr", "openxlsx", "stringr", "scales")) {
  install_if_missing(p)
}

library(ggplot2)
library(dplyr)
library(openxlsx)
library(stringr)
library(scales)

log_msg("============================================================")
log_msg("SCRIPT 104: TABELAS, FIGURAS E TEXTO FINAL")
log_msg("============================================================")

final_prior <- safe_read("04_results/repaired_gse13205_gpl570/Table_103G_integrated_prioritization_GPL570_DE_AUC.csv")
de <- safe_read("04_results/repaired_gse13205_gpl570/Table_103E_GSE13205_differential_expression_GPL570.csv")
auc <- safe_read("04_results/repaired_gse13205_gpl570/Table_103F_GSE13205_univariate_AUC_GPL570.csv")
probe <- safe_read("04_results/repaired_gse13205_gpl570/Table_103C_GSE13205_GPL570_feature_matching_audit.csv")
string_edges <- safe_read("04_results/repaired/string/Table_R02_STRING_interactions_repaired_api.csv")
string_cent <- safe_read("04_results/repaired/string/Table_R03_STRING_centrality_repaired_api.csv")
go <- safe_read("04_results/enrichment/Table_03_GO_terms_selected_genes.csv")
kegg <- safe_read("04_results/enrichment/Table_05_KEGG_pathways_selected_genes.csv")
prox <- safe_read("04_results/networks/Table_07_gene_gene_functional_proximity_GO_KEGG.csv")
expr_long <- safe_read("04_results/repaired_gse13205_gpl570/Table_103D_GSE13205_selected_gene_expression_long_GPL570.csv")

if (nrow(final_prior) == 0) {
  stop("Tabela final de priorização não encontrada ou vazia.")
}

# ============================================================
# 1. LIMPEZA DA PRIORIZAÇÃO
# ============================================================

final_prior$gene_symbol <- as.character(final_prior$gene_symbol)

clean_prior <- final_prior %>%
  filter(
    !is.na(gene_symbol),
    str_trim(gene_symbol) != "",
    !is.na(integrated_score_final_GPL570),
    integrated_score_final_GPL570 > 0
  ) %>%
  arrange(
    desc(integrated_score_final_GPL570),
    GSE13205_GPL570_best_adj_p,
    desc(GSE13205_GPL570_best_auc_directional),
    gene_symbol
  )

manuscript_cols <- c(
  "gene_symbol",
  "evidence_group",
  "analysis_priority",
  "integrated_score_final_GPL570",
  "has_GO_annotation",
  "has_KEGG_pathway",
  "has_STRING_interaction_repaired",
  "has_GSE13205_GPL570_probe_mapping",
  "has_GSE13205_GPL570_DE",
  "has_GSE13205_GPL570_AUC",
  "GSE13205_GPL570_best_adj_p",
  "GSE13205_GPL570_best_auc_directional",
  "GSE13205_GPL570_direction"
)

manuscript_cols <- manuscript_cols[manuscript_cols %in% names(clean_prior)]

manuscript_table <- clean_prior[, manuscript_cols, drop = FALSE]

write.csv(
  clean_prior,
  file.path(tabdir, "Table_Final_01_clean_integrated_gene_prioritization.csv"),
  row.names = FALSE
)

write.csv(
  manuscript_table,
  file.path(tabdir, "Table_Final_02_manuscript_gene_prioritization.csv"),
  row.names = FALSE
)

# ============================================================
# 2. TABELA RESUMIDA DE DE + AUC
# ============================================================

de_auc <- de %>%
  select(
    gene_symbol,
    n_probes,
    matched_probes,
    n_sepsis,
    n_control,
    mean_sepsis,
    mean_control,
    logFC_sepsis_minus_control,
    p_value_wilcoxon,
    adj_p_value_BH,
    direction
  ) %>%
  left_join(
    auc %>%
      select(
        gene_symbol,
        auc_sepsis_control,
        auc_abs_directional
      ),
    by = "gene_symbol"
  ) %>%
  arrange(adj_p_value_BH, desc(auc_abs_directional))

write.csv(
  de_auc,
  file.path(tabdir, "Table_Final_03_GSE13205_DE_AUC_summary.csv"),
  row.names = FALSE
)

top_transcriptomic <- de_auc %>%
  filter(!is.na(adj_p_value_BH)) %>%
  arrange(adj_p_value_BH, desc(auc_abs_directional)) %>%
  head(12)

write.csv(
  top_transcriptomic,
  file.path(tabdir, "Table_Final_04_top_transcriptomic_candidates.csv"),
  row.names = FALSE
)

# ============================================================
# 3. EXCEL FINAL
# ============================================================

xlsx_file <- file.path(tabdir, "Final_sepsis_gene_results_for_manuscript.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "clean_prioritization")
writeData(wb, "clean_prioritization", clean_prior)

addWorksheet(wb, "manuscript_table")
writeData(wb, "manuscript_table", manuscript_table)

addWorksheet(wb, "GSE13205_DE_AUC")
writeData(wb, "GSE13205_DE_AUC", de_auc)

addWorksheet(wb, "top_transcriptomic")
writeData(wb, "top_transcriptomic", top_transcriptomic)

addWorksheet(wb, "probe_mapping")
writeData(wb, "probe_mapping", probe)

addWorksheet(wb, "STRING_edges")
writeData(wb, "STRING_edges", string_edges)

addWorksheet(wb, "STRING_centrality")
writeData(wb, "STRING_centrality", string_cent)

addWorksheet(wb, "GO_terms")
writeData(wb, "GO_terms", go)

addWorksheet(wb, "KEGG_pathways")
writeData(wb, "KEGG_pathways", kegg)

addWorksheet(wb, "functional_proximity")
writeData(wb, "functional_proximity", prox)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# ============================================================
# 4. FIGURAS
# ============================================================

# Figura 1: Top escore integrado
fig1_df <- manuscript_table %>%
  filter(!is.na(gene_symbol)) %>%
  arrange(desc(integrated_score_final_GPL570)) %>%
  head(25)

fig1_df$gene_symbol <- factor(fig1_df$gene_symbol, levels = rev(fig1_df$gene_symbol))

p1 <- ggplot(fig1_df, aes(x = gene_symbol, y = integrated_score_final_GPL570)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Integrated prioritization of sepsis-related gene candidates",
    x = "Gene",
    y = "Integrated score"
  )

ggsave(
  file.path(figdir, "Figure_Final_01_integrated_gene_prioritization.png"),
  p1,
  width = 9,
  height = 7,
  dpi = 300
)

# Figura 2: DE/AUC top genes
fig2_df <- de_auc %>%
  filter(!is.na(adj_p_value_BH), !is.na(auc_abs_directional)) %>%
  arrange(adj_p_value_BH, desc(auc_abs_directional)) %>%
  head(20)

fig2_df$gene_symbol <- factor(fig2_df$gene_symbol, levels = rev(fig2_df$gene_symbol))
fig2_df$neglog10_fdr <- -log10(fig2_df$adj_p_value_BH)

p2 <- ggplot(fig2_df, aes(x = gene_symbol, y = auc_abs_directional, size = neglog10_fdr)) +
  geom_point(alpha = 0.8) +
  coord_flip() +
  theme_bw(base_size = 12) +
  scale_y_continuous(limits = c(0.5, 1.02)) +
  labs(
    title = "Transcriptomic discrimination in GSE13205",
    x = "Gene",
    y = "Directional AUC",
    size = "-log10(FDR)"
  )

ggsave(
  file.path(figdir, "Figure_Final_02_GSE13205_AUC_FDR_bubble.png"),
  p2,
  width = 9,
  height = 7,
  dpi = 300
)

# Figura 3: direção do logFC
fig3_df <- de_auc %>%
  filter(!is.na(logFC_sepsis_minus_control)) %>%
  arrange(desc(abs(logFC_sepsis_minus_control))) %>%
  head(20)

fig3_df$gene_symbol <- factor(fig3_df$gene_symbol, levels = rev(fig3_df$gene_symbol))

p3 <- ggplot(fig3_df, aes(x = gene_symbol, y = logFC_sepsis_minus_control)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Expression differences between sepsis and control samples",
    x = "Gene",
    y = "Mean expression difference: sepsis minus control"
  )

ggsave(
  file.path(figdir, "Figure_Final_03_GSE13205_logFC_direction.png"),
  p3,
  width = 9,
  height = 7,
  dpi = 300
)

# Figura 4: Boxplots dos top genes
if (nrow(expr_long) > 0) {
  top_box <- de_auc %>%
    filter(!is.na(adj_p_value_BH)) %>%
    arrange(adj_p_value_BH, desc(auc_abs_directional)) %>%
    head(12) %>%
    pull(gene_symbol)

  box_df <- expr_long %>%
    filter(gene_symbol %in% top_box)

  box_df$group <- factor(box_df$group, levels = c("control", "sepsis"))

  if (nrow(box_df) > 0) {
    p4 <- ggplot(box_df, aes(x = group, y = expression)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, alpha = 0.75) +
      facet_wrap(~ gene_symbol, scales = "free_y") +
      theme_bw(base_size = 12) +
      labs(
        title = "Top transcriptomic candidates in GSE13205",
        x = "Group",
        y = "Expression"
      )

    ggsave(
      file.path(figdir, "Figure_Final_04_top_gene_expression_boxplots.png"),
      p4,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}

# Figura 5: STRING centrality
if (nrow(string_cent) > 0 && "gene_symbol" %in% names(string_cent)) {
  degree_col <- names(string_cent)[grepl("degree", names(string_cent), ignore.case = TRUE)][1]

  if (!is.na(degree_col)) {
    fig5_df <- string_cent %>%
      arrange(desc(.data[[degree_col]])) %>%
      head(20)

    fig5_df$gene_symbol <- factor(fig5_df$gene_symbol, levels = rev(fig5_df$gene_symbol))

    p5 <- ggplot(fig5_df, aes(x = gene_symbol, y = .data[[degree_col]])) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(
        title = "STRING-based network centrality",
        x = "Gene",
        y = "Degree"
      )

    ggsave(
      file.path(figdir, "Figure_Final_05_STRING_network_centrality.png"),
      p5,
      width = 9,
      height = 7,
      dpi = 300
    )
  }
}

# Figura 6: matriz de evidências
evidence_cols <- c(
  "has_GO_annotation",
  "has_KEGG_pathway",
  "has_STRING_interaction_repaired",
  "has_GSE13205_GPL570_probe_mapping",
  "has_GSE13205_GPL570_DE",
  "has_GSE13205_GPL570_AUC"
)

evidence_cols <- evidence_cols[evidence_cols %in% names(manuscript_table)]

fig6_df <- manuscript_table %>%
  arrange(desc(integrated_score_final_GPL570)) %>%
  head(25) %>%
  select(gene_symbol, all_of(evidence_cols))

if (length(evidence_cols) > 0) {
  long6 <- data.frame()

  for (cc in evidence_cols) {
    tmp <- data.frame(
      gene_symbol = fig6_df$gene_symbol,
      evidence = cc,
      present = as.logical(fig6_df[[cc]]),
      stringsAsFactors = FALSE
    )
    long6 <- rbind(long6, tmp)
  }

  long6$gene_symbol <- factor(long6$gene_symbol, levels = rev(unique(fig6_df$gene_symbol)))
  long6$evidence <- factor(
    long6$evidence,
    levels = evidence_cols,
    labels = c(
      "GO",
      "KEGG",
      "STRING/PPI",
      "GPL570 probe",
      "DE",
      "AUC"
    )[seq_along(evidence_cols)]
  )

  p6 <- ggplot(long6, aes(x = evidence, y = gene_symbol, fill = present)) +
    geom_tile(color = "white") +
    theme_bw(base_size = 12) +
    labs(
      title = "Evidence matrix for top prioritized candidates",
      x = "Evidence layer",
      y = "Gene",
      fill = "Present"
    )

  ggsave(
    file.path(figdir, "Figure_Final_06_evidence_matrix_top_candidates.png"),
    p6,
    width = 9,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# 5. TEXTO DE RESULTADOS
# ============================================================

n_total <- nrow(clean_prior)
n_de <- nrow(de)
n_auc <- nrow(auc)
n_probe <- sum(probe$n_probes_in_expression_matrix > 0, na.rm = TRUE)
n_go <- nrow(go)
n_kegg <- nrow(kegg)
n_prox <- nrow(prox)
n_string_edges <- nrow(string_edges)
n_string_nodes <- nrow(string_cent)

top_integrated <- manuscript_table %>%
  arrange(desc(integrated_score_final_GPL570)) %>%
  head(6) %>%
  pull(gene_symbol)

top_de <- de_auc %>%
  arrange(adj_p_value_BH, desc(auc_abs_directional)) %>%
  head(8) %>%
  pull(gene_symbol)

text_en <- paste0(
"Post-selection functional and transcriptomic prioritization identified a biologically coherent set of sepsis-related candidate genes. ",
"The final analytical framework integrated GO annotation, KEGG pathway mapping, gene-pathway functional proximity, STRING-based protein-protein interaction evidence, GPL570 probe remapping, differential expression, and univariate AUC analysis. ",
"Overall, ", n_total, " gene-level candidates with nonzero evidence were retained after cleaning. ",
"The functional layer yielded ", n_go, " GO records, ", n_kegg, " KEGG pathway records, and ", n_prox, " gene-gene functional proximity pairs, supporting convergence across immune, inflammatory, and stress-response biological processes. ",
"STRING recovery identified ", n_string_edges, " interactions and ", n_string_nodes, " centrality records, adding network-level support to the prioritization. ",
"\n\n",
"For transcriptomic validation, the corrected GSE13205 phenotype comprised 13 septic patients and 8 controls. ",
"After GPL570 probe remapping using the hgu133plus2.db annotation framework, ", n_probe, " genes were recovered in the expression matrix and ", n_de, " genes were evaluated by differential expression and AUC analyses. ",
"The strongest transcriptomic candidates were ", paste(top_de, collapse = ", "), ", which combined low FDR-adjusted p-values and high directional AUC values. ",
"In the integrated ranking, ", paste(top_integrated, collapse = ", "), " achieved the highest overall evidence scores, reflecting convergence across functional annotation, pathway participation, STRING/PPI support, and transcriptomic recovery. ",
"\n\n",
"These findings suggest that the selected genes capture complementary biological dimensions of sepsis, including inflammatory signaling, acute-phase response, tissue remodeling, oxidative or metabolic stress, and immune regulation. ",
"Because GSE13205 is based on vastus lateralis muscle biopsies, transcriptomic findings should be interpreted as tissue-specific evidence rather than direct blood-based diagnostic validation."
)

text_pt <- paste0(
"A priorização funcional e transcriptômica pós-seleção identificou um conjunto biologicamente coerente de genes candidatos relacionados à sepse. ",
"A estrutura analítica final integrou anotação GO, mapeamento em vias KEGG, proximidade funcional gene-via, evidência de interação proteína-proteína pelo STRING, remapeamento de sondas GPL570, expressão diferencial e análise univariada de AUC. ",
"No total, ", n_total, " candidatos gênicos com evidência diferente de zero foram mantidos após a limpeza. ",
"A camada funcional gerou ", n_go, " registros GO, ", n_kegg, " registros KEGG e ", n_prox, " pares de proximidade funcional gene-gene, sustentando convergência em processos biológicos relacionados à resposta imune, inflamação e estresse celular. ",
"A recuperação pelo STRING identificou ", n_string_edges, " interações e ", n_string_nodes, " registros de centralidade, adicionando suporte de rede à priorização. ",
"\n\n",
"Para a validação transcriptômica, o fenótipo corrigido do GSE13205 incluiu 13 pacientes sépticos e 8 controles. ",
"Após o remapeamento das sondas GPL570 com base no pacote hgu133plus2.db, ", n_probe, " genes foram recuperados na matriz de expressão e ", n_de, " genes foram avaliados por expressão diferencial e AUC. ",
"Os candidatos transcriptômicos mais fortes foram ", paste(top_de, collapse = ", "), ", combinando baixos valores de p ajustados por FDR e elevados valores de AUC direcional. ",
"Na classificação integrada, ", paste(top_integrated, collapse = ", "), " alcançaram os maiores escores globais, refletindo convergência entre anotação funcional, participação em vias, suporte STRING/PPI e recuperação transcriptômica. ",
"\n\n",
"Esses achados sugerem que os genes selecionados capturam dimensões biológicas complementares da sepse, incluindo sinalização inflamatória, resposta de fase aguda, remodelamento tecidual, estresse oxidativo ou metabólico e regulação imune. ",
"Como o GSE13205 é baseado em biópsias do músculo vasto lateral, os achados transcriptômicos devem ser interpretados como evidência tecido-específica, e não como validação diagnóstica direta em sangue."
)

writeLines(text_en, file.path(repdir, "Results_text_EN.txt"))
writeLines(text_pt, file.path(repdir, "Results_text_PT.txt"))

# ============================================================
# 6. RELATÓRIO EXECUTIVO
# ============================================================

sink(file.path(repdir, "Final_executive_report.txt"))

cat("============================================================\n")
cat("FINAL EXECUTIVE REPORT - SEPSIS GENE PRIORITIZATION\n")
cat("============================================================\n\n")

cat("Generated at:", as.character(Sys.time()), "\n\n")

cat("Input summary:\n")
cat("Clean gene candidates:", n_total, "\n")
cat("GO records:", n_go, "\n")
cat("KEGG records:", n_kegg, "\n")
cat("Functional proximity pairs:", n_prox, "\n")
cat("STRING interactions:", n_string_edges, "\n")
cat("STRING centrality records:", n_string_nodes, "\n")
cat("GPL570 recovered genes:", n_probe, "\n")
cat("Differential expression records:", n_de, "\n")
cat("AUC records:", n_auc, "\n\n")

cat("Top integrated candidates:\n")
print(manuscript_table %>% head(20))
cat("\n")

cat("Top transcriptomic candidates:\n")
print(top_transcriptomic)
cat("\n")

cat("Generated tables:\n")
cat(file.path(tabdir, "Table_Final_01_clean_integrated_gene_prioritization.csv"), "\n")
cat(file.path(tabdir, "Table_Final_02_manuscript_gene_prioritization.csv"), "\n")
cat(file.path(tabdir, "Table_Final_03_GSE13205_DE_AUC_summary.csv"), "\n")
cat(file.path(tabdir, "Table_Final_04_top_transcriptomic_candidates.csv"), "\n")
cat(xlsx_file, "\n\n")

cat("Generated figures:\n")
cat(file.path(figdir, "Figure_Final_01_integrated_gene_prioritization.png"), "\n")
cat(file.path(figdir, "Figure_Final_02_GSE13205_AUC_FDR_bubble.png"), "\n")
cat(file.path(figdir, "Figure_Final_03_GSE13205_logFC_direction.png"), "\n")
cat(file.path(figdir, "Figure_Final_04_top_gene_expression_boxplots.png"), "\n")
cat(file.path(figdir, "Figure_Final_05_STRING_network_centrality.png"), "\n")
cat(file.path(figdir, "Figure_Final_06_evidence_matrix_top_candidates.png"), "\n\n")

cat("Results text files:\n")
cat(file.path(repdir, "Results_text_EN.txt"), "\n")
cat(file.path(repdir, "Results_text_PT.txt"), "\n")

sink()

log_msg("Script 104 finalizado.")
log_msg("Relatório final:", file.path(repdir, "Final_executive_report.txt"))
RSCRIPT

Rscript scripts/104_make_final_figures_tables_and_results_text.R 2>&1 | tee "$LOG"

echo ""
echo "============================================================"
echo "SCRIPT 104 FINALIZADO"
echo "Relatório:"
echo "04_results/final_manuscript/reports/Final_executive_report.txt"
echo ""
echo "Excel final:"
echo "04_results/final_manuscript/tables/Final_sepsis_gene_results_for_manuscript.xlsx"
echo ""
echo "Figuras:"
echo "04_results/final_manuscript/figures/"
echo ""
echo "Textos:"
echo "04_results/final_manuscript/reports/Results_text_EN.txt"
echo "04_results/final_manuscript/reports/Results_text_PT.txt"
echo ""
echo "Log:"
echo "$LOG"
echo "============================================================"
