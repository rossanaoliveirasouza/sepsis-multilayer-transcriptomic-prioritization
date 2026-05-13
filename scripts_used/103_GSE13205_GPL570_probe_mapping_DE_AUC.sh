#!/usr/bin/env bash
set -u

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT" || exit 1

mkdir -p scripts logs 04_results/repaired_gse13205_gpl570

LOG="logs/103_GSE13205_GPL570_probe_mapping_DE_AUC.log"

cat > scripts/103_GSE13205_GPL570_probe_mapping_DE_AUC.R <<'RSCRIPT'
options(stringsAsFactors = FALSE)
options(timeout = 1200)

root <- "/mnt/f/Marcos/rossana"
setwd(root)

outdir <- "04_results/repaired_gse13205_gpl570"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n")
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      log_msg("Erro lendo:", path, conditionMessage(e))
      data.frame()
    }
  )
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- toupper(trimws(x))
  gsub("[^A-Z0-9]", "", x)
}

auc_rank <- function(y, score) {
  ok <- !is.na(y) & !is.na(score)
  y <- y[ok]
  score <- score[ok]

  if (length(unique(y)) != 2) return(NA_real_)

  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)

  if (n1 < 1 || n0 < 1) return(NA_real_)

  r <- rank(score, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  as.numeric(auc)
}

classify_GSE13205 <- function(pdat) {
  txt <- apply(pdat, 1, function(z) paste(z, collapse = " | "))
  low <- tolower(txt)

  group <- rep(NA_character_, length(low))

  is_sepsis <- grepl("^muscle septic", low) |
    grepl("vastus lateralis, septic patient", low, fixed = TRUE) |
    grepl("septic patient", low, fixed = TRUE)

  is_control <- grepl("^muscle control", low) |
    grepl("vastus lateralis, control subject", low, fixed = TRUE) |
    grepl("control subject", low, fixed = TRUE)

  group[is_sepsis] <- "sepsis"
  group[is_control] <- "control"

  group
}

ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)

  log_msg("Pacote ausente:", pkg)

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }

    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }

  requireNamespace(pkg, quietly = TRUE)
}

log_msg("============================================================")
log_msg("SCRIPT 103: GSE13205 + GPL570/hgu133plus2.db + DE/AUC")
log_msg("============================================================")

for (p in c("GEOquery", "Biobase", "AnnotationDbi", "hgu133plus2.db")) {
  ensure_pkg(p, bioc = TRUE)
}

library(GEOquery)
library(Biobase)
library(AnnotationDbi)
library(hgu133plus2.db)

priority_file <- "04_results/repaired/Table_R08_integrated_gene_candidate_prioritization_repaired.csv"

if (!file.exists(priority_file)) {
  priority_file <- "04_results/tables/Table_22_integrated_gene_candidate_prioritization.csv"
}

priority <- safe_read_csv(priority_file)

if (nrow(priority) == 0) {
  stop("Tabela de priorização não encontrada.")
}

if (!"gene_symbol" %in% names(priority)) {
  stop("Coluna gene_symbol ausente.")
}

if (!"selected_identifier" %in% names(priority)) {
  priority$selected_identifier <- priority$gene_symbol
}

priority$gene_symbol_clean <- clean_gene(priority$gene_symbol)

genes_df <- priority[
  !is.na(priority$gene_symbol_clean) &
    priority$gene_symbol_clean != "",
  ,
  drop = FALSE
]

genes_to_test <- sort(unique(genes_df$gene_symbol_clean))

log_msg("Total de genes válidos para teste:", length(genes_to_test))
log_msg("Genes:", paste(genes_to_test, collapse = ", "))

log_msg("Baixando/carregando GSE13205...")
gsets <- GEOquery::getGEO("GSE13205", GSEMatrix = TRUE, getGPL = FALSE)

if (!is.list(gsets)) gsets <- list(gsets)

# Mapeamento oficial da plataforma GPL570
probe_keys <- keys(hgu133plus2.db, keytype = "PROBEID")

probe_annot <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_keys,
  columns = c("SYMBOL", "ENTREZID", "GENENAME", "REFSEQ"),
  keytype = "PROBEID"
)

probe_annot$SYMBOL_CLEAN <- clean_gene(probe_annot$SYMBOL)
probe_annot <- probe_annot[!is.na(probe_annot$SYMBOL_CLEAN) & probe_annot$SYMBOL_CLEAN != "", ]

write.csv(
  probe_annot,
  file.path(outdir, "Table_103A_GPL570_probe_annotation_hgu133plus2db.csv"),
  row.names = FALSE
)

log_msg("Anotação GPL570 carregada:", nrow(probe_annot), "linhas")

all_pheno <- data.frame()
feature_audit <- data.frame()
expr_long <- data.frame()
de_all <- data.frame()
auc_all <- data.frame()

for (i in seq_along(gsets)) {
  eset <- gsets[[i]]
  series_name <- paste0("GSE13205_series", i)

  expr <- Biobase::exprs(eset)
  pdat <- Biobase::pData(eset)

  group <- classify_GSE13205(pdat)

  pheno <- data.frame(
    dataset = "GSE13205",
    series = series_name,
    sample = colnames(expr),
    group = group,
    title = if ("title" %in% names(pdat)) pdat$title else NA,
    source_name_ch1 = if ("source_name_ch1" %in% names(pdat)) pdat$source_name_ch1 else NA,
    stringsAsFactors = FALSE
  )

  all_pheno <- rbind(all_pheno, pheno)

  log_msg(
    series_name,
    "| features:", nrow(expr),
    "| amostras:", ncol(expr),
    "| sepse:", sum(group == "sepsis", na.rm = TRUE),
    "| controle:", sum(group == "control", na.rm = TRUE),
    "| NA:", sum(is.na(group))
  )

  if (sum(group == "sepsis", na.rm = TRUE) < 2 || sum(group == "control", na.rm = TRUE) < 2) {
    log_msg("Contraste insuficiente:", series_name)
    next
  }

  y <- ifelse(group == "sepsis", 1, ifelse(group == "control", 0, NA))

  expr_probes <- rownames(expr)

  for (gene in genes_to_test) {
    probes <- unique(probe_annot$PROBEID[probe_annot$SYMBOL_CLEAN == gene])
    probes <- probes[probes %in% expr_probes]

    feature_audit <- rbind(
      feature_audit,
      data.frame(
        dataset = "GSE13205",
        series = series_name,
        gene_symbol = gene,
        n_probes_in_GPL570 = sum(probe_annot$SYMBOL_CLEAN == gene, na.rm = TRUE),
        n_probes_in_expression_matrix = length(probes),
        matched_probes = paste(probes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )

    if (length(probes) == 0) next

    mat <- expr[probes, , drop = FALSE]

    if (nrow(mat) > 1) {
      gene_expr <- colMeans(mat, na.rm = TRUE)
    } else {
      gene_expr <- as.numeric(mat[1, ])
    }

    ok <- !is.na(y) & !is.na(gene_expr)

    if (sum(ok) < 4 || length(unique(y[ok])) < 2) next

    x_sepsis <- gene_expr[ok & y == 1]
    x_control <- gene_expr[ok & y == 0]

    if (length(x_sepsis) < 2 || length(x_control) < 2) next

    pval <- tryCatch(
      wilcox.test(x_sepsis, x_control, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    logfc <- mean(x_sepsis, na.rm = TRUE) - mean(x_control, na.rm = TRUE)
    auc <- auc_rank(y[ok], gene_expr[ok])

    de_all <- rbind(
      de_all,
      data.frame(
        dataset = "GSE13205",
        series = series_name,
        gene_symbol = gene,
        n_probes = length(probes),
        matched_probes = paste(probes, collapse = ";"),
        n_sepsis = length(x_sepsis),
        n_control = length(x_control),
        mean_sepsis = mean(x_sepsis, na.rm = TRUE),
        mean_control = mean(x_control, na.rm = TRUE),
        median_sepsis = median(x_sepsis, na.rm = TRUE),
        median_control = median(x_control, na.rm = TRUE),
        logFC_sepsis_minus_control = logfc,
        p_value_wilcoxon = pval,
        direction = ifelse(logfc > 0, "higher_in_sepsis", "lower_in_sepsis"),
        stringsAsFactors = FALSE
      )
    )

    auc_all <- rbind(
      auc_all,
      data.frame(
        dataset = "GSE13205",
        series = series_name,
        gene_symbol = gene,
        n_probes = length(probes),
        matched_probes = paste(probes, collapse = ";"),
        n_sepsis = length(x_sepsis),
        n_control = length(x_control),
        auc_sepsis_control = auc,
        auc_abs_directional = ifelse(is.na(auc), NA_real_, max(auc, 1 - auc)),
        mean_sepsis = mean(x_sepsis, na.rm = TRUE),
        mean_control = mean(x_control, na.rm = TRUE),
        direction = ifelse(logfc > 0, "higher_in_sepsis", "lower_in_sepsis"),
        stringsAsFactors = FALSE
      )
    )

    expr_long <- rbind(
      expr_long,
      data.frame(
        dataset = "GSE13205",
        series = series_name,
        sample = colnames(expr),
        group = group,
        gene_symbol = gene,
        expression = gene_expr,
        n_probes = length(probes),
        stringsAsFactors = FALSE
      )
    )
  }
}

if (nrow(de_all) > 0) {
  de_all$adj_p_value_BH <- p.adjust(de_all$p_value_wilcoxon, method = "BH")
  de_all <- de_all[order(de_all$adj_p_value_BH, de_all$p_value_wilcoxon), ]
}

if (nrow(auc_all) > 0) {
  auc_all <- auc_all[order(-auc_all$auc_abs_directional), ]
}

write.csv(
  all_pheno,
  file.path(outdir, "Table_103B_GSE13205_corrected_phenotype.csv"),
  row.names = FALSE
)

write.csv(
  feature_audit,
  file.path(outdir, "Table_103C_GSE13205_GPL570_feature_matching_audit.csv"),
  row.names = FALSE
)

write.csv(
  expr_long,
  file.path(outdir, "Table_103D_GSE13205_selected_gene_expression_long_GPL570.csv"),
  row.names = FALSE
)

write.csv(
  de_all,
  file.path(outdir, "Table_103E_GSE13205_differential_expression_GPL570.csv"),
  row.names = FALSE
)

write.csv(
  auc_all,
  file.path(outdir, "Table_103F_GSE13205_univariate_AUC_GPL570.csv"),
  row.names = FALSE
)

# Atualizar priorização final
priority2 <- priority

priority2$has_GSE13205_GPL570_probe_mapping <- clean_gene(priority2$gene_symbol) %in%
  feature_audit$gene_symbol[feature_audit$n_probes_in_expression_matrix > 0]

priority2$has_GSE13205_GPL570_DE <- clean_gene(priority2$gene_symbol) %in% clean_gene(de_all$gene_symbol)
priority2$has_GSE13205_GPL570_AUC <- clean_gene(priority2$gene_symbol) %in% clean_gene(auc_all$gene_symbol)

priority2$GSE13205_GPL570_best_adj_p <- NA_real_
priority2$GSE13205_GPL570_best_auc_directional <- NA_real_
priority2$GSE13205_GPL570_direction <- NA_character_

if (nrow(de_all) > 0) {
  for (g in unique(de_all$gene_symbol)) {
    idx <- clean_gene(priority2$gene_symbol) == clean_gene(g)
    sub <- de_all[de_all$gene_symbol == g, , drop = FALSE]
    sub <- sub[order(sub$adj_p_value_BH), , drop = FALSE]

    priority2$GSE13205_GPL570_best_adj_p[idx] <- sub$adj_p_value_BH[1]
    priority2$GSE13205_GPL570_direction[idx] <- sub$direction[1]
  }
}

if (nrow(auc_all) > 0) {
  for (g in unique(auc_all$gene_symbol)) {
    idx <- clean_gene(priority2$gene_symbol) == clean_gene(g)
    sub <- auc_all[auc_all$gene_symbol == g, , drop = FALSE]
    sub <- sub[order(-sub$auc_abs_directional), , drop = FALSE]

    priority2$GSE13205_GPL570_best_auc_directional[idx] <- sub$auc_abs_directional[1]
  }
}

score_cols <- grep("^has_", names(priority2), value = TRUE)

for (cc in score_cols) {
  priority2[[cc]] <- as.logical(priority2[[cc]])
  priority2[[cc]][is.na(priority2[[cc]])] <- FALSE
}

priority2$integrated_score_final_GPL570 <- rowSums(priority2[, score_cols, drop = FALSE], na.rm = TRUE)

priority2 <- priority2[order(
  -priority2$integrated_score_final_GPL570,
  priority2$GSE13205_GPL570_best_adj_p,
  -priority2$GSE13205_GPL570_best_auc_directional,
  priority2$gene_symbol
), ]

write.csv(
  priority2,
  file.path(outdir, "Table_103G_integrated_prioritization_GPL570_DE_AUC.csv"),
  row.names = FALSE
)

# Figuras simples, se houver expressão
if (nrow(expr_long) > 0) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  }

  library(ggplot2)

  top_auc <- head(auc_all$gene_symbol, 12)
  plot_df <- expr_long[expr_long$gene_symbol %in% top_auc, , drop = FALSE]

  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = group, y = expression)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, alpha = 0.7) +
      facet_wrap(~ gene_symbol, scales = "free_y") +
      theme_bw(base_size = 12) +
      labs(
        title = "GSE13205: selected gene expression by corrected phenotype",
        x = "Group",
        y = "Expression"
      )

    ggsave(
      filename = file.path(outdir, "Figure_103_GSE13205_top_gene_expression_boxplots.png"),
      plot = p,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}

sink(file.path(outdir, "report_103_GSE13205_GPL570_DE_AUC.txt"))

cat("============================================================\n")
cat("RELATÓRIO 103: GSE13205 COM MAPEAMENTO GPL570/hgu133plus2.db\n")
cat("============================================================\n\n")

cat("Fenótipo corrigido:\n")
print(table(all_pheno$group, useNA = "ifany"))
cat("\n")

cat("Total de genes testados:", length(genes_to_test), "\n")
cat("Genes testados:\n")
cat(paste(genes_to_test, collapse = ", "), "\n\n")

cat("Genes com probes encontradas na matriz:\n")
print(feature_audit[feature_audit$n_probes_in_expression_matrix > 0, ])
cat("\n")

cat("Expressão diferencial:\n")
cat("Linhas:", nrow(de_all), "\n")
if (nrow(de_all) > 0) print(utils::head(de_all, 50))
cat("\n")

cat("AUC univariada:\n")
cat("Linhas:", nrow(auc_all), "\n")
if (nrow(auc_all) > 0) print(utils::head(auc_all, 50))
cat("\n")

cat("Priorização final GPL570:\n")
print(utils::head(priority2, 50))
cat("\n")

cat("Arquivos gerados:\n")
cat(file.path(outdir, "Table_103A_GPL570_probe_annotation_hgu133plus2db.csv"), "\n")
cat(file.path(outdir, "Table_103B_GSE13205_corrected_phenotype.csv"), "\n")
cat(file.path(outdir, "Table_103C_GSE13205_GPL570_feature_matching_audit.csv"), "\n")
cat(file.path(outdir, "Table_103D_GSE13205_selected_gene_expression_long_GPL570.csv"), "\n")
cat(file.path(outdir, "Table_103E_GSE13205_differential_expression_GPL570.csv"), "\n")
cat(file.path(outdir, "Table_103F_GSE13205_univariate_AUC_GPL570.csv"), "\n")
cat(file.path(outdir, "Table_103G_integrated_prioritization_GPL570_DE_AUC.csv"), "\n")
cat(file.path(outdir, "Figure_103_GSE13205_top_gene_expression_boxplots.png"), "\n")

sink()

log_msg("Script 103 finalizado.")
log_msg("Relatório:", file.path(outdir, "report_103_GSE13205_GPL570_DE_AUC.txt"))
RSCRIPT

Rscript scripts/103_GSE13205_GPL570_probe_mapping_DE_AUC.R 2>&1 | tee "$LOG"

echo ""
echo "============================================================"
echo "SCRIPT 103 FINALIZADO"
echo "Relatório:"
echo "04_results/repaired_gse13205_gpl570/report_103_GSE13205_GPL570_DE_AUC.txt"
echo ""
echo "Tabelas principais:"
echo "04_results/repaired_gse13205_gpl570/Table_103E_GSE13205_differential_expression_GPL570.csv"
echo "04_results/repaired_gse13205_gpl570/Table_103F_GSE13205_univariate_AUC_GPL570.csv"
echo "04_results/repaired_gse13205_gpl570/Table_103G_integrated_prioritization_GPL570_DE_AUC.csv"
echo ""
echo "Log:"
echo "$LOG"
echo "============================================================"
