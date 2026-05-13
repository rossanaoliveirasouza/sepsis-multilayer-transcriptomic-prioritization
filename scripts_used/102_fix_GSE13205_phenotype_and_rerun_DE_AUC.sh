#!/usr/bin/env bash
set -u

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT" || exit 1

mkdir -p scripts logs 04_results/repaired_gse13205

LOG="logs/102_fix_GSE13205_phenotype_and_rerun_DE_AUC.log"

cat > scripts/102_fix_GSE13205_phenotype_and_rerun_DE_AUC.R <<'RSCRIPT'
options(stringsAsFactors = FALSE)
options(timeout = 600)

root <- "/mnt/f/Marcos/rossana"
setwd(root)

dir.create("04_results/repaired_gse13205", recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n")
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- toupper(trimws(x))
  x <- gsub("[^A-Z0-9_\\.-]", "", x)
  x
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- toupper(trimws(x))
  x <- gsub("[^A-Z0-9]", "", x)
  x
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

find_feature_columns <- function(fdat) {
  if (is.null(fdat) || ncol(fdat) == 0) return(character(0))

  nms <- names(fdat)
  hit <- grepl(
    "symbol|gene|refseq|gb_acc|accession|entrez|id|transcript|probe",
    nms,
    ignore.case = TRUE
  )

  unique(c(nms[hit], nms))
}

match_gene_rows <- function(eset, gene_symbol, selected_identifier) {
  expr <- Biobase::exprs(eset)
  fdat <- Biobase::fData(eset)

  rn <- rownames(expr)
  rn_clean <- clean_text(rn)

  gene_clean <- clean_gene(gene_symbol)
  id_clean <- clean_text(selected_identifier)

  hits <- rep(FALSE, length(rn))

  if (!is.na(id_clean) && nchar(id_clean) > 0) {
    hits <- hits | rn_clean == id_clean
  }

  if (!is.na(gene_clean) && nchar(gene_clean) > 0) {
    hits <- hits | clean_gene(rn) == gene_clean
  }

  if (!is.null(fdat) && nrow(fdat) == nrow(expr) && ncol(fdat) > 0) {
    cols <- find_feature_columns(fdat)

    for (cc in cols) {
      vec <- as.character(fdat[[cc]])
      vec_clean_text <- clean_text(vec)
      vec_clean_gene <- clean_gene(vec)

      if (!is.na(id_clean) && nchar(id_clean) > 0) {
        hits <- hits | vec_clean_text == id_clean
        hits <- hits | grepl(id_clean, vec_clean_text, fixed = TRUE)
      }

      if (!is.na(gene_clean) && nchar(gene_clean) > 0) {
        hits <- hits | vec_clean_gene == gene_clean
        hits <- hits | grepl(
          paste0("(^|[^A-Z0-9])", gene_clean, "([^A-Z0-9]|$)"),
          toupper(vec),
          perl = TRUE
        )
      }
    }
  }

  which(hits)
}

log_msg("============================================================")
log_msg("SCRIPT 102: CORREÇÃO FENOTÍPICA GSE13205 + DE + AUC")
log_msg("============================================================")

for (p in c("GEOquery", "Biobase")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Pacote ausente: ", p)
  }
}

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

genes_df <- unique(priority[, c("dataset", "selected_identifier", "gene_symbol", "evidence_group", "analysis_priority")])
genes_df$gene_symbol_clean <- clean_gene(genes_df$gene_symbol)
genes_df$selected_identifier_clean <- clean_text(genes_df$selected_identifier)

genes_df <- genes_df[
  !is.na(genes_df$gene_symbol_clean) &
    genes_df$gene_symbol_clean != "" &
    genes_df$dataset == "GSE13205",
]

if (nrow(genes_df) == 0) {
  log_msg("Nenhum gene específico de GSE13205 encontrado na priorização. Usando todos os genes com símbolo válido.")
  genes_df <- unique(priority[, c("dataset", "selected_identifier", "gene_symbol", "evidence_group", "analysis_priority")])
  genes_df$gene_symbol_clean <- clean_gene(genes_df$gene_symbol)
  genes_df$selected_identifier_clean <- clean_text(genes_df$selected_identifier)
  genes_df <- genes_df[!is.na(genes_df$gene_symbol_clean) & genes_df$gene_symbol_clean != "", ]
}

log_msg("Genes avaliados em GSE13205:", paste(unique(genes_df$gene_symbol_clean), collapse = ", "))

gsets <- GEOquery::getGEO("GSE13205", GSEMatrix = TRUE, getGPL = FALSE)

if (!is.list(gsets)) gsets <- list(gsets)

de_all <- data.frame()
auc_all <- data.frame()
expr_long <- data.frame()
pheno_all <- data.frame()
feature_audit <- data.frame()

for (i in seq_along(gsets)) {
  eset <- gsets[[i]]
  series_name <- paste0("GSE13205_series", i)

  expr <- Biobase::exprs(eset)
  pdat <- Biobase::pData(eset)

  group <- classify_GSE13205(pdat)

  pheno <- data.frame(
    dataset = "GSE13205",
    series = series_name,
    sample = rownames(pdat),
    group = group,
    title = if ("title" %in% names(pdat)) pdat$title else NA,
    source_name_ch1 = if ("source_name_ch1" %in% names(pdat)) pdat$source_name_ch1 else NA,
    characteristics_ch1 = apply(pdat, 1, function(z) {
      paste(z[grepl("characteristics", names(z), ignore.case = TRUE)], collapse = " | ")
    }),
    metadata_text = apply(pdat, 1, function(z) paste(z, collapse = " | ")),
    stringsAsFactors = FALSE
  )

  pheno_all <- rbind(pheno_all, pheno)

  log_msg(series_name, "| amostras:", ncol(expr),
          "| sepse:", sum(group == "sepsis", na.rm = TRUE),
          "| controle:", sum(group == "control", na.rm = TRUE),
          "| NA:", sum(is.na(group)))

  if (sum(group == "sepsis", na.rm = TRUE) < 2 || sum(group == "control", na.rm = TRUE) < 2) {
    log_msg("Contraste insuficiente em", series_name)
    next
  }

  y <- ifelse(group == "sepsis", 1, ifelse(group == "control", 0, NA))

  for (j in seq_len(nrow(genes_df))) {
    gene <- genes_df$gene_symbol_clean[j]
    ident <- genes_df$selected_identifier_clean[j]

    rows <- match_gene_rows(eset, gene_symbol = gene, selected_identifier = ident)

    feature_audit <- rbind(
      feature_audit,
      data.frame(
        dataset = "GSE13205",
        series = series_name,
        gene_symbol = gene,
        selected_identifier = genes_df$selected_identifier[j],
        n_matched_features = length(rows),
        matched_feature_ids = paste(rownames(expr)[rows], collapse = ";"),
        stringsAsFactors = FALSE
      )
    )

    if (length(rows) == 0) next

    mat <- expr[rows, , drop = FALSE]

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
      wilcox.test(x_sepsis, x_control)$p.value,
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
        selected_identifier = genes_df$selected_identifier[j],
        evidence_group = genes_df$evidence_group[j],
        analysis_priority = genes_df$analysis_priority[j],
        n_matched_features = length(rows),
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
        selected_identifier = genes_df$selected_identifier[j],
        evidence_group = genes_df$evidence_group[j],
        analysis_priority = genes_df$analysis_priority[j],
        n_matched_features = length(rows),
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
        selected_identifier = genes_df$selected_identifier[j],
        expression = gene_expr,
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
  pheno_all,
  "04_results/repaired_gse13205/Table_102_GSE13205_corrected_phenotype_map.csv",
  row.names = FALSE
)

write.csv(
  feature_audit,
  "04_results/repaired_gse13205/Table_103_GSE13205_feature_matching_audit.csv",
  row.names = FALSE
)

write.csv(
  expr_long,
  "04_results/repaired_gse13205/Table_104_GSE13205_selected_gene_expression_long.csv",
  row.names = FALSE
)

write.csv(
  de_all,
  "04_results/repaired_gse13205/Table_105_GSE13205_differential_expression_sepsis_vs_control.csv",
  row.names = FALSE
)

write.csv(
  auc_all,
  "04_results/repaired_gse13205/Table_106_GSE13205_univariate_AUC_sepsis_vs_control.csv",
  row.names = FALSE
)

# Atualiza priorização
priority2 <- priority

priority2$has_GSE13205_corrected_DE <- FALSE
priority2$has_GSE13205_corrected_AUC <- FALSE
priority2$GSE13205_best_adj_p <- NA_real_
priority2$GSE13205_best_auc_directional <- NA_real_
priority2$GSE13205_direction <- NA_character_

if (nrow(de_all) > 0) {
  for (g in unique(de_all$gene_symbol)) {
    idx <- clean_gene(priority2$gene_symbol) == clean_gene(g)

    sub_de <- de_all[de_all$gene_symbol == g, , drop = FALSE]
    sub_de <- sub_de[order(sub_de$adj_p_value_BH), , drop = FALSE]

    priority2$has_GSE13205_corrected_DE[idx] <- TRUE
    priority2$GSE13205_best_adj_p[idx] <- sub_de$adj_p_value_BH[1]
    priority2$GSE13205_direction[idx] <- sub_de$direction[1]
  }
}

if (nrow(auc_all) > 0) {
  for (g in unique(auc_all$gene_symbol)) {
    idx <- clean_gene(priority2$gene_symbol) == clean_gene(g)

    sub_auc <- auc_all[auc_all$gene_symbol == g, , drop = FALSE]
    sub_auc <- sub_auc[order(-sub_auc$auc_abs_directional), , drop = FALSE]

    priority2$has_GSE13205_corrected_AUC[idx] <- TRUE
    priority2$GSE13205_best_auc_directional[idx] <- sub_auc$auc_abs_directional[1]
  }
}

score_cols <- grep("^has_", names(priority2), value = TRUE)

for (cc in score_cols) {
  priority2[[cc]] <- as.logical(priority2[[cc]])
  priority2[[cc]][is.na(priority2[[cc]])] <- FALSE
}

priority2$integrated_score_final_GSE13205 <- rowSums(priority2[, score_cols, drop = FALSE], na.rm = TRUE)

priority2 <- priority2[order(
  -priority2$integrated_score_final_GSE13205,
  priority2$GSE13205_best_adj_p,
  -priority2$GSE13205_best_auc_directional,
  priority2$gene_symbol
), ]

write.csv(
  priority2,
  "04_results/repaired_gse13205/Table_107_integrated_prioritization_with_GSE13205_DE_AUC.csv",
  row.names = FALSE
)

sink("04_results/repaired_gse13205/report_102_GSE13205_DE_AUC.txt")

cat("============================================================\n")
cat("RELATÓRIO 102: GSE13205 COM FENÓTIPO CORRIGIDO\n")
cat("============================================================\n\n")

cat("A regra corrigida foi:\n")
cat("Muscle septic / septic patient = sepsis\n")
cat("Muscle control / control subject = control\n\n")

cat("Distribuição dos grupos:\n")
print(table(pheno_all$group, useNA = "ifany"))
cat("\n")

cat("Amostras por grupo:\n")
print(pheno_all[, c("sample", "group", "title", "source_name_ch1")])
cat("\n")

cat("Genes com features mapeadas:\n")
print(feature_audit[feature_audit$n_matched_features > 0, ])
cat("\n")

cat("Expressão diferencial:\n")
cat("Linhas:", nrow(de_all), "\n")
if (nrow(de_all) > 0) print(utils::head(de_all, 50))
cat("\n")

cat("AUC univariada:\n")
cat("Linhas:", nrow(auc_all), "\n")
if (nrow(auc_all) > 0) print(utils::head(auc_all, 50))
cat("\n")

cat("Priorização final com GSE13205:\n")
print(utils::head(priority2, 50))
cat("\n")

cat("Arquivos gerados:\n")
cat("04_results/repaired_gse13205/Table_102_GSE13205_corrected_phenotype_map.csv\n")
cat("04_results/repaired_gse13205/Table_103_GSE13205_feature_matching_audit.csv\n")
cat("04_results/repaired_gse13205/Table_104_GSE13205_selected_gene_expression_long.csv\n")
cat("04_results/repaired_gse13205/Table_105_GSE13205_differential_expression_sepsis_vs_control.csv\n")
cat("04_results/repaired_gse13205/Table_106_GSE13205_univariate_AUC_sepsis_vs_control.csv\n")
cat("04_results/repaired_gse13205/Table_107_integrated_prioritization_with_GSE13205_DE_AUC.csv\n")
cat("04_results/repaired_gse13205/report_102_GSE13205_DE_AUC.txt\n")

sink()

log_msg("Script 102 finalizado.")
log_msg("Relatório:", "04_results/repaired_gse13205/report_102_GSE13205_DE_AUC.txt")
RSCRIPT

Rscript scripts/102_fix_GSE13205_phenotype_and_rerun_DE_AUC.R 2>&1 | tee "$LOG"

echo ""
echo "============================================================"
echo "SCRIPT 102 FINALIZADO"
echo "Relatório:"
echo "04_results/repaired_gse13205/report_102_GSE13205_DE_AUC.txt"
echo ""
echo "Tabelas principais:"
echo "04_results/repaired_gse13205/Table_105_GSE13205_differential_expression_sepsis_vs_control.csv"
echo "04_results/repaired_gse13205/Table_106_GSE13205_univariate_AUC_sepsis_vs_control.csv"
echo "04_results/repaired_gse13205/Table_107_integrated_prioritization_with_GSE13205_DE_AUC.csv"
echo ""
echo "Log:"
echo "$LOG"
echo "============================================================"
