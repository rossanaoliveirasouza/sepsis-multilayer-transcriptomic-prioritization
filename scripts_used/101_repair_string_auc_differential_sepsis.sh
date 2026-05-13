#!/usr/bin/env bash
set -u

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT" || exit 1

mkdir -p scripts logs 04_results/repaired 04_results/repaired/string 04_results/repaired/geo

LOG="logs/101_repair_string_auc_differential_sepsis.log"

cat > scripts/101_repair_string_auc_differential_sepsis.R <<'RSCRIPT'
options(stringsAsFactors = FALSE)
options(timeout = 600)

root <- "/mnt/f/Marcos/rossana"
setwd(root)

dir.create("04_results/repaired", recursive = TRUE, showWarnings = FALSE)
dir.create("04_results/repaired/string", recursive = TRUE, showWarnings = FALSE)
dir.create("04_results/repaired/geo", recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n")
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      log_msg("ERRO lendo:", path, " :: ", conditionMessage(e))
      data.frame()
    }
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

bh_adjust <- function(p) {
  if (length(p) == 0) return(p)
  p.adjust(p, method = "BH")
}

infer_group <- function(pdat) {
  txt <- apply(pdat, 1, function(z) {
    paste(z, collapse = " | ")
  })

  low <- tolower(txt)

  control_regex <- paste(
    c(
      "healthy control",
      "normal control",
      "control",
      "healthy",
      "normal",
      "non sepsis",
      "non-sepsis",
      "nonsepsis",
      "non septic",
      "non-septic",
      "nonseptic",
      "uninfected",
      "no sepsis",
      "without sepsis"
    ),
    collapse = "|"
  )

  sepsis_regex <- paste(
    c(
      "sepsis",
      "septic",
      "septic shock",
      "septicemia",
      "septicaemia",
      "infection",
      "infected",
      "sirs",
      "patient",
      "case"
    ),
    collapse = "|"
  )

  is_control <- grepl(control_regex, low, perl = TRUE)
  is_sepsis <- grepl(sepsis_regex, low, perl = TRUE)

  group <- rep(NA_character_, length(low))

  group[is_control] <- "control"
  group[!is_control & is_sepsis] <- "sepsis"

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
        hits <- hits | grepl(paste0("(^|[^A-Z0-9])", gene_clean, "([^A-Z0-9]|$)"), toupper(vec), perl = TRUE)
      }
    }
  }

  which(hits)
}

log_msg("============================================================")
log_msg("REPARO DIRIGIDO: STRING + AUC + DIFERENCIAL")
log_msg("============================================================")

priority_file <- "04_results/tables/Table_22_integrated_gene_candidate_prioritization.csv"
input_file <- "04_results/tables/Table_01_input_selected_genes_and_identifiers.csv"

priority <- safe_read_csv(priority_file)
input <- safe_read_csv(input_file)

if (nrow(priority) == 0) {
  stop("Tabela principal de priorização não encontrada ou vazia.")
}

if (!"gene_symbol" %in% names(priority)) {
  stop("A tabela principal não possui coluna gene_symbol.")
}

if (!"selected_identifier" %in% names(priority)) {
  priority$selected_identifier <- priority$gene_symbol
}

genes_df <- unique(priority[, c("dataset", "selected_identifier", "gene_symbol", "evidence_group", "analysis_priority")])
genes_df$gene_symbol_clean <- clean_gene(genes_df$gene_symbol)
genes_df$selected_identifier_clean <- clean_text(genes_df$selected_identifier)

genes_df <- genes_df[genes_df$gene_symbol_clean != "" & !is.na(genes_df$gene_symbol_clean), ]

symbols <- unique(genes_df$gene_symbol_clean)
symbols <- symbols[!is.na(symbols) & symbols != ""]

datasets <- unique(genes_df$dataset)
datasets <- datasets[!is.na(datasets) & datasets != ""]

if (length(datasets) == 0) {
  datasets <- c("GSE12624", "GSE13205", "GSE69063")
}

log_msg("Genes únicos para reparo:", paste(symbols, collapse = ", "))
log_msg("Datasets para reparo:", paste(datasets, collapse = ", "))

# ============================================================
# 1) REPARO STRING VIA API
# ============================================================

log_msg("Tentando recuperar mapeamento e interações STRING via API...")

string_mapping <- data.frame()
string_network <- data.frame()

if (length(symbols) > 0) {
  encoded_ids <- utils::URLencode(paste(symbols, collapse = "\r"), reserved = TRUE)

  map_url <- paste0(
    "https://string-db.org/api/tsv/get_string_ids?identifiers=",
    encoded_ids,
    "&species=9606&limit=1&echo_query=1&caller_identity=sepsis_post_gene_pipeline"
  )

  string_mapping <- tryCatch(
    read.delim(map_url, sep = "\t", header = TRUE, quote = "", comment.char = ""),
    error = function(e) {
      log_msg("STRING mapping API falhou:", conditionMessage(e))
      data.frame()
    }
  )

  if (nrow(string_mapping) > 0) {
    write.csv(
      string_mapping,
      "04_results/repaired/string/Table_R01_STRING_mapping_repaired_api.csv",
      row.names = FALSE
    )
    log_msg("STRING mapping recuperado:", nrow(string_mapping), "linhas")
  } else {
    log_msg("STRING mapping continuou vazio.")
  }

  net_url <- paste0(
    "https://string-db.org/api/tsv/network?identifiers=",
    encoded_ids,
    "&species=9606&required_score=400&caller_identity=sepsis_post_gene_pipeline"
  )

  string_network <- tryCatch(
    read.delim(net_url, sep = "\t", header = TRUE, quote = "", comment.char = ""),
    error = function(e) {
      log_msg("STRING network API falhou:", conditionMessage(e))
      data.frame()
    }
  )

  if (nrow(string_network) > 0) {
    write.csv(
      string_network,
      "04_results/repaired/string/Table_R02_STRING_interactions_repaired_api.csv",
      row.names = FALSE
    )

    n1 <- if ("preferredName_A" %in% names(string_network)) string_network$preferredName_A else string_network[, 1]
    n2 <- if ("preferredName_B" %in% names(string_network)) string_network$preferredName_B else string_network[, 2]

    all_nodes <- sort(unique(c(as.character(n1), as.character(n2))))
    degree <- sapply(all_nodes, function(g) sum(n1 == g | n2 == g, na.rm = TRUE))

    centrality <- data.frame(
      gene_symbol = all_nodes,
      string_degree = as.integer(degree),
      stringsAsFactors = FALSE
    )

    if ("score" %in% names(string_network)) {
      centrality$string_mean_score <- sapply(all_nodes, function(g) {
        mean(string_network$score[n1 == g | n2 == g], na.rm = TRUE)
      })
    }

    centrality <- centrality[order(-centrality$string_degree), ]

    write.csv(
      centrality,
      "04_results/repaired/string/Table_R03_STRING_centrality_repaired_api.csv",
      row.names = FALSE
    )

    log_msg("STRING interações recuperadas:", nrow(string_network), "linhas")
    log_msg("STRING centralidade recuperada:", nrow(centrality), "genes")
  } else {
    log_msg("STRING network continuou vazio.")
  }
}

# ============================================================
# 2) REPARO GEO: AUC E DIFERENCIAL
# ============================================================

log_msg("Carregando pacotes GEOquery/Biobase para reparo GEO...")

needed <- c("GEOquery", "Biobase")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Pacote ausente: ", p, ". Instale antes de rodar o reparo.")
  }
}

geo_summary <- data.frame()
de_all <- data.frame()
auc_all <- data.frame()
metadata_audit_all <- data.frame()

for (gse in datasets) {
  log_msg("Processando GEO para reparo:", gse)

  gsets <- tryCatch(
    GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) {
      log_msg("Falha ao baixar/carregar", gse, ":", conditionMessage(e))
      NULL
    }
  )

  if (is.null(gsets)) next

  if (!is.list(gsets)) gsets <- list(gsets)

  for (i in seq_along(gsets)) {
    eset <- gsets[[i]]
    series_name <- paste0(gse, "_series", i)

    expr <- Biobase::exprs(eset)
    pdat <- Biobase::pData(eset)

    group <- infer_group(pdat)

    meta_audit <- data.frame(
      dataset = gse,
      series = series_name,
      sample = rownames(pdat),
      inferred_group = group,
      metadata_text = apply(pdat, 1, function(z) paste(z, collapse = " | ")),
      stringsAsFactors = FALSE
    )

    metadata_audit_all <- rbind(metadata_audit_all, meta_audit)

    n_sepsis <- sum(group == "sepsis", na.rm = TRUE)
    n_control <- sum(group == "control", na.rm = TRUE)
    n_unknown <- sum(is.na(group))

    geo_summary <- rbind(
      geo_summary,
      data.frame(
        dataset = gse,
        series = series_name,
        n_features = nrow(expr),
        n_samples = ncol(expr),
        n_sepsis = n_sepsis,
        n_control = n_control,
        n_unknown = n_unknown,
        can_test = n_sepsis >= 2 && n_control >= 2,
        stringsAsFactors = FALSE
      )
    )

    log_msg(series_name, "amostras:", ncol(expr), "| sepse:", n_sepsis, "| controle:", n_control, "| desconhecido:", n_unknown)

    if (!(n_sepsis >= 2 && n_control >= 2)) {
      log_msg("Sem contraste suficiente em", series_name)
      next
    }

    genes_this <- genes_df[genes_df$dataset == gse | is.na(genes_df$dataset) | genes_df$dataset == "", ]

    if (nrow(genes_this) == 0) {
      genes_this <- genes_df
    }

    for (j in seq_len(nrow(genes_this))) {
      gene <- genes_this$gene_symbol_clean[j]
      ident <- genes_this$selected_identifier_clean[j]

      rows <- match_gene_rows(eset, gene_symbol = gene, selected_identifier = ident)

      if (length(rows) == 0) next

      mat <- expr[rows, , drop = FALSE]

      if (nrow(mat) > 1) {
        gene_expr <- colMeans(mat, na.rm = TRUE)
      } else {
        gene_expr <- as.numeric(mat[1, ])
      }

      y <- ifelse(group == "sepsis", 1, ifelse(group == "control", 0, NA))

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
          dataset = gse,
          series = series_name,
          gene_symbol = gene,
          selected_identifier = genes_this$selected_identifier[j],
          n_matched_features = length(rows),
          n_sepsis = length(x_sepsis),
          n_control = length(x_control),
          mean_sepsis = mean(x_sepsis, na.rm = TRUE),
          mean_control = mean(x_control, na.rm = TRUE),
          logFC_sepsis_minus_control = logfc,
          p_value_wilcoxon = pval,
          direction = ifelse(logfc > 0, "higher_in_sepsis", "lower_in_sepsis"),
          stringsAsFactors = FALSE
        )
      )

      auc_all <- rbind(
        auc_all,
        data.frame(
          dataset = gse,
          series = series_name,
          gene_symbol = gene,
          selected_identifier = genes_this$selected_identifier[j],
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
    }
  }
}

if (nrow(de_all) > 0) {
  de_all$adj_p_value_BH <- bh_adjust(de_all$p_value_wilcoxon)
  de_all <- de_all[order(de_all$adj_p_value_BH, de_all$p_value_wilcoxon), ]
}

if (nrow(auc_all) > 0) {
  auc_all <- auc_all[order(-auc_all$auc_abs_directional), ]
}

write.csv(
  geo_summary,
  "04_results/repaired/geo/Table_R04_GEO_group_inference_summary.csv",
  row.names = FALSE
)

write.csv(
  metadata_audit_all,
  "04_results/repaired/geo/Table_R05_GEO_metadata_group_inference_audit.csv",
  row.names = FALSE
)

write.csv(
  de_all,
  "04_results/repaired/geo/Table_R06_repaired_differential_expression_sepsis_control.csv",
  row.names = FALSE
)

write.csv(
  auc_all,
  "04_results/repaired/geo/Table_R07_repaired_univariate_AUC_sepsis_control.csv",
  row.names = FALSE
)

# ============================================================
# 3) PRIORIZAÇÃO REPARADA
# ============================================================

priority2 <- priority

priority2$has_STRING_mapping_repaired <- FALSE
priority2$has_STRING_interaction_repaired <- FALSE
priority2$has_sepsis_control_AUC_repaired <- FALSE
priority2$has_differential_expression_repaired <- FALSE

if (nrow(string_mapping) > 0) {
  possible_cols <- intersect(c("queryItem", "preferredName", "stringId"), names(string_mapping))
  txt <- apply(string_mapping[, possible_cols, drop = FALSE], 1, paste, collapse = " ")
  mapped_genes <- unique(clean_gene(txt))
  priority2$has_STRING_mapping_repaired <- clean_gene(priority2$gene_symbol) %in% mapped_genes |
    clean_gene(priority2$selected_identifier) %in% mapped_genes
}

if (nrow(string_network) > 0) {
  net_names <- c()
  for (cc in c("preferredName_A", "preferredName_B", "stringId_A", "stringId_B")) {
    if (cc %in% names(string_network)) net_names <- c(net_names, string_network[[cc]])
  }
  net_genes <- unique(clean_gene(net_names))
  priority2$has_STRING_interaction_repaired <- clean_gene(priority2$gene_symbol) %in% net_genes
}

if (nrow(auc_all) > 0) {
  priority2$has_sepsis_control_AUC_repaired <- clean_gene(priority2$gene_symbol) %in% clean_gene(auc_all$gene_symbol)
}

if (nrow(de_all) > 0) {
  priority2$has_differential_expression_repaired <- clean_gene(priority2$gene_symbol) %in% clean_gene(de_all$gene_symbol)
}

score_cols_original <- c(
  "has_orgdb_annotation",
  "has_GO_annotation",
  "has_KEGG_pathway",
  "has_STRING_mapping",
  "has_UniProt_sequence",
  "detected_in_GEO_expression",
  "has_sepsis_control_AUC",
  "has_clinical_numeric_association",
  "has_clinical_binary_association"
)

score_cols_repaired <- c(
  score_cols_original,
  "has_STRING_mapping_repaired",
  "has_STRING_interaction_repaired",
  "has_sepsis_control_AUC_repaired",
  "has_differential_expression_repaired"
)

for (cc in score_cols_repaired) {
  if (!cc %in% names(priority2)) priority2[[cc]] <- FALSE
  priority2[[cc]] <- as.logical(priority2[[cc]])
  priority2[[cc]][is.na(priority2[[cc]])] <- FALSE
}

priority2$integrated_score_repaired <- rowSums(priority2[, score_cols_repaired, drop = FALSE], na.rm = TRUE)

priority2 <- priority2[order(-priority2$integrated_score_repaired, priority2$gene_symbol), ]

write.csv(
  priority2,
  "04_results/repaired/Table_R08_integrated_gene_candidate_prioritization_repaired.csv",
  row.names = FALSE
)

# ============================================================
# 4) RELATÓRIO FINAL
# ============================================================

sink("04_results/repaired/repaired_pipeline_report.txt")

cat("============================================================\n")
cat("RELATÓRIO DO REPARO DIRIGIDO\n")
cat("STRING + AUC + DIFERENCIAL SEPSE-CONTROLE\n")
cat("============================================================\n\n")

cat("Data:", as.character(Sys.time()), "\n\n")

cat("Genes avaliados:\n")
cat(paste(symbols, collapse = ", "), "\n\n")

cat("Datasets avaliados:\n")
cat(paste(datasets, collapse = ", "), "\n\n")

cat("1) STRING\n")
cat("Mapeamentos recuperados:", nrow(string_mapping), "\n")
cat("Interações recuperadas:", nrow(string_network), "\n\n")

cat("2) GEO: inferência de grupos\n")
print(geo_summary)
cat("\n")

cat("3) Expressão diferencial reparada\n")
cat("Linhas:", nrow(de_all), "\n")
if (nrow(de_all) > 0) {
  print(utils::head(de_all, 30))
}
cat("\n")

cat("4) AUC reparada\n")
cat("Linhas:", nrow(auc_all), "\n")
if (nrow(auc_all) > 0) {
  print(utils::head(auc_all, 30))
}
cat("\n")

cat("5) Priorização integrada reparada\n")
cat("Linhas:", nrow(priority2), "\n")
print(utils::head(priority2, 30))
cat("\n")

cat("Arquivos principais gerados:\n")
cat("04_results/repaired/string/Table_R01_STRING_mapping_repaired_api.csv\n")
cat("04_results/repaired/string/Table_R02_STRING_interactions_repaired_api.csv\n")
cat("04_results/repaired/string/Table_R03_STRING_centrality_repaired_api.csv\n")
cat("04_results/repaired/geo/Table_R04_GEO_group_inference_summary.csv\n")
cat("04_results/repaired/geo/Table_R05_GEO_metadata_group_inference_audit.csv\n")
cat("04_results/repaired/geo/Table_R06_repaired_differential_expression_sepsis_control.csv\n")
cat("04_results/repaired/geo/Table_R07_repaired_univariate_AUC_sepsis_control.csv\n")
cat("04_results/repaired/Table_R08_integrated_gene_candidate_prioritization_repaired.csv\n")
cat("04_results/repaired/repaired_pipeline_report.txt\n")

sink()

log_msg("Reparo finalizado.")
log_msg("Relatório:", "04_results/repaired/repaired_pipeline_report.txt")
log_msg("Priorização reparada:", "04_results/repaired/Table_R08_integrated_gene_candidate_prioritization_repaired.csv")
RSCRIPT

Rscript scripts/101_repair_string_auc_differential_sepsis.R 2>&1 | tee "$LOG"

echo ""
echo "============================================================"
echo "REPARO FINALIZADO"
echo "Veja o relatório:"
echo "04_results/repaired/repaired_pipeline_report.txt"
echo ""
echo "Tabela de priorização reparada:"
echo "04_results/repaired/Table_R08_integrated_gene_candidate_prioritization_repaired.csv"
echo ""
echo "Log:"
echo "$LOG"
echo "============================================================"
