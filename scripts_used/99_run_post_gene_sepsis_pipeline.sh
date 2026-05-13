#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/f/Marcos/rossana"
cd "$ROOT"

mkdir -p logs 04_results/logs

echo "============================================================"
echo "RODANDO PIPELINE PÓS-SELEÇÃO DE GENES PARA SEPSE"
echo "============================================================"
echo "Diretório: $ROOT"
echo ""

Rscript scripts/02_post_gene_sepsis_functional_clinical_pipeline.R 2>&1 | tee logs/02_post_gene_sepsis_functional_clinical_pipeline.log

echo ""
echo "============================================================"
echo "FINALIZADO"
echo "Relatório:"
echo "04_results/reports/post_gene_sepsis_pipeline_report.txt"
echo ""
echo "Tabela principal de priorização:"
echo "04_results/tables/Table_22_integrated_gene_candidate_prioritization.csv"
echo ""
echo "Log:"
echo "logs/02_post_gene_sepsis_functional_clinical_pipeline.log"
echo "============================================================"
