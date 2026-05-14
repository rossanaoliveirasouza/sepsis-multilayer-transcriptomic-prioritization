# sepsis-multilayer-transcriptomic-prioritization

Repository associated with the manuscript:

**Multilayer Transcriptomic Prioritization Identifies Candidate Immune, Redox, and Immunometabolic Gene Networks in Sepsis**

This repository contains scripts and outputs used to prioritize candidate genes in sepsis through an integrated transcriptomic and functional evidence framework.

## Overview

This repository contains the scripts and curated outputs associated with a computational transcriptomic study for the identification, prioritization, exploratory validation, and biological characterization of candidate genes associated with sepsis.

The workflow integrates model-based gene selection, dimensionality reduction, transcriptomic validation, functional annotation, network analysis, STRING protein interaction recovery, and KEGG-based route interpretation.

The analysis used three publicly available sepsis-related transcriptomic datasets:

- `GSE12624`: blood samples from polytrauma patients with and without subsequent sepsis;
- `GSE13205`: vastus lateralis muscle biopsies from critically ill septic patients and controls;
- `GSE69063`: RNA-seq profiles from patients with sepsis and healthy controls.

Candidate genes were initially selected independently within each dataset using an in-house modified logistic regression model adapted for high-dimensional transcriptomic data. Selected genes were then integrated across multiple evidence layers, including gene annotation, GO terms, KEGG pathways, STRING interactions, GEO expression support, GPL570 probe remapping, differential expression, and univariate AUC analysis.

The GSE13205 dataset was also used as an exploratory validation layer after corrected GPL570 probe remapping. Because this dataset is based on skeletal muscle biopsies rather than blood, its validation results should be interpreted as tissue-specific exploratory evidence, not as direct blood-based diagnostic validation.

Functional interpretation was performed through GO, KEGG, gene-gene functional proximity, STRING protein interaction evidence, and KEGG route/compound annotation. KEGG-derived compounds and metabolites are interpreted as pathway annotation outputs, not as experimentally measured metabolomic data.

The final repository includes reproducible scripts, logs, ranked candidate gene tables, transcriptomic validation outputs, functional and metabolic annotation tables, network results, manuscript-ready figures, and final packaging files.

## Repository structure

```text
.
├── scripts/
├── logs/
├── 04_results/
│   ├── repaired/
│   ├── repaired_gse13205_gpl570/
│   ├── final_manuscript/
│   └── FINAL_SEPSIS_MANUSCRIPT_PACKAGE/
└── README.md
````

## Main outputs

Final manuscript-ready files are organized in:

```text
04_results/FINAL_SEPSIS_MANUSCRIPT_PACKAGE/
├── tables/
├── figures/
├── reports/
├── scripts_used/
├── logs_used/
├── supplementary/
└── manifest/
```

The final package includes ranked gene candidates, differential expression and AUC summaries, GO/KEGG evidence, STRING network outputs, figures, reports, and a manifest of generated files.

## Main scripts

```text
99_run_post_gene_sepsis_pipeline.sh
101_repair_string_auc_differential_sepsis.sh
102_fix_GSE13205_phenotype_and_rerun_DE_AUC.sh
103_GSE13205_GPL570_probe_mapping_DE_AUC.sh
104_make_final_figures_tables_and_results_text.sh
105_organize_final_package_and_functional_metabolic_figures.sh
```

## Reproducibility

The scripts were designed to run from the project root directory.

Required R packages include:

```text
GEOquery
Biobase
AnnotationDbi
hgu133plus2.db
ggplot2
dplyr
openxlsx
stringr
scales
```

Python scripts require:

```text
pandas
numpy
matplotlib
```

## Citation

If you use this repository, please cite:

Rossana Souza et al. (2026).
**Multilayer Transcriptomic Prioritization Reveals Immune, Redox, and Immunometabolic Candidate Gene Networks in Sepsis**.
Manuscript in preparation.

## License

Code is released under the MIT License.

```
```
