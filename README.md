# QIIME 2 Amplicon Analysis Pipelines (16S & 18S)

This repository contains reproducible pipelines for analyzing **16S and 18S rRNA amplicon sequencing data** using QIIME 2 within a Singularity container. Both pipelines include data import, primer trimming, denoising, taxonomic classification, phylogeny, and diversity analyses.

---

## **Software Environment**

- **QIIME 2**: Tested with version 2023.8 (Singularity container)
- **Singularity**: Used to ensure reproducible environments
- **Dependencies**: Cutadapt, DADA2, MAFFT, FastTree, scikit-learn

---

## **Data Preparation**

1. Place raw FASTQ files in `data/raw/`.
2. Prepare sample metadata files (`sample_metadata*.txt`) in `data/metadata/`.
3. Generate QIIME 2 manifest files in `data/manifests/` (16S and 18S pipelines handle this automatically).

---

## **Pipeline Overview**

### **16S Pipeline**
- Import paired-end reads
- Trim 16S primers (V4 region)
- Denoising using DADA2
- Rarefaction
- Taxonomic classification (SILVA 515-806)
- Phylogenetic tree construction
- Core diversity metrics (alpha, beta diversity)
- Taxonomic summaries (collapsed tables, relative frequencies)

Script: `scripts/16S_pipeline.sh`

---

### **18S Pipeline**
- Import paired-end reads
- Trim 18S primers (V9 region)
- Denoising using DADA2
- Rarefaction
- Taxonomic classification
  - SILVA V9 classifier
  - Optional: PR2 reference classifier
  - Consensus classifier using vsearch
- Phylogenetic tree construction
- Core diversity metrics (alpha, beta diversity)
- Taxonomic summaries (collapsed tables, relative frequencies)

Script: `scripts/18S_pipeline.sh`

---

## **Outputs**

- **.qza**: QIIME 2 artifacts (feature tables, representative sequences, taxonomy)
- **.qzv**: Visualizations (summaries, alpha/beta diversity plots, barplots)
- **Exported tables**: Rarefied tables, relative abundance tables (CSV compatible)

---

## **Usage Instructions**

1. Start Singularity container
2. Navigate to working directory
3. Run the desired pipeline
4. View results using qiime tools view or in QIIME 2 View (https://view.qiime2.org/)

---

## Notes
Sampling depths for rarefaction should be set per dataset.
Primer sequences and parameters can be modified in scripts.
Alpha and beta diversity analyses rely on filtered metadata files.

---

## References

QIIME 2: https://qiime2.org
SILVA rRNA Database: https://www.arb-silva.de/
PR2 rRNA Database: https://pr2-database.org/

