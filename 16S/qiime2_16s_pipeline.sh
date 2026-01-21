###############################################################################
# QIIME 2 16S rRNA amplicon analysis pipeline
# Environment: Singularity + QIIME 2
# Data type: Paired-end Illumina FASTQ (Phred33)
# Region: 16S rRNA (V4, 515F–806R)
###############################################################################

############################
# 1. ENVIRONMENT SETUP
############################

# Enter QIIME 2 Singularity container
cd /srv/scratch/singularity/
singularity shell qiime2.sif
# Alternative with explicit bind:
# singularity shell --bind /srv/scratch/singularity qiime2.sif

# Configure Matplotlib cache (required on shared HPC systems)
export MPLCONFIGDIR=/yourdirectory
mkdir -p /yourdirectory

# Configure temporary directory
export TMPDIR=/yourdirectory

# Move to working directory containing FASTQ files if not already there
cd /yourdirectory


############################
# 2. MANIFEST FILE CREATION
############################
# NOTE:
# - This uses the legacy QIIME 2 manifest format (comma-separated).
# - Newer QIIME 2 versions recommend tab-separated manifests.
# - Sample ID is inferred from filename prefix before first underscore.

# Forward reads
ls *_R1.fastq.gz | \
awk -v P="$(pwd)" -F'_' 'BEGIN{OFS=","}{print $1, P"/"$0, "forward"}' \
> forward_manifest_temp

# Reverse reads
ls *_R2.fastq.gz | \
awk -v P="$(pwd)" -F'_' 'BEGIN{OFS=","}{print $1, P"/"$0, "reverse"}' \
> reverse_manifest_temp

# Combine into paired-end manifest
cat <(echo sample-id,absolute-filepath,direction) \
    forward_manifest_temp \
    reverse_manifest_temp \
> paired-33-manifest


############################
# 3. IMPORT RAW READS
############################

qiime tools import \
  --input-path paired-33-manifest \
  --output-path paired-end-demux-vqPCR.qza \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33

# Generate demultiplexing summary
qiime demux summarize \
  --i-data paired-end-demux-vqPCR.qza \
  --o-visualization paired-end-demux-vqPCR.qzv


############################
# 4. PRIMER TRIMMING (CUTADAPT)
############################
# Primers:
# - Forward: 515F  (TGYCAGCMGCCGCGGTAA)
# - Reverse: 806R  (GACTACNVGGGTWTCTAAT)
# Strict matching (error-rate = 0)
# Reads without primers are discarded

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux-vqPCR.qza \
  --p-front-f TGYCAGCMGCCGCGGTAA \
  --p-front-r GACTACNVGGGTWTCTAAT \
  --p-error-rate 0 \
  --p-cores 52 \
  --p-discard-untrimmed \
  --o-trimmed-sequences trimmed16Sq.qza \
  --verbose

# Summarize trimmed reads
qiime demux summarize \
  --i-data trimmed16Sq.qza \
  --o-visualization trimmed16Sq_summary.qzv

############################
# 5. DENOISING WITH DADA2
############################
# Parameters chosen based on quality profile inspection
# Chimera filtering stringency increased via:
# --p-min-fold-parent-over-abundance 2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed16Sq.qza \
  --p-trunc-len-f 130 \
  --p-trunc-len-r 129 \
  --p-min-overlap 6 \
  --p-min-fold-parent-over-abundance 2 \
  --p-n-reads-learn 500000 \
  --p-n-threads 50 \
  --o-representative-sequences rep-seqs-dada2-16S.qza \
  --o-table table-dada2-16S.qza \
  --o-denoising-stats stats-dada2-16S.qza \
  --verbose

# Visualize denoising statistics
qiime metadata tabulate \
  --m-input-file stats-dada2-16S.qza \
  --o-visualization stats-dada2-16S.qzv

# Rename outputs for clarity
mv rep-seqs-dada2-16S.qza rep-seqs-16S.qza
mv table-dada2-16S.qza table-16S.qza


############################
# 6. FEATURE TABLE SUMMARIES
############################

qiime feature-table summarize \
  --i-table table-16S.qza \
  --o-visualization table-16S.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-16S.qza \
  --o-visualization rep-seqs-16S.qzv


############################
# 7. RAREFACTION
############################
# Sampling depth chosen based on minimum sequencing depth
# (example: conventional threshold ~127,940 reads)

qiime feature-table rarefy \
  --i-table table-16S.qza \
  --p-sampling-depth 490875 \
  --o-rarefied-table rarefied-table-16S.qza

qiime feature-table summarize \
  --i-table rarefied-table-16S.qza \
  --o-visualization rarefied-table-16S.qzv


############################
# 8. TAXONOMIC CLASSIFICATION
############################
# Reference database: SILVA 138 (99%), trimmed to 515–806 region

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs-515-806.qza \
  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
  --o-classifier silva-138-99--515-806-classifier2908.qza \
  --verbose

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99--515-806-classifier2908.qza \
  --i-reads rep-seqs-16S.qza \
  --p-n-jobs 48 \
  --o-classification taxonomy-16S.qza

qiime metadata tabulate \
  --m-input-file taxonomy-16S.qza \
  --o-visualization taxonomy-16S.qzv

# Taxonomic composition barplots
qiime taxa barplot \
  --i-table rarefied-table-16S.qza \
  --i-taxonomy taxonomy-16S.qza \
  --o-visualization taxa-barplot-16S.qzv


############################
# 9. PHYLOGENETIC TREE
############################

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-16S.qza \
  --p-n-threads 48 \
  --o-alignment aligned-rep-seqs-16S.qza \
  --o-masked-alignment masked-aligned-rep-seqs-16S.qza \
  --o-tree unrooted-tree-16S.qza \
  --o-rooted-tree rooted-tree-16S.qza \
  --verbose


############################
# 10. SAMPLE FILTERING
############################
# Keep only samples flagged as "oui" in metadata

qiime feature-table filter-samples \
  --i-table table-16S.qza \
  --m-metadata-file sample_metadata3.txt \
  --p-where "[A garder]='oui'" \
  --o-filtered-table filtered-table-16S.qza


############################
# 11. CORE DIVERSITY METRICS
############################
# Uses filtered but non-rarefied table
# Metadata file is mandatory

qiime diversity core-metrics-phylogenetic \
  --i-table filtered-table-16S.qza \
  --i-phylogeny rooted-tree-16S.qza \
  --m-metadata-file sample_metadata3.txt \
  --p-sampling-depth 127940 \
  --output-dir core-metrics-results-16S \
  --verbose


############################
# 12. ALPHA DIVERSITY
############################

# Rarefaction curves (max depth adjusted to visualize plateau)
qiime diversity alpha-rarefaction \
  --i-table table-16S.qza \
  --i-phylogeny rooted-tree-16S.qza \
  --p-max-depth 1000000 \
  --o-visualization alpha_rarefaction_curvesmax.qzv \
  --verbose

# Additional alpha diversity indices (computed on rarefied table)
for metric in chao1 observed_features simpson_e shannon
do
  qiime diversity alpha \
    --i-table rarefied-table-16S.qza \
    --p-metric ${metric} \
    --o-alpha-diversity alpha_diversity_${metric}.qza

  qiime tools export \
    --input-path alpha_diversity_${metric}.qza \
    --output-path alpha_diversity_${metric}_exported
done


############################
# 13. ALPHA DIVERSITY STATISTICS
############################

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-16S/faith_pd_vector.qza \
  --m-metadata-file sample_metadata3.txt \
  --o-visualization core-metrics-results-16S/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-16S/shannon_vector.qza \
  --m-metadata-file sample_metadata3.txt \
  --o-visualization core-metrics-results-16S/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-16S/evenness_vector.qza \
  --m-metadata-file sample_metadata3.txt \
  --o-visualization core-metrics-results-16S/evenness-group-significance.qzv


############################
# 14. BETA DIVERSITY STATISTICS
############################
# Distance metrics: Bray–Curtis, Jaccard, unweighted UniFrac
# PERMANOVA tests with pairwise comparisons

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-16S/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt \
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-16S/BrayCurtis-construction.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-16S/jaccard_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt\
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-16S/jaccard-construction.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-16S/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt\
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-16S/unWunifrac-construction.qzv \
  --p-pairwise


############################
# 15. EXPORTATION RAREFIED DATA TABLE
############################
qiime taxa collapse \
--i-table core-metrics-results/rarefied_table.qza \
--o-collapsed-table collapseL6table.qza \
--p-level 6 \
--i-taxonomy taxonomy.qza

qiime metadata tabulate \
  --m-input-file collapseL6table.qza \
  --o-visualization collapseL6table.qzv

# Exportation of relative frequences
qiime feature-table relative-frequency \
--i-table collapseL6table.qza \
--o-relative-frequency-table collapse_frequencyL6.table.qza 

qiime metadata tabulate \
  --m-input-file collapse_frequencyL6.table.qza \
  --o-visualization frequencycollapseL6table.qz
