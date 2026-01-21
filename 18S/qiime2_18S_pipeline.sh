###############################################################################
# QIIME 2 18S rRNA amplicon analysis pipeline
# Environment: Singularity + QIIME 2
# Data type: Paired-end Illumina FASTQ (Phred33)
# Target region: 18S rRNA (V9)
###############################################################################

############################
# 1. ENVIRONMENT SETUP
############################

# Enter QIIME 2 Singularity container
cd /yourdirectory/
singularity shell qiime2.sif

# Configure temporary directory
export TMPDIR=/yourdirectory/

# Move to working directory containing FASTQ files
cd //yourdirectory/


############################
# 2. MANIFEST FILE CREATION
############################
# NOTE:
# - Legacy QIIME 2 manifest format (comma-separated)
# - Sample ID inferred from filename prefix before first underscore

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
  --output-path paired-end-demux-v2.qza \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33

# Generate demultiplexing summary
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv


############################
# 4. PRIMER TRIMMING (CUTADAPT)
############################
# Primers (second set, V9 region)
# - Forward: TACACACCGCCCGTC
# - Reverse: TGATCCTTCYGCAGGTTCACCTAC
# Strict matching (error-rate = 0), discard untrimmed reads

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux-vqPCR.qza \
  --p-front-f TACACACCGCCCGTC \
  --p-front-r TGATCCTTCYGCAGGTTCACCTAC \
  --p-error-rate 0 \
  --p-cores 48 \
  --p-discard-untrimmed \
  --o-trimmed-sequences trimmed18SYq.qza \
  --verbose

# Summarize trimmed reads
qiime demux summarize \
  --i-data trimmed18SYq.qza \
  --o-visualization trimmed18SYq_summary.qzv

############################
# 5. DENOISING WITH DADA2
############################
# Parameters chosen based on quality profiles and chimera stringency

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed18SYq.qza \
  --p-trunc-len-f 120 \
  --p-trunc-len-r 90 \
  --p-min-fold-parent-over-abundance 5 \
  --p-n-reads-learn 500000 \
  --p-n-threads 65 \
  --o-representative-sequences rep-seqs-dada2-18S5.qza \
  --o-table table-dada2-18S5.qza \
  --o-denoising-stats stats-dada2-18S5.qza \
  --verbose

# Visualize denoising stats
qiime metadata tabulate \
  --m-input-file stats-dada2-18S5.qza \
  --o-visualization stats-dada2-18S5.qzv

# Rename outputs for clarity
mv rep-seqs-dada2-18S5.qza rep-seqs-18S5.qza
mv table-dada2-18S5.qza table-18S5.qza

# Feature table summaries
qiime feature-table summarize \
  --i-table table-18S5.qza \
  --o-visualization table-18S5.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-18S5.qza \
  --o-visualization rep-seqs-18S5.qzv


############################
# 6. RAREFACTION
############################
# Sampling depth chosen based on minimum sequencing depth in dataset

qiime feature-table rarefy \
  --i-table table-18S5.qza \
  --p-sampling-depth 199417 \
  --o-rarefied-table rarefied-table-18S5.qza

qiime feature-table summarize \
  --i-table rarefied-table-18S5.qza \
  --o-visualization rarefied-table-18S5.qzv


############################
# 7. TAXONOMIC CLASSIFICATION
############################
# SILVA V9 classifier

# Extract V9 region from SILVA
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer TACACACCGCCCGTC \
  --p-r-primer TGATCCTTCYGCAGGTTCACCTAC \
  --p-trunc-len 0 \
  --o-reads silva-138-99-V9-seqs.qza \
  --verbose

# Train Naive Bayes classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-V9-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-classifier-V9-CJ.qza \
  --verbose

# Assign taxonomy with confidence 0.5 to increase assignation
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-classifier-V9-CJ.qza \
  --i-reads rep-seqs-18S5.qza \
  --p-n-jobs 48 \
  --p-confidence 0.5 \
  --o-classification taxonomy-18S-05.qza \
  --verbose

# Visualize taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy-18S-05.qza \
  --o-visualization taxonomy-18S-05.qzv

# Filter to retain only Eukaryotes
qiime taxa filter-table \
  --i-table rarefied-table-18S5.qza \
  --i-taxonomy taxonomy-18S-05.qza \
  --p-include Eukaryota \
  --o-filtered-table filtered-table-18S-05.qza

# Taxonomic barplot
qiime taxa barplot \
  --i-table filtered-table-18S-05.qza \
  --i-taxonomy taxonomy-18S-05.qza \
  --o-visualization rartaxa-bar-plot-05.qzv


############################
# 8. ALTERNATIVE CLASSIFIERS
############################
# Consensus-based classification (vsearch) using SILVA or PR2 references

# SILVA consensus
qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seqs-18S5.qza \
  --i-reference-reads silva-138-99-V9-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --output-dir vsearch \
  --p-threads 48 \
  --verbose

qiime taxa filter-table \
  --i-table rarefied-table-18S5.qza \
  --i-taxonomy vsearch/classification.qza \
  --p-include Eukaryota \
  --o-filtered-table filtered-table-18S-cons.qza

qiime taxa barplot \
  --i-table filtered-table-18S-cons.qza \
  --i-taxonomy vsearch/classification.qza \
  --o-visualization rartaxa-bar-plot-cons.qzv

# PR2 reference classifier
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path pr2_version_4.13.0_18S_mothur.fasta \
  --output-path pr2_v4.13.0.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path pr2_version_4.13.0_18S_mothur.tax \
  --output-path pr2_v4.13.0_tax.qza

# Extract V9 region from PR2
qiime feature-classifier extract-reads \
  --i-sequences pr2_v4.13.0.qza \
  --p-f-primer TACACACCGCCCGTC \
  --p-r-primer TGATCCTTCYGCAGGTTCACCTAC \
  --o-reads pr2_v4.13.0_v9_extracts.qza

# Train PR2 classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads pr2_v4.13.0_v9_extracts.qza \
  --i-reference-taxonomy pr2_v4.13.0_tax.qza \
  --o-classifier pr2_v4.13.0_v9_classifier.qza

# Assign taxonomy using PR2
qiime feature-classifier classify-sklearn \
  --p-n-jobs 48 \
  --i-classifier pr2_v4.13.0_v9_classifier.qza \
  --i-reads rep-seqs-18S5.qza \
  --o-classification taxonomy-18S-PR2.qza

qiime metadata tabulate \
  --m-input-file taxonomy-18S-PR2.qza \
  --o-visualization taxonomy-18S-PR2.qzv

# Filter Eukaryotes and visualize
qiime taxa filter-table \
  --i-table rarefied-table-18S5.qza \
  --i-taxonomy taxonomy-18S-PR2.qza \
  --p-include Eukaryota \
  --o-filtered-table filtered-table-18S-pr2.qza

qiime taxa barplot \
  --i-table filtered-table-18S-pr2.qza \
  --i-taxonomy taxonomy-18S-PR2.qza \
  --o-visualization rartaxa-bar-plot-pr2.qzv


############################
# 9. PHYLOGENY
############################

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-18S5.qza \
  --p-n-threads 48 \
  --o-alignment aligned-rep-seqs-18S.qza \
  --o-masked-alignment masked-aligned-rep-seqs-18S.qza \
  --o-tree unrooted-tree-18S.qza \
  --o-rooted-tree rooted-tree-18S.qza \
  --verbose


############################
# 10. SAMPLE FILTERING & CORE DIVERSITY
############################

qiime feature-table filter-samples \
  --i-table table-18S5.qza \
  --m-metadata-file sample_metadata3.txt \
  --p-where "[A garder]='oui'" \
  --o-filtered-table filtered-table-18S.qza

qiime diversity core-metrics-phylogenetic \
  --i-table filtered-table-18S.qza \
  --i-phylogeny rooted-tree-18S.qza \
  --m-metadata-file sample_metadata3.txt \
  --p-sampling-depth 199417 \
  --output-dir core-metrics-results-18S \
  --verbose


############################
# 11. ALPHA DIVERSITY
############################

# Rarefaction curves
qiime diversity alpha-rarefaction \
  --i-table table-18S5.qza \
  --i-phylogeny rooted-tree-18S.qza \
  --p-max-depth 400000 \
  --o-visualization alpha_rarefaction_curves_max.qzv \
  --verbose

# Additional alpha metrics
for metric in chao1 observed_features simpson_e shannon
do
  qiime diversity alpha \
    --i-table rarefied-table-18S5.qza \
    --p-metric ${metric} \
    --o-alpha-diversity alpha_diversity_${metric}.qza

  qiime tools export \
    --input-path alpha_diversity_${metric}.qza \
    --output-path alpha_diversity_${metric}_exported
done

# Statistical testing
for metric_file in faith_pd_vector shannon_vector evenness_vector
do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results-18S/${metric_file}.qza \
    --m-metadata-file sample_metadata3.txt \
    --o-visualization core-metrics-results-18S/${metric_file}-group-significance.qzv
done


############################
# 12. BETA DIVERSITY
############################

# Bray–Curtis, Jaccard, unweighted UniFrac (pairwise PERMANOVA)

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-18S/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt \
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-18S/BrayCurtis-sample-construction.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-18S/jaccard_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt \
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-18S/jaccard-construction.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-18S/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample_metadata2.txt \
  --m-metadata-column 'Construction' \
  --o-visualization core-metrics-results-18S/unWunifrac-Construction.qzv \
  --p-pairwise

############################
# 13. TAXONOMIC COLLAPSE & RELATIVE FREQUENCIES
############################

# Collapse tables to different taxonomic levels
for lvl in 2 5 6
do
  qiime taxa collapse \
    --i-table core-metrics-results/rarefied-table-18S5.qza \
    --p-level $lvl \
    --i-taxonomy taxonomy.qza \
    --o-collapsed-table collapseL${lvl}table.qza

  qiime metadata tabulate \
    --m-input-file collapseL${lvl}table.qza \
    --o-visualization collapseL${lvl}table.qzv
done

# Relative frequency for collapsed L6 table
qiime feature-table relative-frequency \
  --i-table collapseL6table.qza \
  --o-relative-frequency-table collapse_frequencyL6.table.qza

qiime metadata tabulate \
  --m-input-file collapse_frequencyL6.table.qza \
  --o-visualization frequencycollapseL6table.qzv
