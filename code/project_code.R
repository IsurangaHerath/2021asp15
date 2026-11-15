
# Stage 1: Dataset Inspection

getwd()

files <- list.files("data/raw")
print(files)

file.info(paste0("data/raw/", files))$size

connection <- gzfile(
  "data/raw/GSE97861_Raw_counts_cohort1.tsv.gz",
  "rt"
)
first_lines <- readLines(connection, n = 5)
close(connection)
print(first_lines)

raw_data <- read.table(
  gzfile("data/raw/GSE97861_Raw_counts_cohort1.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

head(raw_data)
dim(raw_data)
colnames(raw_data)


# Stage 2: Basic Data Quality Checking

sum(is.na(raw_data))
sum(duplicated(raw_data$GENES))

sapply(raw_data[-1], is.numeric)

any(as.matrix(raw_data[-1]) < 0)

min(as.matrix(raw_data[-1]))
max(as.matrix(raw_data[-1]))

colSums(raw_data[-1])


# Stage 3: GEO Sample Metadata Inspection

connection <- gzfile(
  "data/raw/GSE97861_series_matrix.txt.gz",
  "rt"
)

metadata_lines <- readLines(connection, n = 30)
close(connection)

cat(metadata_lines, sep = "\n")


# Extract individual sample information

connection <- gzfile(
  "data/raw/GSE97861_series_matrix.txt.gz",
  "rt"
)

all_metadata <- readLines(connection)
close(connection)

sample_metadata <- grep(
  "^!Sample_",
  all_metadata,
  value = TRUE
)

cat(sample_metadata, sep = "\n")


# Create sample information table

get_sample_values <- function(field) {
  line <- grep(
    paste0("^", field, "\t"),
    all_metadata,
    value = TRUE
  )
  
  values <- strsplit(line, "\t")[[1]][-1]
  
  gsub('^"|"$', "", values)
}

sample_table <- data.frame(
  GEO_ID = get_sample_values("!Sample_geo_accession"),
  Sample = get_sample_values("!Sample_title"),
  stringsAsFactors = FALSE
)


# Add group information

sample_table$Group <- ifelse(
  grepl("Negative", sample_table$Sample),
  "Negative",
  "Secondary"
)


# Add cell type information

sample_table$Cell_Type <- sub(
  " .*",
  "",
  sub("^[^_]*_", "", sample_table$Sample)
)


# Add donor information

sample_table$Donor <- sub(
  "_.*",
  "",
  sample_table$Sample
)


# Display sample information

print(sample_table)


# Check sample counts

table(sample_table$Group)
table(sample_table$Cell_Type)

write.csv(
  sample_table,
  "results/sample_information.csv",
  row.names = FALSE
)

# Check sample names in expression data

print(colnames(raw_data)[-1])

# Check sample names in metadata

print(sample_table$Sample)

# Create common sample identifiers

expression_samples <- sub(
  "_RNASeq$",
  "",
  colnames(raw_data)[-1]
)

metadata_samples <- sub(
  " _.*",
  "",
  sample_table$Sample
)

# Check whether sample names match

all(expression_samples == metadata_samples)

# Explore gene expression data

# Summary of gene expression values
summary(as.matrix(raw_data[-1]))

# Calculate average expression for each sample
sample_means <- colMeans(raw_data[-1])

print(sample_means)

# Visualize sample expression distributions
boxplot(
  raw_data[-1],
  las = 2,
  main = "Gene Expression Distribution Across Samples",
  ylab = "Raw Read Counts",
  outline = FALSE
)

# Visualize sequencing depth

sample_totals <- colSums(raw_data[-1])

barplot(
  sample_totals,
  las = 2,
  main = "Sequencing Depth Across Samples",
  ylab = "Total Read Counts",
  col = "lightblue"
)

# Check low-expression genes

gene_counts <- raw_data[-1]

low_expression_genes <- rowSums(gene_counts >= 10) < 3

sum(low_expression_genes)

sum(!low_expression_genes)

# Import TPM data

tpm_data <- read.table(
  gzfile("data/raw/GSE97861_TPM_cohort1.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Inspect TPM data

dim(tpm_data)
head(tpm_data)

# Check TPM data quality

sum(is.na(tpm_data))
sum(duplicated(tpm_data$GENES))

# Check numeric expression values

sapply(tpm_data[-1], is.numeric)

# Check for negative values

any(as.matrix(tpm_data[-1]) < 0)

# Compare genes between datasets

identical(raw_data$GENES, tpm_data$GENES)

# Compare sample names between datasets

identical(colnames(raw_data), colnames(tpm_data))

# Stage 4: TPM Data Quality Checking

# Check missing values
sum(is.na(tpm_data))

# Check duplicate genes
sum(duplicated(tpm_data$GENES))

# Check numeric expression columns
sapply(tpm_data[-1], is.numeric)

# Check negative values
any(as.matrix(tpm_data[-1]) < 0)

# Compare genes between datasets
identical(raw_data$GENES, tpm_data$GENES)

# Compare sample names
identical(colnames(raw_data), colnames(tpm_data))

# Stage 5: TPM Expression Exploration

# Summary of TPM values
summary(as.matrix(tpm_data[-1]))

# Calculate average TPM for each sample
tpm_means <- colMeans(tpm_data[-1])

# Display average TPM values
print(tpm_means)


# Stage 6: TPM Expression Visualization

# Convert TPM values for visualization
log_tpm <- log2(as.matrix(tpm_data[-1]) + 1)

# Create a boxplot
boxplot(
  log_tpm,
  las = 2,
  outline = FALSE,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "log2(TPM + 1)",
  col = "lightblue"
)

# Save TPM expression boxplot

png(
  "results/tpm_expression_boxplot.png",
  width = 1200,
  height = 800,
  res = 120
)

boxplot(
  log_tpm,
  las = 2,
  outline = FALSE,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "log2(TPM + 1)",
  col = "lightblue"
)

dev.off()

# Stage 7: Inspect Sample Groups

# Display sample information columns
colnames(sample_table)

# Display sample information
head(sample_table)

# Count samples in each group
table(sample_table$Group)

# Stage 8: Sample Distribution by Cell Type

# Count samples by cell type and group
table(sample_table$Cell_Type, sample_table$Group)

# Stage 9: Standardize Cell Type Names

# Convert all cell type names to uppercase
sample_table$Cell_Type <- toupper(sample_table$Cell_Type)

# Check the updated distribution
table(sample_table$Cell_Type, sample_table$Group)

# Stage 10: Select TEMRA Samples

# Select only TEMRA samples
temra_samples <- sample_table[
  sample_table$Cell_Type == "TEMRA",
]

# Display selected samples
print(temra_samples)

# Count samples in each group
table(temra_samples$Group)

# Check total number of TEMRA samples
nrow(temra_samples)


# Stage 11: Extract TEMRA Raw Counts

# Create the expression column names
temra_sample_names <- paste0(
  temra_samples$Donor,
  "_TEMRA_RNASeq"
)

# Check whether all samples exist in raw counts
temra_sample_names %in% colnames(raw_data)

# Extract TEMRA expression data
temra_counts <- raw_data[
  ,
  c("GENES", temra_sample_names)
]

# Check extracted data dimensions
dim(temra_counts)

# Display first six rows
head(temra_counts)

# Stage 12: Prepare TEMRA Sample Information

# Create sample information
temra_metadata <- temra_samples[
  match(
    temra_sample_names,
    paste0(temra_samples$Donor, "_TEMRA_RNASeq")
  ),
]

# Select required columns
temra_metadata <- temra_metadata[
  ,
  c("Donor", "Group", "Cell_Type")
]

# Display sample information
print(temra_metadata)

# Check group distribution
table(temra_metadata$Group)

# Stage 13: Final TEMRA Data Validation

# Check sample names and metadata order
identical(
  colnames(temra_counts)[-1],
  paste0(temra_metadata$Donor, "_TEMRA_RNASeq")
)

# Check missing values
sum(is.na(temra_counts))

# Check duplicate genes
sum(duplicated(temra_counts$GENES))

# Check whether all expression values are integers
all(
  as.matrix(temra_counts[-1]) ==
    round(as.matrix(temra_counts[-1]))
)

# Check group distribution
table(temra_metadata$Group)

# Stage 14: Prepare Data for DESeq2

# Create the count matrix
count_matrix <- as.matrix(temra_counts[-1])

# Add gene names as row names
rownames(count_matrix) <- temra_counts$GENES

# Create sample information
col_data <- temra_metadata

# Set sample names as row names
rownames(col_data) <- paste0(
  col_data$Donor,
  "_TEMRA_RNASeq"
)

# Set Negative as the reference group
col_data$Group <- factor(
  col_data$Group,
  levels = c("Negative", "Secondary")
)

# Check count matrix dimensions
dim(count_matrix)

# Check sample information
print(col_data)

# Check sample order
identical(colnames(count_matrix), rownames(col_data))


# Stage 15: Install and Check DESeq2

# Install DESeq2 only if it is not already installed
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  
  # Install BiocManager if necessary
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install DESeq2
  BiocManager::install("DESeq2")
}

# Check DESeq2 installation
requireNamespace("DESeq2", quietly = TRUE)

print(requireNamespace("DESeq2", quietly = TRUE))


# Stage 16: Create DESeq2 Dataset

# Load DESeq2
library(DESeq2)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ Group
)

# Check the DESeq2 dataset
dds

# Check the number of genes and samples
dim(dds)

# Check the group distribution
table(dds$Group)

# Stage 17: Filter Low-Expression Genes

# Count genes before filtering
nrow(dds)

# Keep genes with at least 10 counts in at least 3 samples
keep_genes <- rowSums(counts(dds) >= 10) >= 3

# Apply filtering
dds <- dds[keep_genes, ]

# Count genes after filtering
nrow(dds)

# Calculate removed genes
23368 - nrow(dds)

# Stage 18: Differential Gene Expression Analysis

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract results: Secondary compared with Negative
res <- results(
  dds,
  contrast = c("Group", "Secondary", "Negative")
)

# Display summary of results
summary(res)

# Sort results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Display the first 10 results
head(res_ordered, 10)

# Stage 19: Identify Significant Genes

# Select genes with adjusted p-value below 0.05
significant_genes <- res[
  !is.na(res$padj) & res$padj < 0.05,
]

# Sort significant genes by adjusted p-value
significant_genes <- significant_genes[
  order(significant_genes$padj),
]

# Count significant genes
nrow(significant_genes)

# Display significant genes
print(significant_genes)

# Stage 20: Save Differential Expression Results

# Convert DESeq2 results into data frames
all_results <- as.data.frame(res_ordered)

significant_results <- as.data.frame(significant_genes)

# Save all gene results
write.csv(
  all_results,
  "results/differential_expression_all_genes.csv",
  row.names = TRUE
)

# Save significant gene results
write.csv(
  significant_results,
  "results/significant_genes.csv",
  row.names = TRUE
)

# Confirm saved files
file.exists("results/differential_expression_all_genes.csv")
file.exists("results/significant_genes.csv")

# Stage 21: Create Volcano Plot

# Convert results into a data frame
plot_data <- as.data.frame(res)

# Remove genes with missing adjusted p-values
plot_data <- plot_data[!is.na(plot_data$padj), ]

# Identify significant genes
plot_data$Significance <- ifelse(
  plot_data$padj < 0.05,
  "Significant",
  "Not Significant"
)

# Create volcano plot
plot(
  plot_data$log2FoldChange,
  -log10(plot_data$padj),
  col = ifelse(
    plot_data$Significance == "Significant",
    "red",
    "grey"
  ),
  pch = 20,
  main = "Volcano Plot of Differential Gene Expression",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value"
)

# Add significance threshold
abline(
  h = -log10(0.05),
  col = "blue",
  lty = 2
)

# Add vertical reference line
abline(
  v = 0,
  col = "black",
  lty = 2
)

# Add legend
legend(
  "topright",
  legend = c("Significant", "Not Significant"),
  col = c("red", "grey"),
  pch = 20
)

# Stage 22: Save Volcano Plot

png(
  "results/volcano_plot.png",
  width = 1200,
  height = 900,
  res = 120
)

# Create volcano plot
plot(
  plot_data$log2FoldChange,
  -log10(plot_data$padj),
  col = ifelse(
    plot_data$Significance == "Significant",
    "red",
    "grey"
  ),
  pch = 20,
  main = "Volcano Plot of Differential Gene Expression",
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value"
)

# Add threshold lines
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = 0, col = "black", lty = 2)

# Add legend
legend(
  "topright",
  legend = c("Significant", "Not Significant"),
  col = c("red", "grey"),
  pch = 20
)

# Close the image device
dev.off()

# Check saved file
file.exists("results/volcano_plot.png")

# Stage 23: Visualize Top 10 Significant Genes

# Select the top 10 significant genes
top_genes <- head(significant_genes, 10)

# Create bar plot
barplot(
  top_genes$log2FoldChange,
  names.arg = rownames(top_genes),
  las = 2,
  main = "Top 10 Differentially Expressed Genes",
  xlab = "Genes",
  ylab = "Log2 Fold Change",
  col = "lightblue"
)

# Add reference line at zero
abline(h = 0, col = "red", lty = 2)

# Stage 24: Save Top 10 Gene Plot

png(
  "results/top_10_significant_genes.png",
  width = 1200,
  height = 900,
  res = 120
)

# Recreate the bar plot
barplot(
  top_genes$log2FoldChange,
  names.arg = rownames(top_genes),
  las = 2,
  main = "Top 10 Differentially Expressed Genes",
  xlab = "Genes",
  ylab = "Log2 Fold Change",
  col = "lightblue"
)

abline(h = 0, col = "red", lty = 2)

# Close image device
dev.off()

# Verify saved image
file.exists("results/top_10_significant_genes.png")

# Stage 25: Prepare Significant Gene Summary

# Create a summary table
gene_summary <- data.frame(
  Gene = rownames(significant_genes),
  Log2FoldChange = significant_genes$log2FoldChange,
  Adjusted_P_Value = significant_genes$padj
)

# Identify expression direction
gene_summary$Expression <- ifelse(
  gene_summary$Log2FoldChange > 0,
  "Higher in Secondary",
  "Lower in Secondary"
)

# Display all significant genes
print(gene_summary)

# Count upregulated and downregulated genes
table(gene_summary$Expression)

# Stage 26: Save Significant Gene Summary

write.csv(
  gene_summary,
  "results/significant_gene_summary.csv",
  row.names = FALSE
)

# Verify saved file
file.exists("results/significant_gene_summary.csv")

# Stage 27: Examine CD160 Expression

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Extract CD160 expression
cd160_counts <- normalized_counts["CD160", ]

# Combine expression with sample information
cd160_data <- data.frame(
  Sample = names(cd160_counts),
  Expression = as.numeric(cd160_counts),
  Group = col_data$Group
)

# Display CD160 expression
print(cd160_data)

# Stage 28: Visualise CD160 Expression

# Create log2-transformed expression values
cd160_data$Log2_Expression <- log2(cd160_data$Expression + 1)

# Create a boxplot
boxplot(
  Log2_Expression ~ Group,
  data = cd160_data,
  main = "CD160 Expression in TEMRA Cells",
  xlab = "Group",
  ylab = "Log2 Normalized Expression",
  col = c("lightblue", "lightgreen")
)

# Display individual sample values
stripchart(
  Log2_Expression ~ Group,
  data = cd160_data,
  vertical = TRUE,
  method = "jitter",
  add = TRUE,
  pch = 19
)

# Stage 28: Save CD160 Expression Data

write.csv(
  cd160_data,
  "results/CD160_expression.csv",
  row.names = FALSE
)

# Check whether the file was saved
file.exists("results/CD160_expression.csv")

# Stage 29: Final Analysis Summary

cat("===== Differential Expression Analysis Summary =====\n")

cat("Total samples:", ncol(dds), "\n")
cat("Negative samples:", sum(col_data$Group == "Negative"), "\n")
cat("Secondary samples:", sum(col_data$Group == "Secondary"), "\n")

cat("Genes analysed:", nrow(dds), "\n")
cat("Significant genes (adjusted p-value < 0.05):",
    nrow(significant_genes), "\n")

cat("Higher in Secondary:",
    sum(significant_genes$log2FoldChange > 0), "\n")

cat("Lower in Secondary:",
    sum(significant_genes$log2FoldChange < 0), "\n")


# Stage 30: Principal Component Analysis (PCA)

# Load required package
library(DESeq2)
library(ggplot2)

# Transform count data for visualisation
vsd <- vst(dds, blind = FALSE)

# Create PCA plot
pca_plot <- plotPCA(
  vsd,
  intgroup = "Group"
)

# Display the PCA plot
print(pca_plot)

# Save the PCA plot
ggsave(
  filename = "results/PCA_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Check whether the plot was saved
file.exists("results/PCA_plot.png")
