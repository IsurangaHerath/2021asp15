
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
