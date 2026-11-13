
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