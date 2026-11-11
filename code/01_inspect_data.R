# Stage 1: Dataset Inspection

# 1. Check working directory
getwd()

# 2. Get downloaded file names
files <- list.files("data/raw")

# 3. Display file names
print(files)

# 4. Check file sizes
file.info(paste0("data/raw/", files))$size

# 5. Inspect the first five lines of raw count data
connection <- gzfile(
  "data/raw/GSE97861_Raw_counts_cohort1.tsv.gz",
  "rt"
)
first_lines <- readLines(connection, n = 5)
close(connection)
print(first_lines)

# 6. Import raw gene expression data
raw_data <- read.table(
  gzfile("data/raw/GSE97861_Raw_counts_cohort1.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 7. Display the first six rows
head(raw_data)

# 8. Check dataset dimensions
dim(raw_data)

# 9. Check column names
colnames(raw_data)




# Stage 1.5: Basic Data Quality Checking

# Check missing values
sum(is.na(raw_data))

# Check duplicate gene names
sum(duplicated(raw_data$GENES))

# Check whether expression columns are numeric
sapply(raw_data[-1], is.numeric)

# Check for negative expression counts
any(as.matrix(raw_data[-1]) < 0)

# Check minimum and maximum expression counts
min(as.matrix(raw_data[-1]))
max(as.matrix(raw_data[-1]))

# Check total expression counts in each sample
colSums(raw_data[-1])

