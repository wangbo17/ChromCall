# chromcall full workflow manual test script
# Run manually for integration testing / debugging

rm(list = ls())

library(chromcall)

# --------------------------------------------------------
# STEP 1: Build chromcall sample A from BAM and annotation files
# --------------------------------------------------------
sampleA <- build_chromcall_sample(
  sample_name     = "SampleA",
  experiments     = list(
    H3K27me3 = file.path(system.file("extdata", package = "chromcall"), "h3k27me3_sampleA.bam"),
    H3K4me3  = file.path(system.file("extdata", package = "chromcall"), "h3k4me3_sampleA.bam"),
    Control  = file.path(system.file("extdata", package = "chromcall"), "control_sampleA.bam")
  ),
  control_name    = "Control",
  genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file     = system.file("extdata", "example.bed", package = "chromcall"),
  window_size     = 10000,
  blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall")
)

# Preview sampleA object structure
sampleA
SummarizedExperiment::assayNames(sampleA)
SummarizedExperiment::colData(sampleA)
head(SummarizedExperiment::assay(sampleA, "counts"))

# --------------------------------------------------------
# STEP 2: Perform Poisson-based statistical testing
# --------------------------------------------------------
sampleA_tested <- test_region_counts(sampleA)

# Preview logFC and adjusted p-values
SummarizedExperiment::assayNames(sampleA_tested)
SummarizedExperiment::colData(sampleA_tested)
head(SummarizedExperiment::assay(sampleA_tested, "counts"))
head(SummarizedExperiment::assay(sampleA_tested, "logFC"))
head(SummarizedExperiment::assay(sampleA_tested, "lambda_t"))
head(SummarizedExperiment::assay(sampleA_tested, "p_value"))
head(SummarizedExperiment::assay(sampleA_tested, "p_adj"))

# --------------------------------------------------------
# STEP 3: Export per-region results for H3K4me3
# --------------------------------------------------------
write_experiment_results(sampleA_tested, "H3K4me3", "SampleA_H3K4me3.tsv")

# Preview output file
head(read.delim("SampleA_H3K4me3.tsv"))

# --------------------------------------------------------
# STEP 4: Build a second sample for comparison (here duplicated for demonstration)
# --------------------------------------------------------
sampleB <- build_chromcall_sample(
  sample_name     = "SampleB",
  experiments     = list(
    H3K27me3 = file.path(system.file("extdata", package = "chromcall"), "h3k27me3_sampleB.bam"),
    H3K4me3  = file.path(system.file("extdata", package = "chromcall"), "h3k4me3_sampleB.bam"),
    Control  = file.path(system.file("extdata", package = "chromcall"), "control_sampleB.bam")
  ),
  control_name    = "Control",
  genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file     = system.file("extdata", "example.bed", package = "chromcall"),
  window_size     = 10000,
  blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall")
)

# Test sampleB
sampleB_tested <- test_region_counts(sampleB)

# --------------------------------------------------------
# STEP 5: Compare sampleA and sampleB
# --------------------------------------------------------
comparison <- compare_samples(sampleA_tested, sampleB_tested)

# Preview delta scores
comparison

# --------------------------------------------------------
# STEP 6: Export comparison results
# --------------------------------------------------------
write_comparison_results(comparison, "comparison_SampleA_vs_SampleB.tsv")

# Preview comparison output
head(read.delim("comparison_SampleA_vs_SampleB.tsv"))

# ======== EXPRESSION-AWARE WORKFLOW ========

# STEP 1e: Build chromcall sampleA with expression
sampleA_expr <- build_chromcall_sample(
  sample_name     = "SampleA",
  experiments     = list(
    H3K27me3 = file.path(system.file("extdata", package = "chromcall"), "h3k27me3_sampleA.bam"),
    H3K4me3  = file.path(system.file("extdata", package = "chromcall"), "h3k4me3_sampleA.bam"),
    Control  = file.path(system.file("extdata", package = "chromcall"), "control_sampleA.bam")
  ),
  control_name    = "Control",
  genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file     = system.file("extdata", "example.bed", package = "chromcall"),
  window_size     = 10000,
  blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall"),
  expression_file = system.file("extdata", "expression_sampleA.bed", package = "chromcall")
)

# STEP 2e: Run test_region_counts on sampleA_expr
sampleA_expr_tested <- test_region_counts(sampleA_expr)

# STEP 3e: Build expression-aware sampleB
sampleB_expr <- build_chromcall_sample(
  sample_name     = "SampleB",
  experiments     = list(
    H3K27me3 = file.path(system.file("extdata", package = "chromcall"), "h3k27me3_sampleB.bam"),
    H3K4me3  = file.path(system.file("extdata", package = "chromcall"), "h3k4me3_sampleB.bam"),
    Control  = file.path(system.file("extdata", package = "chromcall"), "control_sampleB.bam")
  ),
  control_name    = "Control",
  genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file     = system.file("extdata", "example.bed", package = "chromcall"),
  window_size     = 10000,
  blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall"),
  expression_file = system.file("extdata", "expression_sampleB.bed", package = "chromcall")
)

# STEP 4e: Run test_region_counts on sampleB_expr
sampleB_expr_tested <- test_region_counts(sampleB_expr)

# STEP 5e: Run expression-aware compare_samples
comparison_expr <- compare_samples(sampleA_expr_tested, sampleB_expr_tested)

# STEP 6e: Export comparison results with expression
write_comparison_results(comparison_expr, "comparison_SampleA_vs_SampleB_expression.tsv")

# STEP 7e: Preview expression-aware output
head(read.delim("comparison_SampleA_vs_SampleB_expression.tsv"))

