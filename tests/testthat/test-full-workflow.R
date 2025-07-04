test_that("Full chromcall workflow runs without errors", {
  expect_error(
    local({
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

      sampleA_tested <- test_region_counts(sampleA)
      write_experiment_results(sampleA_tested, "H3K4me3", tempfile())

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

      sampleB_tested <- test_region_counts(sampleB)

      comparison <- compare_samples(sampleA_tested, sampleB_tested)
      write_comparison_results(comparison, tempfile())
    }),
    NA  # This means: expect NO error
  )
})
