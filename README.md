# ChromCall
*Assigning chromatin status to predefined genomic regions from epigenomic profiling data*

<p align="center">
  <img src="https://img.shields.io/badge/R-4.3.2+-blue" />
  <img src="https://img.shields.io/badge/License-MIT-green" />
  <img src="https://img.shields.io/badge/Version-1.0.0-orange" />
  <img src="https://img.shields.io/github/last-commit/wangbo17/ChromCall" />
</p>

---

## 🔍 Overview

**ChromCall** is an R package for **region-based chromatin enrichment analysis** of epigenomic profiling data, including **ChIP-seq, CUT&RUN, CUT&Tag, and ATAC-seq**. It provides a transparent, statistically principled framework to **quantify enrichment at predefined genomic regions** (e.g. promoters or enhancers), enabling region-matched comparisons across samples and experiments without relying on data-dependent peak boundaries.

---

## 🚀 Key Features

- **Region-centric analysis**  
  Quantifies chromatin enrichment directly within predefined genomic windows, enabling consistent, region-matched comparisons across samples and experiments.

- **Transparent statistical framework**  
  Employs a Poisson-based background model incorporating:
  - experiment-specific genome-wide background estimation  
  - region-specific, control-derived modulation factors

- **Control-aware enrichment testing**  
  Explicitly integrates matched control experiments to account for both technical background and local biological variability.

- **Multiple complementary metrics**  
  For each region and experiment, ChromCall reports:
  - FDR-adjusted p-values  
  - enrichment score (log₂ observed / expected)  
  - Poisson z-score  
  - binary chromatin status (present / absent)

- **Multi-experiment and multi-sample support**  
  Supports joint analysis of multiple chromatin marks and pairwise comparisons between samples within a unified framework.

- **Optional expression integration**  
  Region-level gene expression values (e.g. TSS-associated expression) can be incorporated to enable integrated chromatin–transcription analyses.

---

## 🧠 Statistical Model Overview

ChromCall models read counts as **Poisson-distributed events**, an appropriate approximation for sparse and independent fragment occurrences across fixed genomic windows. After accounting for genome-wide background signal and local, control-derived modulation, this framework enables transparent and analytically tractable region-level inference.

### Background Estimation

For each experiment, a genome-wide background rate \( \lambda_g \) is estimated as the mean read count per non-blacklisted genomic tile:

$$
\lambda_g = \frac{1}{N} \sum_{i=1}^{N} y_i
$$

where \( y_i \) denotes the read count in the *i*th tile and *N* is the total number of non-blacklisted tiles.  
Zero-count tiles are retained by default to avoid upward bias in sparse datasets and to ensure that \( \lambda_g \) reflects global background signal rather than local enrichment.

### Control-based Local Modulation

To account for region-specific biological variability, ChromCall derives a modulation factor from the matched control experiment:

$$
m_i = \max\left(1, \frac{y_i^{(\mathrm{ctrl})}}{\lambda_g^{(\mathrm{ctrl})}}\right)
$$

The expected signal for region *i* in experiment *j* is then defined as:

$$
\lambda_{t,i}^{(j)} = m_i \times \lambda_g^{(j)}
$$

This formulation integrates global sequencing depth with local control variation into a unified and interpretable background model, while preventing deflation in regions with low control signal.

### Statistical Testing and Effect Sizes

ChromCall evaluates region-level enrichment using a one-sided Poisson test:

$$
p_i^{(j)} = P\left(Y \ge y_i^{(j)} \mid Y \sim \mathrm{Pois}(\lambda_{t,i}^{(j)})\right)
$$

Multiple testing correction is applied across all regions using the **Benjamini–Hochberg false discovery rate (FDR)** procedure.

In addition to significance testing, ChromCall reports complementary effect-size metrics:

- **Enrichment score**

$$
s_i^{(j)} = \log_2\left(\frac{y_i^{(j)} + \epsilon}{\lambda_{t,i}^{(j)} + \epsilon}\right)
$$

- **Poisson z-score**

$$
z_i^{(j)} = \frac{y_i^{(j)} - \lambda_{t,i}^{(j)}}{\sqrt{\lambda_{t,i}^{(j)}}}
$$

Together, these metrics provide complementary measures of enrichment strength, effect size, and statistical confidence.

---
## 🧬 Implementation and Data Structures

ChromCall is implemented in **R** and builds upon the **Bioconductor** ecosystem, ensuring interoperability with standard genomic data structures and downstream analysis workflows:

- `GRanges` for representing genomic intervals  
- `SummarizedExperiment` for storing structured assay outputs and metadata  
- `GenomicAlignments` for importing aligned sequencing reads from BAM files  
- `GenomeInfoDb` and `Seqinfo` for genome annotation and consistency checks  

Each processed sample is returned as a `SummarizedExperiment` object containing:

- raw region-level read counts  
- genome-wide and locally adjusted background estimates (\( \lambda_g \), \( \lambda_t \))  
- p-values and FDR-adjusted p-values  
- enrichment scores and Poisson z-scores  

Pairwise sample comparisons generate region-level **Δ enrichment** and **Δ z-score** metrics, enabling direct comparative analysis of chromatin states across biological conditions.

---

## 📦 Installation

ChromCall is available as a development version on GitHub and can be installed using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("GliomaGenomics/ChromCall")
```

---

## 🧪 Basic Workflow

### Build a ChromCall sample

```r
sample <- build_chromcall_sample(
  sample_name   = "sampleA",
  experiments   = list(
    H3K27me3 = "h3k27me3.bam",
    H3K4me3  = "h3k4me3.bam",
    Control  = "control.bam"
  ),
  control_name   = "Control",
  genome_file    = "genome.txt",
  region_file    = "promoters.bed",
  window_size    = 2000,
  blacklist_file = "blacklist.bed",
  expression_file = "expression_tss.bed"
)
```

### Perform region-level enrichment testing

```r
result <- test_region_counts(sample)
```

### Compare two samples

```r
comparison <- compare_samples(resultA, resultB, threshold = 0.05)
```

### Export results

```r
write_experiment_results(result, "H3K4me3", "results.tsv")
write_comparison_results(comparison, "comparison.tsv")
```

---

## 📈 Outputs

| Metric                           | Description                                   |
| -------------------------------- | --------------------------------------------- |
| `counts`                         | Raw read count per region                     |
| `lambda_g`                       | Genome-wide background rate                   |
| `lambda_t`                       | Locally adjusted expected signal              |
| `p_value`, `p_adj`               | Poisson test p-values and FDR-adjusted values |
| `score`                          | log₂(Observed / Expected) enrichment          |
| `z_pois`                         | Poisson z-score                               |
| `DeltaEnrichment`, `DeltaZscore` | Pairwise comparison metrics                   |

---

## 💡 Contact

For questions, issues, or feature requests, please open a 👉 [GitHub issue](https://github.com/wangbo17/ChromCall/issues)

