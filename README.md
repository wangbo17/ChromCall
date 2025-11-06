# ChromCall
*An R package for region-level chromatin enrichment analysis of epigenomic data*


<p align="center">
  <img src="https://img.shields.io/badge/R-4.3.2+-blue" />
  <img src="https://img.shields.io/badge/License-MIT-green" />
  <img src="https://img.shields.io/badge/Status-Beta-orange" />
  <img src="https://img.shields.io/github/last-commit/wangbo17/ChromCall" />
</p>

## 🔍 Overview

**ChromCall** is an R package for *region-level chromatin enrichment analysis* of epigenomic data such as **ChIP-seq**, **CUT&RUN**, **CUT&Tag**, and **ATAC-seq**. It provides a **transparent**, **reproducible**, and **statistically principled** framework to quantify enrichment across predefined genomic regions (e.g. promoters, enhancers, or transcription factor binding sites).  

ChromCall implements a **Poisson-based background model** that integrates genome-wide signal estimation with local modulation derived from matched control data. This unified approach enables consistent and interpretable quantification of chromatin enrichment both within and across experiments.

## 🚀 Key Features

- **Region-centered analysis:** Quantify enrichment within predefined genomic windows (e.g. promoters, enhancers).  
- **Transparent statistical model:** Uses a Poisson-based framework with genome-wide background and control-derived local modulation.  
- **Reproducible quantification:** Provides comparable enrichment estimates across experiments.  
- **Comprehensive outputs:** Effect sizes, FDR-adjusted *p*-values, enrichment scores, and Poisson *z*-scores.  
- **Differential comparison:** Supports pairwise delta-based metrics to assess chromatin pattern changes between samples.  
- **Bioconductor integration:** Fully compatible with GRanges, SummarizedExperiment, and downstream workflows.  
- **Modular and extensible:** Designed for future extension to additional statistical models and multi-sample designs.

## 🧠 Method Overview

ChromCall models read counts as Poisson-distributed events, suitable for sparse and independent fragment occurrences across fixed genomic windows.  
Read counts within each region are assumed to arise from independent sampling with a constant underlying rate, forming a tractable probabilistic basis for background estimation and enrichment testing.

### Background Estimation and Local Modulation

1. **Global background estimation**  

The genome-wide background rate \( \lambda_g \) is estimated as the mean read count per non-blacklisted genomic tile:

$$
\lambda_g = \frac{1}{N} \sum_{i=1}^{N} y_i
$$

Zero-count tiles are retained to avoid upward bias in sparse datasets.

2. **Control-based modulation**  

Local variability is corrected using a control-derived modulation factor:

$$
m_i = \max(1, \frac{y_i^{(ctrl)}}{\lambda_g^{(ctrl)}})
$$

The expected signal for each region *i* in experiment *j* is then:

$$
\lambda_{t,i}^{(j)} = m_i \times \lambda_g^{(j)}
$$

This integrates both global sequencing depth and local control variation into a unified, interpretable background model.

### Statistical Testing and Enrichment Scoring

ChromCall performs one-sided Poisson tests to evaluate whether observed counts significantly exceed the expected signal:

$$
p_i^{(j)} = P(Y \ge y_i^{(j)} \mid Y \sim Pois(\lambda_{t,i}^{(j)}))
$$

Additional quantitative metrics include:

- **Enrichment score:**  

$$
s_i^{(j)} = \log_2\left( \frac{y_i^{(j)} + \epsilon}{\lambda_{t,i}^{(j)} + \epsilon} \right)
$$

- **Poisson z-score:**  

$$
z_i^{(j)} = \frac{y_i^{(j)} - \lambda_{t,i}^{(j)}}{\sqrt{\lambda_{t,i}^{(j)}}}
$$

Multiple testing correction is performed using the **Benjamini–Hochberg (FDR)** procedure.  
Together, these metrics provide complementary measures of enrichment strength, effect size, and statistical confidence.

### Implementation and Data Structure

ChromCall is implemented in **R** and builds upon the **Bioconductor** framework:

- **GRanges** for genomic intervals  
- **SummarizedExperiment** for structured assay outputs  
- **GenomicAlignments** for BAM import  
- **GenomicRanges** and **Seqinfo** for annotation management  

Each processed sample is returned as a `SummarizedExperiment` containing:

- Raw counts  
- Background estimates (\(\lambda_g\), \(\lambda_t\))  
- p-values and FDR-adjusted p-values  
- Enrichment scores and Poisson z-scores  

ChromCall also supports pairwise comparison outputs, including **Δ enrichment** and **Δ z-score** metrics for comparative chromatin analysis.

This modular architecture facilitates future extensions to additional statistical models or multi-sample designs.

## 📦 Installation

```r
# Install development version from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("wangbo17/ChromCall")
```

## 🧬 Example Usage

```r
library(ChromCall)

# Input files
bam_exp <- "H3K27me3_treated.bam"
bam_ctrl <- "H3K27me3_input.bam"
regions  <- "promoters.bed"

# Run ChromCall
result <- run_chromcall(
    bam_exp = bam_exp,
    bam_ctrl = bam_ctrl,
    regions = regions,
    blacklist = "blacklist.bed"
)

# Access results
assay(result, "enrichment_score")
metadata(result)
```

## 📈 Outputs

| Output metric            | Description                                   |
| ------------------------ | --------------------------------------------- |
| `count`                  | Raw read count per region                     |
| `lambda_g`               | Genome-wide background rate                   |
| `lambda_t`               | Locally adjusted expected signal              |
| `pval` / `padj`          | Poisson test p-values and FDR-adjusted values |
| `score`                  | log₂(Observed/Expected) enrichment            |
| `z_pois`                 | Poisson z-score                               |
| `delta_score`, `delta_z` | Pairwise comparison metrics                   |

## 📄 License

MIT License © 2025 Bo Wang

## 💡 Contact

For questions, issues, or feature requests, please open a [GitHub issue](https://github.com/wangbo17/ChromCall/issues).
