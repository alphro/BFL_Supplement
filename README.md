# BFLpkg_Support

Supporting materials for [**Bayesian Federated Cause-of-Death Classification and Quantification Under Distribution Shift
** by Yu Zhu, Jason Teng, and Zehang Richard Li](https://arxiv.org/abs/2505.02257).

This repository includes a standalone reproducible example on real PHMRC data, and the scripts to produce all figures in the paper.

## Repository layout

```
├── BFL_Example/                     standalone reproducible example (real PHMRC, target = AP)
│   ├── BFL_AP_example.Rmd           walkthrough: base / domain / partial / mix across the three
│   │                                label-shift cases, local-self & local-avg baselines, S3 methods
│   ├── helpers.R                    LCVA + PHMRC data glue and per-variant input builders
│   ├── phmrc_clean.csv              bundled PHMRC data
│   └── README.md                    how to run the example
│
├── Results/                         paper reproduction
│   ├── CSV/                         result tables
│   └── Figures/
│       ├── scripts/                 R scripts that generate every figure (PHMRC + CHAMPS)
│       ├── PHMRC/                   PHMRC figures
│       ├── CHAMPS/                  CHAMPS figures
│
```

## Quick start

**Install the package** (needs `LCVA`, not on CRAN):

```r
remotes::install_github("richardli/LCVA")
remotes::install_github("alphro/BFL")
```

**Obtain and process PHMRC data**:

The PHMRC data (csv in the script below) can be obtained from [https://ghdx.healthdata.org/record/ihme-data/population-health-metrics-research-consortium-gold-standard-verbal-autopsy-data-2005-2011]

```r
library(openVA)
PHMRC <- read.csv("IHME_PHMRC_VA_DATA_ADULT_Y2013M09D11_0.csv")
# turn into binary data
binarydata <- ConvertData.phmrc(input = PHMRC, input.test = PHMRC, cause = "gs_text34")
#  Create binary X matrix
X <- binarydata$output[, -2]
tmp <- binarydata$output[, -c(1, 2)]
X0 <- matrix(0, dim(tmp)[1], dim(tmp)[2])
X0[tmp == "Y"] <- 1
X0[tmp == "."] <- NA
dim(X0)
Y <- binarydata$output[,2]
sites <- PHMRC$site
data <- data.frame(X0)
data <- cbind(data, cause = Y, site = sites)
write.csv(data, file = "phmrc_clean.csv", row.names = FALSE)
```

**Run the reproducible example** (from the example folder, so the bundled data resolves):

```r
setwd("BFL_Example")
rmarkdown::render("BFL_AP_example.Rmd")
```

It reproduces the paper's four BFL variants (base / domain / partial / mix) plus
the local-self and local-avg LCVA baselines, across the no-shift, mild, and severe
within-target label-shift settings, and demonstrates the package's S3 methods.

**Regenerate figures** — the scripts in `Results/Figures/scripts/` are CSV-driven;
each script's header lists its input CSV and outputs. Most read the tables in
`Results/CSV/` (and `Others/`) and write PDFs into `Results/Figures/`. Two CHAMPS
scripts (`figure_08`, `figure_10`) recompute from cached model fits and must run
where those caches live; the rest run anywhere with the CSVs.
