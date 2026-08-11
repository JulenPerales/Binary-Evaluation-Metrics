# binmetrics

**An R package for the evaluation of binary classifiers, built around the CSI Framework.**

`binmetrics` is the companion package to the doctoral dissertation
*"The CSI Framework for the Evaluation of Binary Classifications in Rare Events"*
(Julen Perales, UPNA, 2026). It provides a complete, coherent toolkit for evaluating
binary classifiers — from single-threshold metrics to full multi-threshold analysis —
with a particular focus on applications where the event of interest is rare.

---

## Why this package?

Standard evaluation tools scatter metrics across disciplines and often mix up two
fundamentally different concepts:

- **Agreement** — how similar are the predictions to the observations?
- **Skill** — does the classifier outperform a random baseline?

A classifier can have high Agreement simply because the majority class dominates
(e.g. OA is always high when land change is rare). Skill metrics correct for this by
comparing to what random guessing would achieve. `binmetrics` makes this distinction
explicit, names every metric consistently with the dissertation, and lets you evaluate
performance across every possible decision threshold — not just one.

---

## Installation

```r
# From GitHub
remotes::install_github("julenperales/binary-evaluation-metrics")

# Local development install
devtools::install(".")
```

---

## Core concepts

### The contingency table

Every binary evaluation starts here. Given a set of predictions and observations,
four outcomes are possible:

|                    | Observed presence | Observed absence |
|--------------------|:-----------------:|:----------------:|
| **Predicted presence** | **H** (Hits / TP) | **F** (False Alarms / FP) |
| **Predicted absence**  | **M** (Misses / FN) | **C** (Correct Rejections / TN) |

```r
# Build a confusion matrix directly from counts
cm <- confusion_matrix(hits = 28, fa = 72, misses = 23, cr = 2680)
summary(cm)
```

### Agreement vs Skill

| Category | Metrics | What they measure |
|----------|---------|-------------------|
| **Agreement** | OA, PA, SP, FAR, UA, FOR, CSI, F1, F-beta | Similarity between predictions and observations |
| **Skill** | GSS, HSS, PSS, MCC | Improvement over a random baseline |

Agreement is not enough. A model that classifies everything as "no change" in a land
use dataset can score OA > 95% while being completely useless. Skill metrics expose
this by subtracting what random chance would achieve.

### The CSI (Critical Success Index)

CSI — also known as Figure of Merit, Threat Score, or Jaccard Index depending on
the field — is the primary single-threshold metric of the CSI Framework:

```
CSI = H / (H + F + M)
```

It focuses exclusively on the event of interest (presence), ignoring Correct
Rejections. This makes it appropriate when the event is rare and C dominates the
table.

### TOC space

The **Threshold Operating Characteristic** curve is the foundation for multi-threshold
analysis. It plots Hits (H) on the y-axis against the number of predicted positives
`k = H + F` on the x-axis. At every value of `k`, the package computes:

- **Upper bound** — `min(k, P)` — best possible classifier (pure Allocation Agreement)
- **Lower bound** — `max(0, k − Q)` — worst possible classifier
- **Baseline** — `k × prevalence` — random classifier with no discrimination ability

The CSI curve plots `CSI(k)` on the y-axis against `k`, inheriting the same bounds
and baseline. The area between the CSI curve and these bounds gives the integrated
metrics.

---

## Quick start

### Single-threshold analysis

```r
library(binmetrics)

# Finley's tornado forecast data (the classic example)
cm <- confusion_matrix(hits = 28, fa = 72, misses = 23, cr = 2680)

# Agreement metrics (how similar are predictions to observations?)
oa(cm)          # Overall Agreement:        0.973
pa(cm)          # Producer's Accuracy (PA): 0.549  (= Sensitivity)
ua(cm)          # User's Accuracy (UA):     0.280  (= Precision)
sp(cm)          # Specificity (SP):         0.974
csi(cm)         # Critical Success Index:   0.227

# Skill metrics (how much better than random?)
gss(cm)         # Revised Gilbert Skill Score: 0.221
hss(cm)         # Heidke Skill Score:          0.362
pss(cm)         # Peirce Skill Score:          0.523
mcc(cm)         # Matthews Correlation Coeff:  0.381

# All at once — returns a tibble with Agreement/Skill classification.
# Skill metrics include the value of the corresponding agreement metric
# that a random classifier would achieve (baseline_agreement) and the
# name of that metric (corresponding_agreement_metric).
all_metrics(cm)
#> # A tibble: 12 x 5
#>    metric type      value baseline_agreement corresponding_agreement_metric
#>    <chr>  <chr>     <dbl>              <dbl> <chr>
#>  1 OA     Agreement 0.973             NA     NA
#>  2 PA     Agreement 0.549             NA     NA
#>  ...
#>  9 GSS    Skill     0.221            0.0122  CSI
#> 10 HSS    Skill     0.362            0.9474  OA
#> 11 PSS    Skill     0.523            0.0357  PA
#> 12 MCC    Skill     0.381            0.9474  OA
```

### Full multi-threshold analysis

```r
set.seed(42)
scores <- runif(1000)
obs    <- scores > 0.5 + rnorm(1000, sd = 0.25)

clf <- BinaryClassifier$new(scores, obs, name = "My Classifier")
clf
#> BinaryClassifier <My Classifier>
#>   N = 1000 | P = 393 | Q = 607 | prev = 0.3930

clf$summary()   # point metrics at optimal threshold + all integrated metrics
```

### Plots

```r
# Basic plots (continuous classifier curve with bounds and baseline)
clf$plot_toc()          # TOC curve with upper/lower bounds and baseline
clf$plot_csi()          # CSI curve — same structure as TOC
clf$plot_roc()          # ROC curve
clf$plot_all_metrics()  # 2x4 panel: OA / CSI / PA / SP (Agreement) + GSS / HSS / PSS / MCC (Skill)

# Overlay a single-threshold forecast as a reference point.
# For ROC, (FPR, TPR) are normalised — the point is always on [0, 1]^2.
# For CSI and TOC (percent = TRUE), the reference cm's own N and P are used.
clf$plot_toc(percent = TRUE, reference_cm = cm, reference_label = "Finley (1884)")
clf$plot_csi(reference_cm = cm, reference_label = "Finley (1884)")
clf$plot_roc(reference_cm = cm, reference_label = "Finley (1884)")
```

### Interactive exploration

```r
launch_threshold_explorer(clf$toc)   # slide through thresholds, watch metrics update
launch_bounds_explorer()             # explore how prevalence shapes bounds
```

---

## Integrated metrics

Multi-threshold evaluation summarises performance across all possible decision
thresholds. `binmetrics` implements the four integrated metrics introduced in the
dissertation:

| Name | Formula | Description |
|------|---------|-------------|
| **AUCSI** | ∫(CSI − CSI_lower) / ∫(CSI_upper − CSI_lower) | Normalised area under the CSI curve |
| **AUCSI_baseline** | Same integral for the random classifier | Chance-level AUCSI |
| **AUCSIS** | (AUCSI − AUCSI_baseline) / (1 − AUCSI_baseline) | CSI skill score across all thresholds |
| **AUC** | ∫ TPR d(FPR) | Standard area under the ROC curve |
| **AUCS** | 2 × (AUC − 0.5) | AUC skill score |
| **MaxCSI** | max(CSI curve) | Best CSI achievable at optimal threshold |

```r
# Access integrated metrics from the classifier object
clf$aucsi         # AUCSI
clf$aucsis        # AUCSIS (CSI skill across thresholds)
clf$aucs          # AUCS  (AUC skill)
clf$mcsi          # MaxCSI (best single-threshold CSI)

# Or compute directly from a toc_data object
td <- toc_from_scores(scores, obs)
integrated_metrics(td)
```

---

## The Revised Gilbert Skill Score

The dissertation identifies a conceptual inconsistency in the standard GSS formula:
it uses H+F+M (the actual union) for the denominator of CSI_random, but the random
classifier has a different union. The corrected formula is:

```
Revised GSS = (CSI − CSI_random) / (1 − CSI_random)
```

where `CSI_random = H_rand / (H_rand + F_rand + M_rand)` with
`H_rand = k·P/E`, `F_rand = k·Q/E`, `M_rand = P·(E−k)/E`.

The `gss()` function in `binmetrics` implements this corrected version.

---

## Working with toc_data

The `toc_data` object is the package's central data structure. It is a tibble with
one row per threshold `k` and the following key columns:

| Column | Description |
|--------|-------------|
| `k` | Number of predicted positives (0 to N) |
| `hits` | True positives at this threshold |
| `fa`, `misses`, `cr` | False alarms, Misses, Correct Rejections |
| `csi` | CSI value of the classifier |
| `csi_upper`, `csi_lower`, `csi_baseline` | CSI bounds and random baseline |
| `tpr`, `fpr` | True Positive Rate and False Positive Rate (for ROC) |

```r
td <- toc_from_scores(scores, obs)   # build from a score vector
td <- toc_data(hits_vector, P, Q)    # build from pre-computed hit counts

# Evaluate any metric across all thresholds
mc <- toc_metric_curve(td, pss, "PSS")
plot_metric_curve(mc, ylab = "PSS")

# Evaluate at a specific threshold
cm_at_100 <- td[td$k == 100, ]
```

### Building a classifier from pre-computed hits

Useful when you have a published TOC table rather than raw scores:

```r
P <- 5; Q <- 10
hits <- c(0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5)  # perfect classifier
clf  <- BinaryClassifier_from_hits(hits, P, Q, name = "Perfect")
clf$mcsi   # 1
clf$auc    # 1
```

---

## Metric reference

### Agreement metrics

| Function | Name | Formula |
|----------|------|---------|
| `oa()` | Overall Agreement | (H + C) / N |
| `pa()` | Producer's Accuracy | H / (H + M) |
| `ua()` | User's Accuracy | H / (H + F) |
| `sp()` | Specificity | C / (C + F) |
| `far()` | False Alarm Rate | F / (C + F) |
| `for_metric()` | False Omission Rate | M / (C + M) |
| `csi()` | Critical Success Index | H / (H + F + M) |
| `f1_score()` | F1 Score | 2H / (2H + F + M) |
| `f_beta()` | F-beta Score | (1+β²)H / ((1+β²)H + β²M + F) |

`fom()` is kept as an alias for `csi()` for compatibility with land change literature.

### Skill metrics

| Function | Name | Baseline |
|----------|------|----------|
| `gss()` | Revised Gilbert Skill Score | Random classifier (correct union) |
| `hss()` | Heidke Skill Score | Random classifier (same quantity as predicted) |
| `pss()` | Peirce Skill Score | Unbiased random classifier |
| `mcc()` | Matthews Correlation Coefficient | Random classifier |

Note: `kappa_score()` is equivalent to HSS for binary classification. Its use is
discouraged by the dissertation; `hss()` is preferred.

---

## Package structure

```
R/
├── confusion_matrix.R   # S3 class: confusion_matrix()
├── metrics.R            # All single-threshold metrics
├── bounds.R             # Theoretical bounds for hits and CSI
├── toc.R                # toc_data(), toc_from_scores(), toc_metric_curve()
├── integration.R        # AUCSI, AUCSIS, AUC, AUCS, MaxCSI, AUPRC, …
├── BinaryClassifier.R   # R6 class — object-oriented API
├── plots.R              # ggplot2 visualisations
├── shiny_apps.R         # Interactive Shiny explorers
└── utils.R              # Theme and colour palette
```

---

## Citation

If you use `binmetrics` in your research, please cite the dissertation:

> Perales Barriendo, J. (2026). *The CSI Framework for the Evaluation of Binary
> Classifications in Rare Events: Applications to Land Use Change Models and
> Environmental Valuation*. Doctoral dissertation, Public University of Navarra (UPNA).

---

## License

MIT © Julen Perales
