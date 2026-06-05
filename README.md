# binmetrics

**Binary classification evaluation metrics with full threshold analysis.**

`binmetrics` evaluates binary classifiers across the complete range of decision
thresholds, not just at a single operating point. It provides:

- **Agreement metrics**: OA, BA, MCC, Kappa, F1
- **Skill metrics**: FOM/CSI, GSS (= ETS), PSS, HSS
- **Integrated metrics**: AFOM, AUFOM, DFOM, AUC, DAUC, AUCPRC
- **Theoretical bounds**: upper and lower limits for any metric at every threshold
- **Visualisation**: publication-quality ggplot2 plots and interactive Shiny apps

The package is developed alongside the doctoral thesis on binary classification
evaluation by Julen Perales (UPNA).

---

## Installation

```r
# From GitHub
remotes::install_github("julenperales/binary-evaluation-metrics")

# Local development install
devtools::install(".")
```

## Quick start

```r
library(binmetrics)

# --- Single threshold ---
cm <- confusion_matrix(hits = 28, fa = 72, misses = 23, cr = 2680)
summary(cm)
fom(cm)   # Figure of Merit / CSI
gss(cm)   # Gilbert Skill Score

# --- Full threshold analysis from scores ---
set.seed(42)
scores <- runif(1000)
obs    <- scores > 0.5 + rnorm(1000, sd = 0.25)

clf <- BinaryClassifier$new(scores, obs, name = "My Classifier")
clf$summary()

clf$plot_toc()         # TOC curve
clf$plot_fom()         # FOM with bounds + area metrics
clf$plot_roc()         # ROC curve
clf$plot_all_metrics() # 8-panel grid

# --- Interactive apps ---
launch_threshold_explorer(clf$toc)
launch_bounds_explorer()
```

## Key concepts

### TOC space

The Threshold Operating Characteristic space represents classifier performance
as Hits (TP count) versus the number of predicted positives `k`. At every
threshold `k`, the theoretical bounds are:

- **Upper bound**: `min(k, P)` hits — best possible classifier
- **Lower bound**: `max(0, k − Q)` hits — worst-case ranking
- **Random baseline**: `k × prevalence` hits — no-discrimination classifier

### Integrated metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| AFOM   | ∫(FOM − FOM\_lower) / ∫(FOM\_upper − FOM\_lower) | Normalised area under FOM curve |
| AUFOM  | Same integral for the random classifier | Random-classifier AFOM baseline |
| DFOM   | (AFOM − AUFOM) / (1 − AUFOM) | FOM skill vs. random classifier |
| DAUC   | 2 × (AUC − 0.5) | AUC skill vs. random classifier |

## Package structure

```
R/
├── confusion_matrix.R   # S3 class: confusion_matrix()
├── metrics.R            # oa, ba, mcc, kappa_score, f1_score, fom, gss, pss, hss, …
├── bounds.R             # hits_upper/lower/random, fom_upper/lower/random
├── toc.R                # toc_data(), toc_from_scores(), toc_metric_curve()
├── integration.R        # afom, aufom, dfom, auc_roc, dauc, mfom, auc_prc
├── BinaryClassifier.R   # R6 class with OO API
├── plots.R              # plot_toc, plot_roc, plot_fom, plot_metric_curve, …
├── shiny_apps.R         # launch_threshold_explorer, launch_bounds_explorer
└── utils.R              # binm_theme, colour palette
```

## License

MIT © Julen Perales