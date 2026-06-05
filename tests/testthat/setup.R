# When running tests without installing the package, source the R files directly.
# When installed, library(binmetrics) is used via devtools::test().
if (requireNamespace("binmetrics", quietly = TRUE)) {
  library(binmetrics)
} else {
  r_files <- c(
    "../../R/utils.R", "../../R/confusion_matrix.R", "../../R/metrics.R",
    "../../R/bounds.R", "../../R/toc.R", "../../R/integration.R",
    "../../R/BinaryClassifier.R"
  )
  invisible(lapply(r_files, source))
}
