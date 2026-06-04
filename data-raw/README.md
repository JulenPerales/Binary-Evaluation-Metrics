# data-raw

Reference scripts used during thesis development. These scripts contain the
original R code that the `binmetrics` package is based on. They are kept here
for historical reference and reproducibility, but are not part of the package.

| File | Description |
|------|-------------|
| `Para simulaciones.R` | Monte Carlo simulation of synthetic classifiers; produces Excel outputs and comparison plots |
| `Para simulaciones artículo 2.R` | Extended simulations for Article 2 with OA/PSS analysis |
| `Tesis_Figure_Skill.R` | Generates CSI/OA comparison figures for the thesis |
| `Dibujar Cualquier Métrica en multiple threshold.R` | Generic metric-curve plotter; direct inspiration for `toc_metric_curve()` and `plot_metric_curve()` |
| `Figure Simulations pero interactiva.R` | Interactive Shiny app showing AUC/FOM bounds; source for `launch_bounds_explorer()` |
| `Para Mostrar Deslizante con Classificaciones.R` | Interactive Shiny app with classification maps; source for `launch_threshold_explorer()` |
| `Formulas de SKILL SCORE para FINLEY.R` | Reference implementations of ETS, HSS, PSS using the Finley (1884) example |
