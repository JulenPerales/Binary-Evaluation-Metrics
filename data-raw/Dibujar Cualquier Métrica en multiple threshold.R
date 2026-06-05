library(ROCR)
library(PRROC)
E<-100
Pr<-0.4
P<-E*Pr
Q<-E-P
# l<-50/(65*0.45^2)
# f1<-function(x){(-l*x^2+l*0.9*x)*0.65}
# f2<-function(x){2.46*exp(-4.4*x)+0.9*0.3}
Upper<-function(x){min(x,P)}
Lower<-function(x){max(0,(x-(E-P)))}
Random<-function(x){Pr*x}
a <- readRDS(file = "C:/Users/julen.perales/OneDrive - UPNA/Tesis/Tabla_clasificador_Skill.rds")
Curve<-a$FOM
THEME<-list(
  legend.key.width = unit(4, "line"),
  legend.key.height = unit(4, "line"),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
  plot.background = element_rect(fill = "white", color = "white"),
  panel.background = element_rect(fill = "white", color = "white"),
  plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 4),
  plot.subtitle = element_text(size = 18, hjust = 0.5),
  plot.caption = element_text(size = 10, hjust = 0.5),
  axis.text.x = element_text(size = 20, vjust=-3),
  axis.text.y = element_text(size = 20, hjust=-0),
  axis.title.x = element_text(size = 28, face = "bold", vjust = -4),
  axis.title.y = element_text(size = 28, face = "bold", vjust = 3),
  axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
  panel.border = element_rect(color = "black", fill = NA),
  panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = 18, hjust = 0, vjust = 5),
  legend.key = element_rect(color = "white"),
  legend.background = element_rect(fill = "white", color = "white"),
  legend.text = element_text(size = 28),
  legend.position = "right")



# 1) Build Contingency once
Contingency <- data.frame(
  Threshold = a$`Hits+FalseAlarms`,
  Curve     = a$Acc*Pr/0.3,
  Upper     = sapply(x1, Upper),
  Lower     = sapply(x1, Lower),
  Random    = sapply(x1, Random)
)

# 2) Metric functions: signature (x, Threshold, P, E)

OA  <- function(x, Threshold, P, E) (x + E - P - (Threshold - x)) / E

BA <- function(x, Threshold, P, E) {
  0.5 * (x / P + (E - Threshold - P + x) / (E - P))
}

MCC <- function(x, Threshold, P, E) {
  TP <- x
  FP <- Threshold - x
  FN <- P - x
  TN <- E - Threshold - P + x
  
  num <- TP * TN - FP * FN
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  
  ifelse(den == 0, 0, num / den)
}

Pierce <- function(x, Threshold, P, E) {
  x / P - (Threshold - x) / (E - P)
}

Kappa <- function(x, Threshold, P, E) {
  TP <- x
  TN <- E - Threshold - P + x
  
  po <- (TP + TN) / E
  pe <- (Threshold * P + (E - Threshold) * (E - P)) / E^2
  
  ifelse(1 - pe == 0, NA_real_, (po - pe) / (1 - pe))
}

F1  <- function(x, Threshold, P, E) {
  denom <- Threshold + P
  ifelse(denom == 0, NA_real_, 2 * x / denom)
}
F2 <- function(x, Threshold, P, E) {
  num <- 5 * x
  den <- 5 * x + 4 * (P - x) + (Threshold - x)
  ifelse(den == 0, NA_real_, num / den)
}


FOM <- function(x, Threshold, P, E) x / (P + (Threshold - x))

GSS <- function(x, Threshold, P, E) {
  Hr <- (Threshold * P) / E
  num <- x - Hr
  den <- (x + (Threshold - x) + (P - x)) - Hr  # TP + FP + FN - Hr
  
  ifelse(den == 0, NA_real_, num / den)
}


# 3) One function that: builds df + plots it
plot_metric_from_contingency <- function(Contingency, metric_fun, ylab,
                                         P, E, threshold_scale = 100) {
  
  T <- Contingency$Threshold
  
  df <- data.frame(
    Threshold = T / threshold_scale,
    Curve     = metric_fun(Contingency$Curve,  T, P, E),
    Upper     = metric_fun(Contingency$Upper,  T, P, E),
    Lower     = metric_fun(Contingency$Lower,  T, P, E),
    Random    = metric_fun(Contingency$Random, T, P, E)
  )
  if(sum(df$Random)<1){type<-"Skill"}else{type<-"Agreement"}
  
  # ---- Area ratio: area(Curve-Lower) / area(Upper-Lower) ----
  x <- df$Threshold
  if(type=="Agreement"){
  area_num <- trapz(x, pmax(0, df$Curve - df$Lower))
  area_den <- trapz(x, pmax(0, df$Upper - df$Lower))}
  else{
  area_num <- trapz(x, pmax(0, df$Curve - df$Random))
  area_den <- trapz(x, pmax(0, df$Upper - df$Random))}
  area_ratio <- ifelse(area_den == 0, NA_real_, area_num / area_den)
  
  ggplot(df, aes(x = Threshold)) +
    geom_line(aes(y = Curve,  color = "Classifier",        linetype = "Classifier"),        size = 1.25) +
    geom_line(aes(y = Upper,  color = "Upper Bound",       linetype = "Upper Bound"),       size = 1.25) +
    geom_line(aes(y = Lower,  color = "Lower Bound",       linetype = "Lower Bound"),       size = 1.25) +
    geom_line(aes(y = Random, color = "Random Classifier", linetype = "Random Classifier"), size = 1) +
    labs(x = "Hits + False Alarms (% of Extent)", y = ylab) +
    scale_color_manual(
      name = "",
      values = c("Classifier" = "#298d8d", "Upper Bound" = "grey20", "Lower Bound" = "grey65",
                 "Random Classifier" = "#950000", "Naive Model" = "purple3"),
      limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Classifier" = "solid", "Upper Bound" = "solid", "Lower Bound" = "solid",
                 "Random Classifier" = "dotted", "Naive Model" = "twodash"),
      limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = c(0, 0), labels = percent_format()) +
    scale_y_continuous(limits = c(0, 1.0001), breaks = seq(0, 1, by = 0.1), expand = c(0, 0), labels = percent_format()) +
    theme(!!!THEME) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(reverse = TRUE),
           linetype = guide_legend(reverse = TRUE)) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = sprintf("Area ratio = %.3f", area_ratio),
      hjust = 1.1, vjust = 1.4,
      size = 10
    )
}

AUCPRC<-trapz(
  Contingency$Curve/P,
  ifelse(Contingency$Threshold == 0, 1, Contingency$Curve/Contingency$Threshold)
)

ggplot(Contingency, aes(x = Threshold)) +
  geom_line(aes(y = Curve,  color = "Classifier",        linetype = "Classifier"),        size = 1.25) +
  geom_line(aes(y = Random, color = "Random Classifier", linetype = "Random Classifier"), size = 1) +
  labs(x = "Hits + False Alarms (% of Extent)", y = "Hits") +
  scale_color_manual(
    name = "",
    values = c("Classifier" = "#298d8d", "Upper Bound" = "grey20", "Lower Bound" = "grey65",
               "Random Classifier" = "#950000", "Naive Model" = "purple3"),
    limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Classifier" = "solid", "Upper Bound" = "solid", "Lower Bound" = "solid",
               "Random Classifier" = "dotted", "Naive Model" = "twodash"),
    limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
  ) +
  scale_x_continuous(limits = c(0, E), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, P), expand = c(0, 0)) +  theme(!!!THEME) +
  guides(color = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE))+
  coord_fixed(ratio = 1/Pr)

ggplot(Contingency, aes(x = Curve/P)) +
  geom_line(aes(y = Curve/Threshold,  color = "Classifier",        linetype = "Classifier"),        size = 1.25) +
  # geom_line(aes(y = Upper/Threshold, x = Upper/P,  color = "Upper Bound",       linetype = "Upper Bound"),       size = 1.25) +
  # geom_line(aes(y = Lower,  color = "Lower Bound",       linetype = "Lower Bound"),       size = 1.25) +
  geom_line(aes(y = Pr, color = "Random Classifier", linetype = "Random Classifier"), size = 1) +
  labs(x = "Recall", y = "Precision") +
  scale_color_manual(
    name = "",
    values = c("Classifier" = "#298d8d", "Upper Bound" = "grey20", "Lower Bound" = "grey65",
               "Random Classifier" = "#950000", "Naive Model" = "purple3"),
    limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Classifier" = "solid", "Upper Bound" = "solid", "Lower Bound" = "solid",
               "Random Classifier" = "dotted", "Naive Model" = "twodash"),
    limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound", "Lower Bound")
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = c(0, 0), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1.0001), breaks = seq(0, 1, by = 0.1), expand = c(0, 0), labels = percent_format()) +
  theme(!!!THEME) +
  guides(color = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE)) +
  coord_fixed(ratio = 1)+
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("AUCPRC = %.3f", AUCPRC),
    hjust = 1.1, vjust = 1.4,
    size = 10
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("Skill = %.3f", (AUCPRC - Pr)/(1 - Pr)),
    hjust = 1.1, vjust = 2.8,
    size = 10
  )
  
  
 
# 4) Call it once per metric
plot_metric_from_contingency(Contingency, OA,    "OA",    P = P, E = E)
plot_metric_from_contingency(Contingency, BA,    "BA",    P = P, E = E)
plot_metric_from_contingency(Contingency, MCC,   "MCC",   P = P, E = E)
plot_metric_from_contingency(Contingency, Pierce,"Pierce",P = P, E = E)
plot_metric_from_contingency(Contingency, Kappa, "Kappa", P = P, E = E)

plot_metric_from_contingency(Contingency, FOM,   "FOM",   P = P, E = E)
plot_metric_from_contingency(Contingency, F1,    "F1",    P = P, E = E)
plot_metric_from_contingency(Contingency, F2,    "F2",    P = P, E = E)
plot_metric_from_contingency(Contingency, GSS, "GSS", P = P, E = E)

