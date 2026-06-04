gc()
# Install and load required packages if not already installed
if (!requireNamespace("readxl", quietly = TRUE)) {
install.packages("readxl")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
install.packages("dplyr")
}
library(readxl)
library(dplyr)
library(pracma)
library(ggplot2)
library(SciViews)
library(patchwork)
library(cowplot)
library(pbapply)
library(tibble)
# Set the directory containing Excel files
folder_path <- "C:/Users/julen.perales/Desktop/GG/new"
# List all Excel files in the directory
excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
# Initialize an empty list to store data frames
dfs <- list()
# Read each Excel file and store data frames in the list
for (file in excel_files) {
df <- read_excel(file)
dfs[[length(dfs) + 1]] <- df
}
# Merge all data frames in the list
merged_df <- bind_rows(dfs)
# print the merged data frame
tabla<-merged_df[!is.na(merged_df$AFOM), ]
# tabla<-merged_df
my_theme <- theme(
  plot.background = element_rect(fill = "white", color = "white"),
  panel.background = element_rect(fill = "white", color = "white"),
  plot.title = element_text(size = 22, face = "bold", hjust = 0.5, vjust = 1),
  plot.subtitle = element_text(size = 18, hjust = 0.5),
  plot.caption = element_text(size = 10, hjust = 0.5),
  axis.text = element_text(size = 24),
  axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
  axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
  axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
  legend.text = element_text(size = 24),
  legend.title = element_text(size = 24, hjust = 0, vjust = 5),
  legend.position = "right",
  legend.key = element_rect(color = "white"),
  legend.background = element_rect(fill = "white", color = "white"),
  legend.key.width = unit(2, "cm"),  # ??? makes lines in the legend longer
  panel.border = element_rect(color = "black", fill = NA),
  panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
  panel.grid.minor = element_blank(),
  plot.margin = margin(t = 20, r = 15, b = 20, l = 15)
)


FOM<-function(a, p, pop){
  b<-a
  c<-a
  d<-a

  #calculamos FOM para AUC=1
  b["Hits"]<-a$`Hits+FalseAlarms`
  b["Hits"][b["Hits"]>p]<-p
  b["misses"]<-p-b$Hits
  b["FOM"]<-b$Hits/(b$`Hits+FalseAlarms`+b$misses)
  #calculamos FOM para AUC=0.5
  c["Hits"]<-a$`Hits+FalseAlarms`*p/pop
  c["misses"]<-p-c$Hits
  c["FOM"]<-c$Hits/(c$`Hits+FalseAlarms`+c$misses)
  #calculamos FOM para AUC=0.0
  # d["Hits"]<-pmax(0, a$`Hits+FalseAlarms`-p)
  d$Hits <- 0
  d$Hits[d$`Hits+FalseAlarms`>(pop-p)] <- d$`Hits+FalseAlarms`[d$`Hits+FalseAlarms`>(pop-p)]-(pop-p)
  d["misses"]<-p-d$Hits
  d["FOM"]<-d$Hits/(d$`Hits+FalseAlarms`+d$misses)
  # Remove the last row from the dataframe
  # a <- a[-nrow(a), ]
  # b <- b[-nrow(b), ]
  # c <- c[-nrow(c), ]
  # d <- d[-nrow(c), ]
  library(pracma)
  x <- a$`Hits+FalseAlarms`
  y1 <- a$FOM
  y2 <- b$FOM
  y3 <- c$FOM
  y4 <- d$FOM
  area_over <- trapz(x, y2-y1)/trapz(x, y2-y4)
  # area_over <- trapz(x, rep(0,length(a$FOM)))/trapz(x, y2)
  area_under <- trapz(x, y1-y4)/trapz(x, y2-y4)
  area_under_random <- trapz(x, y3-y4)/trapz(x, y2-y4)
  return(list(area_under,area_under_random, b, c, d))
}
pr<-0.1
pop<-10000
p<-pop*pr
E<-pop

f1 <- function(x) {
  result <- -yt*x*exp(x*xt)/(1-exp(x*xt))-pr
  
  # Check for NaN values
  if (is.nan(result)) {
    return(1e10)  # or any large value or sentinel value that makes sense in your context
  }
  
  return(result)
}
f2 <- function(x) {
  result <- (p-yt+pr/x)/(pr/(x*exp(x*xt)))-exp(x*E)
  
  # Check for NaN values
  if (is.nan(result)) {
    return(1e10)  # or any large value or sentinel value that makes sense in your context
  }
  
  return(round(result,4))
}

Curva<-function(xt,yt){
  # b1<-optimize(f = f1, lower = -0.1, upper = 0, tol = 1e-2)
  # b1<-b1$minimum
  interval <- c(-1, 0)
  # Define precision criterion
  precision <- 1e-8
  while(interval[2] - interval[1] > precision) {
    # Use optimize to find the minimum (closest to zero) of the absolute value of f2(x) within the current interval
    root <- optimize(f = f1, interval = interval)
    
    # Update the interval based on the root found
    if (f1(root$minimum) < 0) {
      interval <- c(root$minimum, interval[2])
    } else {
      interval <- c(interval[1], root$minimum)
    }
  }
  b1<-root$minimum
  a1<-yt/(1-exp(b1*xt))
  k1<--a1
  y1<-function(x){min(a1+k1*exp(b1*x),x)}
  x1<-seq(0,xt, length.out=1000)
  x1<-seq(0,xt-1)
  # plot(x,sapply(x,y1))
  interval <- c(-1, 0)
  while(interval[2] - interval[1] > precision) {
    # Use optimize to find the minimum (closest to zero) of the absolute value of f2(x) within the current interval
    root <- optimize(f = f2, interval = interval)
    
    # Update the interval based on the root found
    if (f2(root$minimum) < 0) {
      interval <- c(root$minimum, interval[2])
    } else {
      interval <- c(interval[1], root$minimum)
    }
  }
  b2<-root$minimum
  # b2<-optimize(f = f2, lower = -0.00001, upper = 0, tol = 1e-2)$minimum
  # b2<-b2$minimum
  a2<-yt-pr/b2
  k2<-pr/(b2*exp(b2*xt))
  y2<-function(x){min(a2+k2*exp(b2*x),p)}
  # plot(x,sapply(x,y2))
  x2<-seq(xt, E, length.out=1000)
  x2<-seq(xt, E)
  x<-c(x1,x2)
  df<-data.frame(
    Sim=x,
    Upper=pmin(x,p),
    Lower=pmax(0,x-(E-p)),
    Curve=c(sapply(x1,y1),sapply(x2,y2))
  )
  if (any(is.na(df$Curve))) {df$Curve[is.na(df$Curve)]<-p}
  # print(ggplot(df, aes(x=Sim))+
  #   geom_line(aes(y=Upper))+
  #   geom_line(aes(y=Lower))+
  #   geom_line(aes(y=Curve)))
  # plot(c(x1,x2),c(sapply(x1,y1),sapply(x2,y2)))
  return(c(sapply(x1,y1),sapply(x2,y2)))
}

for (pr in c(0.01, 0.05,0.1,0.25)) {
  pop<-10000
  p<-pop*pr
  E<-pop
  tabla<-as.data.frame(matrix(0,1,ncol = 14))
  names(tabla)<-c("xt","yt","AFOM","AUFOM","AUC","DAUC", "DFOM", "MFOM","L","bMFOM","ThresholdMFOM","GSS","Julen", "ThresholdMSKILL")
  pb <- txtProgressBar(min = 1, max = E, style = 3)
  
  for(xt in c(1,seq(10,E, by=10))){
    # print((xt/E)*100)
    setTxtProgressBar(pb, xt)
    for(yt in c(1,seq(p/100,p,by=p/100))){
      if (yt<=xt&yt>pr*xt){
      gc()
      a<-as.data.frame(matrix(0,nrow=pop+1,ncol = 2))
      names(a)<-c("Hits+FalseAlarms", "FOM")
      a[1]<-c(0:pop)
      # a["rate"]<-r1+a[1]/e
      a["rate"]<-1
      a["Hits"]<-0
      a["Acc"]<-Curva(xt,yt)
      a["Hits"]<-c(0,diff(a$Acc))
      a["misses"]<-p-a$Acc
  
      a[2]<-a[5]/(a[1]+a[6])
  
      AUC<-0
      b<-a["Hits+FalseAlarms"]-a["Acc"]
      for (t in 1:(nrow(a)-1)){
        AUC<-AUC+((a[["Acc"]][t+1]+a[["Acc"]][t])/p)*(b[[t+1,1]]-b[[t,1]])/(pop-p)/2
      }
      # area<-area-P^2/(2*E)
      # AUC<-trapz(x=c(a$`Hits+FalseAlarms`),y=c(a$`Acc`))/(p*(pop)-p^2)
      lista<-FOM(a,p,pop)
      # abline(a = 0, b = 1, col = "red")
      a_df <- data.frame(x = a[[1]], y = a[[2]] * 100)
      # Additional calculations
      AFOM <- lista[[1]] * 100
      AUFOM <- lista[[2]] * 100
      AUC <- AUC * 100
      DAUC <- (AUC / 100 - 0.5) * 2 * 100
      DFOM <- ((lista[[1]] - lista[[2]]) / (1 - lista[[2]])) * 100
      MFOM <- max(a_df$y)
      thr<-which(a_df$y==max(a_df$y))
      bMFOM<-(thr*pr)/((1-pr)*thr+p)*100
      if (any(is.na(a))) {
        tabla[nrow(tabla)+1,]<-c(xt,yt,0,0,0,0,0,0,0, 0, 0, 0, 0, 0)
        
        next  # Skip this iteration if any NA is found in the row
      }
  
      a$F<-a$`Hits+FalseAlarms`-a$Acc
      a$C<-E-(a$misses+a$Acc+a$F)
      a$E<-a$C+a$misses+a$Acc+a$F
      a$Hr<-pr*a$`Hits+FalseAlarms`
      a$Fr<-(1-pr)*a$`Hits+FalseAlarms`
      a$Mr<-p-a$Hr
      a$Cr<-E-a$Hr-a$Fr-a$Mr
      a$GSS<-(a$Acc/(a$`Hits+FalseAlarms`+a$misses)-a$Hr/(a$`Hits+FalseAlarms`+a$misses))/(1-a$Hr/(a$`Hits+FalseAlarms`+a$misses))*100
      a$Julen<-(a$Acc/(a$`Hits+FalseAlarms`+a$misses)-a$Hr/(a$Hr +a$Fr+a$Mr))/(1-a$Hr/(a$Hr +a$Fr+a$Mr))*100
      thrSkill<-which(a$Julen==max(a$Julen))
      
      #de regalo, te calculo también OA, para la tesis
      OA<-max((a$Acc+a$C)/(a$E),na.rm=TRUE)
      PSS<-max((a$Acc / (a$Acc + a$misses)) - (a$F / (a$C+ a$F)),na.rm=TRUE)
      max(a$Julen,na.rm=TRUE)
      max(a$FOM,na.rm=TRUE)
      # plot(a$`Hits+FalseAlarms`, (a$Acc+a$C)/(a$E))
      # plot(a$`Hits+FalseAlarms`, (a$Acc / (a$Acc + a$misses)) - (a$F / (a$C+ a$F)))
      # plot(a$`Hits+FalseAlarms`, a$Acc)
      
      tabla[nrow(tabla)+1,]<-c(xt,yt,AFOM,AUFOM,AUC,DAUC,DFOM,MFOM,which(a["FOM"]==max(a["FOM"]))-p, bMFOM, thr, a$GSS[thr], a$Julen[thr], thrSkill)
     
    }
  }
  }
  #   write_xlsx(tabla,"C:/Users/julen.perales/OneDrive - UPNA/Tesis/Simulaciones.xlsx")
  # tabla <- tabla[tabla$MFOM >= pr*100, ]
  tabla <- tabla[tabla$DAUC >=0, ]
  
  writexl::write_xlsx(tabla, paste0("C:/Users/julen.perales/OneDrive - UPNA/Articulo/ASCI/Simulations/simulations_",pr*100, ".xlsx"))
}

#PRIMERO PLOTEAMOS LAS CURVAS ROC Y CSI PARA LA FIGURE 3
tabla<- read_xlsx("C:/Users/julen.perales/OneDrive - UPNA/Articulo/ASCI/Simulations/simulations_5.xlsx")
pr<-0.05
pop<-10000
p<-pop*pr
E<-pop
AUC_A<-90.2624332834417
AUC_B<-89.5284917717438
#este pertenece a xt=780 y yt=350
A_index<-which(round(tabla$AUC,5)==round(AUC_A,5))
#este pertenece a xt=2960 y yt=460
B_index<-which(round(tabla$AUC,5)==round(AUC_B,5))

# Get xt and yt for A_index and B_index
xt_A <- tabla[[A_index, 1]]
yt_A <- tabla[[A_index, 2]]

xt_B <- tabla[[B_index, 1]]
yt_B <- tabla[[B_index, 2]]

xt<-xt_A
yt<-yt_A

a <- as.data.frame(matrix(0, nrow = pop + 1, ncol = 2))
names(a) <- c("Hits+FalseAlarms", "FOM")
a[["Hits+FalseAlarms"]] <- 0:pop
a[["rate"]] <- 1
a[["Hits"]] <- Curva(xt, yt)
a[["misses"]] <- pop * pr - a$Hits

roc_df_A <-  data.frame(
  FPR = (a$`Hits+FalseAlarms` - a$Hits) / (E * (1 - pr)),
  TPR = a$Hits / p,
  FOM = a[["Hits"]]/(a[["Hits+FalseAlarms"]]+a[["misses"]]),
  HF = a[["Hits+FalseAlarms"]]
)

xt<-xt_B
yt<-yt_B

b <- as.data.frame(matrix(0, nrow = pop + 1, ncol = 2))
names(b) <- c("Hits+FalseAlarms", "FOM")
b[["Hits+FalseAlarms"]] <- 0:pop
b[["rate"]] <- 1
b[["Hits"]] <- Curva(xt, yt)
b[["misses"]] <- pop * pr - b$Hits

roc_df_B <-  data.frame(
  FPR = (b$`Hits+FalseAlarms` - b$Hits) / (E * (1 - pr)),
  TPR = b$Hits / p,
  FOM = b[["Hits"]]/(b[["Hits+FalseAlarms"]]+b[["misses"]]),
  HF = b[["Hits+FalseAlarms"]]
)

# Add labels for legend
roc_df_A$model <- "A_index"
roc_df_B$model <- "B_index"

# Combine both data frames
roc_combined <- rbind(roc_df_A, roc_df_B)


ggplot() +
  geom_line(data = roc_combined, aes(x = FPR, y = TPR, color = model, linetype = model), linewidth = 1.5) +
  geom_line(data = data.frame(FPR = c(0, 1), TPR = c(0, 1), model = "Uniform Line"), 
            aes(x = FPR, y = TPR, color = model, linetype = model), linewidth = 1) +
  scale_color_manual(name="",values = c("A_index" = "#1f77b4", "B_index" = "#bb3115", "Uniform Line" = "gray70")) +
  scale_linetype_manual(name="",values = c("A_index" = "solid", "B_index" = "solid", "Uniform Line" = "dashed")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  coord_equal() + my_theme



baselineFOM<-function(pr,pop){
  x<-seq(0,pop)
  hits<-x*pr
  fa<-x*(1-pr)
  return(hits/(fa+p))
}
upperFOM<-function(pr,pop){
  x<-seq(0,pop)
  hits<-pmin(pr*pop, x)
  fa<-pmax(0,x-pr*pop)
  return(hits/(fa+p))
}
lowerFOM<-function(pr,pop){
  x<-seq(0,pop)
  hits<-pmax(0, x-(pop-pr*pop))
  fa<-pmax(x,(1-pr)*pop)
  return(hits/(fa+p))
}
x_vals <- seq(0, pop)

baseline_df <- data.frame(HF = x_vals,
  FOM = baselineFOM(pr, pop) * 100,
  type = "Uniform Line")

upper_df <- data.frame(HF = x_vals,
  FOM = upperFOM(pr, pop) * 100,
  type = "Upper Bound")

lower_df <- data.frame(HF = x_vals,
  FOM = lowerFOM(pr, pop) * 100,
  type = "Lower Bound")

fom_bounds <- rbind(baseline_df, upper_df, lower_df) 
max(roc_df_A$FOM*100-baseline_df$FOM)
max(roc_df_B$FOM*100-baseline_df$FOM)
# Plot both ROC curves
ggplot() +
  geom_line(data = roc_combined, aes(x = HF, y = FOM * 100, color = model, linetype = model), linewidth = 1.5) +
  geom_line(data = fom_bounds, aes(x = HF, y = FOM, color = type, linetype = type), linewidth = 1) +
  scale_color_manual(name = "", values = c(
    "A_index" = "#1f77b4", "B_index" = "#bb3115", "Uniform Line" = "gray70",
    "Upper Bound" = "grey1", "Lower Bound" = "grey40", name = NULL
  ),
  breaks = c("A_index", "B_index", "Uniform Line", "Upper Bound", "Lower Bound")) +
  scale_linetype_manual(name = "", values = c(
    "A_index" = "solid", "B_index" = "solid", "Uniform Line" = "dashed",
    "Upper Bound" = "twodash", "Lower Bound" = "longdash"
  , name = NULL),
  breaks = c("A_index", "B_index", "Uniform Line", "Upper Bound", "Lower Bound")) +
  scale_x_continuous(expand = c(0, 0), name = "False Positive Rate") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "True Positive Rate") +
  labs(legend = "") +
  coord_fixed(ratio = 10000 / 100) +
  my_theme

# 
# # Function to extend and compute FOM curves
# extend_FOM_curve <- function(tocd, model_label) {
#   a <- tocd
#   lastHF <- a$`Hits+FalseAlarms`[nrow(a) - 1]
#   lastH <- a$Hits[nrow(a) - 1]
#   
#   new_rows <- data.frame(matrix(0, nrow = 100, ncol = ncol(a)))
#   colnames(new_rows) <- colnames(a)
#   new_rows$`Hits+FalseAlarms` <- seq(lastHF, pop, length.out = 100)
#   new_rows$Hits <- seq(lastH, p, length.out = 100)
#   
#   a <- rbind(a[1:(nrow(a) - 1), ], new_rows, a[nrow(a), ])
#   a$Hits <- pmin(a$Hits, p)
#   a$Misses <- p - a$Hits
#   a$FOM <- a$Hits / (a$`Hits+FalseAlarms` + a$Misses)
#   a$model <- model_label
#   
#   return(a)
# }
# 
# # Generate A_index FOM curve
# xt <- xt_A
# yt <- yt_A
# a <- as.data.frame(matrix(0, nrow = pop + 1, ncol = 2))
# names(a) <- c("Hits+FalseAlarms", "FOM")
# a$`Hits+FalseAlarms` <- 0:pop
# a$Hits <- Curva(xt, yt)
# a_A <- extend_FOM_curve(a, "A_index")
# 
# # Generate B_index FOM curve
# xt <- xt_B
# yt <- yt_B
# a <- as.data.frame(matrix(0, nrow = pop + 1, ncol = 2))
# names(a) <- c("Hits+FalseAlarms", "FOM")
# a$`Hits+FalseAlarms` <- 0:pop
# a$Hits <- Curva(xt, yt)
# a_B <- extend_FOM_curve(a, "B_index")
# 
# # Combine both
# combined_fom <- rbind(a_A, a_B)
# combined_fom$baseline <- combined_fom$`Hits+FalseAlarms` * p / pop
# combined_fom$baselineFOM <- combined_fom$baseline / (combined_fom$`Hits+FalseAlarms` - combined_fom$baseline + p)
# combined_fom$upper <- pmin(combined_fom$`Hits+FalseAlarms`, p)
# combined_fom$lower <- pmax(0, combined_fom$`Hits+FalseAlarms` - (pop - p))
# combined_fom$upperFOM <- combined_fom$upper / (combined_fom$`Hits+FalseAlarms` - combined_fom$upper + p)
# combined_fom$lowerFOM <- combined_fom$lower / (combined_fom$`Hits+FalseAlarms` - combined_fom$lower + p)
# 
# # Plot both FOM curves
# ggplot(combined_fom, aes(x = `Hits+FalseAlarms` / 1e2, y = FOM * 100, color = model)) +
#   geom_line(size = 1.5) +
#   geom_line(
#     data = data.frame(HitsFalseAlarms = seq(p + 1, pop, length.out = 1000)),
#     aes(x = HitsFalseAlarms / 1e2, y = (p / HitsFalseAlarms) * 100),
#     linetype = "twodash", color = "gray50"
#   ) +
#   geom_line(aes(x = `Hits+FalseAlarms` / 1e2, y = baselineFOM * 100, color = "Baseline", linetype = "Baseline"),
#     size = 1
#   ) +
#   geom_line(aes(x = `Hits+FalseAlarms` / 1e2, y = upperFOM * 100, color = "Upper Bound", linetype = "Upper Bound"),
#     size = 2
#   ) +
#   geom_line(aes(x = `Hits+FalseAlarms` / 1e2, y = lowerFOM * 100, color = "Lower Bound", linetype = "Lower Bound"),
#     size = 2
#   ) +
#   labs(
#     x = "Hits + False Alarms (thousands)",
#     y = "CSI (FOM)",
#     title = "FOM Curves for A_index and B_index",
#     color = "Model"
#   ) +
#   scale_color_manual(
#     values = c("A_index" = "#1f77b4", "B_index" = "#bb3115", "Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), 
#       "Lower Bound" = "gray30"
#     )
#   ) +
#   scale_linetype_manual(
#     values = c("Baseline" = "solid","Upper Bound" = "dotted","Lower Bound" = "dashed", "A_index" = "solid","B_index" = "solid"
#     )
#   ) +
#   scale_x_continuous(
#     expand = c(0, 0), 
#     breaks = seq(0, pop / 1e2, by = pop / 5 / 1e2)
#   ) +
#   scale_y_continuous(
#     expand = c(0, 0), 
#     limits = c(0, 100), 
#     breaks = seq(0, 100, by = 20)
#   ) +
#   coord_fixed(ratio = (pop / 1e2) / 100) +
#   my_theme

for (pr in c(0.01, 0.05,0.1,0.25)) {
  xt_A=780
  yt_A=350
  
  xt_B=2960
  yt_B=460
  
  tabla<-read_xlsx(paste0("C:/Users/julen.perales/OneDrive - UPNA/Articulo/ASCI/Simulations/simulations_",pr*100, ".xlsx"))
  
  p<-pop*pr
  # Create a color palette for positive values (greens)
  # Your color palettes
  greens_palette <- colorRampPalette(c("grey85", "darkgreen"))(n = 1000)
  reds_palette <- colorRampPalette(c("grey85", "darkred"))(n = 1000)
  # # Replace values >= 1000 with 1000 in a new column
  # tabla$adjusted_L <- ifelse(tabla$xt >= (p*5), 0, ifelse(tabla$xt >= p, p*5,tabla$xt))
  tabla$adjusted_L_FOM <- ifelse(tabla$L+p >= (p*2), 2, (tabla$L+p)/p)
  tabla$adjusted_L_TOC <- ifelse(tabla$xt >= (p*2), 2, ifelse(tabla$xt >= (p*2),2,tabla$xt/p))
  # # Assign colors based on the specified conditions
  # tabla$color <- ifelse(tabla$L >= 1000, tail(greens_palette, 1),
  #   ifelse(tabla$L >= 0, reds_palette[cut(tabla$L/1000*1000, breaks = 1000)],
  #   ifelse(tabla$L >= -1000, greens_palette[cut(tabla$L/(min(tabla$L))*1000, breaks = 1000)], "darkgreen")))
  tabla$color_FOM <- ifelse(tabla$adjusted_L_FOM >= 2, tail(reds_palette, 1),
                        ifelse(tabla$adjusted_L_FOM >= 1, reds_palette[cut(c(1:2), breaks = 1000)],
                               ifelse(tabla$adjusted_L_FOM >= 0, greens_palette[cut(c(0.1,1), breaks = 1000)], "darkgreen")))
  tabla$color_TOC <- ifelse(tabla$adjusted_L_TOC >= 2, tail(reds_palette, 1),
                        ifelse(tabla$adjusted_L_TOC >= 1, reds_palette[cut(c(1,2), breaks = 1000)],
                               ifelse(tabla$adjusted_L_TOC >= 0, greens_palette[cut(c(0.1,1), breaks = 1000)], "darkgreen")))
  # # Create a dataframe for the legend
  # legend_data <- data.frame(value = c(min(tabla$adjusted_L), 0, 1000),
  #                           color = c("darkgreen", "grey85", tail(reds_palette, 1)))
  legend_data_FOM <- data.frame(value = c(0, 1, 2),
                            color = c("darkgreen", "grey85", tail(reds_palette, 1)))
  legend_data_TOC <- data.frame(value = c(0, 1, 2),
                            color = c("darkgreen", "grey85", tail(reds_palette, 1)))
  # # Adjust labels for the legend
  # # legend_labels <- c(p+min(tabla$adjusted_L), "P", ">=2*P")
  legend_labels_FOM <- c(0, "P",">2*P")
  legend_labels_TOC <- legend_labels_FOM
  
  # fDAUC<-function(FOM){(2-pr/FOM-pr*FOM)/(2*(1-pr))}
  
  #consideramos dos triangulos diferentes por que cambian los supuestos de MFOM
  maxAUC<-function(FOM){
    a<-FOM*p
    if (FOM>=pr){
    result<-(2*p*(E-p)-(p^2/a-p)*(p-a))/(2*p*(E-p))}
    else{
      result<-(2*p*(E-p)-(E-p)*(p-a))/(2*p*(E-p))
    }
    
    # result<-(p/2*(2*E-p/FOM-FOM*p))/(p*(E-p))
    return(result)
    }
  minAUC<-function(FOM){
    
    result<-(1+FOM)/2
    return(result)
  }
  tabla_backup<-tabla
  tabla$min<-sapply(tabla$MFOM/100,minAUC)*100
  tabla <- tabla[tabla$AUC >= tabla$min, ]
  
  
  #yt=666	xt=806
  #yt=918	xt=3338
  
  # Generate dataframes for MaxDAUC and MinDAUC
  sequencia<-seq(pr,1,length.out=300) 
  Range<-data.frame(
    MFOM=sequencia*100,
    Max=sapply(sequencia,maxAUC)*100,
    Min=pmax(sapply(sequencia,minAUC),pr)*100
  )
  Range$Max[1]<-Range$Min[1]<-50
  #EMPEZAMOS CON AGREEMENT
  #GRAFICO PARA AUC Y MAXFOM
  pl<-ggplot(tabla, aes(x = MFOM, y = AUC, color = adjusted_L_TOC)) +
  # ggplot(tabla, aes(x = MFOM, y = AUC)) +
    # geom_point(alpha = 0.5, color="blue") +
    geom_point(alpha = 0.8,size=2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
    # geom_line(data = df, aes(x = V1*100, y = V2*100), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MFOM, y = Max), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MFOM, y = Min), color = "black", size = 1.5) +
    # geom_point(aes(x = MFOM[A_index], y = AUC[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # geom_point(aes(x = MFOM[B_index], y = AUC[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    labs(x = "Maximun CSI (%)", y = "AUC (%)") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    annotate(geom = "text", x = tabla$MFOM[A_index], y = tabla$AUC[A_index], label = "Classifier A", vjust = 2, hjust = -0, size = 10, fontface = "bold") +
    annotate(geom = "text", tabla$MFOM[B_index], y = tabla$AUC[B_index], label = "Classifier B", vjust = -1, hjust = 0.75, size = 10, fontface = "bold") +
    
    scale_color_gradientn(
      colors = c("darkgreen", "grey85", tail(reds_palette, 1)),
      values = scales::rescale(c(0, 1, max(tabla$adjusted_L_FOM, tabla$adjusted_L_TOC))),
      breaks = legend_data_TOC$value,
      labels = legend_labels_TOC,
      limits = c(0, 2))+
      my_theme +
      coord_fixed() +
      guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))+
    theme(plot.margin = margin(t = 20, r = 20, b = 30, l = 20))
  
  print(pl)
  
  maxAFOM<-function(FOM){

    Tot<-p*(1/2-ln(pr)-pr/2)

    numerator<- -p*(pr/2-ln(FOM/pr)-1+FOM^2/2)
    return(numerator/Tot)
  }
  Fmin<-function(FOM){
    h<-function(x){pmin(x, (p-FOM*p)/(E-FOM*p)*(x-FOM*p)+FOM*p)}
    f<-function(x){h(x)/(x-h(x)+p)}
    numerator<-integrate(f, lower=0, upper=E)$value-p*pr/2
    Tot<-p*(1/2-ln(pr)-pr/2)
    return(numerator/Tot)
  }
  tabla$min<-sapply(tabla$MFOM/100,Fmin)*100
  tabla <- tabla[tabla$AFOM >= tabla$min, ]
  # Generate dataframes for MaxDAUC and MinDAUC
  sequencia<-seq(pr,1,length.out=300)
  Range<-data.frame(
    MFOM=sequencia*100,
    Max=pmax(sapply(sequencia,maxAFOM),pr)*100,
    Min=sapply(sequencia, Fmin)*100
  )
 
  A_index<-which(round(tabla$AUC,6)==round(89.5541931175204,6))
  B_index<-which(round(tabla$AUC,6)==round(89.5140094111686,6))
  
  #GRAFICO PARA AFOM Y MAXFOM
  pl<-ggplot(tabla, aes(x = MFOM, y = AFOM, color = adjusted_L_FOM)) +
  # ggplot(tabla, aes(x = MFOM, y = AFOM)) +
    geom_point(alpha = 1, size=2) +
    # geom_point(alpha = 0.5, color="blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
    # geom_line(data = df, aes(x = V1*100, y = V2*100), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MFOM, y = Max), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MFOM, y = Min), color = "black", size = 1.5) +
    # geom_point(aes(x = MFOM[A_index], y = AFOM[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # geom_point(aes(x = MFOM[B_index], y = AFOM[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # labs(x = "DMFOM (% points)", y = "DAUC (%)", title = "Difference with baseline of AUC and Maximum FOM") +
    labs(x = "Maximun CSI (%)", y = "AUCSI (%)") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    scale_color_gradientn(
      colors = c("darkgreen", "grey85", tail(reds_palette, 1)),
      values = scales::rescale(c(0, 1, max(tabla$adjusted_L_FOM, tabla$adjusted_L_TOC))),
      breaks = legend_data_TOC$value,
      labels = legend_labels_TOC,
      limits = c(0, 2)
    ) +
    annotate(geom = "text", x = tabla$MFOM[A_index], y = tabla$AFOM[A_index], label = "Classifier A", vjust = 2, hjust = -0, size = 10, fontface = "bold") +
    annotate(geom = "text", tabla$MFOM[B_index], y = tabla$AFOM[B_index], label = "Classifier B", vjust = -1, hjust = 0.9, size = 10, fontface = "bold") +
  
    my_theme +
    coord_fixed() +
    guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))+
    theme(plot.margin = margin(t = 20, r = 20, b = 30, l = 20))
  
  print(pl)
  UpperAUC<-function(FOM){
    x_A<-FOM*p
    y_A<-FOM
    S__A<-y_A-x_A*pr/((1-pr)*x_A+p)
    XB<-function(x){
      p/x-y_A-pr*x/((1-pr)*x+p)
    }
    result<-p
    # Define precision criterion
    precision <- 1e-8
    while(round(XB(result),8) > precision) {
      if (XB(result) > 0) {
        result<-result+XB(result)
      }
    }
    x_B<-result
    y_B<-p/x_B
    AUCS<-(p*(E-p)-(p-x_A)*(x_B-p)/2)/(p*(E-p))*2-1
  
    return(AUCS)
  }
  UpperAFOM<-function(FOM){
    x_A<-FOM*p
    y_A<-FOM
    S__A<-y_A-x_A*pr/((1-pr)*x_A+p)
    XB<-function(x){
      p/x-y_A-pr*x/((1-pr)*x+p)
    }
    result<-p
    # Define precision criterion
    precision <- 1e-8
    while(round(XB(result),8) > precision) {
      if (XB(result) > 0) {
        result<-result+XB(result)
      }
    }
    x_B<-result
    y_B<-p/x_B
  
    h<-function(x){pmin(x, (p-x_A)/(x_B-x_A)*(x-x_B)+p,p)}
    u<-function(x){pr*x/((1-pr)*x+p)}
    f<-function(x){h(x)/(x-h(x)+p)}
    numerator<-integrate(f, lower=0, upper=E)$value-p*pr/2
    Tot<-p*(1/2-ln(pr)-pr/2)
    AFOM<-numerator/Tot
    AUFOM<-(integrate(u, lower=0, upper=E)$value-p*pr/2)/Tot
    
    return(((AFOM-AUFOM)/(1-AUFOM)))
  }
  #SEGUIMOS CON SKILL
  #AUCS Y MAXJULEN
  tabla$Skill<-(tabla$MFOM-tabla$bMFOM)/(100-tabla$bMFOM)*100
  tabla$Skill<-tabla$Julen
  
  #Ahora "sequencia" hace referencia a Skill, por lo que hay que traducirla
  sequenciaSkill<-seq(0,1,length.out=300)
  sequenciaFOM<-sapply(sequenciaSkill, function(x){x+(pr*x)/((1-pr)*x+p)})
  
  Range<-data.frame(
    MSKILL=sequenciaSkill*100,
    Max=(sapply(sequenciaFOM,UpperAUC))*100,
    Min=(pmax(sapply(sequenciaFOM,minAUC),pr)*2-1)*100
  )
  
  pl<- ggplot(tabla, aes(x = Skill, y = DAUC, color=adjusted_L_TOC)) +
  # ggplot(tabla, aes(x = Skill, y = DAUC)) +
    # geom_point(alpha = 0.5, color="blue") +
    geom_point(alpha = 1, size=2)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
    geom_line(data = Range, aes(x = MSKILL, y = Max), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MSKILL, y = Min), color = "black", size = 1.5) +
    # geom_point(aes(x = Skill[A_index], y = DAUC[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # geom_point(aes(x = Skill[B_index], y = DAUC[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    labs(x = "CSI Skill (%)", y = "AUC Skill (%)") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    annotate(geom = "text", x = tabla$Skill[A_index], y = tabla$DAUC[A_index], label = "Classifier A", vjust = 2, hjust = -0, size = 10, fontface = "bold") +
    annotate(geom = "text", tabla$Skill[B_index], y = tabla$DAUC[B_index], label = "Classifier B", vjust = -1, hjust = 0.9, size = 10, fontface = "bold") +  
    scale_color_gradientn(
      colors = c("darkgreen", "grey85", tail(reds_palette, 1)),
      values = scales::rescale(c(0, 1, max(tabla$adjusted_L_FOM, tabla$adjusted_L_TOC))),
      breaks = legend_data_TOC$value,
      labels = legend_labels_TOC,
      limits = c(0, 2))+
    my_theme +
    coord_fixed() +
    guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))+
    theme(plot.margin = margin(t = 20, r = 20, b = 30, l = 20))
  print(pl)
  Range<-data.frame(
    MSKILL=sequenciaSkill*100,
    Max=(sapply(sequenciaFOM,UpperAFOM))*100,
    Min=(sapply(sequenciaFOM, Fmin)-max(tabla$AUFOM)/100)/(1-max(tabla$AUFOM)/100)*100
  )
  
  #AFOMS Y MAXJULEN
  pl<-ggplot(tabla, aes(x = Skill, y = DFOM, color = adjusted_L_FOM)) +
    # ggplot(tabla, aes(x = MFOM, y = AFOM)) +
    geom_point(alpha = 1, size=2) +
    # geom_point(alpha = 0.5, color="blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
    geom_line(data = Range, aes(x = MSKILL, y = Max), color = "black", size = 1.5) +
    geom_line(data = Range, aes(x = MSKILL, y = Min), color = "black", size = 1.5) +
    # geom_point(aes(x = Skill[A_index], y = DFOM[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # geom_point(aes(x = Skill[B_index], y = DFOM[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
    # labs(x = "DMFOM (% points)", y = "DAUC (%)", title = "Difference with baseline of AUC and Maximum FOM") +
    labs(x = "CSI Skill (% points)", y = "AUCSI Skill(%)") +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
  
    annotate(geom = "text", x = tabla$Skill[A_index], y = tabla$DFOM[A_index], label = "Classifier A", vjust = 2, hjust = -0, size = 10, fontface = "bold") +
    annotate(geom = "text", tabla$Skill[B_index], y = tabla$DFOM[B_index], label = "Classifier B", vjust = -1, hjust = 0.9, size = 10, fontface = "bold") +  
    scale_color_gradientn(
      colors = c("darkgreen", "grey85", tail(reds_palette, 1)),
      values = scales::rescale(c(0, 1, max(tabla$adjusted_L_FOM, tabla$adjusted_L_TOC))),
      breaks = legend_data_TOC$value,
      labels = legend_labels_TOC,
      limits = c(0, 2))+  
    my_theme +
    coord_fixed() +
    guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))+
    theme(plot.margin = margin(t = 20, r = 20, b = 30, l = 20))
print(pl)
  }


#CARGO LOS PLOTS DIRECTAMENTE PARA LA FIGURE DE LA TESIS. LOS CLASIFICADORES SON DIFERENTES A LOS ORIGINALES, PERO ME ESTOY VOLVIENDO LOCO
A<-readRDS(paste(folder_path,"/ClassifierA.rds",sep=""))
B<-readRDS(paste(folder_path,"/ClassifierB.rds",sep=""))

nuevo_dataframe <- data.frame(
  `Hits+FalseAlarms` = A$`Hits+FalseAlarms`, # Columna de A
  FOM_A = A$FOM,                             # Columna FOM de A
  FOM_B = B$FOM,                             # Columna FOM de B
  Upper = sapply(A$`Hits+FalseAlarms`, function(x) min(x / p, p / x)), # Calculated Upper column
  Lower = sapply(A$`Hits+FalseAlarms`, function(x) max(0, x / p - (pop - p))), # Calculated Lower column
  Uniform = sapply(A$`Hits+FalseAlarms`, function(val) (val * pr) / ((1 - pr) * val + (E * pr))),
  check.names = FALSE
)

ggplot(nuevo_dataframe, aes(x = `Hits+FalseAlarms`, y = FOM_A*100)) +
  geom_line(data=vline, aes(x=x, y=y, color = "P", linetype ="P"), size = 1.25) + 
  geom_line(aes(color = "Classifier A", linetype ="Classifier A"), size = 2) +
  geom_line(aes(x = `Hits+FalseAlarms`, y = FOM_B*100, color = "Classifier B", linetype ="Classifier B"), size = 2) +
  labs(x = "Hits+FalseAlarms", y = "CSI (as %)") +
  theme_minimal() +
  geom_line(aes(x = `Hits+FalseAlarms`, y = Upper * 100, colour = "Upper Bound", linetype ="Upper Bound"), size = 1.5) +
  geom_line( aes(x = `Hits+FalseAlarms`, y = Uniform * 100, colour = "Uniform Line",linetype ="Uniform Line"), size = 1.5) +
  geom_line(aes(x = `Hits+FalseAlarms`, y = Lower * 100, colour = "Lower Bound",linetype ="Lower Bound"), size = 1.5) +
  scale_x_continuous(limits = c(0, pop), expand = c(0, 0), breaks = seq(0,pop,by=2000)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=20)) +
  scale_color_manual(name = "", values = c("Classifier A" = "#BB3115","Classifier B" = "#074F69","P"="grey40","Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), "Uniform Line"="grey", "Lower Bound" = rgb(79, 129, 189, maxColorValue = 255)),
                     breaks = c("Classifier A", "Classifier B", "P", "Uniform Line", "Upper Bound", "Lower Bound")) +
  scale_linetype_manual(name = "", values = c("Classifier A" = "solid","Classifier B" = "solid", "P"="longdash","Upper Bound" = "twodash", "Uniform Line" = "dashed", "Lower Bound" = "longdash"),
                        breaks = c("Classifier A", "Classifier B", "P", "Uniform Line", "Upper Bound", "Lower Bound")) +
  
  my_theme +
  coord_fixed(ratio=pop/100)+
  guides(fill = guide_legend(byrow = TRUE))


nuevo_dataframe <- data.frame(
  `Hits+FalseAlarms` = A$`Hits+FalseAlarms`, # Columna de A
  TOC_A = A$Acc,                             # Columna FOM de A
  TOC_B = B$Acc,                             # Columna FOM de B
  Upper=sapply(A$`Hits+FalseAlarms`,function(val) min(val, p)),
  Lower=sapply(A$`Hits+FalseAlarms`,function(val) max(0, val-(E-p) )),
  Uniform=sapply(A$`Hits+FalseAlarms`,function(val) {(val * pr)}),
  check.names = FALSE
)


ggplot(nuevo_dataframe, aes(x = `Hits+FalseAlarms`, y = TOC_A)) +
  geom_line(aes(color = "Classifier A", linetype ="Classifier A"), size = 2) +
  geom_line(aes(x = `Hits+FalseAlarms`, y = TOC_B, color = "Classifier B", linetype ="Classifier B"), size = 2) +
  labs(x = "Hits+FalseAlarms", y = "Hits") +
  theme_minimal() +
  geom_line(data=vline, aes(x=x, y=y*100, color = "P", linetype ="P"), size = 1.25) + 
  geom_line(aes(x = `Hits+FalseAlarms`, y = Upper, colour = "Upper Bound", linetype ="Upper Bound"), size = 1.5) +
  geom_line( aes(x = `Hits+FalseAlarms`, y = Uniform, colour = "Uniform Line",linetype ="Uniform Line"), size = 1.5) +
  geom_line(aes(x = `Hits+FalseAlarms`, y = Lower, colour = "Lower Bound",linetype ="Lower Bound"), size = 1.5) +
  scale_x_continuous(limits = c(0, pop), expand = c(0, 0), breaks = seq(0,pop,by=2000)) +
  scale_y_continuous(limits = c(0, p), expand = c(0, 0), breaks = seq(0,p,by=200)) +
  scale_color_manual(name = "", values = c("Classifier A" = "#BB3115","Classifier B" = "#074F69","P"="grey40","Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), "Uniform Line"="grey", "Lower Bound" = rgb(79, 129, 189, maxColorValue = 255)),
                     breaks = c("Classifier A", "Classifier B", "P", "Uniform Line", "Upper Bound", "Lower Bound")) +
  scale_linetype_manual(name = "", values = c("Classifier A" = "solid","Classifier B" = "solid", "P"="longdash","Upper Bound" = "twodash", "Uniform Line" = "dashed", "Lower Bound" = "longdash"),
                        breaks = c("Classifier A", "Classifier B", "P", "Uniform Line", "Upper Bound", "Lower Bound")) +
  
  my_theme +
  coord_fixed(ratio=pop/p)+
  guides(fill = guide_legend(byrow = TRUE))

