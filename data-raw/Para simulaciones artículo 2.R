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
  # print(ggplot(df, aes(x=Sim))+
  #   geom_line(aes(y=Upper))+
  #   geom_line(aes(y=Lower))+
  #   geom_line(aes(y=Curve)))
  # plot(c(x1,x2),c(sapply(x1,y1),sapply(x2,y2)))
  return(c(sapply(x1,y1),sapply(x2,y2)))
}


tabla<-as.data.frame(matrix(0,1,ncol = 14))
names(tabla)<-c("xt","yt","AFOM","AUFOM","AUC","DAUC", "DFOM", "MFOM","L","bMFOM","ThresholdMFOM","GSS","Julen", "ThresholdMSKILL")
pb <- txtProgressBar(min = 1, max = E, style = 3)

for(xt in seq(1,E, by=20)){
  # print((xt/E)*100)
  setTxtProgressBar(pb, xt)
  for(yt in seq(1,p,by=2)){
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
    plot(a$`Hits+FalseAlarms`, a$Acc)
    
    tabla[nrow(tabla)+1,]<-c(xt,yt,AFOM,AUFOM,AUC,DAUC,DFOM,MFOM,which(a["FOM"]==max(a["FOM"]))-p, bMFOM, thr, a$GSS[thr], a$Julen[thr], thrSkill)
   
  }
}
}
#   write_xlsx(tabla,"C:/Users/julen.perales/OneDrive - UPNA/Tesis/Simulaciones.xlsx")
# tabla <- tabla[tabla$MFOM >= pr*100, ]
tabla <- tabla[tabla$DAUC >=0, ]
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
# legend_labels_TOC <- legend_labels_FOM

# fDAUC<-function(FOM){(2-pr/FOM-pr*FOM)/(2*(1-pr))}

#consideramos dos triangulos diferentes por que cambian los supuestos de MFOM
maxAUC<-function(FOM){
  a<-FOM*p
  if (FOM>=pr){
  result<-(2*p*q-(p^2/a-p)*(p-a))/(2*p*q)}
  else{
    result<-(2*p*q-(E-p)*(p-a))/(2*p*q)
  }
  
  # result<-(p/2*(2*E-p/FOM-FOM*p))/(p*q)
  return(result)
  }
minAUC<-function(FOM){
  
  result<-(1+FOM)/2
  return(result)
}
tabla_backup<-tabla
tabla$min<-sapply(tabla$MFOM/100,minAUC)*100
tabla <- tabla[tabla$AUC >= tabla$min, ]
#este pertenece a xt=2472 y yt=820
A_index<-which(round(tabla$AUC,5)==round(89.5541931175204,5))
#este pertenece a xt=3434 y yt=930
B_index<-which(round(tabla$AUC,4)==round(89.5140094111686,4))

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
ggplot(tabla, aes(x = MFOM, y = AUC, color = adjusted_L_TOC)) +
# ggplot(tabla, aes(x = MFOM, y = AUC)) +
  # geom_point(alpha = 0.5, color="blue") +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
  # geom_line(data = df, aes(x = V1*100, y = V2*100), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MFOM, y = Max), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MFOM, y = Min), color = "black", size = 1.5) +
  # geom_point(aes(x = MFOM[A_index], y = AUC[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # geom_point(aes(x = MFOM[B_index], y = AUC[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  labs(x = "MFOM (% points)", y = "AUC (%)") +
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
  theme(
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
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white")
  ) +
  coord_fixed() +
  guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))

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
ggplot(tabla, aes(x = MFOM, y = AFOM, color = adjusted_L_FOM)) +
# ggplot(tabla, aes(x = MFOM, y = AFOM)) +
  geom_point(alpha = 0.5) +
  # geom_point(alpha = 0.5, color="blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
  # geom_line(data = df, aes(x = V1*100, y = V2*100), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MFOM, y = Max), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MFOM, y = Min), color = "black", size = 1.5) +
  geom_point(aes(x = MFOM[A_index], y = AFOM[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  geom_point(aes(x = MFOM[B_index], y = AFOM[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # labs(x = "DMFOM (% points)", y = "DAUC (%)", title = "Difference with baseline of AUC and Maximum FOM") +
  labs(x = "MFOM (% points)", y = "AFOM (%)") +
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

  theme(
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
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white")
  ) +
  coord_fixed() +
  guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))

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
  AUCS<-(p*q-(p-x_A)*(x_B-p)/2)/(p*q)*2-1

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
# 
# sequenciaFOM<-seq(0,1,length.out=300)
# sequenciaSkill<-(sapply(sequenciaFOM, function(x){x-(pr*x*p)/((1-pr)*x*p+p)}))/
#                 (sapply(sequenciaFOM, function(x){1-(pr*x*p)/((1-pr)*x*p+p)}))
# 
# sequenciaFOM1<-sapply(seq(0,p,length.out=300), x/p)
# sequenciaFOM2<-sapply(seq(p,E,length.out=300), p/x)
# sequenciaSkill1<-(sapply(sequenciaFOM1, function(x){x-(pr*x*p)/((1-pr)*x*p+p)}))/
#   (sapply(sequenciaFOM, function(x){1-(pr*x*p)/((1-pr)*x*p+p)}))
# sequenciaSkill2<-(sapply(sequenciaFOM2, function(x){x-(pr*p/x)/((1-pr)*x*p+p)}))/
#   (sapply(sequenciaFOM, function(x){1-(pr*x*p)/((1-pr)*p/x+p)}))

Range<-data.frame(
  MSKILL=sequenciaSkill*100,
  Max=(sapply(sequenciaFOM,UpperAUC))*100,
  Min=(pmax(sapply(sequenciaFOM,minAUC),pr)*2-1)*100
)

ggplot(tabla, aes(x = Skill, y = DAUC, color=adjusted_L_TOC)) +
# ggplot(tabla, aes(x = Skill, y = DAUC)) +
  # geom_point(alpha = 0.5, color="blue") +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
  geom_line(data = Range, aes(x = MSKILL, y = Max), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MSKILL, y = Min), color = "black", size = 1.5) +
  # geom_point(aes(x = Skill[A_index], y = DAUC[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # geom_point(aes(x = Skill[B_index], y = DAUC[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  labs(x = "FOM Skill (%)", y = "AUC Skill (%)") +
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
  theme(
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
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white")
  ) +
  coord_fixed() +
  guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))

Range<-data.frame(
  MSKILL=sequenciaSkill*100,
  Max=(sapply(sequenciaFOM,UpperAFOM))*100,
  Min=(sapply(sequenciaFOM, Fmin)-max(tabla$AUFOM)/100)/(1-max(tabla$AUFOM)/100)*100
)

#AFOMS Y MAXJULEN
ggplot(tabla, aes(x = Skill, y = DFOM, color = adjusted_L_FOM)) +
  # ggplot(tabla, aes(x = MFOM, y = AFOM)) +
  geom_point(alpha = 0.5) +
  # geom_point(alpha = 0.5, color="blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
  geom_line(data = Range, aes(x = MSKILL, y = Max), color = "black", size = 1.5) +
  geom_line(data = Range, aes(x = MSKILL, y = Min), color = "black", size = 1.5) +
  # geom_point(aes(x = Skill[A_index], y = DFOM[A_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # geom_point(aes(x = Skill[B_index], y = DFOM[B_index]), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # labs(x = "DMFOM (% points)", y = "DAUC (%)", title = "Difference with baseline of AUC and Maximum FOM") +
  labs(x = "FOM SKill (% points)", y = "AFOM Skill(%)") +
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
  theme(
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
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white")
  ) +
  coord_fixed() +
  guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))


t<-seq(0,pop,length=max(pop,1000))
U<-trapz(t,pr*t/((1-pr)*t+p))
TOT<-trapz(t,pmin(t/p,p/t))

# MaxAFOM<-function(FOM){max((2-pr-FOM^2+2*ln(FOM/pr))/(1-pr-2*ln(pr)),U/TOT)}
# # MaxAFOM<-function(FOM){max((2-pr-FOM^2+2*ln(FOM/pr))/(1-pr-2*ln(pr)),U)}
# # tabla$MaxDAFOM<-(sapply(tabla$MFOM/100,MaxAFOM)-tabla$AUFOM/100)/(1-tabla$AUFOM/100)
# # tabla$MaxDAFOM<-pmin(tabla$MaxDAFOM,tabla$DFOM)
# # tabla$MaxDAFOM<-sapply(tabla$MFOM/100,MaxAFOM)
# 
H<-function(FOM){
  if (FOM>pr){
  a<-integral(function(x){pmax(pr*x,pmin(FOM*x,p))/(x-pmax(pr*x,pmin(FOM*x,p))+p)},0,pop)
  }else{a<-U}
  return(a)}
S<-function(FOM){
  # if (FOM>pr){
    a<-integral(function(x){(pmin(x,FOM*p)+pmax((p-FOM*p)/(pop-FOM*p)*(x-FOM*p),0))/(x-(pmin(x,FOM*p)+pmax((p-FOM*p)/(pop-FOM*p)*(x-FOM*p),0))+p)},0,pop)
  # }else{a<-U}
  return(a)}
#

h <- pbsapply(sequencia, H)
s <- pbsapply(sequencia, S)
# Fmax_df <- data.frame(MFOM = sequencia*100, MaxDAUC = (sapply(sequencia,MaxAFOM)-U/TOT)/(1-U/TOT))
# Fmin_df <- data.frame(MFOM = sequencia*100, MinDAUC = (pmin(h,s)-U)*100/(TOT-U))
# rm(x)
sequencia<-seq(0,1,length.out=100) 
FOM<-function(x){min(x/p,B(x)+DFOM,p/x)}
H<-function(x){(B(x)+DFOM)*(x+p)/(1+B(x)+DFOM)}
se<-seq(0,1-pr/(2-pr),length.out=100)

Fmax_df <- data.frame(MFOM = se*100, MaxDAUC = 0)
Fmin_df <- data.frame(MFOM = (sequencia-sapply(sequencia*p,B))*100, MinDAUC = (s-U)*100/(TOT-U))
i<-1
for (DFOM in seq(0,1-pr/(2-pr),length.out=100)) {
  FOM<-function(x){min(x/p,B(x)+DFOM,p/x)}
  Fmax_df$MaxDAUC[[i]]<-(integral(Vectorize(FOM), 0,pop)-U)/(TOT-U)
  i<-i+1
}

# tabla$MinDAFOM<-pmin(tabla$H,tabla$S)
# tabla$MinDAFOM<-(tabla$MinDAFOM-U)/(TOT-U)

plot5 <- ggplot(tabla, aes(x = MFOM, y = DFOM, color = adjusted_L_FOM)) +
  geom_point(alpha=0.01) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", size = 1) +
  geom_line(data = Fmax_df, aes(x = MFOM, y = MaxDAUC * 100), color = "black", size = 1.5) +
  geom_line(data = Fmin_df, aes(x = MFOM, y = MinDAUC), color = "black", size = 1.5) +
  geom_point(data = data.frame(x = 56, y = 70), aes(x = x, y = y), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # geom_text(label="classifier A") +  # Add text without aes()
  geom_point(data = data.frame(x = 27, y = 58), aes(x = x, y = y), shape = 21, colour = "black", fill = "white", size = 5, stroke = 3) +  # Triangle shape
  # labs(x = "DMFOM (% points)", y = "DAFOM (%)", title = "Difference with baseline of AFOM and Maximun FOM") +
  labs(x = "DMFOM (% points)", y = "DAFOM (%)") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10))+
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    # sec.axis = sec_axis(~., name = "AFOM (%)", labels = scales::percent_format())
  ) +
  scale_color_gradientn(
    colors = c("darkgreen", "grey85", tail(reds_palette, 1)),
    values = scales::rescale(c(0, 1, max(tabla$adjusted_L_FOM, tabla$adjusted_L_TOC))),
    breaks = legend_data_FOM$value,
    labels = legend_labels_FOM,
    limits = c(0, 2)
  ) +
  annotate(geom = "text", x = 54, y = 68, label = "Classifier A", vjust = -1, hjust = -0.5, size = 10, fontface = "bold") +
  annotate(geom = "text", x = 25, y = 56, label = "Classifier B", vjust = -1, hjust = -0.5, size = 10, fontface = "bold") +
  theme(
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
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    # Add the guide_colorbar() function here for consistent formatting
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),  # Removes border around legend
    legend.background = element_rect(fill = "white", color = "white"),  # Sets legend background color
    plot.margin = unit(c(10,3,3,3), "lines")
  ) +
  coord_fixed()+
  guides(color = guide_colorbar(title = "Tipping Point \n(in times of P)", title.position = "top"))

# Combine the two plots into a single layout
plot4 <- plot4 + theme(legend.position = "none")
combined_plot <- plot4 + plot5 + plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a" ,tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 28,face="bold"))
# combined_plot <- plot4 + plot5
# Display the combined plot
combined_plot

pl <- pl + theme(legend.position = "none")
combined_plot <- pl + plot1 + plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a" ,tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 28,face="bold"))
combined_plot

plot_toc<-function(a,prevalence, population){
  tocd_a<-as.data.frame(matrix(0,ncol = 3,nrow=nrow(a)))
  tocd_a[2]<-a[1]
  tocd_a[3]<-a[5]
  tocd_a[1]<-seq(1,0,length=nrow(tocd_a))
  names(tocd_a)<-c("Threshold","Hits+FalseAlarms","Hits")
  tocd_a$`Hits+FalseAlarms`<-as.numeric(tocd_a$`Hits+FalseAlarms`)
  tocd_c<-as.data.frame(matrix(0,ncol = 3,nrow=nrow(c)))
  tocd_c[2]<-c[1]
  tocd_c[3]<-c[5]
  tocd_c[1]<-seq(1,0,length=nrow(tocd_c))
  names(tocd_c)<-c("Threshold","Hits+FalseAlarms","Hits")
  nticks<-pop
  digits<-2
  tocd_c$`Hits+FalseAlarms`<-as.numeric(tocd_c$`Hits+FalseAlarms`)
  
  pl <- ggplot(tocd_a, aes(x = `Hits+FalseAlarms`, y = Hits)) +
    # Scatter plot
    geom_line(aes(x = `Hits+FalseAlarms`, y = Hits, color = "Classifier A", linetype = "Classifier A"), size = 2) +
    geom_line(data=tocd_c, aes(x = `Hits+FalseAlarms`, y = Hits, color = "Classifier B", linetype = "Classifier B"), size = 2) +
    
    # Maximum line
    geom_line(data = data.frame(x = c(0, population, prevalence * population),
                                y = c(0, prevalence * population, prevalence * population)),
              aes(x = x, y = y, colour = "Upper Bound", linetype ="Upper Bound"), size = 2) +
    # Minimum line
    geom_line(data = data.frame(x = c(0, population, (1 - prevalence) * population),
                                y = c(0, prevalence * population, 0)),
              aes(x = x, y = y, colour = "Lower Bound", linetype ="Lower Bound"), size = 2) +
    # # Hits+Misses line
    # geom_line(aes(x = `Hits+FalseAlarms`, y = rep(prevalence * population, nrow(tocd)), colour = "Upper Bound", linetype ="Upper Bound"),
    #           size = 5, col = rgb(146, 208, 80, maxColorValue = 255)) +
    # Uniform line
    geom_line(data = data.frame(x = c(0, population),
                                y = c(0, prevalence * population)),
              aes(x = x, y = y,  colour = "Uniform Line", linetype ="Uniform Line"), size = 2) +
    geom_line(data=data.frame(x = p, y = seq(0,p)), aes(x=x, y=y, color = "P", linetype ="P"), size = 1.25) + 
    
    # Set x and y axis limits without expanding
    scale_x_continuous(limits = c(0, population), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, prevalence * population), expand = c(0, 0)) +
    scale_color_manual(name = "", values = c("Classifier A" = "darkgreen","Classifier B" = "blue","P"="red","Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), "Uniform Line"="grey", "Lower Bound" = rgb(79, 129, 189, maxColorValue = 255))) +
    scale_linetype_manual(name = "", values = c("Classifier A" = "solid","Classifier B" = "solid", "P"="longdash","Upper Bound" = "twodash", "Uniform Line" = "dashed", "Lower Bound" = "twodash")) +  
    
    # Customize the theme
    theme_minimal() +
    # Add labels
    labs(x = "Hits+False Alarms", y = "Hits") +
    # Format the theme
    theme(
      legend.key.width = unit(3, "line"),
      # legend.spacing.y = unit(10, 'cm'),
      # legend.box.spacing = unit(10, "cm"),
      legend.key.height = unit(4, "line"), # Adjust the distance between lines
      
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
      plot.subtitle = element_text(size = 18, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 20),
      axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
      axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
      axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
      
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 24, hjust = 0, vjust = 5),
      legend.key = element_rect(color = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.text = element_text(size = 28),
      legend.position = "right"
    )+
    coord_fixed(ratio=1/prevalence)
    # guides(color = guide_colorbar(title = "Location of \nthe peak \nin times of P", title.position = "top"))
  
  # Display the plot
  print(pl)
 }

library(ggplot2)
library(cowplot) # for plot_grid function

pl <- pl + theme(legend.position = "none")
combined_plot <- pl + plot1 + plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a" ,tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 28,face="bold"))
combined_plot

# Display the combined plot
print(combined_plot)

b<-plot+geom_line(data = a_df, aes(x = x, y = y), color = "orange", size = 1.5)
c<-pl+geom_line(data = tocd, aes(x = `Hits+FalseAlarms`, y = Hits), color = "orange", size = 1.5)

#para buscar curvas para el articulo
# Assuming 'tabla' is your data frame

tabla$diff<-tabla$DAUC-tabla$MFOM
# Filter values of tabla$DAUC rounded with no decimals equal to 90
filtered_values <- tabla$DAUC[round(tabla$DAUC) == 70]
filtered_diff <- tabla$diff[round(tabla$DAUC)== 70]
r<-round(length(filtered_values)*0.05)
sorted_diff <- sort(filtered_diff)
sorted_diff2 <- sort(filtered_diff, decreasing = TRUE)

# Get the 10th maximum value of tabla$diff
tenth_max_value <- sorted_diff2[r]

# Get the 10th minimum value of tabla$diff
tenth_min_value <- sorted_diff[r]
mean_diff <- mean(filtered_diff)

# Find the index of the smallest value of tabla$diff corresponding to filtered_values
index1 <- which.min(tabla$diff[tabla$DAUC %in% filtered_values])
index1 <- which(round(filtered_diff) == round(tenth_min_value))[1]
index2 <- which.max(tabla$diff[tabla$DAUC %in% filtered_values])
index2 <- which(round(filtered_diff) == round(tenth_max_value))[1]

# Get the corresponding value of tabla$DAUC
selected_DAUC1 <- filtered_values[index1]
selected_DAUC2 <- filtered_values[index2]

# Get the corresponding tabla$yt and tabla$xt values
yt1 <- tabla$yt[tabla$DAUC == selected_DAUC1][1]
xt1 <- tabla$xt[tabla$DAUC == selected_DAUC1][1]
yt2 <- tabla$yt[tabla$DAUC == selected_DAUC2][1]
xt2 <- tabla$xt[tabla$DAUC == selected_DAUC2][1]

values<-list(c(xt1,yt1),c(xt2,yt2))
pl<-list()
prints<-list()
for (i in 1:length(values)) {
  xt<-values[[i]][[1]]
  yt<-values[[i]][[2]]
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
    if (result<0) {
      return(-1)
    }else{
      return(round(result,6))
    }
  }
  
  
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
  # 
  interval <- c(-1, 0)
  # Define precision criterion
  precision <- 1e-16
  while(interval[2] - interval[1] > precision) {
    # Use optimize to find the minimum (closest to zero) of the absolute value of f2(x) within the current interval
    root <- optimize(f = f2, interval = interval)
    # print(root)
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
  
  # gc()
  a<-as.data.frame(matrix(0,nrow=pop+1,ncol = 2))
  names(a)<-c("Hits+FalseAlarms", "FOM")
  a[1]<-c(0:pop)
  a["Acc"]<-0
  a$Acc[1:xt]<-sapply(a$`Hits+FalseAlarms`[1:xt],y1)
  a$Acc[xt:pop+1]<-sapply(a$`Hits+FalseAlarms`[xt:pop+1],y2)
  a$Acc<-pmin(a$`Hits+FalseAlarms`,a$Acc)
  a$FOM<-a$Acc/(a$`Hits+FalseAlarms`-a$Acc+p)
  
  a["baseline"]<-(a$`Hits+FalseAlarms`*pr)/((1-pr)*a$`Hits+FalseAlarms`+p)
  a["baseline_change"]<-c(0,diff(a$"baseline"))
  a["change"]<-c(0,diff(a$"FOM"))
  
  L_FOM<-which(a$change - a$baseline_change < 0)[1]+1
  L_TOC<-xt
  AUC<-0
  b<-a["Hits+FalseAlarms"]-a["Acc"]
  
  for (t in 1:(nrow(a)-1)){
    AUC<-AUC+((a[["Acc"]][t+1]+a[["Acc"]][t])/p)*(b[[t+1,1]]-b[[t,1]])/(pop-p)/2
  }
  
  lista<-FOM(a,p,pop)
  
  AFOM =lista[[1]]*100
  AUFOM =lista[[2]]*100
  AUC =AUC*100
  DAUC =(AUC/100-0.5)*2*100
  DFOM =((lista[[1]]-lista[[2]])/(1-lista[[2]]))*100
  MFOM = (max(a["FOM"]))*100
  pl[[i]]<-a
 
  prints<-append(prints,print(paste("AFOM =",round(AFOM))))
  prints<-append(prints,print(paste("DAFOM =",round(DFOM))))
  prints<-append(prints,print(paste("AUC =",round(AUC))))
  prints<-append(prints,print(paste("FOM =",round(MFOM))))
  prints<-append(prints,print(paste("MAX =", which(a$change < 0)[1]+1)))
  
}
vline <- data.frame(x = p, y = seq(0,100))

# plot1<-ggplot(a_df, aes(x = x, y = y)) +
plot1<-ggplot(pl[[1]], aes(x = `Hits+FalseAlarms`, y = FOM*100)) +
  geom_line(aes(color = "Classifier A", linetype ="Classifier A"), size = 2) +
  geom_line(data=pl[[2]],aes(x = `Hits+FalseAlarms`, y = FOM*100, color = "Classifier B", linetype ="Classifier B"), size = 2) +
  labs(x = "Hits+FalseAlarms", y = "FOM (as %)") +
  theme_minimal() +
  geom_line(data=vline, aes(x=x, y=y, color = "P", linetype ="P"), size = 1.25) + 
  geom_line(data = lista[[3]], aes(x = `Hits+FalseAlarms`, y = FOM * 100, colour = "Upper Bound", linetype ="Upper Bound"), size = 1.5) +
  geom_line(data = lista[[4]], aes(x = `Hits+FalseAlarms`, y = FOM * 100, colour = "Uniform Line",linetype ="Uniform Line"), size = 1.5) +
  geom_line(data = lista[[5]], aes(x = `Hits+FalseAlarms`, y = FOM * 100, colour = "Lower Bound",linetype ="Lower Bound"), size = 1.5) +
  scale_x_continuous(limits = c(0, pop), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0,100,by=10)) +
  scale_color_manual(name = "", values = c("Classifier A" = "darkgreen","Classifier B" = "blue","P"="red","Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), "Uniform Line"="grey", "Lower Bound" = rgb(79, 129, 189, maxColorValue = 255))) +
  scale_linetype_manual(name = "", values = c("Classifier A" = "solid","Classifier B" = "solid", "P"="longdash","Upper Bound" = "twodash", "Uniform Line" = "dashed", "Lower Bound" = "twodash")) +
  
  theme(
    legend.key.width = unit(3, "line"),
    # legend.spacing.y = unit(10, 'cm'),
    # legend.box.spacing = unit(10, "cm"),
    legend.key.height = unit(4, "line"), # Adjust the distance between lines
    
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
    axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
    axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
    
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = 28),
    legend.position = "right"
  )+
  coord_fixed(ratio=pop/100)+
  guides(fill = guide_legend(byrow = TRUE))

prevalence<-pr
population<-pop
a<-pl[[1]]
c<-pl[[2]]
tocd_a<-as.data.frame(matrix(0,ncol = 3,nrow=nrow(a)))
tocd_a[2]<-a[1]
tocd_a[3]<-a[3]
tocd_a[1]<-seq(1,0,length=nrow(tocd_a))
names(tocd_a)<-c("Threshold","Hits+FalseAlarms","Hits")
tocd_a$`Hits+FalseAlarms`<-as.numeric(tocd_a$`Hits+FalseAlarms`)
tocd_c<-as.data.frame(matrix(0,ncol = 3,nrow=nrow(c)))
tocd_c[2]<-c[1]
tocd_c[3]<-c[3]
tocd_c[1]<-seq(1,0,length=nrow(tocd_c))
names(tocd_c)<-c("Threshold","Hits+FalseAlarms","Hits")
nticks<-pop
digits<-2
tocd_c$`Hits+FalseAlarms`<-as.numeric(tocd_c$`Hits+FalseAlarms`)

plot2<-ggplot(tocd_a, aes(x = `Hits+FalseAlarms`, y = Hits)) +
  # Scatter plot
  geom_line(aes(x = `Hits+FalseAlarms`, y = Hits, color = "Classifier A", linetype = "Classifier A"), size = 2) +
  geom_line(data=tocd_c, aes(x = `Hits+FalseAlarms`, y = Hits, color = "Classifier B", linetype = "Classifier B"), size = 2) +
  
  # Maximum line
  geom_line(data = data.frame(x = c(0, population, prevalence * population),
                              y = c(0, prevalence * population, prevalence * population)),
            aes(x = x, y = y, colour = " Upper Bound", linetype =" Upper Bound"), size = 2) +
  # Minimum line
  geom_line(data = data.frame(x = c(0, population, (1 - prevalence) * population),
                              y = c(0, prevalence * population, 0)),
            aes(x = x, y = y, colour = " Lower Bound", linetype =" Lower Bound"), size = 2) +
  # # Hits+Misses line
  # geom_line(aes(x = `Hits+FalseAlarms`, y = rep(prevalence * population, nrow(tocd)), colour = "Upper Bound", linetype ="Upper Bound"),
  #           size = 5, col = rgb(146, 208, 80, maxColorValue = 255)) +
  # Uniform line
  geom_line(data = data.frame(x = c(0, population),
                              y = c(0, prevalence * population)),
            aes(x = x, y = y,  colour = "Uniform Line", linetype ="Uniform Line"), size = 2) +
  geom_line(data=data.frame(x = p, y = seq(0,p)), aes(x=x, y=y, color = "P", linetype ="P"), size = 1.25) + 
  # Set x and y axis limits without expanding
  scale_x_continuous(limits = c(0, population), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, prevalence * population), expand = c(0, 0)) + 
  scale_color_manual(name = "", values = c("Classifier A" = "darkgreen","Classifier B" = "blue","P"="red"," Upper Bound" = rgb(79, 129, 189, maxColorValue = 255), "Uniform Line"="grey", " Lower Bound" = rgb(79, 129, 189, maxColorValue = 255))) +
  scale_linetype_manual(name = "", values = c("Classifier A" = "solid","Classifier B" = "solid", "P"="longdash", " Upper Bound" = "twodash", "Uniform Line" = "dashed", " Lower Bound" = "twodash")) +  
  
  # Customize the theme
  theme_minimal() +
  # Add labels
  labs(x = "Hits+False Alarms", y = "Hits") +
  # Format the theme
  theme(
    legend.key.width = unit(3, "line"),
    # legend.spacing.y = unit(10, 'cm'),
    # legend.box.spacing = unit(10, "cm"),
    legend.key.height = unit(4, "line"), # Adjust the distance between lines
    
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
    axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
    axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
    
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = 28),
    legend.position = "right"
  )+
  coord_fixed(ratio=1/prevalence)

plot1 <- plot1 + theme(legend.position = "none")
combined_plot <- plot2 + plot1 + plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a" ,tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 28,face="bold"))
  
# combined_plot <- plot4 + plot5
# Display the combined plot
combined_plot
prints


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
  
  theme(
    legend.key.width = unit(3, "line"),
    # legend.spacing.y = unit(10, 'cm'),
    # legend.box.spacing = unit(10, "cm"),
    legend.key.height = unit(4, "line"), # Adjust the distance between lines
    
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
    axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
    axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
    
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = 28),
    legend.position = "right"
  )+
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
  
  theme(
    legend.key.width = unit(3, "line"),
    # legend.spacing.y = unit(10, 'cm'),
    # legend.box.spacing = unit(10, "cm"),
    legend.key.height = unit(4, "line"), # Adjust the distance between lines
    
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"inches"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 28, face = "bold", vjust = -1),
    axis.title.y = element_text(size = 28, face = "bold", vjust = 1),
    axis.title.y.right = element_text(size = 18, face = "bold", vjust = 1),
    
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 24, hjust = 0, vjust = 5),
    legend.key = element_rect(color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = 28),
    legend.position = "right"
  )+
  coord_fixed(ratio=pop/p)+
  guides(fill = guide_legend(byrow = TRUE))

