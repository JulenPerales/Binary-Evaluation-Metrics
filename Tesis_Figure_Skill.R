E<-100
Pr<-0.3
P<-E*Pr
Q<-E-P
# l<-50/(65*0.45^2)
# f1<-function(x){(-l*x^2+l*0.9*x)*0.65}
# f2<-function(x){2.46*exp(-4.4*x)+0.9*0.3}
Upper<-function(x){min(x/P,P/x)}
Lower<-function(x){max(0,(x-(E-P))/E)}
Random<-function(x){Pr*x/((1-Pr)*x+P)}
a <- readRDS(file = "C:/Users/julen.perales/OneDrive - UPNA/Tesis/Tabla_clasificador_Skill.rds")
Curve<-a$FOM

# x1<-seq(0,0.6,length.out=100)
# x2<-seq(0.6,1,length.out=100)
x1<-c(0:E)
power_factor <- 0.75  # Choose a value between 0 (flat) and 1 (original-like)
a$Acc_flat <- rank(a$Acc) / length(a$Acc)  # Scale to [0, 1]
a$Acc_flat <- a$Acc_flat^power_factor      # Apply power adjustment
a$Acc_flat <- a$Acc_flat * 30  
b <- data.frame(Acc = pmin(a$Acc,a$Acc_flat))
b["Hits+FalseAlarms"]<-a$`Hits+FalseAlarms`
b["F"]<-b$`Hits+FalseAlarms`-b$Acc
b["misses"]<-max(b$Acc)-b$Acc
b["CorrectR"]<-100-(b$Acc+b$misses+b$F)
b["FOM"]<-b$Acc/(b$Acc+b$misses+b$F)
  

df<-data.frame(Threshold=x1/100,
               Curve=Curve,
               Upper=sapply(x1,Upper),
               Lower=sapply(x1,Lower),
               Random=sapply(x1,Random),
               Flat=b$FOM
               )
ggplot(df, aes(x = Threshold))+
  geom_line(aes(y = Curve, color="Classifier", linetype="Classifier"), size=1.25)+
  geom_line(aes(y = Upper, color="Upper Bound", linetype="Upper Bound"), size=1.25)+
  geom_line(aes(y = Lower, color="Lower Bound", linetype="Lower Bound"), size=1.25)+
  geom_line(aes(y = Random, color="Random Classifier", linetype="Random Classifier"), size=1)+
  geom_line(aes(y = Flat, color="Naive Model", linetype="Naive Model"), size=1.25)+
  
  # geom_segment(aes(x = 0.5, y = Upper[51], xend = 0.5, yend = Random[51]), linetype = "twodash", color = "#E97132", size=1.25) +  # AB
  geom_segment(aes(x = 0.5, y = Upper[51], xend = 0.5, yend = Flat[51]), linetype = "twodash", color = "#E97132", size=1.25) +  # AB
  geom_segment(aes(x = 0.5, y = Curve[51], xend = 0.5, yend = 1), linetype = "twodash", color = "#00B050", size=1.25) +  # CA
  # geom_segment(aes(x = 0.3, y = Upper[31], xend = 0.3, yend = Random[31]), linetype = "twodash", color = "#00B0F0", size=1.25) +  # ED
  geom_segment(aes(x = 0.3, y = Upper[31], xend = 0.3, yend = Flat[31]), linetype = "twodash", color = "#00B0F0", size=1.25) +  # ED
  
  geom_segment(aes(x = 0.5, y = Random[51], xend = 0.5, yend = 0), linetype = "dashed", color = "grey", size = 1) +  # CA
  geom_segment(aes(x = 0.3, y = Random[31], xend = 0.3, yend = 0), linetype = "dashed", color = "grey", size = 1) +  # ED
  
  # Points with labels
  geom_point(aes(x=0.5, y=Curve[51]), size=6) +
  geom_text(aes(x=0.5, y=Curve[51], label="A"), vjust=-0.3, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.5, y=Random[51]), size=6) +
  geom_text(aes(x=0.5, y=Random[51], label="B"), vjust=1.2, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.5, y=Upper[51]), size=6) +
  geom_text(aes(x=0.5, y=Upper[51], label="C"), vjust=-1.2, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.3, y=Random[31]), size=6) +
  geom_text(aes(x=0.3, y=Random[31], label="D"), vjust=-0.8, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.3, y=Upper[31]), size=6) +
  geom_text(aes(x=0.3, y=Upper[31], label="E"), vjust=1.5, hjust=-1.5, size=10) + 
  
  geom_point(aes(x=0.5, y=1), size=6) +
  geom_text(aes(x=0.5, y=1, label="F"), vjust=1, hjust=-1.5, size=10) + 
  
  # geom_point(aes(x=0.5, y=15/55), size=6, color="red") +
  # geom_text(aes(x=0.5, y=15/55, label="G"), vjust=0.3, hjust=2, size=10) + 
  
  geom_point(aes(x=0.5, y=Flat[51]), size=6) +
  geom_text(aes(x=0.5, y=Flat[51], label="H"), vjust=-1, hjust=-1.2, size=10) + 
  
  geom_point(aes(x=0.3, y=Flat[31]), size=6) +
  geom_text(aes(x=0.3, y=Flat[31], label="I"), vjust=-1, hjust=-1.2, size=10) + 
    
  labs(x="Hits + False Alarms (% of Extent)", y="CSI") +
  
  scale_color_manual(name = "", values = c("Classifier" = "#298d8d", "Upper Bound" = "grey20","Lower Bound" = "grey65", "Random Classifier" = "#950000", "Naive Model"="purple3"),limits = c("Classifier", "Naive Model","Random Classifier", "Upper Bound","Lower Bound")) +
  scale_linetype_manual(name = "", values = c("Classifier" = "solid", "Upper Bound" = "solid","Lower Bound" = "solid", "Random Classifier" = "dotted", "Naive Model"="twodash"),limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound","Lower Bound"))+
  

  scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,by=0.1), expand = c(0, 0), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1.01), breaks=seq(0,1,by=0.1), expand = c(0, 0), labels = percent_format()) +
  
  theme(
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
    legend.position = "right"
  ) +
  coord_fixed(ratio=1) +
  guides(fill = guide_legend(reverse = TRUE))

a$F<-a$`Hits+FalseAlarms`-a$Acc
a$CorrectR<-E-(a$misses+a$Acc+a$F)
a$E<-a$CorrectR+a$misses+a$Acc+a$F

OA<-(a$Acc+a$CorrectR)/E

Upper2<-function(x){min(x+Q,E+P-x)/E}
Lower2<-function(x){max(Q-x,x-Q)/E}
Random2<-function(x){(Pr*x+(1-Pr)*(E-x))/E}


df<-data.frame(Threshold=x1/100,
               OA=OA,
               Upper=sapply(x1,Upper2),
               Lower=sapply(x1,Lower2),
               Random=sapply(x1,Random2),
               Flat=(b$Acc+b$CorrectR)/E
)

ggplot(df, aes(x = Threshold))+
  geom_line(aes(y = OA, color="Classifier", linetype="Classifier"), size=1.25)+
  geom_line(aes(y = Upper, color="Upper Bound", linetype="Upper Bound"), size=1.25)+
  geom_line(aes(y = Lower, color="Lower Bound", linetype="Lower Bound"), size=1.25)+
  geom_line(aes(y = Random, color="Random Classifier", linetype="Random Classifier"), size=1.25)+
  geom_line(aes(y = Flat, color="Naive Model", linetype="Naive Model"), size=1.25)+
  
  # geom_segment(aes(x = 0.5, y = Upper[51], xend = 0.5, yend = Random[51]), linetype = "twodash", color = "#E97132", size=1.25) +  # AB
  geom_segment(aes(x = 0.5, y = Upper[51], xend = 0.5, yend = Flat[51]), linetype = "twodash", color = "#E97132", size=1.25) +  # AB
  geom_segment(aes(x = 0.5, y = OA[51], xend = 0.5, yend = 1), linetype = "twodash", color = "#00B050", size=1.25) +  # CA
  # geom_segment(aes(x = 0.3, y = Upper[31], xend = 0.3, yend = Random[31]), linetype = "twodash", color = "#00B0F0", size=1.25) +  # ED
  geom_segment(aes(x = 0.3, y = Upper[31], xend = 0.3, yend = Flat[31]), linetype = "twodash", color = "#00B0F0", size=1.25) +  # ED
  
  geom_segment(aes(x = 0.5, y = Random[51], xend = 0.5, yend = 0), linetype = "dashed", color = "grey", size = 1) +  # CA
  geom_segment(aes(x = 0.3, y = Random[31], xend = 0.3, yend = 0), linetype = "dashed", color = "grey", size = 1) +  # ED
  
  # Points with labels
  geom_point(aes(x=0.5, y=OA[51]), size=6) +
  geom_text(aes(x=0.5, y=OA[51], label="A"), vjust=1.5, hjust=1.5, size=10) + 
  
  geom_point(aes(x=0.5, y=Random[51]), size=6) +
  geom_text(aes(x=0.5, y=Random[51], label="B"), vjust=1.2, hjust=1.5, size=10) + 
  
  geom_point(aes(x=0.5, y=Upper[51]), size=6) +
  geom_text(aes(x=0.5, y=Upper[51], label="C"), vjust=-0.3, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.3, y=Random[31]), size=6) +
  geom_text(aes(x=0.3, y=Random[31], label="D"), vjust=1.2, hjust=1.5, size=10) + 
  
  geom_point(aes(x=0.3, y=Upper[31]), size=6) +
  geom_text(aes(x=0.3, y=Upper[31], label="E"), vjust=1, hjust=-1.5, size=10) + 
  
  geom_point(aes(x=0.5, y=Upper[31]), size=6) +
  geom_text(aes(x=0.5, y=1, label="F"), vjust=1, hjust=-1.5, size=10) + 
  
  # geom_point(aes(x=0.5, y=15/55), size=6, color="red") +
  # geom_text(aes(x=0.5, y=15/55, label="F"), vjust=0.3, hjust=2, size=10) + 
  # 
  geom_point(aes(x=0.5, y=Flat[51]), size=6) +
  geom_text(aes(x=0.5, y=Flat[51], label="H"),  vjust=-0.3, hjust=-0.8, size=10) + 
  
  geom_point(aes(x=0.3, y=Flat[31]), size=6) +
  geom_text(aes(x=0.3, y=Flat[31], label="I"), vjust=-0.4, hjust=-1.2, size=10) + 
  labs(x="Hits + False Alarms (% of Extent)", y="Overall Agreement") +
  
  scale_color_manual(name = "", values = c("Classifier" = "#298d8d", "Upper Bound" = "grey20","Lower Bound" = "grey65", "Random Classifier" = "#950000", "Naive Model"="purple3"),limits = c("Classifier", "Naive Model","Random Classifier", "Upper Bound","Lower Bound")) +
  scale_linetype_manual(name = "", values = c("Classifier" = "solid", "Upper Bound" = "solid","Lower Bound" = "solid", "Random Classifier" = "dotted", "Naive Model"="twodash"),limits = c("Classifier", "Naive Model", "Random Classifier", "Upper Bound","Lower Bound"))+
  
  
  scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,by=0.1), expand = c(0, 0), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1.01), breaks=seq(0,1,by=0.1), expand = c(0, 0), labels = percent_format()) +
  
  theme(
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
    legend.position = "right"
  ) +
  coord_fixed(ratio=1) +
  guides(fill = guide_legend(reverse = TRUE))

