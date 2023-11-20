

remove(list = ls())
setwd("~/project/mcIdentify/data/")

library(data.table)
library(dplyr)

##method
model_result <- fread("./model_built/pathway_score_213/model_result/method.csv")

data1 <- melt(model_result)

f1 <- data1 %>% filter(variable == "f1")
f1$method <- factor(f1$method,levels = c("mcIdentify","ikarus","SCINA","scMRMA"))

library(ggplot2)
library(ggpubr)
library(ggprism)

p1 <- ggplot(f1,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+ 
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="F1 score", x  = "",title = "F1 score of different methods")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p1  




accuracy <- data1 %>% filter(variable == "accuracy")
accuracy$method <- factor(accuracy$method,levels = c("mcIdentify","ikarus","SCINA","scMRMA"))
p2 <- ggplot(accuracy,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Accuracy", x  = "",title = "Accuracy of different methods")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p2




recall <- data1 %>% filter(variable == "recall")
recall$method <- factor(recall$method,levels = c("mcIdentify","ikarus","SCINA","scMRMA"))
p3 <- ggplot(recall,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+ 
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Recall", x  = "",title = "Recall of different methods")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p3



precisoon <- data1 %>% filter(variable == "precisoon")
precisoon$method <- factor(precisoon$method,levels = c("mcIdentify","ikarus","SCINA","scMRMA"))
p4 <- ggplot(precisoon,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+ 
  scale_color_manual(values=c("black"))+
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Precison", x  = "",title = "Precison of different methods")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p4


library(cowplot)


prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C","D"),
  hjust = -1,
  nrow = 2
)
prow
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(prow, legend, rel_widths = c(3, .4))




##model framwork
model_result <- fread("./model_built/pathway_score_213/model_result/model.csv")

data1 <- melt(model_result)
colnames(data1)[2] <- "method"

f1 <- data1 %>% filter(variable == "f1")
f1$method <- factor(f1$method,levels = c("DNN","FR","LR","SVM","XGBOOST"))

library(ggplot2)
library(ggpubr)
library(ggprism)

p1 <- ggplot(f1,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+ 
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="F1 score", x  = "",title = "F1 score of different model framowrks")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p1  




accuracy <- data1 %>% filter(variable == "accuracy")
accuracy$method <- factor(accuracy$method,levels = c("DNN","FR","LR","SVM","XGBOOST"))
p2 <- ggplot(accuracy,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+ 
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Accuracy", x  = "",title = "Accuracy of different model framowrks")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p2




recall <- data1 %>% filter(variable == "recall")
recall$method <- factor(recall$method,levels = c("DNN","FR","LR","SVM","XGBOOST"))
p3 <- ggplot(recall,aes(x=method,y=value))+ 
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+ 
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Recall", x  = "",title = "Recall of different model framowrks")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p3



precisoon <- data1 %>% filter(variable == "precisoon")
precisoon$method <- factor(precisoon$method,levels = c("DNN","FR","LR","SVM","XGBOOST"))
p4 <- ggplot(precisoon,aes(x=method,y=value))+
  stat_boxplot(geom = "errorbar",width=0.15)+ 
  geom_boxplot(size=0.5,fill="#E8E8E8",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=datasets),width =0.05,shape = 21,size=3)+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442","red","blue"))+  
  scale_color_manual(values=c("black"))+
  ylim(0,1.01)+
  theme_prism(border = T)+
  labs(y="Precison", x  = "",title = "Precison of different model framowrks")+
  theme(legend.text = element_text(size = 13,family = "sans"))
p4


library(cowplot)


prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C","D"),
  hjust = -1,
  nrow = 2
)
prow
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(prow, legend, rel_widths = c(3, .4))

