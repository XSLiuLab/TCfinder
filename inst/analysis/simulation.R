
## pathway_simulation
setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)
data <- fread("./model_built/pathway_score_213/model_result/mcIdentify_simulation_pathway/10_pathway.txt")
colnames(data) <- c("simulation","measure","5%")
for (number in c(20,40,60,80,100,120,140,160,180)) {
  filename <- paste0("./model_built/pathway_score_213/model_result/mcIdentify_simulation_pathway/",number,"_pathway.txt")
  data1 <- fread(filename)
  cname <- paste0(number/2,"%")
  data[,cname] <- data1$V3
  
}


data[,"0%"] <- rep(c(0.98,0.98,0.99,0.98),100)
data2 <- melt(data)
data2$variable <- factor(data2$variable,levels = c("0%","5%","10%","20%","30%","40%",
                                                   "50%","60%","70%","80%","90%"))

draw_data <- data2 %>% filter(measure == "precision:")

library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=draw_data,aes(x=variable,y=value,fill = variable))+
  geom_boxplot(size=1, draw_quantiles = c(0.5))+
  theme_prism(border = T)+theme(legend.position = 'none')+
  labs(y="Precision",title = " ")+
  theme(axis.title.x = element_blank())+
  geom_hline(aes(yintercept=0.9),linetype=5,col="red")+
  geom_hline(aes(yintercept=0.8),linetype=5,col="red")+
  scale_y_continuous(breaks=c(0,0.5,0.6,0.7,0.8,0.9,1))+
  scale_fill_manual(values = c("#F27970", "#BB9727","#54B345","#32B897",
                               "#05B9E2", "#8983BF","#C76DA2","#F27970",
                               "#BB9727","#54B345","#32B897"))+
  ylim(0.5,1)







## simulation gene 
setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)
mcIdentify <- fread("./model_built/pathway_score_213/model_result/simulation_gene/mcIdentify.txt")
SCINA <- fread("./model_built/pathway_score_213/model_result/simulation_gene/SCINA.txt")
scMRMA <- fread("./model_built/pathway_score_213/model_result/simulation_gene/scMRMA.txt")
ikarus <- fread("./model_built/pathway_score_213/model_result/simulation_gene/ikarus_pred.txt")


colnames(SCINA) <- c("ID","measure","SCINA")
colnames(mcIdentify) <- c("ID","measure","mcIdentify")
colnames(ikarus) <- c("ID","measure","ikarus")
colnames(scMRMA) <- c("ID","measure","scMRMA")


data <- left_join(mcIdentify,ikarus) %>% left_join(.,SCINA) %>% left_join(.,scMRMA) %>% arrange(ID)

data$infor <- c(rep("1000gene",40),rep("1500gene",40),rep("2000gene",40),rep("2500gene",40),rep("500gene",40))

data1 <- melt(data)

data1$variable <- factor(data1$variable,levels = c("mcIdentify","scMRMA","SCINA","ikarus"))
data1$infor <- factor(data1$infor,levels = c("500gene","1000gene","1500gene","2000gene","2500gene"))

draw_data <- data1 %>% filter(measure == "F1:") #  & variable == "mcIdentify"
draw_data <- data1 %>% filter(measure == "recall:") #  & variable == "mcIdentify"


a <- draw_data %>% group_by(variable, infor) %>% summarise(mean(value))
colnames(a) <- c("method","gene","value")
V1 <- c("500gene","1000gene","1500gene","2000gene","2500gene")
library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(a,aes(x=gene,y=value,fill=method))+
  geom_bar(position="dodge",stat="identity")+
  labs(x="Number of random genes",y="Recall")+
  theme_prism(border = F)+
  geom_hline(aes(yintercept=0.9),linetype=5,col="red")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,0.9))+
  scale_x_discrete(breaks=V1, labels=c("500","1000","1500","2000","2500"))+
  scale_fill_manual(values = c("#6E9ECE", "#CCCCCC","#E6928F","#8FBC8F"),
                    breaks=c("mcIdentify","scMRMA","SCINA","ikarus"),
                    labels=c("mcIdentify","scMRMA","SCINA","ikarus"))




## 4pathway simulation
setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

data1 <- fread("./model_built/pathway_score_213/model_result/simulation_pathway2/36_without4_pathway.txt")
data1$pathway <- "withoutpathway"

data2 <- fread("./model_built/pathway_score_213/model_result/simulation_pathway2/40_with4_pathway.txt")
data2$pathway <- "withpathway"

data3 <- rbind(data2,data1)
data3$pathway <- factor(data3$pathway,levels = c("withpathway","withoutpathway"))

draw_data <- data3 %>% filter(V2 == "F1:")

my_lists <- list(c("withpathway","withoutpathway"))

ggplot(data=draw_data,aes(x=pathway,y=V3,fill = pathway))+
  geom_boxplot(size=1, draw_quantiles = c(0.5))+
  theme_prism(border = F)+theme(legend.position = 'none')+
  labs(y="F1 score",title = " ")+
  stat_compare_means(method = "wilcox.test")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(breaks=c(0,0.5,0.6,0.7,0.8,0.9,1))+
  scale_fill_manual(values = c("#91D1C2FF", "#FDAF91FF"))+
  ylim(0.5,1)






##simulation gene 2
library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/model_result/simulation_gene.csv")
data1$infor <- rep(c("100gene","200gene","300gene","400gene","500gene"),4)

data2 <- melt(data1)


data2$method <- factor(data2$method,levels = c("TCfinder","scMRMA","SCINA","ikraus"))
data2$infor <- factor(data2$infor,levels = c("100gene","200gene","300gene","400gene","500gene"))

a <- data2 %>% dplyr::filter(variable == "f1")

V1 <- c("100gene","200gene","300gene","400gene","500gene")
library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(a,aes(x=infor,y=value,fill=method))+
  geom_bar(position="dodge",stat="identity")+
  labs(x="Number of randomly inactivate genes",y="F1 Score")+
  theme_prism(border = F)+
  geom_hline(aes(yintercept=0.95),linetype=5,col="red")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,0.95,1))+
  scale_x_discrete(breaks=V1, labels=c("100","200","300","400","500"))+
  scale_fill_manual(values = c("#6E9ECE", "#CCCCCC","#E6928F","#8FBC8F"),
                    breaks=c("TCfinder","scMRMA","SCINA","ikraus"),
                    labels=c("TCfinder","scMRMA","SCINA","ikarus"))






