
setwd("~/project/mcIdentify/data/model_built/pathway_score_213/predict_result/")
remove(list = ls())
library(umap)
library(ggprism)
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")

data <- fread("./GSE309_predict.csv", data.table=F)
colnames(data)[2] <- "predict"
data1 <- data %>% mutate(true = case_when(type == 0 ~ "malignant",type == 1 ~ "normal")) %>%
  mutate(predict = case_when(predict == 0 ~ "malignant",predict == 1 ~ "normal")) %>% select(true,predict,GSE673_diff_path213$hsa)

umap1 <- umap::umap(data1[,3:215])
umap2 <- umap1$layout
df1<-data.frame(umap2,data1$true)
df1$data1.true<-as.factor(df1$data1.true)

p1<-ggplot(data = df1,aes(x=X1,y=X2,color=data1.true))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 True")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p1


df2<-data.frame(umap2,data1$predict)
df2$data1.predict<-as.factor(df2$data1.predict)

p2<-ggplot(data = df2,aes(x=X1,y=X2,color=data1.predict))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 Predict")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p2


library(cowplot)
prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  labels = c(),
  hjust = -1,
  nrow = 1
)
prow
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 5, 0, 5),legend.text = element_text(size = 15,family = "sans"))
)

plot_grid(prow, legend, rel_widths = c(4, .7))







setwd("~/project/mcIdentify/data/model_built/pathway_score_213/predict_result/")
remove(list = ls())
library(umap)
library(ggprism)
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")

data <- fread("./GOSH_predict.csv", data.table=F)
colnames(data)[2] <- "predict"
data1 <- data %>% mutate(true = case_when(type == 0 ~ "malignant",type == 1 ~ "normal")) %>%
  mutate(predict = case_when(predict == 0 ~ "malignant",predict == 1 ~ "normal")) %>% select(true,predict,GSE673_diff_path213$hsa)

umap1 <- umap::umap(data1[,3:215])
umap2 <- umap1$layout

predict_data <- fread("~/project/mcIdentify/data/model_built/pathway_score_213/framwork_umap/XGBoost_309_predict.csv")
predict_data <- predict_data[-1,]
colnames(predict_data) <- "predict"

df2<-data.frame(umap2,predict_data$predict)

df2$predict_data.predict<-as.factor(df2$predict_data.predict)


library(ggplot2)
p3<-ggplot(data = df2,aes(x=X1,y=X2,color=predict_data.predict))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 LR Predict")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p3


p4<-ggplot(data = df2,aes(x=X1,y=X2,color=predict_data.predict))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 RF Predict")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p4


p5<-ggplot(data = df2,aes(x=X1,y=X2,color=predict_data.predict))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 SVM Predict")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p5


p6<-ggplot(data = df2,aes(x=X1,y=X2,color=predict_data.predict))+
  geom_point(size = 0.5)+labs(x="UMAP1",y="UMAP2",color="")+
  guides(fill="none")+theme_classic()+scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  scale_colour_manual(values = c("#F8766D","#00BFC4"))+theme_prism(border = T)+ggtitle("GSE309 XGBoost Predict")+
  theme(axis.text = element_blank(),axis.ticks=element_blank())
p6



library(cowplot)
prow <- plot_grid(
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  p5 + theme(legend.position="none"),
  p6 + theme(legend.position="none"),
  align = 'vh',
  labels = c(),
  hjust = -1,
  nrow = 2
)
prow
legend <- get_legend(
  p3 + theme(legend.box.margin = margin(0, 5, 0, 5),legend.text = element_text(size = 15,family = "sans"))
)

plot_grid(prow, legend, rel_widths = c(4, .7))

