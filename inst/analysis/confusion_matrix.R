

setwd("~/project/mcIdentify/data/")
remove(list = ls())

library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/model_result/confusion/confusion_GOSH.csv")
data2 <- as.data.frame(data1[,-1])
rownames(data2) <- data1$V1



data3 <- round(data2 / rowSums(data2),2)
data3$real <- rownames(data3)
a <- melt(data3)
a$real <- factor(a$real, levels = c("normal","malignant"))
a$variable <- factor(a$variable,levels = c("malignant","normal"))

library(ggplot2)
ggplot(a, aes(real,variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(value))) +
  scale_fill_gradient(low = "#F0F0F0", high = "#3575b5") +
  labs(x = "True", y = "Guess", title = "GOSH") +
  theme_prism(border = T)+
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")



