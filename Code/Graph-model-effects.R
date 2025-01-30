#=== Covariate effects graph for Figure 2 ===#

library(tidyverse)

#Data read in

Coefs<-read.csv("D:/Repositories/GWRSM-SDM/Model-outputs/Output-summaries/Model_coefficients.csv",header=TRUE)

coef.plot <- ggplot(data=subset(Coefs, Covariate != "PO intercept" & Covariate != "PA intercept" 
                   & Covariate != "Intercept"),
                aes(x = Covariate, y = Effect.size, group_by=Model,color=Model), 
                position = position_dodge(0.5)) +
                geom_pointrange(stat="identity", aes(ymin = Lower.95.CrI, ymax = Upper.95.CrI), 
                  position = position_dodge(0.5)) +
                geom_hline(yintercept = 0) +
                #ggtitle(label="A") +
                ylab("Effect size") +
                coord_flip() +
                theme_classic() +
                theme(legend.position="bottom",legend.title=element_blank())

ggsave(coef.plot,
       filename=paste0("Model-outputs/Fig2a.png"),
       device="png",
       height=8, width = 16, units = "cm",dpi="print")


  