#=== Covariate effects graph ===#
# Written by: S Herbert
# Written for R version 4.3.1

library(tidyverse)

#Data read in

Coefs<-read.csv("GWRSM-SDM/Model-outputs/Model_coefficients.csv",header=TRUE) #Or wherever you've cloned the repository

ggplot(data=subset(Coefs, Covariate != "PO intercept" & Covariate != "PA intercept" 
                   & Covariate != "Intercept"),
       aes(x = Covariate, y = Effect.size, group_by=Model,color=Model), 
       position = position_dodge(0.3)) +
  geom_pointrange(stat="identity", aes(ymin = Lower.95.CrI, ymax = Upper.95.CrI), 
                  position = position_dodge(0.3)) +
  geom_hline(yintercept = 0) +
  ggtitle(label="A") +
  ylab("Effect size") +
  coord_flip() +
  theme_bw()



  