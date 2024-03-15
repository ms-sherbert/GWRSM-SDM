#=== Plot for Fig 5 - summary of restorable refugia ===#
# Written by: S Herbert
# Written for R version 4.3.1

library(tidyverse)

refugia <- read.csv("model-outputs/Output-summaries/Restorable_areas.csv", header=TRUE)

refugia["Scenario"][refugia["Scenario"] == "Whole region"] <- "All suitably waterlogged soils"
refugia["Scenario"][refugia["Scenario"] == "Accesible and infection risk < 0.7"] <- "Accessible & infection risk < 0.7"
refugia["Scenario"][refugia["Scenario"] == "Accesible"] <- "Accessible"

#Re-order the factor levels for graphical display
refugia$Scenario = factor(refugia$Scenario,levels=c("Accessible & infection risk < 0.7",
                                                "Infection risk < 0.7",
                                                "Accessible",
                                                "All suitably waterlogged soils")) 

ggplot(data = subset(refugia, Scenario != "Total"), aes(x=Scenario,y=Area_km,fill = Suitability)) +
  geom_bar(stat = "identity",position = "stack",color="black",width = 0.9) +
  scale_fill_manual(values = c("#50f744","white")) +
  labs(y = "Area (sq. km)",
       x = "",
       fill = "Predicted relative abundance") +
  theme_bw() +
  coord_flip() +
  geom_text(aes(x=Scenario,y=c(2839,2469,736,847,1006,1080,219,380),label=Area_km), 
            stat="identity",
            position="stack",hjust="right",size = 3.5) +
  theme(legend.position = c(0.81, 0.25),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
