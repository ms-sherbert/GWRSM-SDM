#=== Plot for Fig 5 - summary of restorable refugia ===#

library(tidyverse)

refugia <- read.csv("Model-outputs/Area_suitable-2025-01-27.csv", header=TRUE)

refugia["Scenario"][refugia["Scenario"] == "All"] <- "All suitably waterlogged soils"
refugia["Scenario"][refugia["Scenario"] == "IR1"] <- "Low infection risk (Scenario 1)"
refugia["Scenario"][refugia["Scenario"] == "IR2"] <- "Low infection risk (Scenario 2)"
refugia["Scenario"][refugia["Scenario"] == "Accessible_IR1"] <- "Accessible & low infection risk (Scenario 1)"
refugia["Scenario"][refugia["Scenario"] == "Accessible_IR2"] <- "Accessible & low infection risk (Scenario 2)"

refugia["Intensity"][refugia["Intensity"] == "High"] <- "Most suitable"
refugia["Intensity"][refugia["Intensity"] == "Low"] <- "Least suitable"

#Re-order the factor levels for graphical display
refugia$Scenario = factor(refugia$Scenario,levels=c("Accessible & low infection risk (Scenario 2)",
                                                    "Accessible & low infection risk (Scenario 1)",
                                                    "Low infection risk (Scenario 2)",    
                                                    "Low infection risk (Scenario 1)",
                                                    "Accessible",
                                                    "All suitably waterlogged soils")) 

BPl <- ggplot(data = subset(refugia, Scenario != "Total"), 
              aes(x=Scenario,y=Area_km2,fill = Intensity)) +
        geom_bar(stat = "identity",position = "stack",color="black",width = 0.9,show.legend=FALSE) +
        facet_wrap(~factor(Intensity)) +      
        scale_fill_manual(values = c("#440154","#2a788e")) +
        labs(y = expression(paste("Area (",km^2,")")),
            x = "",
            fill = "Predicted relative abundance") +
        lims(y=c(0,700)) +
        theme_bw() +
        coord_flip() +
        geom_text(aes(x=Scenario,y=Area_km2+20,label=round(Area_km2,0)), 
            stat="identity",
            position="stack",hjust="left",size = 3,color="black")
       #theme(legend.position = c(0.7, 0.2))
            #axis.text.x = element_text(size = 11),
            #axis.text.y = element_text(size = 11),
            #legend.text = element_text(size = 11),
            #legend.title = element_text(size = 12))


ggsave(BPl,
       filename=paste0("Model-outputs/Figure5.png"),
       device="png",
       height=9, width = 16, units = "cm",dpi="print")