#=== Plot for Fig 5 - summary of restorable refugia ===#

library(tidyverse)

refugia <- read.csv("Model-outputs/Output-summaries/Restorable_areas_100m.csv", header=TRUE)

refugia["Scenario"][refugia["Scenario"] == "Whole region"] <- "All suitably waterlogged soils"
refugia["Scenario"][refugia["Scenario"] == "Accesible and infection risk < 0.7"] <- "Accessible & infection risk < 0.7"
refugia["Scenario"][refugia["Scenario"] == "Accesible"] <- "Accessible"

#Re-order the factor levels for graphical display
refugia$Scenario = factor(refugia$Scenario,levels=c("Accessible & infection risk < 0.7",
                                                "Infection risk < 0.7",
                                                "Accessible",
                                                "All suitably waterlogged soils")) 

BPl <- ggplot(data = subset(refugia, Scenario != "Total"), 
              aes(x=Scenario,y=Area_km,fill = Suitability)) +
        geom_bar(stat = "identity",position = "stack",color="black",width = 0.9,show.legend=FALSE) +
        facet_wrap(~factor(Suitability)) +      
        scale_fill_manual(values = c("#2a788e","#440154")) +
        labs(y = expression(paste("Area (",km^2,")")),
            x = "",
            fill = "Predicted relative abundance") +
        theme_bw() +
        coord_flip() +
        geom_text(aes(x=Scenario,y=Area_km-20,label=Area_km), 
            stat="identity",
            position="stack",hjust="right",size = 3.0,color="white") 
       # theme(legend.position = c(0.7, 0.2))
            #axis.text.x = element_text(size = 11),
            #axis.text.y = element_text(size = 11),
            #legend.text = element_text(size = 11),
            #legend.title = element_text(size = 12))


ggsave(BPl,
       filename=paste0("Model-outputs/Figure5_100a.png"),
       device="png",
       height=9, width = 16, units = "cm",dpi="print")