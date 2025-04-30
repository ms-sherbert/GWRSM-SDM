### Process for downloading myrtle rust risk data from NIWA website 
### (https://myrtlerust.com/disease/plant/niwa-maps)

# Code originally written by Stephanie Tomscha, June 2022
# Updated by Sarah Herbert in January 2025 for R version 4.3.1

#install packages as needed
install.packages("httr")
install.packages("XML")

# importing packages
library(httr)
library(XML)
library(terra)

###Getting the file names, which are organized by date###
a<-seq(as.Date("2020-08-01"), as.Date("2025-01-11"), by="weeks") #New download on 16/01/2025
a<-gsub("-", "", a)
a <- paste("wk_myrinf_", a, sep="")
names_file<-paste(a,".asc", sep="")
names_file2<-paste("https://niwa-myrtle-rust.s3.ap-southeast-2.amazonaws.com/unprocessed/", names_file, sep="")


destfile <- "Z:/BioProtectAotearoa/New_meanIR_download/" #Or wherever you want these files downloaded
destination<-paste(destfile, names_file, sep="")


for(i in seq_along(names_file2)){
  download.file(names_file2[i], destination[i], mode="wb")
}