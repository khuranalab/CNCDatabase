library(xlsx)
library(plyr)
library(data.table)

date<-"Sep_23_2020"

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_collection",sep="")
fileName<-"cancer_collection_curated.csv"
#fileName<-"cancer_collection_merged.csv"
fileName<-file.path(filePath,fileName)

cancer_collection_dat<-fread(file=fileName,sep=",",data.table=FALSE)

#####

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_summary",sep="")
fileName<-"cancer_summary_v6.xlsx"
#fileName<-"cancer_collection_merged.csv"
fileName<-file.path(filePath,fileName)

cancer_summary_dat<-read.xlsx(file=fileName,sheetName = "Sheet1")
colnames(cancer_summary_dat)[1]<-"cancer_summary_id"

#######

cancer_summary_merged_dat<-merge(cancer_collection_dat,cancer_summary_dat,all.x=TRUE,by="cancer_summary_id")

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_collection",sep="")
fileName<-"cancer_collection_w_cancer_summary_v7.txt"
#fileName<-"cancer_collection_merged.csv"
fileName<-file.path(filePath,fileName)
  
write.table(cancer_summary_merged_dat,fileName,sep="\t",row.names = FALSE,col.names = TRUE)
    
