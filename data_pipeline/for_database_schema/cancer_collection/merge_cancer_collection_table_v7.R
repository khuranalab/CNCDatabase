library(xlsx)
library(plyr)

date<-"Sep_23_2020"

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_collection",sep="")
fileName<-"list.txt"
fileName<-file.path(filePath,fileName)

fileList<-read.table(fileName,stringsAsFactors = FALSE)
fileList<-fileList$V1

datWhole<-{}

for(idx in 1:length(fileList)){

  file<-fileList[idx]
  fileName<-file.path(filePath,file)  
  #dat<-read.xlsx(fileName,sheetName = "Sheet3" )
  cat(sprintf("Reading %s\n",fileName))
  
  dat<-read.xlsx(fileName,sheetIndex = 1 ,stringsAsFactors=FALSE)
  #str(dat)
  #str<-colnames(dat)
  #cat(sprintf("%s\t",str))
  #cat(sprintf("\n"))
  #pmidList<-unique(dat$pmid)
  #if(length(pmidList)!=1){cat(sprintf("More than one number in the pmid field\n"))}
 
  # check colenames
  #cat(sprintf("%s\n",colnames(dat)))
  
  #evidenceTypeList<-unique(dat$evidencetype)
  #cat(sprintf("There are %s type(s) of evidence(s) in this publication.\n",length(evidenceTypeList)))
  #cat(sprintf("%s\n",evidenceTypeList))
  
  #dat2<-split(dat,dat$element)
  
  datWhole[[idx]]<-dat
  
  
  
}


datWhole<-rbind.fill(datWhole)
datWhole$id<-seq(1,nrow(datWhole))

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_collection",sep="")
fileName<-"cancer_collection_curated.csv"
#fileName<-"cancer_collection_merged.csv"
fileName<-file.path(filePath,fileName)

write.table(datWhole,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)


