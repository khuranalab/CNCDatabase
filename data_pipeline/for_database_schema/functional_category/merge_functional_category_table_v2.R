library(xlsx)
library(plyr)


date<-"Sep_23_2020"

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_category",sep="")
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

datWhole$id<-seq(1,nrow(datWhole),1)


filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_category",sep="")
fileName<-"functional_category_curated.txt"
#fileName<-"functional_category_merged_01_25_2019.csv"
fileName<-file.path(filePath,fileName)

write.table(datWhole,file=fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)


