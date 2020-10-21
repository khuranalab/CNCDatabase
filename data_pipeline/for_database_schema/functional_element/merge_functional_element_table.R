library(xlsx)
library(plyr)
library(data.table)


filePath<-"~/work/Ekta_lab/cncDatabase_dev/db_data/Jan_16_2019/functional_element_with_snv_by_cancer"
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
  str<-colnames(dat)
  cat(sprintf("%s\t",str))
  cat(sprintf("\n"))
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
#datWhole$id<-seq(1,nrow(datWhole))

######

datWhole$index<-paste(datWhole$gene_symbol,datWhole$functional_element_type,datWhole$study_id,sep="#")

unique_functional_element<-unique(datWhole$index)
functional_element_id<-seq(1,length(unique_functional_element),1)

functional_element_tmp<-data.frame(functional_element_id,unique_functional_element,stringsAsFactors = FALSE)
colnames(functional_element_tmp)<-c("functional_element_id","index")

datWhole_modified<-merge(datWhole,functional_element_tmp,by="index")

######
# make functional_element table and functional_element_gene_association table
######
tmp<-strsplit(functional_element_tmp$index,"#")
tmp2<-data.frame(do.call(rbind,tmp),stringsAsFactors = FALSE)
colnames(tmp2)<-c("gene_symbol","functional_element_type","study_id")

functional_element_tmp<-cbind(functional_element_tmp,tmp2)

functional_element_table<-functional_element_tmp[,c("functional_element_id","study_id")]
colnames(functional_element_table)<-c("id","study_id")

######
# make functional_element_gene_association table
######
functional_element_gene_association_table<-functional_element_tmp[,c("functional_element_id","functional_element_type","gene_symbol")]
colnames(functional_element_gene_association_table)<-c("id","functional_element_type","gene_symbol")


filePath<-"~/work/Ekta_lab/cncDatabase_dev/db_data/Jan_16_2019/gene_summary"
fileName<-"hgnc_complete_set.txt"
fileName<-file.path(filePath,fileName)

hgnc_table<-fread(fileName,data.table=FALSE,stringsAsFactors = FALSE)
hgnc_table$gene_summary_id<-seq(1,nrow(hgnc_table),1)
colnames(hgnc_table)[2]<-c("gene_symbol")
hgnc_table_concise<-hgnc_table[,c("gene_summary_id","gene_symbol")]

cc<-merge(functional_element_gene_association_table,hgnc_table_concise,by="gene_symbol",all.x=TRUE)


######

filePath<-"~/work/Ekta_lab/cncDatabase_dev/db_data/Jan_16_2019/functional_element"

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"functional_element_table_merged.csv"
fileName<-file.path(filePath,fileName)

write.table(functional_element_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)

#######

filePath<-"~/work/Ekta_lab/cncDatabase_dev/db_data/Jan_16_2019/functional_element_gene_association"

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"functional_element__gene_association_table_merged.csv"
fileName<-file.path(filePath,fileName)

write.table(functional_element_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)



