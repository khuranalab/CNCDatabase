library(xlsx)
library(plyr)
library(data.table)

date<-"Sep_23_2020"

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")
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
datWhole$idx<-seq(1,nrow(datWhole))

######
# fix emsembl id
######
filePath<-"~/work/Ekta_lab/compositeDriver_data/ensembl_hg37"
fileName<-"hg37_ensembl_table_v2.txt"
fileName<-file.path(filePath,fileName)

geneID<-fread(file=fileName,sep="\t",header=TRUE,data.table=FALSE,stringsAsFactors = FALSE)
#colnames(geneID)[2]<-"gene"

######
tmpDat<-datWhole
#dat3<-merge(dat2,geneID,by.x="gene")

for(idx in 1:nrow(tmpDat)){
  
  #cat(sprintf("idx: %s\n",idx))     
  geneName<-tmpDat[idx,]$gene_symbol
  
  # ensembl_gene_id may have multiple return value, it need to fix 
  ensembl_gene_id<-unique(geneID[geneID$external_gene_name %in% geneName,]$ensembl_gene_id)[1]
  tmpDat[idx,]$ensemblid<-ensembl_gene_id
  
  if(is.na(ensembl_gene_id)){
    cat(sprintf("idx:%s,geneName:%s without ensembl id\n",idx,geneName))
  }
  
}

#datWhole<-tmpDat

#######
# fix gene name without ensembl id
fix_list<-c("AC007431.1","C19orf22","FLJ41941",
            "RFPL3-AS1","C11orf10","AP001465.5.1",
            "AE000659.41.1","SAMMSON","lncRNA-ATB",
            "RP11.1101K5.1","KIAA1737","MIR205",
            "G029190","G025135","G029190","hsa-mir-142","Ala.TGC","Ala.TGC","Met.CAT","Gly.GCC","Ala.TGC",
            "G062818","G079632")

find_ensembl_id<-c("ENSG00000166329","ENSG00000198858","ENSG00000235295",
                  "ENSG00000205853","ENSG00000134825","none",
                  "none","ENSG00000240405","ENSG00000244306",
                  "ENSG00000253434","ENSG00000198894","ENSG00000284485",
                  "none","none","none","ENSG00000284353","none","none","none","none","none",
                  "none","none")

fix_missing_ensembl_id<-data.frame(fix_list,find_ensembl_id,stringsAsFactors = FALSE)
colnames(fix_missing_ensembl_id)<-c("gene_symbol","ensemblid")

#tmpDat[tmpDat$gene_symbol %in% fix_list,]
#geneID[geneID$hgnc_symbol %in% fix_list,]
#geneID[geneID$ensembl_gene_id %in% find_ensembl_id,]

tmpDat[tmpDat$gene_symbol %in% fix_missing_ensembl_id$gene_symbol,]$ensemblid<-fix_missing_ensembl_id$ensemblid

tmpDat[tmpDat$gene_symbol %in% fix_missing_ensembl_id$gene_symbol,]

datWhole<-tmpDat

# 1594 entries in total after removing duplidated entries from study 3

#######


filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_collection",sep="")
fileName<-"cancer_collection_w_cancer_summary_v7.txt"
fileName<-file.path(filePath,fileName)

cancer_summary_merged<-read.table(fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE)
colnames(cancer_summary_merged)[2]<-"cancer_collection_id"
cancer_summary_simple<-cancer_summary_merged[,c(2,10:18)]


cc<-merge(datWhole,cancer_summary_simple,all.x=TRUE,by="cancer_collection_id")

# 18 genes with gene expression association
unique(cc[cc$evidence_type %in% "gene expression association",]$gene_symbol)

# 21 genes with functional validation
unique(cc[cc$evidence_type %in% "experimental validation",]$gene_symbol)

#######

aggreated_cancer_collections<-c("PanCancer",
                               "Pancan-no-skin-melanoma-lymph",
                               "Carcinoma","Hematopoietic system",
                               "Adenocarcinoma","Lymphatic system",
                               "Digestive tract","Female reprodutive system",
                               "Squamous","Breast",
                               "Lung","CNS",
                               "Giloma","Sarcoma",
                               "Kidney","Myeloid")

c1<-cc[!(cc$cancer_full_name %in% aggreated_cancer_collections),]
c1<-c1[c1$evidence_type %in% "computational prediction",]

single_cancer_type_prediction<-unique(c1$gene_symbol)

c2<-cc[(cc$cancer_full_name %in% aggreated_cancer_collections),]
c2<-c2[c1$evidence_type %in% "computational prediction",]

pancancer_type_prediction<-unique(c2$gene_symbol)

length(intersect(single_cancer_type_prediction,pancancer_type_prediction))


#######
if(FALSE){
  
  
filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")
fileName<-paste("functional_element_with_snv_by_cancer_merged_",date,".txt",sep="")
fileName<-file.path(filePath,fileName)
write.table(cc,fileName,sep="\t",row.names = FALSE,col.names = TRUE)

#######

unique(datWhole[datWhole$evidence_type %in% "computational prediction",]$gene_symbol)

datWhole[grepl("breast",datWhole$cancer_name_from_study),]
}

#######

#datWhole$index<-paste(datWhole$gene_symbol,datWhole$functional_element_type,datWhole$study_id,sep="#")
datWhole$index<-paste(datWhole$gene_symbol,datWhole$functional_element_type,datWhole$study_id,datWhole$ensemblid, sep="#")

unique_functional_element<-datWhole[,c("index","idx")]
unique_functional_element<-unique_functional_element[!(duplicated(unique_functional_element$index)),]
unique_functional_element$functional_element_id<-seq(1,nrow(unique_functional_element),1)

functional_element_tmp<-unique_functional_element[,c("functional_element_id","index")]
colnames(functional_element_tmp)<-c("functional_element_id","index")

datWhole_modified<-merge(datWhole,functional_element_tmp,by="index",all.x=TRUE)

datWhole_modified<-datWhole_modified[,2:ncol(datWhole_modified)]
datWhole_modified<-datWhole_modified[order(datWhole_modified$idx),]

datWhole_modified$index2<-paste(datWhole_modified$gene_symbol,
                                datWhole_modified$functional_element_type,
                                datWhole_modified$cancer_collection_id,
                                datWhole_modified$study_id,
                                datWhole_modified$evidence_type,sep="#")

unique_functional_element_with_snv_by_cancer<-datWhole_modified[,c("index2","idx")]
unique_functional_element_with_snv_by_cancer<-unique_functional_element_with_snv_by_cancer[!(duplicated(unique_functional_element_with_snv_by_cancer$index2)),]
unique_functional_element_with_snv_by_cancer$functional_element_with_snv_by_cancer_id<-seq(1,nrow(unique_functional_element_with_snv_by_cancer),1)

functional_element_with_snv_by_cancer_tmp<-unique_functional_element_with_snv_by_cancer[,c("functional_element_with_snv_by_cancer_id","index2")]
colnames(functional_element_with_snv_by_cancer_tmp)<-c("functional_element_with_snv_by_cancer_id","index2")

datWhole_modified2<-merge(datWhole_modified,functional_element_with_snv_by_cancer_tmp,by="index2",all.x=TRUE)

datWhole_modified2<-datWhole_modified2[order(datWhole_modified2$idx),]

######
#  handle cancer_driver_evidence table from evidence_summary
######

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/evidence_summary",sep="")
fileName<-"evidence_summary.xlsx"
fileName<-file.path(filePath,fileName)

evidence_summary<-read.xlsx(fileName,sheetIndex = 1 ,stringsAsFactors=FALSE)
colnames(evidence_summary)[1]<-c("evidence_summary_id")

evidence_summary$evidence_summary_idx<-paste(evidence_summary$evidence_type,evidence_summary$evidence_method,sep="#")

# to include "evidence_description"
evidence_summary_short<-evidence_summary[,c("evidence_summary_id","evidence_summary_idx","evidence_description")]

datWhole_modified2$evidence_summary_idx<-paste(datWhole_modified2$evidence_type,datWhole_modified2$evidence_method,sep="#")

cc<-merge(datWhole_modified2,evidence_summary_short,by="evidence_summary_idx",all.x=TRUE)

# to include "evidence_description"                                     
datWhole_modified2<-cc

tmpTable<-cc[,c("functional_element_with_snv_by_cancer_id","evidence_summary_id")]
tmpTable$index<-paste(tmpTable$functional_element_with_snv_by_cancer_id,tmpTable$evidence_summary_id,sep="#")
tmpTable<-tmpTable[!(duplicated(tmpTable$index)),]
tmpTable<-tmpTable[order(tmpTable$functional_element_with_snv_by_cancer_id),]
tmpTable$id<-seq(1,nrow(tmpTable),1)           

cancer_driver_evidence_table<-tmpTable[,c("id","functional_element_with_snv_by_cancer_id","evidence_summary_id")]


filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/cancer_driver_evidence",sep="")

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"cancer_driver_evidence_curated.csv"
fileName<-file.path(filePath,fileName)

write.table(cancer_driver_evidence_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)



                                                 
######
# make functional_element table 
######
tmp<-strsplit(functional_element_tmp$index,"#")
tmp2<-data.frame(do.call(rbind,tmp),stringsAsFactors = FALSE)
colnames(tmp2)<-c("gene_symbol","functional_element_type","study_id","ensemblid")

functional_element_tmp<-cbind(functional_element_tmp,tmp2)

functional_element_table<-functional_element_tmp[,c("functional_element_id","study_id")]
colnames(functional_element_table)<-c("id","study_id")


filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element",sep="")

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"functional_element_curated.csv"
fileName<-file.path(filePath,fileName)

write.table(functional_element_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)

########
# make gene_summary table
########

## test MIR250 in this file?
filePath<-"~/work/Ekta_lab/compositeDriver_data/ensembl_hg37"
fileName<-"hg37_ensembl_table_v2.txt"
fileName<-file.path(filePath,fileName)

gene_summary_table<-fread(fileName,data.table=FALSE)
id<-seq(1,nrow(gene_summary_table),1)
gene_summary_table<-cbind(id,gene_summary_table)
colnames(gene_summary_table)<-c("id","external_gene_name","ensemblid",
                                "chromosome","start_position","end_position",
                                "description","external_gene_source","entrezgene",
                                "gene_type","hgncid","hgnc_symbol")


filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/gene_summary",sep="")

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"gene_summary_curated.txt"
fileName<-file.path(filePath,fileName)

write.table(gene_summary_table,file=fileName,sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

######
# make functional_element_gene_association table
######

functional_element_gene_summary_tmp<-merge(functional_element_tmp,gene_summary_table,all.x = TRUE, by="ensemblid")
functional_element_gene_summary_tmp<-functional_element_gene_summary_tmp[order(functional_element_gene_summary_tmp$functional_element_id),]

functional_element_gene_association_table<-functional_element_gene_summary_tmp[,c("functional_element_id","functional_element_type","gene_symbol","ensemblid","id")]
colnames(functional_element_gene_association_table)<-c("functional_element_id","functional_element_type","gene_symbol","ensemblid","gene_summary_id")
functional_element_gene_association_table$id<-seq(1,nrow(functional_element_gene_association_table),1)
functional_element_gene_association_table<-functional_element_gene_association_table[,c("id","functional_element_id","functional_element_type","gene_symbol","ensemblid","gene_summary_id")]

####
#functional_element_gene_association_table<-functional_element_tmp[,c("functional_element_id","functional_element_type","gene_symbol")]
#colnames(functional_element_gene_association_table)<-c("functional_element_id","functional_element_type","gene_symbol")
#functional_element_gene_association_table$id<-seq(1,nrow(functional_element_gene_association_table),1)
#functional_element_gene_association_table<-functional_element_gene_association_table[,c("id","functional_element_id","functional_element_type","gene_symbol")]

if(FALSE){
filePath<-"~/work/Ekta_lab/cncDatabase_dev/db_data/Jan_16_2019/gene_summary"
fileName<-"hgnc_complete_set.txt"
fileName<-file.path(filePath,fileName)

hgnc_table<-fread(fileName,data.table=FALSE,stringsAsFactors = FALSE)
hgnc_table$gene_summary_id<-seq(1,nrow(hgnc_table),1)
colnames(hgnc_table)[2]<-c("gene_symbol")
hgnc_table_concise<-hgnc_table[,c("gene_summary_id","gene_symbol")]

cc<-merge(functional_element_gene_association_table,hgnc_table_concise,by="gene_symbol",all.x=TRUE)
}

#######

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_gene_association",sep="")

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"functional_element_gene_association_curated.csv"
fileName<-file.path(filePath,fileName)

write.table(functional_element_gene_association_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)

#######
# make functional_element_with_snv_by_cancer table
#######

functional_element_with_snv_by_cancer_table<-datWhole_modified2[,c("functional_element_with_snv_by_cancer_id",
                                                                   "functional_element_id",
                                                                   "cancer_collection_id",
                                                                   "num_of_mutations",
                                                                   "num_of_mutated_samples",
                                                                   "somatic",
                                                                   "germline")]

colnames(functional_element_with_snv_by_cancer_table)[1]<-c("id")
functional_element_with_snv_by_cancer_table<-functional_element_with_snv_by_cancer_table[!(duplicated(functional_element_with_snv_by_cancer_table$id)),]

                                                                                         
filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")

if( !file.exists(paste(filePath,sep="/")) ){
  dir.create(paste(filePath,sep=""),recursive=TRUE)
}


fileName<-"functional_element_with_snv_by_cancer_curated.csv"
fileName<-file.path(filePath,fileName)

write.table(functional_element_with_snv_by_cancer_table,file=fileName,sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)


#######
# for my own reference
#######
if(TRUE){
datWhole_modified2$evidence_summary_idx2<-paste(datWhole_modified2$gene_symbol,datWhole_modified2$functional_element_type,datWhole_modified2$evidence_summary_idx,sep="#")
datWhole_modified3<-datWhole_modified2[!duplicated(datWhole_modified2$evidence_summary_idx2),]

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")
fileName<-"functional_element_with_snv_by_cancer_merged_Feb_20_2020.txt"
fileName<-file.path(filePath,fileName)
write.table(datWhole_modified3,fileName,sep="\t",row.names = FALSE,col.names = TRUE)
}

########

######
#  load study table from study
######

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/study",sep="")
fileName<-"study.xlsx"
fileName<-file.path(filePath,fileName)

study<-read.xlsx(fileName,sheetIndex = 1 ,stringsAsFactors=FALSE)
colnames(study)[1]<-c("study_id")
study2<-study[,c(1,3:9)]


#ww<-datWhole_modified2[datWhole_modified2$pmid %in% 0,]

#
# merge dataWhole, cancer_summary_simple, study tables
#

# old solution
#datWhole_part1<-datWhole_modified2[,c(2:18)]

# solution to add "evidence_description"
datWhole_part1<-datWhole_modified2[,c(3:18,23,19)]

datWhole_merged_with_cancer_summary<-merge(datWhole_part1,cancer_summary_simple,by="cancer_collection_id")

datWhole_merged_with_cancer_summary_with_study<-merge(datWhole_merged_with_cancer_summary,study2,by="study_id")
#colnames(datWhole_merged_with_cancer_summary_with_study)[6]<-"cancer_class"


#######
# add element description from "functional_category" table
#######

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_category",sep="")
fileName<-"functional_category_curated.txt"
fileName<-file.path(filePath,fileName)

# tab delimited format
functional_category<-fread(fileName,sep="\t",header=TRUE,stringsAsFactors = FALSE,data.table=FALSE)
functional_category$functional_category_idx<-paste(functional_category$study_id,functional_category$functional_element_type,sep="#")
functional_category_short<-functional_category[,c(4:6)]

datWhole_merged_with_cancer_summary_with_study$functional_category_idx<-paste(
  datWhole_merged_with_cancer_summary_with_study$study_id,
  datWhole_merged_with_cancer_summary_with_study$functional_element_type,
  sep="#"
)

datWhole_merged_with_cancer_summary_with_study<-merge(datWhole_merged_with_cancer_summary_with_study,functional_category_short,by="functional_category_idx",all.x=TRUE)

#######

# old solution
#datWhole_merged_with_cancer_summary_with_study2<-datWhole_merged_with_cancer_summary_with_study[,c(17,3:16,18:30)]

# solution to add "evidence_description"
datWhole_merged_with_cancer_summary_with_study2<-datWhole_merged_with_cancer_summary_with_study[,c(19,4:11,36,37,12:18,20:35)]

datWhole_merged_with_cancer_summary_with_study2<-datWhole_merged_with_cancer_summary_with_study2[order(datWhole_merged_with_cancer_summary_with_study2$idx),]

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")
fileName<-"noncoding_cancer_driver_Sep_23_2020.txt"
fileName<-file.path(filePath,fileName)
write.table(datWhole_merged_with_cancer_summary_with_study2,fileName,sep="\t",row.names = FALSE,col.names = TRUE)



#######

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,sep="")
fileName<-"COSMIC_v80_603_cancer_genes.xlsx"
fileName<-file.path(filePath,fileName)

cancerGeneList<-read.xlsx(fileName,sheetIndex = 1 ,stringsAsFactors=FALSE)

datWhole_merged_with_cancer_summary_with_study2$cosmic_cgc<-"No"
datWhole_merged_with_cancer_summary_with_study2[datWhole_merged_with_cancer_summary_with_study2$gene_symbol %in% cancerGeneList$COSMIC,]$cosmic_cgc<-"Yes"

filePath<-paste("~/work/Ekta_lab/cncDatabase_dev/db_data/",date,"/functional_element_with_snv_by_cancer",sep="")
fileName<-"noncoding_cancer_driver_Sep_23_2020_w_COSMIC_cancer_gene_annotation.txt"
fileName<-file.path(filePath,fileName)
write.table(datWhole_merged_with_cancer_summary_with_study2,fileName,sep="\t",row.names = FALSE,col.names = TRUE)

