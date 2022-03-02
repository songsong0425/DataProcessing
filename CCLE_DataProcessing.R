#Change directory
dir <-"Z:/users/songyeon/DrugTargetIdentification/DepMap(21Q4)/CCLE" ###경로바꾸기
if(!is.dir(dir))
  dir.create(dir, recursive=TRUE)
setwd(dir)

######### 파일 불러오기 #########
ex<-read.csv("CCLE_expression.csv")
cn<-read.csv("CCLE_gene_cn.csv")
mut<-read.csv("CCLE_mutations.csv")
mut1<-read.csv("CCLE_mutations_PrimarySite.csv")
info<-read.csv("sample_info_PrimarySite.csv")
crispr<-read.csv("CRISPR_gene_effect.csv")
dependency<-read.csv("CRISPR_gene_dependency.csv")

######### 파일 전처리 #########
#Set key for merge()
colnames(ex)[1]<-"DepMap_ID"
colnames(cn)[1]<-"DepMap_ID"

#Gene symbol 추출
names(ex)<-gsub(x=names(ex), pattern="\\..\\d+\\.", replacement="") 
names(cn)<-gsub(x=names(cn), pattern="\\..\\d+\\.", replacement="") 
names(crispr)<-gsub(x=names(crispr), pattern="\\..\\d+\\.", replacement="") 
names(dependency)<-gsub(x=names(dependency), pattern="\\..\\d+\\.", replacement="") 

#Gene ID 추출
names(ex) <- sub('.*\\.(\\d+)\\.$', '\\1', names(ex))
names(cn_protein_gene) <- sub('.*\\.(\\d+)\\.$', '\\1', names(cn_protein_gene))

ex<-merge(info, ex, key="DepMap_ID")
cn<-merge(info, cn, key="DepMap_ID")
cn_protein_gene<-merge(info, cn_protein_gene, key="DepMap_ID")
mut<-merge(info, mut, key="DepMap_ID")
crispr<-merge(info, crispr, key="DepMap_ID")
dependency<-merge(info, dependency, key="DepMap_ID")

######### 파일 저장 #########
write.csv(ex, "CCLE_expression_ID_PrimarySite.csv", row.names=FALSE)
write.csv(cn, "CCLE_gene_cn_Primary_20220126.csv", row.names=FALSE)
write.csv(cn_protein_gene, "CCLE_gene_cn_ID_PrimarySite_ProteinCodingGenes.csv", row.names=FALSE)
write.csv(mut, "CCLE_mutations_PrimarySite.csv", row.names=FALSE)
write.csv(crispr, "CRISPR_gene_effect_Primary_20220126.csv", row.names=FALSE)
write.csv(dependency, "CRISPR_gene_dependency_Primary_20220126.csv", row.names=FALSE)

######### Protein-coding gene만 추출하기 #########
library(tidyverse)
cn_pcg<-cn
ex_pcg<-ex
pcg_list<-list(cn_pcg, ex_pcg)
cols <- reduce(pcg_list, .f = ~ intersect(colnames(.x), colnames(.y)))
pcg<-map_dfr(pcg_list, ~.[cols])
cn_protein_gene<-pcg[1:1750,]
write.csv(cn_protein_gene, "CCLE_gene_cn_Primary_ProteinCodingGenes.csv", row.names=FALSE)

######### Mutation dataset 정리하기 #########
library(reshape2)

length(unique(mut$Entrez_Gene_Id))
mut_without_entrez0<-mut
mut_without_entrez0<-mut_without_entrez0[mut_without_entrez0$Entrez_Gene_Id != 0,]
length(unique(mut_without_entrez0$Entrez_Gene_Id))
length(unique(mut_without_entrez0$Hugo_Symbol))

mut_without_entrez0<-mut_without_entrez0[mut_without_entrez0$Variant_Classification != 'Silent',]

sum(is.na(mut_without_entrez0$Hugo_Symbol))

mut_subset<-mut_without_entrez0[,c(1,2,4)]
mut_subset_count_2<-dcast(mut_subset, DepMap_ID~Entrez_Gene_Id) #핵심
mut_subset_count_2<-merge(info, mut_subset_count_2, key="DepMap_ID")
write.csv(mut_subset_count_2, "CCLE_mutations_ID_PrimarySite_Counting.csv", row.names=FALSE)

ex_primary<-read.csv("CCLE_expression_Symbol_PrimarySite.csv")
crispr<-read.csv("Z:/users/songyeon/DrugTargetIdentification/DepMap(21Q4)/CRISPR/CRISPR_gene_dependency_Primary_20220126.csv")
bb<-as.data.frame(table(ex_primary$primary_disease))
table(cn_protein_gene$primary_disease)
table(mut_subset_count_1$primary_disease)
table(crispr$primary_disease)
