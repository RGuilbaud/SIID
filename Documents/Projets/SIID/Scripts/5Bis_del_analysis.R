################ Look at the del ###############
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(forcats)

### Files and repertories 

## ref gen genes
setwd("~/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
infos<-read_xlsx("emergen_romane_lineage.xlsx")

## muations signature of variants
mut_20I<-read.table(file="Mut_20I_Alpha_V1.txt",sep=":",header = F)
mut_20I$V2<-gsub("-","del",mut_20I$V2)
mut_21J<-read.table(file="Mut_21J_Delta.txt",sep=":",header = F)
mut_21J$V2<-gsub("-","del",mut_21J$V2)
mut_21K<-read.table(file="Mut_21K_Omicron.txt",sep=":",header = F)
mut_21K$V2<-gsub("-","del",mut_21K$V2)
mut_21L<-read.table(file="Mut_21L_Omicron.txt",sep=":",header = F)
mut_21L$V2<-gsub("-","del",mut_21L$V2)
mut_22B<-read.table(file="Mut_22B_Omicron.txt",sep=":",header = F)
mut_22B$V2<-gsub("-","del",mut_22B$V2)


### import mutations from the sequencing runs from bichat
setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_1-2/")
tab12<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_2-3/")
tab23<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_4/")
tab4<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_5/")
tab5<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_6/")
tab6<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_7/")
tab7<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)

## remove samples redone
tab12<-tab12[which((tab12$Id_sample%in%tab5$Id_sample)==F),]
tab23<-tab23[which((tab23$Id_sample%in%tab5$Id_sample)==F),]
tab4<-tab4[which((tab4$Id_sample%in%tab5$Id_sample)==F),]
tab4<-tab4[which((tab4$Id_sample%in%tab6$Id_sample)==F),]
tab6<-tab6[which((tab6$Id_sample%in%tab23$Id_sample)==F),]

### import mutations from pitié
setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/")
tabP<-read.table(file = "Mutations J0 .txt",sep = "\t",header = T)
tabP<-tabP[which(tabP$Id_sample%in%infos$sample[which(infos$origin=="pitie")]),]


### merge the tables
## only bichat
tab<-rbind(tab12,rbind(tab23,rbind(tab4,rbind(tab5, rbind(tab6,tab7)))))
## TO UNCOMMENT IF PITIE + BICHAT (UNCOMMENT LINES 153-159)
#tab<-rbind(tab12,rbind(tab23,rbind(rbind(tab4,rbind(tab5, rbind(tab6,tab7))),tabP)))


### sel del with freq >= 70
tab<-tab[which(str_detect(tab$mutation,"del")),]
tab<-tab[which(tab$freq>=70),]

length(unique(tab$mutation))

tab$aa_pos<-str_sub(tab$mutation,2, -4)
tab$gene_mutation<-paste(tab$gene,tab$mutation,sep="_")

### nb del per sample
nb_del<-tab %>% group_by(Id_sample) %>% summarise(n=n())

ggplot(data=nb_del, aes(x=as.factor(Id_sample), y=n)) +
  geom_bar(stat="identity")+
  ylab("number of del mutations")+
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())


### nb samples per del
nb_sample<- tab %>% group_by(gene_mutation) %>% summarise(n=n())
nb_sample$AA_pos<-as.numeric(tab$aa_pos[match(nb_sample$gene_mutation,tab$gene_mutation)])
nb_sample$gene<-tab$gene[match(nb_sample$gene_mutation,tab$gene_mutation)]
nb_sample$gene_start<-genes$Start[match(nb_sample$gene,genes$GeneId)]
nb_sample$gene_order<-nb_sample$gene_start * 10000
nb_sample$order<-nb_sample$gene_order + nb_sample$AA_pos

ggplot(data=nb_sample, aes(x=fct_reorder(gene_mutation, order), y=n)) +
  geom_bar(stat="identity")+
  ylab("number of samples")+
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())


### del specific of a variant 
tab$mut_spe<-NA
for (m in 1:nrow(tab)){
  if(tab$mutation[m] %in% mut_20I$V2[which(mut_20I$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"20I_Alpha_V1"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "20I_Alpha_V1", sep=";")
    }
  }
  if(tab$mutation[m] %in% mut_21J$V2[which(mut_21J$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"21J_Delta"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "21J_Delta", sep=";")
    }
  }
  if(tab$mutation[m] %in% mut_21K$V2[which(mut_21K$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"21K_Omicron"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "21K_Omicron", sep=";")
    }
  }
  if(tab$mutation[m] %in% mut_21L$V2[which(mut_21L$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"21L_Omicron"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "21L_Omicron", sep=";")
    }
  }
  if(tab$mutation[m] %in% mut_22B$V2[which(mut_22B$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"22B_Omicron"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "22B_Omicron", sep=";")
    }
  }
}



nb_sample$mut_spe<-tab$mut_spe[match(nb_sample$gene_mutation,tab$gene_mutation)]

#### plot the number of samples carrying the dels
### Bichat only
setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
png(filename = "nb_samples_variant_of_dels_seq1-2-3-4-5-6-7.png", width=15, height=8,units="in",res = 200)
ggplot(data=nb_sample, aes(x=fct_reorder(gene_mutation, order), y=n, fill=mut_spe)) +
  geom_bar(stat="identity")+
  ylab("number of samples")+
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())
dev.off()

### TO UNCOMMENT IF PITIE + BICHAT
# setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
# png(filename = "nb_samples_variant_of_dels_seq1-2-3-4-5-6-7-pitie.png", width=30, height=8,units="in",res = 200)
# ggplot(data=nb_sample, aes(x=fct_reorder(gene_mutation, order), y=n, fill=mut_spe)) +
#   geom_bar(stat="identity")+
#   ylab("number of samples")+
#   theme(axis.text.x = element_text(angle = 90, size = 4), axis.title.x = element_blank())
# dev.off()


tab$variant<-infos$clade[match(tab$Id_sample,infos$sample)]
tab$individual<-infos$individual[match(tab$Id_sample,infos$sample)]


########### Look at ORF7 and ORF8 --> PITIE + BICHAT (nothing to see in bichat)
tab_orf78<-tab[which(tab$gene=="ORF7a" | tab$gene=="ORF7b" | tab$gene=="ORF8"),]
infos_del78<-infos[which(infos$sample%in%tab_orf78$Id_sample),]
infos_del78$nb_del<-NA
for (i in 1:nrow(infos_del78)){
  infos_del78$nb_del[i]<-length(which(tab_orf78$Id_sample==infos_del78$sample[i]))
}
infos_del78p<-infos_del78[which(infos_del78$origin=="pitie"),]

## keep the samples that have at least 10bp del in OFR7 and 8
infos_del78_supp10aa<-infos_del78[which(infos_del78$nb_del>10),]
setwd("~/Projets/SIID/2_Filter_trim/Accurate")
write.table(infos_del78_supp10aa,file="samples_del_supp10aa_ORF78.txt",sep="\t",row.names = F,quote = F)
# ggplot(data=infos_del78, aes(x=nb_del)) +
#   geom_bar(stat="count")

nb_sample78<-nb_sample[which(nb_sample$gene=="ORF7a" | nb_sample$gene=="ORF7b" | nb_sample$gene=="ORF8"),]
nb_sample78<-nb_sample78[which(nb_sample78$n<50),]

setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
png(filename = "nb_samples_variant_of_dels_seq1-2-3-4-pitie_ORF7-8.png", width=15, height=6,units="in",res = 200)
ggplot(data=nb_sample78, aes(x=fct_reorder(gene_mutation, order), y=n)) +
  geom_bar(stat="identity")+
  ylab("number of samples")+
  theme(axis.text.x = element_text(angle = 90, size = 5), axis.title.x = element_blank())+
  scale_y_continuous(breaks = seq(0,max(nb_sample78$n)+ 1, by=1))
dev.off()