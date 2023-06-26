########################## look at the potential coinfections ###################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)

## Ggplot2 theme importation
setwd("~/Projets/SIID/Scripts/")
source(file = "Function_theme_Publication.R")

## import mutation table
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
tab<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie_min5-20pct .txt",sep="\t",header = T)

## Remove filtered samples and mutations
samples_remove<-read.table(file="samples_more_20pctNA.txt",sep="\t",header = F)
tab<-tab[which((tab$Id_sample%in%samples_remove$V1)==F),]

mut_remove<-read.table(file="mut_more50pctNA_bichat_or_pitie.txt",sep="\t", header = F)
tab<-tab[which((tab$mutation%in%mut_remove$V1)==F),]


## import infos on genes 
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")

## import infos on samples
infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos<-infos[infos$individual!="patient_suivi",]

## import variants of the samples
setwd("~/Projets/SIID/2_Filter_trim/Accurate")
nextclade_16<-read.csv(file="nextclade_1-2-3-4-5-6-7.csv",sep=";")
nextclade_pitie<-read.csv(file="nextclade_pitie.csv",sep=";")
nextclade<-rbind(nextclade_16,nextclade_pitie)

infos$virus_nexclade<-nextclade$clade_legacy[match(infos$sample,as.character(nextclade$seqName))]
infos$virus_nexclade<- gsub(" ", "_",
                            gsub("\\(","",
                                 gsub("\\)","",
                                      gsub(",","", infos$virus_nexclade))))






## import mutations specific of each variant (extracted from covariants)
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
mut_20I<-read.table(file="Mut_20I_Alpha_V1.txt",sep=":",header = F)
mut_21J<-read.table(file="Mut_21J_Delta.txt",sep=":",header = F)
mut_21K<-read.table(file="Mut_21K_Omicron.txt",sep=":",header = F)
mut_21L<-read.table(file="Mut_21L_Omicron.txt",sep=":",header = F)
mut_22B<-read.table(file="Mut_22B_Omicron.txt",sep=":",header = F)
mut_22E<-read.table(file="Mut_22E_Omicron.txt",sep=":",header = F)
mut_22F<-read.table(file="Mut_22F_Omicron.txt",sep=":",header = F)


## list of variants sorted by time of apparition
var_samples<-unique(infos$virus_nexclade[which((is.na(infos$virus_nexclade)==F) & infos$virus_nexclade!="recombinant")])
var_samples<-sort(var_samples)
var_samples<-c(var_samples,"22F_Omicron")
variants<-data.frame(variant=var_samples, order=seq(1,length(var_samples),by=1))


### identify the mutations characteristic of variants in our data (tab)
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
  if(tab$mutation[m] %in% mut_22E$V2[which(mut_22E$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"22E_Omicron"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "22E_Omicron", sep=";")
    }
  }
  if(tab$mutation[m] %in% mut_22F$V2[which(mut_22F$V1==tab$gene[m])]){
    if(is.na(tab$mut_spe[m])){
      tab$mut_spe[m]<-"22F_Omicron"
    } else {
      tab$mut_spe[m]<-paste(tab$mut_spe[m], "22F_Omicron", sep=";")
    }
  }
}

#### Freq median of major mutants
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")

tab_maj<-tab[which(tab$obs=="major"),]
tab$variant<-infos$virus_nexclade[match(tab$Id_sample,infos$sample)]
med_maj_freq<-tab_maj %>% group_by(Id_sample) %>% summarize( med_freq = median(freq))

tab_coinf<-tab[which(tab$Id_sample%in%med_maj_freq$Id_sample[which(med_maj_freq$med_freq<80)]),]
samples_coinf<-unique(tab_coinf$Id_sample)
infos_coinf<-infos[which(infos$sample%in%samples_coinf),]

#### frequency graph of the possible coinfections
for (i in 1:length(samples_coinf)){
  tab_tmp<-tab[which(tab$Id_sample==samples_coinf[i]),]
  
  tab_tmp$gene_mutation<-paste(tab_tmp$gene,tab_tmp$mutation)
  tab_tmp$aa_pos<-as.numeric(str_sub(tab_tmp$mutation,2, -2))
  tab_tmp$gene_start<-genes$Start[match(tab_tmp$gene,genes$GeneId)]
  tab_tmp$gene_order<-tab_tmp$gene_start * 10000
  tab_tmp$order<-tab_tmp$gene_order + tab_tmp$aa_pos
  
  ### If we know the variant we can extract the previous and future variants
  if(tab_tmp$variant[1]!="recombinant" | is.na(tab_tmp$variant[1])==F){
    tab_tmp$variant_mut_spe<-"Other"
    ## define current, previous and future variants
    tab_tmp$variant_mut_spe[which(tab_tmp$mut_spe==tab_tmp$variant[1])]<-"Current variant"
    prev_var<-variants$variant[which(variants$order<variants$order[which(variants$variant==tab_tmp$variant[1])])]
    next_var<-variants$variant[which(variants$order>variants$order[which(variants$variant==tab_tmp$variant[1])])]
    
    tab_tmp$nb_prev<-NA
    tab_tmp$nb_current<-NA
    tab_tmp$nb_next<-NA
    
    ## count the number of mutations in each category
    for(l in 1:nrow(tab_tmp)){
      tab_tmp$nb_prev[l]<-length(which( sapply (prev_var, grep, tab_tmp$mut_spe[l]) != 0))
      tab_tmp$nb_current[l]<-length(which( sapply (tab_tmp$variant[1], grep, tab_tmp$mut_spe[l]) != 0))
      tab_tmp$nb_next[l]<-length(which( sapply (next_var, grep, tab_tmp$mut_spe[l]) != 0))
    }
    
    ## classifies the mutations
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev>0 & tab_tmp$nb_current==0 & tab_tmp$nb_next==0)]<-"Prev variant"
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev>0 & tab_tmp$nb_current>0 & tab_tmp$nb_next==0)]<-"Prev and current variant"
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev>0 & tab_tmp$nb_current>0 & tab_tmp$nb_next>0)]<-"Prev, current and future variant"
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev>0 & tab_tmp$nb_current==0 & tab_tmp$nb_next>0)]<-"Prev and future variant"
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev==0 & tab_tmp$nb_current>0 & tab_tmp$nb_next>0)]<-"Current and future variant"
    tab_tmp$variant_mut_spe[which(tab_tmp$nb_prev==0 & tab_tmp$nb_current==0 & tab_tmp$nb_next>0)]<-"Future variant"
    
    ## attributes a color per category 
    tab_tmp$color<-"#B3B3B3"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Current variant")]<-"#B3DE69"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Prev variant")]<-"#FB8072"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Prev and current variant")]<-"#FDB462" 
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Prev, current and future variant")]<-"#8DD3C7"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Prev and future variant")]<-"#80B1D3"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Current and future variant")]<-"#8470ff"
    tab_tmp$color[which(tab_tmp$variant_mut_spe=="Future variant")]<-"#BC80BD"
    
    tab_tmp<-tab_tmp[order(tab_tmp$variant_mut_spe),]
    
    
    color_set<-unique(tab_tmp$color)
    
    ## plot mutations frequencies with color by category
    png(filename=paste("Mutations_frequencies_",samples_coinf[i],".png",sep=""),  width=17, height=7,units="in",res = 200) 
    p0<-ggplot(tab_tmp, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=variant_mut_spe)) +
      #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
      geom_point(aes(color=variant_mut_spe))+
      #scale_shape_manual(values=c(15, 16, 17,18))+
      scale_color_manual(values=color_set)+
      ggtitle(paste(samples_coinf[i],tab_tmp$variant[1]))+
      theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
    print(p0)
    dev.off()
    
  }
  
}

###################### Manual case studies : delta omicron potencial coinfections ##########################
genes<-genes[order(genes$Start),]
#genes$color<-c("#cc4a6d", "#c95733", "#c9a443", "#a4db52", "#5d8d37", "#61d89d", "#5c8ed9", "#605ec0", "#7f47d2", "#d082d0", "#923a85", "#d746c2")

#### 112111060593 delta + omicron?

tab_112111060593<-tab[which(tab$Id_sample=="112111060593"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112111060593$mut_delta_omicron<-"Other"
tab_112111060593$mut_delta_omicron[which((str_detect(tab_112111060593$mut_spe,"Alpha")==F) & (str_detect(tab_112111060593$mut_spe,"21J_Delta")==T) & (str_detect(tab_112111060593$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112111060593$mut_delta_omicron[which((str_detect(tab_112111060593$mut_spe,"Alpha")==F) & (str_detect(tab_112111060593$mut_spe,"21J_Delta")==T) & (str_detect(tab_112111060593$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112111060593$mut_delta_omicron[which((str_detect(tab_112111060593$mut_spe,"Alpha")==F) & (str_detect(tab_112111060593$mut_spe,"21J_Delta")==F) & (str_detect(tab_112111060593$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

## plot mutations frequencies with color by variant
tab_112111060593$gene_mutation<-paste(tab_112111060593$gene,tab_112111060593$mutation)
tab_112111060593$aa_pos<-as.numeric(str_sub(tab_112111060593$mutation,2, -2))
tab_112111060593$gene_start<-genes$Start[match(tab_112111060593$gene,genes$GeneId)]
tab_112111060593$gene_order<-tab_112111060593$gene_start * 10000
tab_112111060593$order<-tab_112111060593$gene_order + tab_112111060593$aa_pos

#col_set<-genes$color[which(genes$GeneId%in%tab_112111060593$gene)]

tab_112111060593$gene<-gsub("ORF","",tab_112111060593$gene)

#png(filename="Mutations_spe_Delta_Omicron_112111060593.png",  width=7, height=17,units="in",res = 200) 
pdf(file="Mutations_spe_Delta_Omicron_112111060593.pdf",  width=7, height=17)
p1<-ggplot(tab_112111060593, aes(x=freq, y=reorder(as.factor(mutation),order), group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron),size=3)+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  #ggtitle("112111060593_21J_Delta")+
  facet_grid(reorder(gene,gene_start) ~ .,scales = "free_y", space ="free_y" ) +
  theme_Publication()+
  theme(legend.position="top", axis.text.y = element_text( size=10),axis.title.y = element_blank(),legend.title = element_blank(),
        strip.text.y = element_text(size = 10))
print(p1)
dev.off()






##### 112112068328 delta + omicron?

tab_112112068328<-tab[which(tab$Id_sample=="112112068328"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112112068328$mut_delta_omicron<-"Other"
tab_112112068328$mut_delta_omicron[which((str_detect(tab_112112068328$mut_spe,"Alpha")==F) & (str_detect(tab_112112068328$mut_spe,"21J_Delta")==T) & (str_detect(tab_112112068328$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112112068328$mut_delta_omicron[which((str_detect(tab_112112068328$mut_spe,"Alpha")==F) & (str_detect(tab_112112068328$mut_spe,"21J_Delta")==T) & (str_detect(tab_112112068328$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112112068328$mut_delta_omicron[which((str_detect(tab_112112068328$mut_spe,"Alpha")==F) & (str_detect(tab_112112068328$mut_spe,"21J_Delta")==F) & (str_detect(tab_112112068328$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

tab_112112068328$gene_mutation<-paste(tab_112112068328$gene,tab_112112068328$mutation)
tab_112112068328$aa_pos<-as.numeric(str_sub(tab_112112068328$mutation,2, -2))
tab_112112068328$gene_start<-genes$Start[match(tab_112112068328$gene,genes$GeneId)]
tab_112112068328$gene_order<-tab_112112068328$gene_start * 10000
tab_112112068328$order<-tab_112112068328$gene_order + tab_112112068328$aa_pos

tab_112112068328$gene<-gsub("ORF","",tab_112112068328$gene)

## plot mutations frequencies with color by variant
#png(filename="Mutations_spe_Delta_Omicron_112112068328.png",  width=17, height=7,units="in",res = 200) 
pdf(file="Mutations_spe_Delta_Omicron_112112068328.pdf",  width=7, height=17)
p2<-ggplot(tab_112112068328, aes(x=freq, y=reorder(as.factor(mutation),order), group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  facet_grid(reorder(gene,gene_start) ~ .,scales = "free_y", space ="free_y" ) +
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron),size=3)+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  #ggtitle("112111060593_21J_Delta")+
  theme_Publication()+
  theme(legend.position="top", axis.text.y = element_text( size=10),axis.title.y = element_blank(),legend.title = element_blank(),
        strip.text.y = element_text(size = 10))
print(p2)
dev.off()


png(filename=paste("Mutations_spe_Delta_Omicron_112111060593_112112068328.png", sep=""),  width=15, height=10,units="in",res = 200) 
grid.arrange(p1,p2, ncol=1, nrow=2) 
dev.off()


##### 112112059092 Delta: not a coinfection case for reference 

tab_112112059092<-tab[which(tab$Id_sample=="112112059092"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112112059092$mut_delta_omicron<-"Other"
tab_112112059092$mut_delta_omicron[which((str_detect(tab_112112059092$mut_spe,"Alpha")==F) & (str_detect(tab_112112059092$mut_spe,"21J_Delta")==T) & (str_detect(tab_112112059092$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112112059092$mut_delta_omicron[which((str_detect(tab_112112059092$mut_spe,"Alpha")==F) & (str_detect(tab_112112059092$mut_spe,"21J_Delta")==T) & (str_detect(tab_112112059092$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112112059092$mut_delta_omicron[which((str_detect(tab_112112059092$mut_spe,"Alpha")==F) & (str_detect(tab_112112059092$mut_spe,"21J_Delta")==F) & (str_detect(tab_112112059092$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

tab_112112059092$gene_mutation<-paste(tab_112112059092$gene,tab_112112059092$mutation)
tab_112112059092$aa_pos<-as.numeric(str_sub(tab_112112059092$mutation,2, -2))
tab_112112059092$gene_start<-genes$Start[match(tab_112112059092$gene,genes$GeneId)]
tab_112112059092$gene_order<-tab_112112059092$gene_start * 10000
tab_112112059092$order<-tab_112112059092$gene_order + tab_112112059092$aa_pos

## plot mutations frequencies with color by variant
png(filename="Mutations_spe_Delta_Omicron_112112059092.png",  width=17, height=7,units="in",res = 200) 
p3<-ggplot(tab_112112059092, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112112059092_21J_Delta")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p3)
dev.off()



#### 112201080111 Omicron: not a coinfection case for reference

tab_112201080111<-tab[which(tab$Id_sample=="112201080111"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112201080111$mut_delta_omicron<-"Other"
tab_112201080111$mut_delta_omicron[which((str_detect(tab_112201080111$mut_spe,"Alpha")==F) & (str_detect(tab_112201080111$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201080111$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112201080111$mut_delta_omicron[which((str_detect(tab_112201080111$mut_spe,"Alpha")==F) & (str_detect(tab_112201080111$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201080111$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112201080111$mut_delta_omicron[which((str_detect(tab_112201080111$mut_spe,"Alpha")==F) & (str_detect(tab_112201080111$mut_spe,"21J_Delta")==F) & (str_detect(tab_112201080111$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

tab_112201080111$gene_mutation<-paste(tab_112201080111$gene,tab_112201080111$mutation)
tab_112201080111$aa_pos<-as.numeric(str_sub(tab_112201080111$mutation,2, -2))
tab_112201080111$gene_start<-genes$Start[match(tab_112201080111$gene,genes$GeneId)]
tab_112201080111$gene_order<-tab_112201080111$gene_start * 10000
tab_112201080111$order<-tab_112201080111$gene_order + tab_112201080111$aa_pos

## plot mutations frequencies with color by variant
png(filename="Mutations_spe_Delta_Omicron_112201080111.png",  width=17, height=7,units="in",res = 200) 
p4<-ggplot(tab_112201080111, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(16, 17,18))+
  scale_color_manual(values=c("purple","red","grey"))+
  ggtitle("112201080111_21K_Omicron")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p4)
dev.off()


########################## 112111060074 delta + omicron?

tab_112111060074<-tab[which(tab$Id_sample=="112111060074"),]

tab_112111060074$mut_delta_omicron<-"Other"
tab_112111060074$mut_delta_omicron[which((str_detect(tab_112111060074$mut_spe,"Alpha")==F) & (str_detect(tab_112111060074$mut_spe,"21J_Delta")==T) & (str_detect(tab_112111060074$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112111060074$mut_delta_omicron[which((str_detect(tab_112111060074$mut_spe,"Alpha")==F) & (str_detect(tab_112111060074$mut_spe,"21J_Delta")==T) & (str_detect(tab_112111060074$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112111060074$mut_delta_omicron[which((str_detect(tab_112111060074$mut_spe,"Alpha")==F) & (str_detect(tab_112111060074$mut_spe,"21J_Delta")==F) & (str_detect(tab_112111060074$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

tab_112111060074$gene_mutation<-paste(tab_112111060074$gene,tab_112111060074$mutation)
tab_112111060074$aa_pos<-as.numeric(str_sub(tab_112111060074$mutation,2, -2))
tab_112111060074$gene_start<-genes$Start[match(tab_112111060074$gene,genes$GeneId)]
tab_112111060074$gene_order<-tab_112111060074$gene_start * 10000
tab_112111060074$order<-tab_112111060074$gene_order + tab_112111060074$aa_pos

png(filename="Mutations_spe_Delta_Omicron_112111060074.png",  width=17, height=7,units="in",res = 200) 
p5<-ggplot(tab_112111060074, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112111060074_21K_Omicron")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p5)
dev.off()

########### portent une série de mut omicron

#### 112108020456 delta + omicron?

tab_112108020456<-tab[which(tab$Id_sample=="112108020456"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112108020456$mut_delta_omicron<-"Other"
tab_112108020456$mut_delta_omicron[which((str_detect(tab_112108020456$mut_spe,"Alpha")==F) & (str_detect(tab_112108020456$mut_spe,"21J_Delta")==T) & (str_detect(tab_112108020456$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112108020456$mut_delta_omicron[which((str_detect(tab_112108020456$mut_spe,"Alpha")==F) & (str_detect(tab_112108020456$mut_spe,"21J_Delta")==T) & (str_detect(tab_112108020456$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112108020456$mut_delta_omicron[which((str_detect(tab_112108020456$mut_spe,"Alpha")==F) & (str_detect(tab_112108020456$mut_spe,"21J_Delta")==F) & (str_detect(tab_112108020456$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

## plot mutations frequencies with color by variant
tab_112108020456$gene_mutation<-paste(tab_112108020456$gene,tab_112108020456$mutation)
tab_112108020456$aa_pos<-as.numeric(str_sub(tab_112108020456$mutation,2, -2))
tab_112108020456$gene_start<-genes$Start[match(tab_112108020456$gene,genes$GeneId)]
tab_112108020456$gene_order<-tab_112108020456$gene_start * 10000
tab_112108020456$order<-tab_112108020456$gene_order + tab_112108020456$aa_pos

png(filename="Mutations_spe_Delta_Omicron_112108020456.png",  width=17, height=7,units="in",res = 200) 
p6<-ggplot(tab_112108020456, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112108020456_21J_Delta")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p6)
dev.off()



#### 112108020460 delta + omicron?

tab_112108020460<-tab[which(tab$Id_sample=="112108020460"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112108020460$mut_delta_omicron<-"Other"
tab_112108020460$mut_delta_omicron[which((str_detect(tab_112108020460$mut_spe,"Alpha")==F) & (str_detect(tab_112108020460$mut_spe,"21J_Delta")==T) & (str_detect(tab_112108020460$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112108020460$mut_delta_omicron[which((str_detect(tab_112108020460$mut_spe,"Alpha")==F) & (str_detect(tab_112108020460$mut_spe,"21J_Delta")==T) & (str_detect(tab_112108020460$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112108020460$mut_delta_omicron[which((str_detect(tab_112108020460$mut_spe,"Alpha")==F) & (str_detect(tab_112108020460$mut_spe,"21J_Delta")==F) & (str_detect(tab_112108020460$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

## plot mutations frequencies with color by variant
tab_112108020460$gene_mutation<-paste(tab_112108020460$gene,tab_112108020460$mutation)
tab_112108020460$aa_pos<-as.numeric(str_sub(tab_112108020460$mutation,2, -2))
tab_112108020460$gene_start<-genes$Start[match(tab_112108020460$gene,genes$GeneId)]
tab_112108020460$gene_order<-tab_112108020460$gene_start * 10000
tab_112108020460$order<-tab_112108020460$gene_order + tab_112108020460$aa_pos

png(filename="Mutations_spe_Delta_Omicron_112108020460.png",  width=17, height=7,units="in",res = 200) 
p7<-ggplot(tab_112108020460, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112108020460_21J_Delta")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p7)
dev.off()





#### 112201016115 delta + omicron?

tab_112201016115<-tab[which(tab$Id_sample=="112201016115"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112201016115$mut_delta_omicron<-"Other"
tab_112201016115$mut_delta_omicron[which((str_detect(tab_112201016115$mut_spe,"Alpha")==F) & (str_detect(tab_112201016115$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201016115$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112201016115$mut_delta_omicron[which((str_detect(tab_112201016115$mut_spe,"Alpha")==F) & (str_detect(tab_112201016115$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201016115$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112201016115$mut_delta_omicron[which((str_detect(tab_112201016115$mut_spe,"Alpha")==F) & (str_detect(tab_112201016115$mut_spe,"21J_Delta")==F) & (str_detect(tab_112201016115$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

## plot mutations frequencies with color by variant
tab_112201016115$gene_mutation<-paste(tab_112201016115$gene,tab_112201016115$mutation)
tab_112201016115$aa_pos<-as.numeric(str_sub(tab_112201016115$mutation,2, -2))
tab_112201016115$gene_start<-genes$Start[match(tab_112201016115$gene,genes$GeneId)]
tab_112201016115$gene_order<-tab_112201016115$gene_start * 10000
tab_112201016115$order<-tab_112201016115$gene_order + tab_112201016115$aa_pos

png(filename="Mutations_spe_Delta_Omicron_112201016115.png",  width=17, height=7,units="in",res = 200) 
p8<-ggplot(tab_112201016115, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112201016115_21J_Delta")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p8)
dev.off()



#### 112201046704 delta + omicron?

tab_112201046704<-tab[which(tab$Id_sample=="112201046704"),]

## of which variant(s) (all Omicron regrouped in one category) are the mutations characteristic of? 
tab_112201046704$mut_delta_omicron<-"Other"
tab_112201046704$mut_delta_omicron[which((str_detect(tab_112201046704$mut_spe,"Alpha")==F) & (str_detect(tab_112201046704$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201046704$mut_spe,"Omicron")==F))]<-"Mut_Delta"
tab_112201046704$mut_delta_omicron[which((str_detect(tab_112201046704$mut_spe,"Alpha")==F) & (str_detect(tab_112201046704$mut_spe,"21J_Delta")==T) & (str_detect(tab_112201046704$mut_spe,"Omicron")==T))]<-"Mut_Delta_Omicron"
tab_112201046704$mut_delta_omicron[which((str_detect(tab_112201046704$mut_spe,"Alpha")==F) & (str_detect(tab_112201046704$mut_spe,"21J_Delta")==F) & (str_detect(tab_112201046704$mut_spe,"Omicron")==T))]<-"Mut_Omicron"

## plot mutations frequencies with color by variant
tab_112201046704$gene_mutation<-paste(tab_112201046704$gene,tab_112201046704$mutation)
tab_112201046704$aa_pos<-as.numeric(str_sub(tab_112201046704$mutation,2, -2))
tab_112201046704$gene_start<-genes$Start[match(tab_112201046704$gene,genes$GeneId)]
tab_112201046704$gene_order<-tab_112201046704$gene_start * 10000
tab_112201046704$order<-tab_112201046704$gene_order + tab_112201046704$aa_pos

png(filename="Mutations_spe_Delta_Omicron_112201046704.png",  width=17, height=7,units="in",res = 200) 
p8<-ggplot(tab_112201046704, aes(x=reorder(as.factor(gene_mutation),order), y=freq, group=mut_delta_omicron)) +
  #geom_line(aes(linetype=mut_delta_omicron, color=mut_delta_omicron))+
  geom_point(aes(color=mut_delta_omicron,  shape=mut_delta_omicron))+
  scale_shape_manual(values=c(15, 16, 17,18))+
  scale_color_manual(values=c("blue","purple","red","grey"))+
  ggtitle("112201046704_21J_Delta")+
  theme(legend.position="top", axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),legend.title = element_blank())
print(p8)
dev.off()