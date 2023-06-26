######## Comparaison patients immuno deprimes vs personel 
detach(package:dplyr)
library(plyr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggforce)
library(ggrepel)
library(gridExtra)

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

## Get only the samples we use in our analysis
infos<-infos[which(infos$sample%in%nextclade$seqName),]

## count number of samples per variant
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
tab_effectifs<-infos %>% group_by(virus_nexclade) %>% count(individual)

tab_effectifs_sort<-arrange(tab_effectifs, virus_nexclade, desc(individual))

tab_effectifs_cumsum <- ddply(tab_effectifs_sort, "virus_nexclade",
                   transform, label_ypos=cumsum(n) - 0.5*n)

png(filename="Effectifs_variants.png",  width=12, height=8,units="in",res = 200) 
ggplot(data=tab_effectifs_cumsum, aes(x=virus_nexclade, y=n, fill=individual)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos +5, label=n), vjust=1.6, 
            color="black", size=3)+
  ylab("Number of samples")+
  scale_fill_discrete(labels=c('Patient', 'Control'))+
  theme_Publication()+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90),legend.title = element_blank())
dev.off()

pdf(file="Effectifs_variants.pdf",  width=8, height=6) 
ggplot(data=tab_effectifs_cumsum, aes(x=virus_nexclade, y=n, fill=individual)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos +5, label=n), vjust=1.6, 
            color="black", size=3.5)+
  ylab("Number of samples")+
  scale_fill_brewer(labels=c('Patient', 'Control'),palette = "Set2")+
  theme_Publication()+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90),legend.title = element_blank())
dev.off()

write.csv(tab_effectifs,file = "effectifs_patients_control_per_virus.csv",row.names = F)

tab_effectifs_pitie<-infos %>% filter(origin=="pitie") %>% group_by(virus_nexclade) %>% count(individual)
write.csv(tab_effectifs_pitie,file = "effectifs_patients_control_per_virus_Pitie.csv",row.names = F)

tab_effectifs_bichat<-infos %>% filter(origin=="bichat") %>% group_by(virus_nexclade) %>% count(individual)
write.csv(tab_effectifs_bichat,file = "effectifs_patients_control_per_virus_Bichat.csv",row.names = F)


## add variant info to mutation table
tab$variant<-infos$virus_nexclade[match(tab$Id_sample,infos$sample)]


## import mutations specific of each variant (extracted from covariants)
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
mut_20I<-read.table(file="Mut_20I_Alpha_V1.txt",sep=":",header = F)
mut_21J<-read.table(file="Mut_21J_Delta.txt",sep=":",header = F)
mut_21K<-read.table(file="Mut_21K_Omicron.txt",sep=":",header = F)
mut_21L<-read.table(file="Mut_21L_Omicron.txt",sep=":",header = F)
mut_22B<-read.table(file="Mut_22B_Omicron.txt",sep=":",header = F)
mut_22E<-read.table(file="Mut_22E_Omicron.txt",sep=":",header = F)
mut_22F<-read.table(file="Mut_22F_Omicron.txt",sep=":",header = F)



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

# 
# # ### TEST
# set<-infos
# mutations<-tab
# Gene<-"S"
# genetab<-genes
# variant<-"21J_Delta"

# variant =
# "20I_Alpha_V1" "21J_Delta" "21K_Omicron"  "21L_Omicron"  "22B_Omicron"        
# "21I_Delta" "21A_Delta"  "21D_Eta"  "20B"     "19B"   "22E_Omicron"  "20C"   "20D"  "22C_Omicron"  "recombinant"  "22A_Omicron"  "20E_EU1"     "20A"          "20H_Beta_V2" 


### function to count the number of patients and personnel that have each mutation (per gene and per variant)
## set= "info" tab to subset the samples per variant
## mutations = "tab" mutation table
## Gene = gene we want to study
## genetab = "genes" list of all the genes 
## variant = variant we want to study
mutation_comp<- function(set, mutations, Gene,genetab, variant){
  # selection of the gene 
  if (Gene != "all"){
    genetab<-genetab[genetab$GeneId==Gene,]
  }
  
  # subset of samples of the chosen variant and gene
  set<-set[which(set$virus_nexclade==variant),]
  mutations<- mutations[which(mutations$gene == genetab$GeneId & mutations$Id_sample %in% set$sample),]
  mutations$group<-set$individual[match(mutations$Id_sample,set$sample)] # group = patient or personnel
  
  # if there are mutations for the chosen gene and variant 
  if(nrow(mutations)>0){
    # creates a data frame containing if the mutation is major, minor or Abs (absent)
    # the samples per row and the mutations per column (+ 2 columns for sample name and group)
    tab_mut<-as.data.frame(matrix("Abs",length(unique(mutations$Id_sample)),length(unique(mutations$mutation))+2))
    colnames(tab_mut)<-c("Id_sample","group",unique(mutations$mutation))
    tab_mut$Id_sample<-unique(mutations$Id_sample)
    tab_mut$group<-mutations$group[match(tab_mut$Id_sample,mutations$Id_sample)]
    
    # count of the number of major and minor samples per column (= per mutation)
    for(i in 3:ncol(tab_mut)){
      tmp<-mutations[which(mutations$mutation==colnames(tab_mut[i])),]
      tmp_maj<-tmp[which(tmp$obs=="major"),]
      tmp_min<-tmp[which(tmp$obs=="minor"),]
      tab_mut[which(tab_mut$Id_sample%in%tmp_maj$Id_sample),i]<-"major"
      tab_mut[which(tab_mut$Id_sample%in%tmp_min$Id_sample),i]<-"minor"
      
      count_tab<-tab_mut %>% group_by(group) %>% count(tab_mut[,i])
      # data frame to store the count of major, minor and total patient and personnel samples par mutation
      mut_count<-data.frame(mutation=NA,mut_spe=NA,nb_major_Patient=0,nb_minor_Patient=0,nb_mut_Patient=0,nb_total_Patient=0,nb_major_Personnel=0,nb_minor_Personnel=0,nb_mut_Personnel=0,nb_total_Personnel=0)
      mut_count$mutation<-colnames(tab_mut[i])
      try(mut_count$nb_major_Patient<-count_tab$n[count_tab$group=="patient_J0" & count_tab$`tab_mut[, i]`=="major"]) # use "try" to skip if there are no mutations in this case
      try(mut_count$nb_minor_Patient<-count_tab$n[count_tab$group=="patient_J0" & count_tab$`tab_mut[, i]`=="minor"])
      mut_count$nb_total_Patient<-length(which(tab_mut$group=="patient_J0"))
      mut_count$nb_mut_Patient<- mut_count$nb_major_Patient + mut_count$nb_minor_Patient
      try(mut_count$nb_major_Personnel<-count_tab$n[count_tab$group=="personnel" & count_tab$`tab_mut[, i]`=="major"])
      try(mut_count$nb_minor_Personnel<-count_tab$n[count_tab$group=="personnel" & count_tab$`tab_mut[, i]`=="minor"])
      mut_count$nb_total_Personnel<-length(which(tab_mut$group=="personnel"))
      mut_count$nb_mut_Personnel<- mut_count$nb_major_Personnel + mut_count$nb_minor_Personnel
      
      # concatenate together the counts for each mutation
      if (i == 3){
        recap_count<-mut_count
      } else {
        recap_count<- rbind(recap_count,mut_count)
      }
      # add the info if the mutation is characteristic of a variant 
      recap_count$mut_spe<-mutations$mut_spe[match(recap_count$mutation,mutations$mutation)]
      
    }
    
    recap_count$gene<-Gene
    
    return(recap_count)
  } else { # if there are no mutantations for this gene and variant we return an empty table 
    recap_count<-as.data.frame(matrix(nrow = 0,ncol = 11))
    colnames(recap_count)<-c("mutation", "mut_spe" ,"nb_major_Patient" ,"nb_minor_Patient",
                             "nb_mut_Patient","nb_total_Patient","nb_major_Personnel","nb_minor_Personnel",
                             "nb_mut_Personnel" ,"nb_total_Personnel","gene" )
    return(recap_count)
  }

}

### SPIKE ONLY

# mut_S_21K<-mutation_comp(infos, tab, "S", genes, "21K_Omicron")
# mut_S_21K_Pa<-mut_S_21K[mut_S_21K$nb_mut_Patient>1 & mut_S_21K$nb_mut_Personnel==0,]
# mut_S_21K_Pa_min<-mut_S_21K_Pa[mut_S_21K_Pa$nb_major_Patient==0,]
# mut_S_21K$variant<-"21K_Omicron"
# 
# 
# mut_S_21L<-mutation_comp(infos, tab, "S", genes, "21L_Omicron")
# mut_S_21L_Pa<-mut_S_21L[mut_S_21L$nb_mut_Patient>1 & mut_S_21L$nb_mut_Personnel==0,]
# mut_S_21L_Pa_min<-mut_S_21L_Pa[mut_S_21L_Pa$nb_major_Patient==0,]
# mut_S_21L$variant<-"21L_Omicron"
# 
# mut_S_22BK<-mutation_comp(infos, tab, "S", genes, "22B_Omicron")
# mut_S_22BK_Pa<-mut_S_22BK[mut_S_22BK$nb_mut_Patient>1 & mut_S_22BK$nb_mut_Personnel==0,]
# mut_S_22BK_Pa_min<-mut_S_22BK_Pa[mut_S_22BK_Pa$nb_major_Patient==0,]
# mut_S_22BK$variant<-"22B_Omicron"
# 
# mut_S_21J<-mutation_comp(infos, tab, "S", genes, "21J_Delta")
# mut_S_21J_Pa<-mut_S_21J[mut_S_21J$nb_mut_Patient>1 & mut_S_21J$nb_mut_Personnel==0,]
# mut_S_21J_Pa_min<-mut_S_21J_Pa[mut_S_21J_Pa$nb_major_Patient==0,]
# mut_S_21J$variant<-"21J_Delta"
# 
# mut_S_20I<-mutation_comp(infos, tab, "S", genes, "20I_Alpha_V1")
# mut_S_20I_Pa<-mut_S_20I[mut_S_20I$nb_mut_Patient>1 & mut_S_20I$nb_mut_Personnel==0,]
# mut_S_20I_Pa_min<-mut_S_20I_Pa[mut_S_20I_Pa$nb_major_Patient==0,]
# mut_S_20I$variant<-"20I_Alpha_V1"
# 
# 
# mut_all<-rbind(mut_S_20I,rbind(mut_S_21J,rbind(mut_S_22BK,rbind(mut_S_21L,mut_S_21K))))


### Patient/personnel mutation count (use of mutation_comp) for each gene and variant (in one table)
for (i in 1:nrow(genes)){
  gene_i<-genes$GeneId[i]
  print(gene_i)
  print("mut_gene_21K")
  mut_gene_21K<-mutation_comp(infos, tab,gene_i , genes, "21K_Omicron")
  mut_gene_21K$variant<-"21K_Omicron"
  
  print("mut_gene_21L")
  mut_gene_21L<-mutation_comp(infos, tab, gene_i, genes, "21L_Omicron")
  mut_gene_21L$variant<-"21L_Omicron"
  
  print("mut_gene_22BK")
  mut_gene_22BK<-mutation_comp(infos, tab, gene_i, genes, "22B_Omicron")
  mut_gene_22BK$variant<-"22B_Omicron"
  
  print("mut_gene_21J")
  mut_gene_21J<-mutation_comp(infos, tab, gene_i, genes, "21J_Delta")
  mut_gene_21J$variant<-"21J_Delta"
  
  print("mut_gene_20I")
  mut_gene_20I<-mutation_comp(infos, tab, gene_i, genes, "20I_Alpha_V1")
  if(nrow(mut_gene_20I)>0){
    mut_gene_20I$variant<-"20I_Alpha_V1"
  }

  
  
  mut_all<-rbind(mut_gene_20I,rbind(mut_gene_21J,rbind(mut_gene_22BK,rbind(mut_gene_21L,mut_gene_21K))))
  
  if (i==1){
      mut_tab<-mut_all
  } else {
    mut_tab<-rbind(mut_tab,mut_all)
  }

}

### selection of the mutations carried only by the patients and that are characteristic of a variant 
mut_tab_Pa<-mut_tab[which(mut_tab$nb_mut_Patient>1 & mut_tab$nb_mut_Personnel==0),]
mut_tab_Pa<-mut_tab_Pa[which(is.na(mut_tab_Pa$mut_spe)==F),]

# keep only the mutations characteristic of a variant that appeared after the samples' variant 
mut_tab_Pa_20I<-mut_tab_Pa[which(mut_tab_Pa$variant=="20I_Alpha_V1" & (str_detect(mut_tab_Pa$mut_spe,"Alpha"))==F),]
mut_tab_Pa_21J<-mut_tab_Pa[which(mut_tab_Pa$variant=="21J_Delta" & (str_detect(mut_tab_Pa$mut_spe,"Alpha"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21J_Delta"))==F),]
mut_tab_Pa_21K<-mut_tab_Pa[which(mut_tab_Pa$variant=="21K_Omicron" & (str_detect(mut_tab_Pa$mut_spe,"Alpha"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21J_Delta"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21K_Omicron"))==F),]
mut_tab_Pa_21L<-mut_tab_Pa[which(mut_tab_Pa$variant=="21L_Omicron" & (str_detect(mut_tab_Pa$mut_spe,"Alpha"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21J_Delta"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21K_Omicron"))==F & (str_detect(mut_tab_Pa$mut_spe,"21L_Omicron"))==F),]
mut_tab_Pa_22B<-mut_tab_Pa[which(mut_tab_Pa$variant=="22B_Omicron" & (str_detect(mut_tab_Pa$mut_spe,"Alpha"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21J_Delta"))==F  & (str_detect(mut_tab_Pa$mut_spe,"21K_Omicron"))==F & (str_detect(mut_tab_Pa$mut_spe,"21L_Omicron"))==F & (str_detect(mut_tab_Pa$mut_spe,"22B_Omicron"))==F),]

mut_tab_Pa2<-rbind(mut_tab_Pa_20I,rbind(mut_tab_Pa_21J,rbind(mut_tab_Pa_21K,rbind(mut_tab_Pa_21L,mut_tab_Pa_22B))))
write.csv(mut_tab_Pa2,file="Mutations_spe_variants_not_yet_present.csv",row.names = F)

# create a gene_mutation variable (T9I exists in gene E and ORF3a)
mut_tab_Pa2$gene_mutation<-paste(mut_tab_Pa2$gene,mut_tab_Pa2$mutation,sep="_")
tab$gene_mutation<-paste(tab$gene,tab$mutation,sep="_")

# get the samples with the mutations characteristic of later variants 
for(j in 1:nrow(mut_tab_Pa2)){
  tmp_tab<-tab[which(tab$gene_mutation==mut_tab_Pa2$gene_mutation[j] & tab$variant==mut_tab_Pa2$variant[j]),]
  if(j==1){
    tab_spe_mut<-tmp_tab
  } else {
    tab_spe_mut<-rbind(tab_spe_mut, tmp_tab)
  }
}
write.table(tab_spe_mut,file="Samples_with_mutations_spe_variants_not_yet_present.txt",row.names = F,sep="\t")

# graph of the samples that carry the mutations 
tab_spe_mut$variant_sample<-paste(tab_spe_mut$variant,tab_spe_mut$Id_sample)
tab_spe_mut$aa_pos<-as.numeric(str_sub(tab_spe_mut$mutation,2, -2))
tab_spe_mut$gene_start<-genes$Start[match(tab_spe_mut$gene,genes$GeneId)]
tab_spe_mut$gene_order<-tab_spe_mut$gene_start * 10000
tab_spe_mut$order<-tab_spe_mut$gene_order + tab_spe_mut$aa_pos

samples<-unique(tab_spe_mut$Id_sample)
tab_spe_mut$Sample<-NA
for(s in 1:length(samples)){
  tab_spe_mut$Sample[which(tab_spe_mut$Id_sample==samples[s])]<-paste("Patient",s,sep="_")
}

tab_spe_mut$mut_spe<-gsub("_Omicron","",tab_spe_mut$mut_spe)
tab_spe_mut$gene_mutation_var<-paste(tab_spe_mut$gene_mutation,tab_spe_mut$mut_spe,sep="_")

variants<-unique(tab_spe_mut$variant)
for(k in 1:length(variants)){
  tab_spe_mut_var<-tab_spe_mut[which(tab_spe_mut$variant==variants[k]),]
  #png(filename=paste("Mutations_spe_variants_after_",variants[k],".png", sep=""),  width=17, height=7,units="in",res = 200) 
  P1<-ggplot(tab_spe_mut_var) +
    geom_bar(aes(x = reorder(as.factor(gene_mutation),order), fill = as.factor(Sample)), position = "fill")+
    theme_Publication()+
    theme(legend.title = element_blank(), legend.text = element_text(size= 12),axis.title.x = element_blank(), axis.text.x = element_text(angle=90,size=12)) + ggtitle(variants[k])+
    scale_fill_manual(values=c("#787aeb","#7fde49", "#b85aeb","#cdcb46","#e24bc7","#65de92","#e94779","#60a441","#cc7dd5","#d28f3c",
                                       "#6994d2", "#e75b38", "#74dfd1","#d76da4","#c0d59a","#a689d3","#94924b","#dab5d9","#609d7e",
                                       "#cf776e","#6db2c9","#cda486","#a57d9d"))
  print(P1)
  #dev.off()
  
  if(k==1){
    Palpha<-P1
  } else if (k==2){
    tab_spe_mut_var<-tab_spe_mut_var[which(tab_spe_mut_var$Id_sample!="112111060593" & tab_spe_mut_var$Id_sample!="112112068328"),]
    P2<-ggplot(tab_spe_mut_var) +
      geom_bar(aes(x = reorder(as.factor(gene_mutation),order), fill = as.factor(Sample)), position = "fill")+
      theme_Publication()+
      theme(legend.title = element_blank(), legend.text = element_text(size= 12),axis.title.x = element_blank(), axis.text.x = element_text(angle=90,size=12)) + 
      ggtitle(variants[k])+
      scale_fill_manual(values=c("#787aeb","#7fde49", "#b85aeb","#cdcb46","#e24bc7","#65de92","#e94779","#60a441","#cc7dd5","#d28f3c",
                                          "#6994d2", "#e75b38", "#74dfd1","#d76da4","#c0d59a","#a689d3","#94924b","#dab5d9","#609d7e",
                                          "#cf776e","#6db2c9","#cda486","#a57d9d"))
    print(P2)
    
    Pdelta1<-P1
    Pdelta2<-P2
  } else if (k==3){
    Pomicron<-P1
  }
  
}

#png(filename=paste("Mutations_spe_variants_later_variants.png", sep=""),  width=10, height=12,units="in",res = 200)
pdf(file=paste("Mutations_spe_variants_later_variants.pdf", sep=""),  width=10, height=12)
grid.arrange(Palpha,Pdelta1,Pomicron, ncol=1, nrow=3) 
dev.off()

#png(filename=paste("Mutations_spe_variants_later_variants_noCoinf.png", sep=""),  width=10, height=12,units="in",res = 200) 
pdf(file=paste("Mutations_spe_variants_later_variants_noCoinf.pdf", sep=""),  width=10, height=12)
grid.arrange(Palpha,Pdelta2,Pomicron, ncol=1, nrow=3) 
dev.off()

# frequencies
tab_spe_mut %>% group_by(variant) %>% summarize(Mean = mean(freq, na.rm=TRUE))
tab_spe_mut %>% group_by(variant) %>% summarize(Med = median(freq, na.rm=TRUE))

### PCA with these variants highlighted 

setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
tab_full<-read.table(file = "Mutations_no_indel_J0_full_seq_1-2-3-4-5-6-7_avec_pitie2_min5-20pct.txt",header=T, sep="\t")

samples_remove<-read.table(file="samples_more_20pctNA.txt",sep="\t",header = F)
mut_remove<-read.table(file = "mut_more50pctNA_bichat_or_pitie.txt",sep="\t",header = F)

tab_full<-tab_full[which((tab_full$Id_sample%in%samples_remove[,1])==F & (tab_full$mutation%in%mut_remove[,1])==F),]


tab_pca<-as.data.frame(matrix(1,length(unique(tab_full$Id_sample)),length(unique(tab_full$mutation))))
colnames(tab_pca)<-unique(tab_full$mutation)
rownames(tab_pca)<-unique(tab_full$Id_sample)


for(i in 1:ncol(tab_pca)){
  tmp<-tab_full[which(tab_full$mutation==colnames(tab_pca[i])),]
  tmp_maj<-tmp[which(tmp$obs=="major"),]
  tmp_min<-tmp[which(tmp$obs=="minor"),]
  tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
  tab_pca[which(rownames(tab_pca)%in%tmp_maj$Id_sample),i]<-2
  tab_pca[which(rownames(tab_pca)%in%tmp_min$Id_sample),i]<-1.5
  tab_pca[which(rownames(tab_pca)%in%tmp_nocov$Id_sample),i]<- NA
  
}

rm(tmp)
rm(tmp_maj)
rm(tmp_min)
rm(tmp_nocov)

## remove non variable columns
tab_pca<-tab_pca[ , which(apply(tab_pca, 2, var,na.rm=T) != 0)]  

## PCA
library(nipals)

pca<-nipals(tab_pca,ncomp=10)


barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))


li<-as.data.frame(pca$score)
li$variant<-infos$clade[match(rownames(li),infos$sample)]
li$variant_N<-infos$virus_nexclade[match(rownames(li),infos$sample)]
li$variant_N[li$variant_N == ''] <- NA
li$sample_type<-infos$individual[match(rownames(li),infos$sample)]


color_set<-c("#B3B3B3",
                      "#FC4E2A",  "#BD0026", "#800026",
                      "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5" ,"#08519C" ,"#08306B" ,
                      "#984EA3","grey")
                      
li2<-li[which(li$variant_N!="19B" & li$variant_N!="20A" & li$variant_N!="20B" & li$variant_N!="20C" & li$variant_N!="20D" & li$variant_N!="20H_Beta_V2" & li$variant_N!="20E_EU1" & li$variant_N!="21D_Eta" ),]
li_mut_spe<-li2[which(rownames(li2)%in%unique(tab_spe_mut$Id_sample)),]


png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_Nextclade_J0_better_zoom_mut_spe_samples.png", width=16, height=12,units="in",res = 400)
ggplot(li2, aes(x=PC1 , y=PC2, color=variant_N, shape=sample_type))+
  geom_point(size=2, aes(shape=sample_type))+
  scale_color_manual(values = color_set)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())+
  geom_text_repel(data=li_mut_spe,aes(x=PC1,y=PC2,label=rownames(li_mut_spe)),max.overlaps = 100, min.segment.length = 0.01)
  
dev.off()



li3<-li[which(li$variant_N=="21K_Omicron"| li$variant_N=="21L_Omicron"|li$variant_N=="22B_Omicron"),]
li_mut_spe2<-li3[which(rownames(li3)%in%unique(tab_spe_mut$Id_sample)),]
color_set2<-c( "#6BAED6", "#2171B5" ,"#08306B" )

png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_Nextclade_J0_better_zoom_mut_spe_samples_omicron.png", width=16, height=12,units="in",res = 400)
ggplot(li3, aes(x=PC1 , y=PC2, color=variant_N, shape=sample_type))+
  geom_point(size=2, aes(shape=sample_type))+
  scale_color_manual(values = color_set2)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())+
  geom_text_repel(data=li_mut_spe2,aes(x=PC1,y=PC2,label=rownames(li_mut_spe2)),max.overlaps = 100, box.padding = 1)

dev.off()



##########################################################################################################################################
#### Boxplot nb mutation per gene patients vs personnels per variant 
boxplot_mut<-function(variant){
  mut_tab_var<-mut_tab[which(mut_tab$variant==variant),]
  # turn counts into percentages 
  mut_tab_var$pct_mut_Patient<-(mut_tab_var$nb_mut_Patient*100) / mut_tab_var$nb_total_Patient
  mut_tab_var$pct_major_Patient<-(mut_tab_var$nb_major_Patient*100) / mut_tab_var$nb_total_Patient
  mut_tab_var$pct_minor_Patient<-(mut_tab_var$nb_minor_Patient*100) / mut_tab_var$nb_total_Patient
  
  mut_tab_var$pct_mut_Personnel<-(mut_tab_var$nb_mut_Personnel*100) / mut_tab_var$nb_total_Personnel
  mut_tab_var$pct_major_Personnel<-(mut_tab_var$nb_major_Personnel*100) / mut_tab_var$nb_total_Personnel
  mut_tab_var$pct_minor_Personnel<-(mut_tab_var$nb_minor_Personnel*100) / mut_tab_var$nb_total_Personnel

  
  setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
  
  ## all mutations
  mut_tab2<-pivot_longer(mut_tab_var,cols = c(13,16), names_to = "patient_type",values_to = "pct_mut")
  mut_tab2$gene_start<-genes$Start[match(mut_tab2$gene,genes$GeneId)]
  
  # png(filename = paste("Nombre_mut_par_gene",variant,"patients_vs_personnels.png"), width=14, height=8,units="in",res = 200)
  # p1<-ggplot(mut_tab2, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
  #   geom_boxplot()+
  #   stat_compare_means(aes(group = patient_type), label = "p.format")+
  #   xlab("Genes") + ylab("% sample avec mutation") + ggtitle(paste(variant, "toutes mutations",sep=" "))+ 
  #   geom_point(position=position_jitterdodge())
  # print(p1)
  # dev.off()
  
  med1 <- ddply(mut_tab2, .(reorder(gene,gene_start), patient_type), summarize, med = median(pct_mut))
  med1$med<-round(med1$med,3)
  
  png(filename = paste("Nombre_mut_par_gene",variant,"patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
  p1<-ggplot(mut_tab2, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
    geom_boxplot()+
    stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
    xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste(variant, "all mutations",sep=" "))+
    geom_point(position=position_jitterdodge())+
    scale_fill_discrete(labels=c('Patients', 'Controls'))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text( size=10,),
          axis.title.x = element_blank(),
          legend.position="top")+
    geom_text(data = med1, aes(x = med1[,1], y = -5, label =med), 
              size = 3, vjust = -1.5, position = position_dodge(width = 1))
  print(p1)
  dev.off()
  
  ## mut maj
  mut_tab_var_maj<-mut_tab_var[which(mut_tab_var$nb_major_Patient != 0 | mut_tab_var$nb_major_Personnel != 0 ),]
  mut_tab_maj<-pivot_longer(mut_tab_var_maj,cols = c(14,17), names_to = "patient_type",values_to = "pct_mut")
  mut_tab_maj$gene_start<-genes$Start[match(mut_tab_maj$gene,genes$GeneId)]
  
  # png(filename = paste("Nombre_mut_Maj_par_gene",variant,"patients_vs_personnels.png"), width=14, height=8,units="in",res = 200)
  # p2<-ggplot(mut_tab_maj, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
  #   geom_boxplot()+
  #   stat_compare_means(aes(group = patient_type), label = "p.format")+
  #   xlab("Genes") + ylab("% sample avec mutation") + ggtitle(paste(variant, "mutations majeures",sep=" "))+
  #   geom_point(position=position_jitterdodge())
  # print(p2)
  # dev.off()
  
  ## mut min
  mut_tab_var_min<-mut_tab_var[which(mut_tab_var$nb_minor_Patient != 0 | mut_tab_var$nb_minor_Personnel != 0 ),]
  mut_tab_min<-pivot_longer(mut_tab_var_min,cols = c(15,18), names_to = "patient_type",values_to = "pct_mut")
  mut_tab_min$gene_start<-genes$Start[match(mut_tab_min$gene,genes$GeneId)]
  
  # png(filename = paste("Nombre_mut_Min_par_gene",variant,"patients_vs_personnels.png"), width=14, height=8,units="in",res = 200)
  # p3<-ggplot(mut_tab_min, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
  #   geom_boxplot()+
  #   stat_compare_means(aes(group = patient_type), label = "p.format")+
  #   xlab("Genes") + ylab("% sample avec mutation") + ggtitle(paste(variant, "mutations mineures",sep=" "))+ 
  #   geom_point(position=position_jitterdodge())
  # print(p3)
  # dev.off()
  # 
  
  med2 <- ddply(mut_tab_min, .(reorder(gene,gene_start), patient_type), summarize, med = median(pct_mut))
  med2$med<-round(med2$med,3)
  
  #png(filename = paste("Nombre_mut_Min_par_gene",variant,"patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
  pdf(file = paste("Nombre_mut_Min_par_gene",variant,"patients_vs_personnels_Clean.pdf"), width=14, height=8)
    p4<-ggplot(mut_tab_min, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
    geom_boxplot(outlier.shape = NA)+
    stat_compare_means(aes(group = patient_type), label = "p.signif",label.y = 105)+
    xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste(variant, "minority mutations",sep=" "))+
    geom_point(position=position_jitterdodge(),alpha=0.2)+
    scale_fill_brewer(labels=c('Patient', 'Control'),palette = "Set2")+
    theme_Publication()+
    theme(legend.title = element_blank(),
          # legend.text = element_text(size=10),
          # axis.text.x = element_text(size=10),
          # axis.text.y = element_text( size=10,),
          axis.title.x = element_blank(),
          legend.position="top")+
    geom_text(data = med2, aes(x = med2[,1], y = -5, label =med), 
              size = 3, vjust = -1.5, position = position_dodge(width = 1))
  print(p4)
  dev.off()

  
  # png(filename = paste("Nombre_mut_Min_par_gene",variant,"patients_vs_personnels_Poster.png"), width=14, height=8,units="in",res = 200)
  # p4<-ggplot(mut_tab_min, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
  #   geom_boxplot()+
  #   stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105, size=5)+
  #   xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste(variant, "minority mutations",sep=" "))+ 
  #   geom_point(position=position_jitterdodge())+
  #   scale_fill_discrete(labels=c('Patients', 'Controls'))+ 
  #   theme(legend.title = element_blank(),
  #         legend.text = element_text(size=16),
  #         axis.text.x = element_text(size=16),
  #         axis.text.y = element_text(size=16,),
  #         axis.title.y = element_text(size=16,),
  #         axis.title.x = element_blank(),
  #         legend.position="top")
  # print(p4)
  # dev.off()
  
  
  #### Only Spike
  mut_tab2_S<-mut_tab2[which(mut_tab2$gene=="S"),]
  
  med1_S <- ddply(mut_tab2_S, .(reorder(gene,gene_start), patient_type), summarize, med = median(pct_mut))
  med1_S$med<-round(med1_S$med,3)
  
  png(filename = paste("Nombre_mut_Spike",variant,"patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
  p1<-ggplot(mut_tab2_S, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
    geom_boxplot()+
    stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
    xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste(variant, "all mutations",sep=" "))+
    geom_point(position=position_jitterdodge())+
    scale_fill_discrete(labels=c('Patients', 'Controls'))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text( size=10,),
          axis.title.x = element_blank(),
          legend.position="top")+
    geom_text(data = med1_S, aes(x = med1_S[,1], y = -5, label =med), 
              size = 3, vjust = -1.5, position = position_dodge(width = 1))
  print(p1)
  dev.off()
  

  ## mut min
  mut_tab_min_S<-mut_tab_min[which(mut_tab_min$gene=="S"),]
  
  med2_S <- ddply(mut_tab_min_S, .(reorder(gene,gene_start), patient_type), summarize, med = median(pct_mut))
  med2_S$med<-round(med2_S$med,3)
  
  png(filename = paste("Nombre_mut_Min_Spike",variant,"patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
  p4<-ggplot(mut_tab_min_S, aes(x=reorder(gene,gene_start), y=pct_mut, fill=patient_type)) +
    geom_boxplot()+
    stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
    xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste(variant, "minority mutations",sep=" "))+
    geom_point(position=position_jitterdodge())+
    scale_fill_discrete(labels=c('Patients', 'Controls'))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text( size=10,),
          axis.title.x = element_blank(),
          legend.position="top")+
    geom_text(data = med2_S, aes(x = med2_S[,1], y = -5, label =med), 
              size = 3, vjust = -1.5, position = position_dodge(width = 1))
  print(p4)
  dev.off()
  
  
  
  
}


boxplot_mut("21K_Omicron")
boxplot_mut("21L_Omicron")
boxplot_mut("22B_Omicron")
boxplot_mut("21J_Delta")
boxplot_mut("20I_Alpha_V1")



#### Boxplot nb mutation per gene patients vs personnels all variants

mut_tab$pct_mut_Patient<-(mut_tab$nb_mut_Patient*100) / mut_tab$nb_total_Patient
mut_tab$pct_major_Patient<-(mut_tab$nb_major_Patient*100) / mut_tab$nb_total_Patient
mut_tab$pct_minor_Patient<-(mut_tab$nb_minor_Patient*100) / mut_tab$nb_total_Patient

mut_tab$pct_mut_Personnel<-(mut_tab$nb_mut_Personnel*100) / mut_tab$nb_total_Personnel
mut_tab$pct_major_Personnel<-(mut_tab$nb_major_Personnel*100) / mut_tab$nb_total_Personnel
mut_tab$pct_minor_Personnel<-(mut_tab$nb_minor_Personnel*100) / mut_tab$nb_total_Personnel


setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")

## all mutations
mut_tab2<-pivot_longer(mut_tab,cols = c(13,16), names_to = "patient_type",values_to = "pct_mut")
mut_tab2$gene_start<-genes$Start[match(mut_tab2$gene,genes$GeneId)]


med1 <- ddply(mut_tab2, .(variant, patient_type), summarize, med = median(pct_mut, na.rm = T))
med1$med<-round(med1$med,3)

png(filename = paste("Nombre_mut_par_gene_patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
p1<-ggplot(mut_tab2, aes(x=variant, y=pct_mut, fill=patient_type)) +
  geom_boxplot()+
  stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
  xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste( "all mutations",sep=" "))+
  geom_point(position=position_jitterdodge())+
  scale_fill_discrete(labels=c('Patients', 'Controls'))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text( size=10,),
        axis.title.x = element_blank(),
        legend.position="top")+
  geom_text(data = med1, aes(x = med1[,1], y = -5, label =med), 
            size = 3, vjust = -1.5, position = position_dodge(width = 1))
print(p1)
dev.off()

## mut min
mut_tab_min<-mut_tab[which(mut_tab$nb_minor_Patient != 0 | mut_tab$nb_minor_Personnel != 0 ),]
mut_tab_min<-pivot_longer(mut_tab_min,cols = c(15,18), names_to = "patient_type",values_to = "pct_mut")
mut_tab_min$gene_start<-genes$Start[match(mut_tab_min$gene,genes$GeneId)]


med2 <- ddply(mut_tab_min, .(variant, patient_type), summarize, med = median(pct_mut, na.rm = T))
med2$med<-round(med2$med,3)

png(filename = paste("Nombre_mut_Min_par_gene_patients_vs_personnels_Clean.png"), width=14, height=8,units="in",res = 200)
p4<-ggplot(mut_tab_min, aes(x=variant, y=pct_mut, fill=patient_type)) +
  geom_boxplot()+
  stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
  xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste("minority mutations",sep=" "))+
  geom_point(position=position_jitterdodge())+
  scale_fill_discrete(labels=c('Patients', 'Controls'))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text( size=10,),
        axis.title.x = element_blank(),
        legend.position="top")+
  geom_text(data = med2, aes(x = med2[,1], y = -5, label =med), 
            size = 3, vjust = -1.5, position = position_dodge(width = 1))
print(p4)
dev.off()


#### Without Alpha variant
mut_tab<-mut_tab[which(mut_tab$variant!="20I_Alpha_V1"),]

setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")

## all mutations
mut_tab2<-pivot_longer(mut_tab,cols = c(13,16), names_to = "patient_type",values_to = "pct_mut")
mut_tab2$gene_start<-genes$Start[match(mut_tab2$gene,genes$GeneId)]


med1 <- ddply(mut_tab2, .(variant, patient_type), summarize, med = median(pct_mut, na.rm = T))
med1$med<-round(med1$med,3)

png(filename = paste("Nombre_mut_par_gene_patients_vs_personnels_Delta_Omicron_Clean.png"), width=14, height=8,units="in",res = 200)
p1<-ggplot(mut_tab2, aes(x=variant, y=pct_mut, fill=patient_type)) +
  geom_boxplot()+
  stat_compare_means(aes(group = patient_type), label = "p.format",label.y = 105)+
  xlab("Genes") + ylab("% samples with mutation") + ggtitle(paste( "all mutations",sep=" "))+
  geom_point(position=position_jitterdodge())+
  scale_fill_discrete(labels=c('Patients', 'Controls'))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text( size=10,),
        axis.title.x = element_blank(),
        legend.position="top")+
  geom_text(data = med1, aes(x = med1[,1], y = -5, label =med), 
            size = 3, vjust = -1.5, position = position_dodge(width = 1))
print(p1)
dev.off()

## mut min
mut_tab_min<-mut_tab[which(mut_tab$nb_minor_Patient != 0 | mut_tab$nb_minor_Personnel != 0 ),]
mut_tab_min<-pivot_longer(mut_tab_min,cols = c(15,18), names_to = "patient_type",values_to = "pct_mut")
mut_tab_min$gene_start<-genes$Start[match(mut_tab_min$gene,genes$GeneId)]


med2 <- ddply(mut_tab_min, .(variant, patient_type), summarize, med = median(pct_mut, na.rm = T))
med2$med<-round(med2$med,3)

#png(filename = paste("Nombre_mut_Min_par_gene_patients_vs_personnels_Delta_Omicron_Clean_publi.png"), width=14, height=8,units="in",res = 200)
pdf(file = paste("Nombre_mut_Min_par_gene_patients_vs_personnels_Delta_Omicron_Clean_publi.pdf"), width=10, height=8)
p4<-ggplot(mut_tab_min, aes(x=variant, y=pct_mut, fill=patient_type)) +
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(aes(group = patient_type), label = "p.signif",label.y = 105)+
  xlab("Genes") + ylab("% samples with mutation") + 
  #ggtitle(paste("minority mutations",sep=" "))+
  geom_point(position=position_jitterdodge(),alpha=0.2)+
  scale_fill_brewer(labels=c('Patient', 'Control'),palette = "Set2")+
  theme_Publication()+
  theme(legend.title = element_blank(),
        # legend.text = element_text(size=10),
        # axis.text.x = element_text(size=10),
        # axis.text.y = element_text( size=10,),
        axis.title.x = element_blank(), legend.key.size = unit(0.5, "cm"))+
  geom_text(data = med2, aes(x = med2[,1], y = -5, label =med), 
            size = 3, vjust = -1.5, position = position_dodge(width = 1))
print(p4)
dev.off()

