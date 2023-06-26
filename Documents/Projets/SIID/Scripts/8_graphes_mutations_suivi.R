###### Mutations graphs and figures
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/Romane/Documents/Projets/SID/data")

genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
genes$length<-genes$End - genes$Start

infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos$virus<-infos$clade
infos$virus<- gsub(" ", "_",
                   gsub("\\(","",
                        gsub("\\)","",
                             gsub(",","", infos$virus))))


setwd("C:/Users/Romane/Documents/Projets/SID")

tab_recap<-read.table(file="Recap_mutations no_indel_suivi .txt",sep = "\t", header = T)
tab_recap$Gene_start<-genes$Start[match(tab_recap$Gene,genes$GeneId)]

tab<-read.table(file="Mutations no_indel_suivi .txt",sep="\t",header = T)
tab$gene_start<-genes$Start[match(tab$gene,genes$GeneId)]


setwd("C:/Users/Romane/Documents/Projets/SID/Suivi/")
#nbe mutations par gene
# 
# png(filename = "Boxplot_mut_per_gene_suivi.png",width = 800,height = 600)
# ggplot(tab_recap, aes(x=reorder(Gene,Gene_start) , y=nb_total))+
#   geom_boxplot() +
#   xlab("Genes") + ylab("nb mutations per sample") + ggtitle("Mutations per gene")
# dev.off()
# 

# nbe major and minor 
tab_recap2<-pivot_longer(tab_recap,cols = 2:3, names_to = "mut_type",values_to = "nb_mut")
tab_recap2$nb_mut_pond<-tab_recap2$nb_mut/genes$length[match(tab_recap2$Gene,genes$GeneId)]

png(filename = "Boxplot_mut_per_gene_major_minor_suivi.png",width = 800,height = 600)
ggplot(tab_recap2, aes(x=reorder(Gene,Gene_start) , y=nb_mut, fill=mut_type))+
  geom_boxplot(position = position_dodge()) +
  xlab("Genes") + ylab("nb mutations per sample") + ggtitle("Mutations per gene")
dev.off()

png(filename = "Boxplot_mut_per_gene_major_minor_suivi_pond.png",width = 800,height = 600)
ggplot(tab_recap2, aes(x=reorder(Gene,Gene_start) , y=nb_mut_pond, fill=mut_type))+
  geom_boxplot(position = position_dodge()) + 
  #geom_point(position=position_jitterdodge())+
  xlab("Genes") + ylab("nb mutations per sample / length of the ORF or gene") + ggtitle("Mutations per gene")
dev.off()


#nbe mutations par sample
png(filename = "Nb_mutation_per_sample_suivi.png", width=20, height=6,units="in",res = 200)
ggplot(tab,aes(x=as.factor(Id_sample)))+
  geom_bar(stat = "count")+
  theme(axis.text.x = element_text(angle=90,size = 8),axis.title.x = element_blank())
dev.off()

# genome coverage
png(filename = "Gene_coverage_suivi.png",width = 600,height = 600)
ggplot(tab_recap,aes(x=gene_cov))+
  geom_bar(stat = "bin", binwidth = 1)
dev.off()

partial_gene_cov<-tab_recap[which(tab_recap$gene_cov!=100),]
length(unique(partial_gene_cov$Id_sample))
length(unique(partial_gene_cov$Gene))


# freq mutation 
png(filename = "Boxplot_freq_per_gene_suivi.png",width = 800,height = 600)
ggplot(tab, aes(x=reorder(gene,gene_start) , y=freq, fill=obs))+
  geom_boxplot() +
  xlab("Genes") + ylab("frequence") 
dev.off()



# Nbe samples per mutation

tab_mut<- tab %>% group_by(gene) %>% count(mutation)
tab_mut %>% ungroup() %>% arrange(desc(n)) %>% slice(1:10)

samples_per_mutation<-function(input_tab,name){
  tab2<- input_tab %>% group_by(gene) %>% count(mutation)
  tab2$gene_start<-genes$Start[match(tab2$gene,genes$GeneId)]
  # 
  # png(filename = paste("Nb_samples_per_mutation_J0_",name,".png",sep=""),width = 800,height = 600)
  # P1<-ggplot(tab2, aes(x = n, fill= reorder(gene,gene_start))) +
  #   geom_histogram(binwidth=10)+
  #   xlab("Nb samples with mutation") + ylab("Nb mutations")+
  #   theme(legend.title = element_blank())
  # print(P1)
  # dev.off()

  # on major mutations
  tab_maj<-input_tab[input_tab$obs=="major",]
  
  tab_maj2<- tab_maj %>% group_by(gene) %>% count(mutation)
  tab_maj2$gene_start<-genes$Start[match(tab_maj2$gene,genes$GeneId)]
  tab_maj2$obs<-"major"
  
  
  tab_min<-input_tab[input_tab$obs=="minor",]
  
  tab_min2<- tab_min %>% group_by(gene) %>% count(mutation)
  tab_min2$gene_start<-genes$Start[match(tab_min2$gene,genes$GeneId)]
  tab_min2$obs<-"minor"
  
  tab_maj_min<-rbind(tab_maj2,tab_min2)

  png(filename = paste("Nb_samples_per_mutation_major_vs_minor_suivi_",name,".png",sep=""),width = 800,height = 600)
  P2<-ggplot(tab_maj_min, aes(x = n, fill= obs)) +
    geom_histogram(binwidth=1)+
    xlab("Nb samples with mutation") + ylab("Nb mutations") + ggtitle(name)
  theme(legend.title = element_blank())
  print(P2)
  dev.off()

  # png(filename = paste("Nb_samples_per_mutation_major_vs_minor_zoom_J0_",name,".png",sep=""),width = 800,height = 600)
  # P3<-ggplot(tab_maj_min, aes(x = n, fill= obs)) +
  #   geom_histogram(binwidth=1)+
  #   xlab("Nb samples with mutation") + ylab("Nb mutations")+ xlim(0,50)+ ggtitle(name)
  # theme(legend.title = element_blank())
  # print(P3)
  # dev.off()

  tab2 %>% ungroup() %>% arrange(desc(n)) %>% head(20)
  tab2 %>% ungroup() %>% arrange(desc(n)) %>% write.table(file = paste("Number_samples_per_mutation_",name,".txt",sep=""),row.names = F,quote=F,sep="\t")
}

tab_20I<-tab[which(tab$Id_sample%in%infos$sample[which(infos$virus=="20I_Alpha_V1")]),]
samples_per_mutation(tab_20I,"20I_Alpha_V1")

tab_21J<-tab[which(tab$Id_sample%in%infos$sample[which(infos$virus=="21J_Delta")]),]
samples_per_mutation(tab_21J,"21J_Delta")


tab_21K<-tab[which(tab$Id_sample%in%infos$sample[which(infos$virus=="21K_Omicron")]),]
samples_per_mutation(tab_21K,"21K_Omicron")

tab_21L<-tab[which(tab$Id_sample%in%infos$sample[which(infos$virus=="21L_Omicron")]),]
samples_per_mutation(tab_21L,"21L_Omicron")

tab_22B<-tab[which(tab$Id_sample%in%infos$sample[which(infos$virus=="22B_Omicron")]),]
samples_per_mutation(tab_22B,"22B_Omicron")


# ACP 
#list_samples<-c("1121110126472", "112201016115", "112201067922")
setwd("C:/Users/Romane/Documents/Projets/SID")
tab_full<-read.table(file = "Mutations_no_indel_Suivi_full.txt",header=T, sep="\t")

tab_pca<-as.data.frame(matrix(0,length(unique(tab_full$Id_sample)),length(unique(tab_full$mutation))))
colnames(tab_pca)<-unique(tab_full$mutation)
rownames(tab_pca)<-unique(tab_full$Id_sample)

# tab_pca_mino<-as.data.frame(matrix(0,length(unique(tab$Id_sample)),length(unique(tab$mutation))))
# colnames(tab_pca_mino)<-unique(tab$mutation)
# rownames(tab_pca_mino)<-unique(tab$Id_sample)

for(i in 1:ncol(tab_pca)){
  tmp<-tab_full[which(tab_full$mutation==colnames(tab_pca[i])),]
  tmp_maj<-tmp[which(tmp$obs=="major"),]
  tmp_min<-tmp[which(tmp$obs=="minor"),]
  tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
  tab_pca[which(rownames(tab_pca)%in%tmp_maj$Id_sample),i]<-1
  tab_pca[which(rownames(tab_pca)%in%tmp_min$Id_sample),i]<-0.5
  tab_pca[which(rownames(tab_pca)%in%tmp_nocov$Id_sample),i]<- NA
  
  # tab_pca_mino[which(rownames(tab_pca_mino)%in%tmp_min$Id_sample),i]<-1
}

rm(tmp)
rm(tmp_maj)
rm(tmp_min)
rm(tmp_nocov)

tab_pca<-tab_pca[ , which(apply(tab_pca, 2, var,na.rm=T) != 0)]

#library(ade4)

#pca<-dudi.pca(tab_pca,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

library(nipals)

pca<-nipals(tab_pca,ncomp=3)

barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))


li<-as.data.frame(pca$scores)
li$variant<-infos$clade[match(rownames(li),infos$sample)]
li$sample_type<-infos$individual[match(rownames(li),infos$sample)]

color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                      "#FC4E2A", "#E7298A", "#BD0026", "#800026",
                      "#88B2BEFF", "#2F8FA8FF", "#4A879CFF", "#316283FF", "#345A78FF" ,"#0F3B6CFF" ,"#081F39FF" ,
                      "#984EA3")
                      
png(filename = "PCA_all_genes_Suivi_Axes12.png",width = 800,height = 800)
ggplot(li, aes(x=PC1 , y=PC2, color=variant , label= rownames(li)))+
  geom_point(size=2)+
  scale_color_manual(values = color_set)
dev.off()

png(filename = "PCA_all_genes_Suivi_Axes23.png",width = 800,height = 800)
ggplot(li, aes(x=PC3 , y=PC2, color=variant , label= rownames(li)))+
  geom_point(size=2)+
  scale_color_manual(values = color_set)
dev.off()

# 
# # no outliers
# li2<-li[which(li$Axis1>-50),]
# png(filename = "PCA_all_genes_J0_zoom.png",width = 800,height = 800)
# ggplot(li2, aes(x=Axis1 , y=Axis2, color=variant))+
#   geom_point(size=2)+
#   scale_color_manual(values = color_set)
# dev.off()
# 

# ACP SPike

tab_S<-tab_full[tab_full$gene=="S",]
tab_pca_S<-as.data.frame(matrix(0,length(unique(tab_S$Id_sample)),length(unique(tab_S$mutation))))
colnames(tab_pca_S)<-unique(tab_S$mutation)
rownames(tab_pca_S)<-unique(tab_S$Id_sample)


for(i in 1:ncol(tab_pca_S)){
  tmp<-tab_full[which(tab_full$mutation==colnames(tab_pca_S[i])),]
  tmp_maj<-tmp[which(tmp$obs=="major"),]
  tmp_min<-tmp[which(tmp$obs=="minor"),]
  tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
  tab_pca_S[which(rownames(tab_pca_S)%in%tmp_maj$Id_sample),i]<-1
  tab_pca_S[which(rownames(tab_pca_S)%in%tmp_min$Id_sample),i]<-0.5
  tab_pca_S[which(rownames(tab_pca_S)%in%tmp_nocov$Id_sample),i]<- NA
  
}

rm(tmp)
rm(tmp_maj)
rm(tmp_min)
rm(tmp_nocov)

tab_pca_S<-tab_pca_S[ , which(apply(tab_pca_S, 2, var,na.rm=T) != 0)]

pca_S<-nipals(tab_pca_S,ncomp=3)
#pca_S<-dudi.pca(tab_pca_S,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

li_S<-as.data.frame(pca_S$scores)
li_S$variant<-infos$clade[match(rownames(li_S),infos$sample)]
li_S$sample_type<-infos$individual[match(rownames(li_S),infos$sample)]

color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                      "#FC4E2A", "#E7298A", "#BD0026", "#800026",
                      "#88B2BEFF", "#2F8FA8FF", "#4A879CFF", "#316283FF", "#345A78FF" ,"#0F3B6CFF" ,"#081F39FF" ,
                      "#984EA3")
                      
png(filename = "PCA_spike_J0.png",width = 800,height = 800)
ggplot(li_S, aes(x=PC1 , y=PC2, color=variant))+
  geom_point(size=2)+
  scale_color_manual(values = color_set)
dev.off()






## Données manquantes

tab_NA<-tab_full[which(tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"),]
nb_mut<-length(unique(paste(tab_full$mutation, tab_full$gene)))
nb_sample<-length(unique(tab_full$Id_sample))
# par échantillon
NA_sample<-tab_full %>% group_by(Id_sample) %>% count
NA_sample$Nb<-NA
for(s in 1:nrow(NA_sample)){
  NA_sample$Nb[s]<-length(tab_full$Id_sample[which(tab_full$Id_sample==NA_sample$Id_sample[s]& (tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"))])
}
NA_sample$pct_NA<-(NA_sample$Nb * 100) / nb_mut

png(filename = "Missing_data_pers_sample_Suivi.png", width=20, height=8,units="in",res = 400)
ggplot(data=NA_sample, aes(x=as.factor(Id_sample), y=pct_NA)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle=90,size = 9))
dev.off()


infos_NA<-infos[which(infos$sample%in%NA_sample$Id_sample[which(NA_sample$pct_NA>=10)]),]

# par mutation
NA_mut<-tab_full %>% group_by(mutation) %>% count()
NA_mut$Nb<-NA
for(s in 1:nrow(NA_mut)){
  NA_mut$Nb[s]<-length(tab_full$mutation[which(tab_full$mutation==NA_mut$mutation[s]& (tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"))])
}
NA_mut$pct_NA<-(NA_mut$Nb * 100) / nb_sample

png(filename = "Missing_data_pers_mut_Suivi.png", width=8, height=8,units="in",res = 400)
ggplot(NA_mut, aes(x=pct_NA)) + 
  geom_histogram(binwidth=1)
dev.off()





# 
# # ACP sans SPike
# 
# tab_nS<-tab[tab$gene!="S",]
# tab_pca_nS<-as.data.frame(matrix(0,length(unique(tab_S$Id_sample)),length(unique(tab_nS$mutation))))
# colnames(tab_pca_nS)<-unique(tab_nS$mutation)
# rownames(tab_pca_nS)<-unique(tab_nS$Id_sample)
# 
# 
# for(i in 1:ncol(tab_pca_nS)){
#   tmp<-tab[which(tab$mutation==colnames(tab_pca_nS[i])),]
#   tmp_maj<-tmp[which(tmp$obs=="major"),]
#   tmp_min<-tmp[which(tmp$obs=="minor"),]
#   tab_pca_nS[which(rownames(tab_pca_nS)%in%tmp_maj$Id_sample),i]<-1
#   tab_pca_nS[which(rownames(tab_pca_nS)%in%tmp_min$Id_sample),i]<-0.5
#   
# }
# 
# rm(tmp)
# rm(tmp_maj)
# rm(tmp_min)
# 
# pca_nS<-dudi.pca(tab_pca_nS)
# 
# li_nS<-pca_nS$li
# li_nS$variant<-infos$clade[match(rownames(li_nS),infos$sample)]
# 
# color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
#              "#FC4E2A", "#E7298A", "#BD0026", "#800026",
#              "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5" ,"#08519C" ,"#08306B" ,
#              "#984EA3")
# 
# png(filename = "PCA_non_spike_J0_Axes12.png",width = 800,height = 800)
# ggplot(li_nS, aes(x=Axis1 , y=Axis2, color=variant))+
#   geom_point(size=2)+
#   scale_color_manual(values = color_set)
# dev.off()
# 
# png(filename = "PCA_non_spike_J0_Axes23.png",width = 800,height = 800)
# ggplot(li_nS, aes(x=Axis3 , y=Axis2, color=variant))+
#   geom_point(size=2)+
#   scale_color_manual(values = color_set)
# dev.off()
# 
# li_nS2<-li_nS[which(li_nS$Axis1<50),]
# 
# png(filename = "PCA_non_spike_J0_zoom.png",width = 800,height = 800)
# ggplot(li_nS2, aes(x=Axis1 , y=Axis2, color=variant))+
#   geom_point(size=2)+
#   scale_color_manual(values = color_set)
# dev.off()
