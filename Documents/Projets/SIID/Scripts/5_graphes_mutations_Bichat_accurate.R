###### Mutations graphs and figures
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(forcats)

## Ggplot2 theme importation
setwd("~/Projets/SIID/Scripts/")
source(file = "Function_theme_Publication.R")

## Work repertory and files 
work_rep<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/"

setwd("C:/Users/Romane/Documents/Projets/SIID/data")

genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
genes$length<-genes$End - genes$Start

infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos$virus<-infos$clade
infos$virus<- gsub(" ", "_",
                   gsub("\\(","",
                        gsub("\\)","",
                             gsub(",","", infos$virus))))

metadata<-read_xlsx("Metadata_VF.xlsx",sheet="Metadata_copie")

infos$type<-metadata$Type[match(infos$dossier,metadata$STAT_DOSSIER)]

setwd(work_rep)

tab_recap<-read.table(file="Recap_mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie_min5-20pct .txt",sep = "\t", header = T)
tab_recap$Gene_start<-genes$Start[match(tab_recap$Gene,genes$GeneId)]

tab<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie_min5-20pct .txt",sep="\t",header = T)
tab$gene_start<-genes$Start[match(tab$gene,genes$GeneId)]


#setwd(paste(work_rep,"J0/",sep=""))

#### nbe mutations par gene
png(filename = "Boxplot_mut_per_gene_J0.png", width=10, height=8,units="in",res = 400)
ggplot(tab_recap, aes(x=reorder(Gene,Gene_start) , y=nb_total))+
  geom_boxplot() +
  xlab("Genes") + ylab("nb mutations per sample") + ggtitle("Mutations per gene")
dev.off()


# nbe major and minor 
tab_recap2<-pivot_longer(tab_recap,cols = 2:3, names_to = "mut_type",values_to = "nb_mut")
tab_recap2$nb_mut_pond<-tab_recap2$nb_mut/genes$length[match(tab_recap2$Gene,genes$GeneId)]

png(filename = "Boxplot_mut_per_gene_major_minor_J0.png", width=10, height=8,units="in",res = 400)
ggplot(tab_recap2, aes(x=reorder(Gene,Gene_start) , y=nb_mut, fill=mut_type))+
  geom_boxplot(position = position_dodge()) +
  xlab("Genes") + ylab("nb mutations per sample") + ggtitle("Mutations per gene")
dev.off()


#png(filename = "Boxplot_mut_per_gene_major_minor_J0_pond.png", width=10, height=8,units="in",res = 400)
pdf(file = "Boxplot_mut_per_gene_major_minor_J0_pond.pdf", width=10, height=8)
ggplot(tab_recap2, aes(x=reorder(Gene,Gene_start) , y=nb_mut_pond, fill=mut_type))+
  geom_boxplot(position = position_dodge()) +
  scale_fill_brewer(palette = "Set2")+
  xlab("Genes") + ylab("nb mutations per sample / length of the ORF or gene") + ggtitle("Mutations per gene")+
  theme_Publication()
dev.off()



pdf(file = "Dotplot_mut_per_gene_major_minor_J0_pond.pdf", width=10, height=8)
ggplot(tab_recap2, aes(x=reorder(Gene,Gene_start) , y=nb_mut_pond, color=mut_type))+
  geom_point(position=position_jitterdodge())+
  scale_color_brewer(palette = "Set2")+
  xlab("Genes") + ylab("nb mutations per sample / length of the ORF or gene") + ggtitle("Mutations per gene")+
  theme_Publication()

dev.off()


#### genome coverage
png(filename = "Gene_coverage_J0.png", width=8, height=8,units="in",res = 400)
ggplot(tab_recap,aes(x=gene_cov))+
  geom_bar(stat = "bin", binwidth = 1)
dev.off()

partial_gene_cov<-tab_recap[which(tab_recap$gene_cov!=100),]
length(unique(partial_gene_cov$Id_sample))
length(unique(partial_gene_cov$Gene))


#### freq mutation 
png(filename = "Boxplot_freq_per_gene_J0.png", width=10, height=8,units="in",res = 400)
ggplot(tab, aes(x=reorder(gene,gene_start) , y=freq, fill=obs))+
  geom_boxplot() +
  xlab("Genes") + ylab("frequence") 
dev.off()



#### Nbe samples per mutation

tab_mut<- tab %>% group_by(gene) %>% count(mutation)
tab_mut %>% ungroup() %>% arrange(desc(n)) %>% slice(1:10)

samples_per_mutation<-function(input_tab,name){
  tab2<- input_tab %>% group_by(gene) %>% count(mutation)
  tab2$gene_start<-genes$Start[match(tab2$gene,genes$GeneId)]
  
  png(filename = paste("Nb_samples_per_mutation_J0_",name,".png",sep=""), width=8, height=8,units="in",res = 400)
  P1<-ggplot(tab2, aes(x = n, fill= reorder(gene,gene_start))) +
    geom_histogram(binwidth=10)+
    xlab("Nb samples with mutation") + ylab("Nb mutations")+
    theme(legend.title = element_blank())
  print(P1)
  dev.off()

  # on major mutations
  tab_maj<-input_tab[input_tab$obs=="major",]
  
  tab_maj2<- tab_maj %>% group_by(gene) %>% count(mutation)
  tab_maj2$gene_start<-genes$Start[match(tab_maj2$gene,genes$GeneId)]
  tab_maj2$obs<-"major"
  
  # on minor mutations
  tab_min<-input_tab[input_tab$obs=="minor",]
  
  tab_min2<- tab_min %>% group_by(gene) %>% count(mutation)
  tab_min2$gene_start<-genes$Start[match(tab_min2$gene,genes$GeneId)]
  tab_min2$obs<-"minor"
  
  tab_maj_min<-rbind(tab_maj2,tab_min2)

  png(filename = paste("Nb_samples_per_mutation_major_vs_minor_J0_",name,".png",sep=""), width=8, height=8,units="in",res = 400)
  P2<-ggplot(tab_maj_min, aes(x = n, fill= obs)) +
    geom_histogram(binwidth=10)+
    xlab("Nb samples with mutation") + ylab("Nb mutations") + ggtitle(name)
  theme(legend.title = element_blank())
  print(P2)
  dev.off()

  png(filename = paste("Nb_samples_per_mutation_major_vs_minor_zoom_J0_",name,".png",sep=""), width=8, height=8,units="in",res = 400)
  P3<-ggplot(tab_maj_min, aes(x = n, fill= obs)) +
    geom_histogram(binwidth=1)+
    xlab("Nb samples with mutation") + ylab("Nb mutations")+ xlim(0,50)+ ggtitle(name)
  theme(legend.title = element_blank())
  print(P3)
  dev.off()

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


############################################################################################################
#### Missing data, PCA and heatmap are now done separately in : 5_Missing_data_PCA_heatmap.R 
############################################################################################################



# ###### ACP 
# setwd(work_rep)
# tab_full<-read.table(file = "Mutations_no_indel_J0_full_avec_pitie.txt",header=T, sep="\t")
# 
# 
# setwd(paste(work_rep,"J0/",sep=""))
# tab_pca<-as.data.frame(matrix(1,length(unique(tab_full$Id_sample)),length(unique(tab_full$mutation))))
# colnames(tab_pca)<-unique(tab_full$mutation)
# rownames(tab_pca)<-unique(tab_full$Id_sample)
# 
# 
# for(i in 1:ncol(tab_pca)){
#   tmp<-tab_full[which(tab_full$mutation==colnames(tab_pca[i])),]
#   tmp_maj<-tmp[which(tmp$obs=="major"),]
#   tmp_min<-tmp[which(tmp$obs=="minor"),]
#   tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
#   tab_pca[which(rownames(tab_pca)%in%tmp_maj$Id_sample),i]<-2
#   tab_pca[which(rownames(tab_pca)%in%tmp_min$Id_sample),i]<-1.5
#   tab_pca[which(rownames(tab_pca)%in%tmp_nocov$Id_sample),i]<- NA
#   
# }
# 
# rm(tmp)
# rm(tmp_maj)
# rm(tmp_min)
# rm(tmp_nocov)
# 
# 
# #library(ade4)
# 
# #pca<-dudi.pca(tab_pca,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
# 
# library(nipals)
# 
# pca<-nipals(tab_pca,ncomp=10)
# 
# 
# barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# 
# 
# li<-as.data.frame(pca$score)
# li$variant<-infos$clade[match(rownames(li),infos$sample)]
# li$sample_type<-infos$individual[match(rownames(li),infos$sample)]
# 
# color_set<-c("#66C2A5", "#FC8D62", "#E78AC3", "#A6D854", "#FFD92F", "#B3B3B3",
#              "#FC4E2A",  "#BD0026", "#800026",
#              "#88B2BEFF", "#2F8FA8FF", "#4A879CFF", "#316283FF", "#345A78FF" ,"#0F3B6CFF" ,"#081F39FF" ,
#              "#984EA3")
# 
# 
# png(filename = "PCA_all_genes_J0_Axes12_patients_personnels.png", width=8, height=8,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=variant , label= rownames(li), shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+ 
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_better_zoom.png", width=16, height=12,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=variant, shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())+
#   facet_zoom(xlim = c(-0.15, 0.03),ylim = c(-0.005, 0.002))
# dev.off()
# 
# 
# ### distinction site
# li$origin<-infos$origin[match(rownames(li),infos$sample)]
# 
# png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_par_site.png", width=8, height=8,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=origin, shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_par_site_better_zoom.png", width=16, height=12,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=origin, shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())+
#   facet_zoom(xlim = c(-0.15, 0.03),ylim = c(-0.005, 0.002))
# dev.off()
# 
# # # no outliers
# # li2<-li[which(li$PC1>-0.1),]
# # 
# # li2$origin<-infos$origin[match(rownames(li2),infos$sample)]
# # 
# # png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_zoom.png", width=8, height=8,units="in",res = 400)
# # ggplot(li2, aes(x=PC1 , y=PC2, color=variant, shape=sample_type))+
# #   geom_point(size=2, aes(shape=sample_type))+
# #   scale_color_manual(values = color_set)+
# #   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
# #   labs(color = "Variant", shape= element_blank())
# # dev.off()
# 
# 
# 
# 
# 
# 
# #### ACP SPike
# 
# tab_S<-tab_full[tab_full$gene=="S",]
# tab_pca_S<-as.data.frame(matrix(0,length(unique(tab_S$Id_sample)),length(unique(tab_S$mutation))))
# colnames(tab_pca_S)<-unique(tab_S$mutation)
# rownames(tab_pca_S)<-unique(tab_S$Id_sample)
# 
# 
# for(i in 1:ncol(tab_pca_S)){
#   tmp<-tab_full[which(tab_full$mutation==colnames(tab_pca_S[i])),]
#   tmp_maj<-tmp[which(tmp$obs=="major"),]
#   tmp_min<-tmp[which(tmp$obs=="minor"),]
#   tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
#   tab_pca_S[which(rownames(tab_pca_S)%in%tmp_maj$Id_sample),i]<-1
#   tab_pca_S[which(rownames(tab_pca_S)%in%tmp_min$Id_sample),i]<-0.5
#   tab_pca_S[which(rownames(tab_pca_S)%in%tmp_nocov$Id_sample),i]<- NA
#   
#   }
# 
# rm(tmp)
# rm(tmp_maj)
# rm(tmp_min)
# rm(tmp_nocov)
# 
# pca_S<-nipals(tab_pca_S,ncomp=3)
# barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# 
# li_S<-as.data.frame(pca_S$scores)
# li_S$variant<-infos$clade[match(rownames(li_S),infos$sample)]
# li_S$sample_type<-infos$individual[match(rownames(li_S),infos$sample)]
# 
# 
# 
# png(filename = "PCA_spike_J0_patients_personnels.png", width=8, height=8,units="in",res = 400)
# ggplot(li_S, aes(x=PC1 , y=PC2, color=variant))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+ 
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# 
# 
# #### ACP Omicrons 
# 
# omicron_samples<-infos$sample[grep("Omicron",infos$virus)]
# tab_pca_O<-tab_pca[which(rownames(tab_pca)%in%omicron_samples),]
# 
# tab_pca_O<-tab_pca_O[ , which(apply(tab_pca_O, 2, var,na.rm=T) != 0)]
# 
# pca_O<-nipals(tab_pca_O,ncomp=3)
# 
# 
# barplot(pca_O$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# 
# 
# li_O<-as.data.frame(pca_O$score)
# li_O$variant<-infos$clade[match(rownames(li_O),infos$sample)]
# li_O$sample_type<-infos$individual[match(rownames(li_O),infos$sample)]
# 
# color_set_O<-c("#88B2BEFF", "#2F8FA8FF", "#4A879CFF", "#316283FF", "#345A78FF" ,"#0F3B6CFF" ,"#081F39FF" )
#                       
# 
# png(filename = "PCA_all_genes_J0_Axes12_patients_personnels_omicron.png", width=8, height=8,units="in",res = 400)
# ggplot(li_O, aes(x=PC1 , y=PC2, color=variant , label= rownames(li_O), shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set_O)+ 
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# ### distinction site
# li_O$origin<-infos$origin[match(rownames(li_O),infos$sample)]
# 
# png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_omicron_par_site.png", width=8, height=8,units="in",res = 400)
# ggplot(li_O, aes(x=PC1 , y=PC2, color=origin, shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# 
# 
# ## per immuno-depression cause
# pca_per_var<-function(variant){
#   samples_var<-infos$sample[which(infos$virus==variant)]
#   tab_var<-tab_full[which(tab_full$Id_sample%in%samples_var),]
#   tab_var$mutation_gene<-paste(tab_var$mutation, tab_var$gene,sep="_")
#   mut_var<-unique(tab_var$mutation_gene[which(tab_var$obs=="major" | tab_var$obs=="minor")])
#   tab_var<-tab_var[which(tab_var$mutation_gene%in%mut_var),]
#   
#   tab_pca_var<-as.data.frame(matrix(0,length(unique(tab_var$Id_sample)),length(unique(tab_var$mutation))))
#   colnames(tab_pca_var)<-unique(tab_var$mutation)
#   rownames(tab_pca_var)<-unique(tab_var$Id_sample)
#   
#   
#   for(i in 1:ncol(tab_pca_var)){
#     tmp<-tab_var[which(tab_var$mutation==colnames(tab_pca_var[i])),]
#     tmp_maj<-tmp[which(tmp$obs=="major"),]
#     tmp_min<-tmp[which(tmp$obs=="minor"),]
#     tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
#     tab_pca_var[which(rownames(tab_pca_var)%in%tmp_maj$Id_sample),i]<-1
#     tab_pca_var[which(rownames(tab_pca_var)%in%tmp_min$Id_sample),i]<-0.5
#     tab_pca_var[which(rownames(tab_pca_var)%in%tmp_nocov$Id_sample),i]<- NA
#     
#   }
#   
#   rm(tmp)
#   rm(tmp_maj)
#   rm(tmp_min)
#   rm(tmp_nocov)
#   
#   ## remove non variable columns
#   tab_pca_var<-tab_pca_var[ , which(apply(tab_pca_var, 2, var,na.rm=T) != 0)]  
# 
#   library(nipals)
# 
#   pca_var<-nipals(tab_pca_var, ncomp=3)
#   
#   
#   barplot(pca_var$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
#   
#   
#   li_var<-as.data.frame(pca_var$score)
#   li_var$variant<-infos$clade[match(rownames(li_var),infos$sample)]
#   li_var$sample_type<-infos$individual[match(rownames(li_var),infos$sample)]
#   li_var$type<-infos$type[match(rownames(li_var),infos$sample)]
#   
#   color_set<-c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
#                         
#   
#   
#   png(filename = paste("PCA_all_genes_J0_Axes12_patients_personnels_",variant,".png",sep=""), width=8, height=8,units="in",res = 400)
#   p1<-ggplot(li_var, aes(x=PC1 , y=PC2, color=type , shape=sample_type))+
#     geom_point(size=2, aes(shape=sample_type))+
#     scale_color_manual(values = color_set)+ 
#     scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))#+
#     #labs(color = "Variant", shape= element_blank())
#   print(p1)
#   dev.off()
#   
# }
# 
# pca_per_var("21K_Omicron")
# pca_per_var("21L_Omicron")
# pca_per_var("22B_Omicron")
# 
# 
# 
# 
# ##### Données manquantes
# setwd(work_rep)
# tab_full<-read.table(file = "Mutations_no_indel_J0_full.txt",header=T, sep="\t")
# 
# 
# tab_NA<-tab_full[which(tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"),]
# nb_mut<-length(unique(paste(tab_full$mutation, tab_full$gene)))
# nb_sample<-length(unique(tab_full$Id_sample))
# # par échantillon
# NA_sample<-tab_full %>% group_by(Id_sample) %>% count
# NA_sample$Nb<-NA
# for(s in 1:nrow(NA_sample)){
#   NA_sample$Nb[s]<-length(tab_full$Id_sample[which(tab_full$Id_sample==NA_sample$Id_sample[s]& (tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"))])
# }
# NA_sample$pct_NA<-(NA_sample$Nb * 100) / nb_mut
# 
# setwd(paste(work_rep,"J0/",sep=""))
# png(filename = "Missing_data_pers_sample_J0.png", width=25, height=8,units="in",res = 400)
# ggplot(data=NA_sample, aes(x=as.factor(Id_sample), y=pct_NA)) +
#   geom_bar(stat="identity")+
#   theme(axis.text.x = element_text(angle=90,size = 8),axis.title.x = element_blank())
# dev.off()
# 
# 
# infos_NA<-infos[which(infos$sample%in%NA_sample$Id_sample[which(NA_sample$pct_NA>=10)]),]
# 
# # par mutation
# NA_mut<-tab_full %>% group_by(mutation) %>% count()
# NA_mut$Nb<-NA
# for(s in 1:nrow(NA_mut)){
#   NA_mut$Nb[s]<-length(tab_full$mutation[which(tab_full$mutation==NA_mut$mutation[s]& (tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"))])
# }
# NA_mut$pct_NA<-(NA_mut$Nb * 100) / nb_sample
# 
# 
# png(filename = "Missing_data_pers_mut_J0.png", width=8, height=8,units="in",res = 400)
# ggplot(NA_mut, aes(x=pct_NA)) + 
#   geom_histogram(binwidth=1)
# dev.off()
# 
# 
# 
# 
# ##### Heatmap des variations
# 
# mut_tab_graph<-tab_full
# mut_tab_graph$mutation_gene<-paste(mut_tab_graph$mutation,mut_tab_graph$gene,sep="_")
# mut_tab_graph$gene_start<-genes$Start[match(mut_tab_graph$gene,genes$GeneId)]
# mut_tab_graph$gene_order<-mut_tab_graph$gene_start * 10000
# mut_tab_graph$order<-mut_tab_graph$gene_order + mut_tab_graph$AA_pos
# mut_tab_graph$origin<-infos$origin[match(mut_tab_graph$Id_sample,infos$sample)]
# mut_tab_graph$virus<-infos$virus[match(mut_tab_graph$Id_sample,infos$sample)]
# mut_tab_graph$origin_virus<-paste(mut_tab_graph$origin,mut_tab_graph$virus,sep="_")
# 
# mut_tab_graph$virus_sample<-paste(mut_tab_graph$virus,mut_tab_graph$Id_sample,sep="_")
# 
# png(filename = paste("heatmap_freq_mut.png"), width=20, height=10,units="in",res = 400)
# p1<-ggplot(data = mut_tab_graph, aes(x=fct_reorder(mutation_gene, order), y=as.factor(virus_sample), fill=freq)) + 
#   geom_tile(color = "black") +
#   ylab("jour")+
#   scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray")+
#   xlab("Mutation")+
#   theme(axis.text.x = element_text(angle=90,size = 2),axis.text.y = element_text(size = 12))
# print(p1)
# dev.off()
# 


