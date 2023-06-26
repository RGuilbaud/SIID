#################### Missing data, PCA and heatmap #########################
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(forcats)
library(stringr)

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

setwd("~/Projets/SIID/2_Filter_trim/Accurate")
nextclade_16<-read.csv(file="nextclade_1-2-3-4-5-6-7.csv",sep=";")
nextclade_pitie<-read.csv(file="nextclade_pitie.csv",sep=";")
nextclade<-rbind(nextclade_16,nextclade_pitie)

infos$virus_nexclade<-nextclade$clade_legacy[match(infos$sample,as.character(nextclade$seqName))]
infos$virus_nexclade[which(infos$sample=="fastq_total_662208034197")]<-NA
infos$pango_nexclade<-nextclade$Nextclade_pango[match(infos$sample,as.character(nextclade$seqName))]

setwd(work_rep)
tab_full<-read.table(file = "Mutations_no_indel_J0_full_seq_1-2-3-4-5-6-7_avec_pitie2_min5-20pct.txt",header=T, sep="\t")



##### Données manquantes

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

#setwd(paste(work_rep,"J0/",sep=""))
png(filename = "Missing_data_pers_sample_J0_seq1-2-3-4-5-6-7_pitie.png", width=25, height=8,units="in",res = 400)
ggplot(data=NA_sample, aes(x=as.factor(Id_sample), y=pct_NA)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle=90,size = 4),axis.title.x = element_blank())
dev.off()

## Samples with more than 20% NA to remove
samples_remove<-NA_sample$Id_sample[which(NA_sample$pct_NA>=20)]
tab_full<-tab_full[which((tab_full$Id_sample%in%samples_remove)==F),]
write.table(samples_remove,file="samples_more_20pctNA.txt",col.names = F, row.names = F, quote = F)

#infos_NA<-infos[which(infos$sample%in%NA_sample$Id_sample[which(NA_sample$pct_NA>=10)]),]

# par mutation
NA_mut<-tab_full %>% group_by(mutation) %>% count()
NA_mut$Nb<-NA
NA_mut$Nb_b<-NA
NA_mut$Nb_p<-NA
tab_full_b<-tab_full[which(tab_full$Id_sample%in%infos$sample[which(infos$origin=="bichat")]),]
tab_full_p<-tab_full[which(tab_full$Id_sample%in%infos$sample[which(infos$origin=="pitie")]),]

for(s in 1:nrow(NA_mut)){
  NA_mut$Nb[s]<-length(tab_full$mutation[which(tab_full$mutation==NA_mut$mutation[s]& (tab_full$obs=="discarded" | tab_full$obs=="discarded_minor"))])
  NA_mut$Nb_b[s]<-length(tab_full_b$mutation[which(tab_full_b$mutation==NA_mut$mutation[s]& (tab_full_b$obs=="discarded" | tab_full_b$obs=="discarded_minor"))])
  NA_mut$Nb_p[s]<-length(tab_full_p$mutation[which(tab_full_p$mutation==NA_mut$mutation[s]& (tab_full_p$obs=="discarded" | tab_full_p$obs=="discarded_minor"))])
}
NA_mut$pct_NA<-(NA_mut$Nb * 100) / length(unique(tab_full$Id_sample))
NA_mut$pct_NA_b<-(NA_mut$Nb_b * 100) / length(unique(tab_full_b$Id_sample))
NA_mut$pct_NA_p<-(NA_mut$Nb_p * 100) / length(unique(tab_full_p$Id_sample))

png(filename = "Missing_data_pers_mut_J0_seq1-2-3-4-5-6-7_pitie.png", width=8, height=8,units="in",res = 400)
ggplot(NA_mut, aes(x=pct_NA)) + 
  geom_histogram(binwidth=1)
dev.off()

## remove mutations with more than 50% bichat ou pitie missing 
mut_remove<-NA_mut$mutation[which(NA_mut$pct_NA_b>=50 | NA_mut$pct_NA_p>=50)]
tab_full<-tab_full[which((tab_full$mutation%in%mut_remove)==F),]
write.table(mut_remove,file="mut_more50pctNA_bichat_or_pitie.txt",col.names = F, row.names = F, quote = F)


###### ACP 

#setwd(paste(work_rep,"J0/",sep=""))
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
li$sample_type<-infos$individual[match(rownames(li),infos$sample)]


color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                      "#FC4E2A", "#E7298A", "#BD0026", "#800026",
                      "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5" ,"#08519C" ,"#08306B" ,
                      "#984EA3","grey")
                      
# 
# png(filename = "PCA_all_genes_J0_Axes12_seq1-2-3-4-5-6-7_pitie.png", width=8, height=8,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=variant , label= rownames(li), shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+ 
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())
# dev.off()
# 
# 
# png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6_pitie_J0_better_zoom.png", width=16, height=12,units="in",res = 400)
# ggplot(li, aes(x=PC1 , y=PC2, color=variant, shape=sample_type))+
#   geom_point(size=2, aes(shape=sample_type))+
#   scale_color_manual(values = color_set)+
#   scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
#   labs(color = "Variant", shape= element_blank())+
#   facet_zoom(xlim = c(-0.02, 0.02),ylim = c(-0.1, 0.12))
# dev.off()


png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_Nextclade_J0_better_zoom_publi.png", width=16, height=12,units="in",res = 400)
ggplot(li, aes(x=PC1 , y=PC2, color=variant_N, shape=sample_type))+
  geom_point(size=2, aes(shape=sample_type))+
  scale_color_manual(values = color_set)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())+
  facet_zoom(xlim = c(-0.02, 0.02),ylim = c(-0.1, 0.12))+
  theme_Publication()
dev.off()


# zoom on delta and omicron
li2<-li[which(li$PC1>-0.02 & li$PC1<0.02 & li$PC2>-0.1 & li$PC2<0.12 ),]

li2$origin<-infos$origin[match(rownames(li2),infos$sample)]

color_set2<-c( "#FFD92F", "#B3B3B3",
                      "#FC4E2A", "#E7298A", "#BD0026", "#800026",
                      "#C6DBEF", "#9ECAE1", "#4292C6", "#2171B5" ,"#08519C" ,"#08306B" ,
                      "#984EA3","grey")

#png(filename = "PCA_all_genes_Axes12_patients_personnels_J0_zoom.png", width=8, height=8,units="in",res = 400)
pdf(file = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_Nextclade_J0_zoom.pdf", width=12, height=12)
ggplot(li2, aes(x=PC1 , y=PC2, color=variant, shape=sample_type))+
  geom_point(size=3, aes(shape=sample_type))+
  scale_color_manual(values = color_set2)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())+
  theme_Publication()
dev.off()


### distinction site
li$origin<-infos$origin[match(rownames(li),infos$sample)]

png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_J0_par_site.png", width=8, height=8,units="in",res = 400)
ggplot(li, aes(x=PC1 , y=PC2, color=origin, shape=sample_type))+
  geom_point(size=2, aes(shape=sample_type))+
  scale_color_manual(values = color_set)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())
dev.off()


png(filename = "PCA_all_genes_Axes12_seq1-2-3-4-5-6-7_pitie_J0_par_site_better_zoom.png", width=16, height=12,units="in",res = 400)
ggplot(li, aes(x=PC1 , y=PC2, color=origin, shape=sample_type))+
  geom_point(size=2, aes(shape=sample_type))+
  scale_color_manual(values = color_set)+
  scale_shape_manual(values = c(1, 17),labels=c('Patients', 'Controls'))+
  labs(color = "Variant", shape= element_blank())+
  facet_zoom(xlim = c(-0.02, 0.02),ylim = c(-0.1, 0.12))
dev.off()


### look at omicron and delta outliers 
li_subset<-li[which(str_detect(li$variant,"Delta") | str_detect(li$variant,"Omicron")),]



##### Heatmap des variations

mut_tab_graph<-tab_full
mut_tab_graph$mutation_gene<-paste(mut_tab_graph$mutation,mut_tab_graph$gene,sep="_")
mut_tab_graph$gene_start<-genes$Start[match(mut_tab_graph$gene,genes$GeneId)]
mut_tab_graph$gene_order<-mut_tab_graph$gene_start * 10000
mut_tab_graph$order<-mut_tab_graph$gene_order + mut_tab_graph$AA_pos
mut_tab_graph$origin<-infos$origin[match(mut_tab_graph$Id_sample,infos$sample)]
mut_tab_graph$virus_nexclade<-infos$virus_nexclade[match(mut_tab_graph$Id_sample,infos$sample)]
mut_tab_graph$origin_virus<-paste(mut_tab_graph$origin,mut_tab_graph$virus_nexclade,sep="_")

mut_tab_graph$virus_sample<-paste(mut_tab_graph$virus_nexclade,mut_tab_graph$Id_sample,sep="_")

png(filename = paste("heatmap_freq_mut_seq1-2-3-4-5-6-7_pitie.png"), width=60, height=30,units="in",res = 400)
p1<-ggplot(data = mut_tab_graph, aes(x=fct_reorder(mutation_gene, order), y=as.factor(virus_sample), fill=freq)) + 
  geom_tile(color = "black") +
  ylab("jour")+
  scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray")+
  xlab("Mutation")+
  theme(axis.text.x = element_text(angle=90,size = 2),axis.text.y = element_text(size = 5))
print(p1)
dev.off()


### samples bichat
mut_tab_graph_b<-tab_full_b
mut_tab_graph_b$mutation_gene<-paste(mut_tab_graph_b$mutation,mut_tab_graph_b$gene,sep="_")
mut_tab_graph_b$gene_start<-genes$Start[match(mut_tab_graph_b$gene,genes$GeneId)]
mut_tab_graph_b$gene_order<-mut_tab_graph_b$gene_start * 10000
mut_tab_graph_b$order<-mut_tab_graph_b$gene_order + mut_tab_graph_b$AA_pos
mut_tab_graph_b$origin<-infos$origin[match(mut_tab_graph_b$Id_sample,infos$sample)]
mut_tab_graph_b$virus_nexclade<-infos$virus_nexclade[match(mut_tab_graph_b$Id_sample,infos$sample)]
mut_tab_graph_b$origin_virus<-paste(mut_tab_graph_b$origin,mut_tab_graph_b$virus_nexclade,sep="_")

mut_tab_graph_b$virus_sample<-paste(mut_tab_graph_b$virus_nexclade,mut_tab_graph_b$Id_sample,sep="_")

png(filename = paste("heatmap_freq_mut_seq1-2-3-4-5-6-7.png"), width=20, height=10,units="in",res = 400)
p1<-ggplot(data = mut_tab_graph_b, aes(x=fct_reorder(mutation_gene, order), y=as.factor(virus_sample), fill=freq)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray")+
  xlab("Mutation")+
  theme(axis.text.x = element_text(angle=90,size = 3),axis.text.y = element_text(size = 10),axis.title.y = element_blank(),axis.title.x = element_blank())
print(p1)
dev.off()




#### personels del pitié 

samples_del<-c("fastq_total_662112062454", "fastq_total_662111089143","fastq_total_662111074221", "fastq_total_662112007007", "fastq_total_662112007055","fastq_total_662112007078","fastq_total_662112022517")
mut_tab_graph_del<-mut_tab_graph[which(mut_tab_graph$Id_sample%in%samples_del),]

mut_del<-unique(mut_tab_graph_del$mutation_gene)

for(j in 1:length(mut_del)){
  tmp<-mut_tab_graph_del[which(mut_tab_graph_del$mutation_gene==mut_del[j]),]
  if (length(which(tmp$obs=="major" | tmp$obs=="minor"))>0){
    if(exists("tab_graph_del")==F){
      tab_graph_del<-tmp
    } else {
      tab_graph_del<-rbind(tab_graph_del,tmp)
    }
  }
}

png(filename = paste("heatmap_freq_mut_cluster_personnel_CF.png"), width=20, height=10,units="in",res = 400)
p1<-ggplot(data = tab_graph_del, aes(x=fct_reorder(mutation_gene, order), y=as.factor(virus_sample), fill=freq)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray")+
  xlab("Mutation")+
  theme(axis.text.x = element_text(angle=90,size = 8),axis.text.y = element_text(size = 10),axis.title.y = element_blank(),axis.title.x = element_blank())
print(p1)
dev.off()
