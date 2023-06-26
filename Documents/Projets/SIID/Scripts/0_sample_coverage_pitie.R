########## Sample coverage 
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(paletteer)
library(grid)
library(gtable)
library(readxl)

### Files to use
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/tsv_files/")
my_files<-dir(pattern="*.tsv", recursive = T)
my_files<-my_files[grep("fastq_total",my_files)]

setwd("C:/Users/Romane/Documents/Projets/SIID/data")
samples<-read_xlsx("Metadata_pitie.xlsx",sheet="Metadata")
samples$Dossier<-paste("fastq_total_",samples$Dossier,sep="")


## calc coverage as sum of allelic counts (ATCG + indel)
for( i in 1:length(my_files)){
  file_name = sub("\\.tsv", "", my_files[i])
  if (file_name%in%samples$Dossier){
    setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/tsv_files/")
    tab<-read_tsv(my_files[i], col_names = TRUE)
    if(nrow(tab)>0){
      tab$cov<-apply(tab[3:7],1,sum)
      tab$Sample<-file_name
      
      
      if(i==1){
        cov_tab<-tab
      } else {
        cov_tab<-rbind(cov_tab,tab)
      }
    }
  }
}


### plot average coverage
mean_cov<-cov_tab %>% group_by(Sample) %>% summarise(mean(cov))
colnames(mean_cov)<-c("sample","cov_moy")
ggplot(mean_cov, aes(x=cov_moy)) + 
  geom_histogram(binwidth = 1)+
  scale_x_continuous(breaks = seq(0,max(mean_cov$cov_moy),by=10))


### plot median coverage
median_cov<-cov_tab %>% group_by(Sample) %>% summarise(median(cov))
colnames(median_cov)<-c("sample","cov_med")
summary(median_cov$cov_med)
samples$med_cov<-median_cov$cov_med[match(samples$Dossier,median_cov$sample)]

Pcovmed<-ggplot(data=median_cov, aes(x=sample, y=cov_med, group=1)) +
  geom_line()+
  geom_point()+
  scale_y_continuous(breaks = seq(0,max(median_cov$cov_med)+ 100,by=100))+
  theme(axis.text.x = element_text(angle=90))



####### Number of reads per sample


setwd(paste("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/"))
nb_reads_tmp<-read.table(file = "nb_reads.txt",sep="\t",header = T)
nb_reads_tmp<-nb_reads_tmp[which(nb_reads_tmp$Sample%in%samples$Dossier[which(samples$seq==seqs[s])]),]
nb_reads_tmp$Nb_raw_reads<-nb_reads_tmp$Nb_raw_reads/4
nb_reads_tmp$nb_reads_filtered<-nb_reads_tmp$nb_reads_filtered/4
  

samples$nb_reads<-nb_reads$Nb_raw_reads[match(samples$Dossier,nb_reads$Sample)]
summary(samples$nb_reads)

graph_tab<-pivot_longer(nb_reads,cols = 2:3, names_to = "Filter",values_to = "nb_reads")
# png(filename = "Number_reads_before_after_filters.png", width=14, height=8,units="in",res = 400)
# ggplot(data=graph_tab, aes(x=as.factor(Sample), y=nb_reads, fill=Filter)) +
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_y_continuous(breaks = seq(0,max(graph_tab$nb_reads)+ 500,by=10000))+
#   theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank(),legend.position = "top")
# dev.off()

png(filename = "Number_raw_reads.png", width=14, height=8,units="in",res = 400)
ggplot(data=nb_reads, aes(x=as.factor(Sample), y=Nb_raw_reads)) +
  geom_bar(stat="identity")+
  scale_y_continuous(breaks = seq(0,max(nb_reads$Nb_raw_reads)+ 500,by=10000))+
  theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank())
dev.off()




###### global graph

samples <- samples[order(samples[,1]),]

Pct<-ggplot(data=samples, aes(x=as.factor(Dossier), y=Ct, group=1)) +
  geom_line()+
  geom_point()+
  scale_y_continuous(breaks = seq(0,max(samples$Ct, na.rm = T)+ 10,by=1))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())


Preads<-ggplot(data=samples, aes(x=as.factor(Dossier), y=nb_reads)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_y_continuous(breaks = seq(0,max(samples$nb_reads)+ 500,by=50000))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "top")


Pcovmed<-ggplot(data=samples, aes(x=as.factor(Dossier), y=med_cov, group=1)) +
  geom_line()+
  geom_point()+
  scale_y_continuous(breaks = seq(0,max(samples$med_cov)+ 100,by=200))+
  theme(axis.text.x = element_text(angle=90, size=5), axis.title.x = element_blank())

summary(samples$Ct)
summary(samples$nb_reads)
summary(samples$med_cov)

setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Bichat_accurate")
png(filename = "Recap_Ct_nbReads_Cov.png", width=20, height=8,units="in",res = 400)
g2 <- ggplotGrob(Pct)
g3 <- ggplotGrob(Preads)
g4 <- ggplotGrob(Pcovmed)
g <- rbind(g2, g3, g4, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths, g4$widths)
grid.newpage()
grid.draw(g)
dev.off()                                       