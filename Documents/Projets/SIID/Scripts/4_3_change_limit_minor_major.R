#################### Change the threshold for minor/major mutants from 50% to 20% ######################
##### turns the minor variants with freq > 20% into major variants 
library(ggplot2)

### Work rep
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/")

### load mutation table, mutation recap table and mutation full table (= has all the mutations for each sample)
tab<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt",sep="\t",header = T)
recap<-read.table(file = "Recap_mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt",sep="\t",header = T)
tab_full<-read.table(file = "Mutations_no_indel_J0_full_seq_1-2-3-4-5-6-7_avec_pitie2.txt",header=T, sep="\t")

### Change the minors > 20% into majors 
tab$obs[which(tab$freq>=20 & tab$obs=="minor")]<-"major"

ggplot(tab, aes(x=gene , y=freq, fill=obs))+
  geom_boxplot() +
  xlab("Genes") + ylab("frequence") 

tab_full$obs[which(tab_full$freq>=20 & tab_full$obs=="minor")]<-"major"


### Update the number of minor and major mutations count
recap$nb_major<-NA
recap$nb_minor<-NA
recap$nb_total<-NA

for(i in 1:nrow(recap)){
  tmp_tab<-tab[which(tab$Id_sample==recap$Id_sample[i] & tab$gene==recap$Gene[i]),]
  recap$nb_major[i]<-length(which(tmp_tab$obs=="major"))
  recap$nb_minor[i]<-length(which(tmp_tab$obs=="minor"))
  recap$nb_total[i]<-recap$nb_major[i] + recap$nb_minor[i]
}


### write the new files
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")

write.table(tab,file = "Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie_min5-20pct .txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(tab_full,file = "Mutations_no_indel_J0_full_seq_1-2-3-4-5-6-7_avec_pitie2_min5-20pct.txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(recap,file = "Recap_mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie_min5-20pct .txt",sep="\t",col.names = T,row.names = F,quote = F)




##################### just bichat #########################
### Work rep
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/")

### load mutation table, mutation recap table and mutation J0 and suivi table 
tab<-read.table(file="Mutations no_indel_J0_seq1-2-3-4 .txt",sep="\t",header = T)
recap<-read.table(file = "Recap_mutations no_indel_J0_seq1-2-3-4 .txt",sep="\t",header = T)
tab_full_J0_suivi<-read.table(file = "Mutations_no_indel_J0_Suivi_full_seq_1-2-3-4-5-6-7_avec_pitie2.txt",header=T, sep="\t")

### Change the minors > 20% into majors 
tab$obs[which(tab$freq>=20 & tab$obs=="minor")]<-"major"

ggplot(tab, aes(x=gene , y=freq, fill=obs))+
  geom_boxplot() +
  xlab("Genes") + ylab("frequence") 

tab_full_J0_suivi$obs[which(tab_full_J0_suivi$freq>=20 & tab_full_J0_suivi$obs=="minor")]<-"major"


### Update the number of minor and major mutations count
recap$nb_major<-NA
recap$nb_minor<-NA
recap$nb_total<-NA

for(i in 1:nrow(recap)){
  tmp_tab<-tab[which(tab$Id_sample==recap$Id_sample[i] & tab$gene==recap$Gene[i]),]
  recap$nb_major[i]<-length(which(tmp_tab$obs=="major"))
  recap$nb_minor[i]<-length(which(tmp_tab$obs=="minor"))
  recap$nb_total[i]<-recap$nb_major[i] + recap$nb_minor[i]
}


### write the new files
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")

write.table(tab,file = "Mutations no_indel_J0_seq1-2-3-4_min5-20pct .txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(recap,file = "Recap_mutations no_indel_J0_seq1-2-3-4_min5-20pct .txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(tab_full_J0_suivi,file = "Mutations_no_indel_J0_Suivi_full_seq_1-2-3-4-5-6-7_min5-20pct.txt",sep="\t",col.names = T,row.names = F,quote = F)
