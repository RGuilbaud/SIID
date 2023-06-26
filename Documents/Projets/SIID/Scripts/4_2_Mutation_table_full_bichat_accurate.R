########### Matrice avec données manquantes quand non couvert ##########

library(dplyr)
library(readxl)
library(stringr)

### Files and repertories 
rep<-"~/Projets/SIID/2_Filter_trim"

setwd("~/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
infos<-read_xlsx("emergen_romane_lineage.xlsx")

setwd(rep)

### Function to create a mutations table including missing data and non mutations
mut_tab_with_no_cov<-function(pop){
  ## Choose the dataset to use
  if(pop == "J0"){
    setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
    #infos<-infos[which(infos$individual=="patient_J0" | infos$individual=="personnel"),]
    #infos<-infos[which(infos$origin=="bichat"),]
    tab<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt",sep="\t",header = T)
    
    samples12<-read.csv(file="samplelist_seq1-2.csv",sep=";",header = F)
    samples23<-read.csv(file="samplelist_seq2-3.csv",sep=";",header = F)
    samples4<-read.csv(file="samplelist_seq4.csv",sep=";",header = F)
    samples5<-read.csv(file="samplelist_seq5.csv",sep=";",header = F)
    samples6<-read.csv(file="samplelist_seq6.csv",sep=";",header = F)
    samples7<-read.csv(file="samplelist_seq7.csv",sep=";",header = F)
    
    # samples redone: remove the first seq that did not work (or the best if both worked)
    samples12<-samples12[which((samples12$V1%in%samples5$V1)==F),]
    samples23<-samples23[which((samples23$V1%in%samples5$V1)==F),]
    samples4<-samples4[which((samples4$V1%in%samples5$V1)==F),]
    
    samples4<-samples4[which((samples4$V1%in%samples6$V1)==F),]
    samples6<-samples6[which((samples6$V1%in%samples23$V1)==F),]
    
    
  } else if (pop == "Suivi"){
    setwd("~/Projets/SIID/2_Filter_trim/Accurate/")

    tab<-read.table(file="Mutations no_indel_Suivi_seq1-2-3-4-5-6-7_pitie .txt",sep="\t",header = T)
    
    samples12<-read.csv(file="samplelist_seq1-2.csv",sep=";",header = F)
    samples23<-read.csv(file="samplelist_seq2-3.csv",sep=";",header = F)
    samples4<-read.csv(file="samplelist_seq4.csv",sep=";",header = F)
    samples5<-read.csv(file="samplelist_seq5.csv",sep=";",header = F)
    samples6<-read.csv(file="samplelist_seq6.csv",sep=";",header = F)
    samples7<-read.csv(file="samplelist_seq7.csv",sep=";",header = F)
    
    # samples redone: remove the first seq that did not work (or the best if both worked)
    samples12<-samples12[which((samples12$V1%in%samples5$V1)==F),]
    samples23<-samples23[which((samples23$V1%in%samples5$V1)==F),]
    samples4<-samples4[which((samples4$V1%in%samples5$V1)==F),]
    
    samples4<-samples4[which((samples4$V1%in%samples6$V1)==F),]
    samples6<-samples6[which((samples6$V1%in%samples23$V1)==F),]
    
  }else if (pop == "J0_Suivi"){
    setwd("~/Projets/SIID/2_Filter_trim/Accurate/")

    tabJ0<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7 .txt",sep="\t",header = T)
    tabJS<-read.table(file="Mutations no_indel_Suivi_seq1-2-3-4-5-6-7 .txt",sep="\t",header = T)
    
    samples12<-read.csv(file="samplelist_seq1-2.csv",sep=";",header = F)
    samples23<-read.csv(file="samplelist_seq2-3.csv",sep=";",header = F)
    samples4<-read.csv(file="samplelist_seq4.csv",sep=";",header = F)
    samples5<-read.csv(file="samplelist_seq5.csv",sep=";",header = F)
    samples6<-read.csv(file="samplelist_seq6.csv",sep=";",header = F)
    samples7<-read.csv(file="samplelist_seq7.csv",sep=";",header = F)
    
    # samples redone: remove the first seq that did not work (or the best if both worked)
    samples12<-samples12[which((samples12$V1%in%samples5$V1)==F),]
    samples23<-samples23[which((samples23$V1%in%samples5$V1)==F),]
    samples4<-samples4[which((samples4$V1%in%samples5$V1)==F),]
    
    samples4<-samples4[which((samples4$V1%in%samples6$V1)==F),]
    samples6<-samples6[which((samples6$V1%in%samples23$V1)==F),]
    
    tab<-rbind(tabJ0, tabJS)
  }
  
  infos<-infos[which(infos$sample%in%tab$Id_sample),]
  
  ## for each gene and each sample
  for (i in 1:nrow(genes)){
    list_mut<-unique(tab$mutation[which(tab$gene==genes$GeneId[i])])
    for(j in 1:nrow(infos)){
      ## table containing all the mutations detected for the panel 
      mut_tab<-data.frame(Id_sample=infos$sample[j],depth=NA,mutation=list_mut,obs=NA,freq=NA,gene=genes$GeneId[i])
      mut_tab$AA_pos<-as.numeric(str_sub(mut_tab$mutation,2, -2))
      
      ## mutations carried by the sample
      tmp_tab<-tab[which(tab$Id_sample==infos$sample[j] & tab$gene==genes$GeneId[i]),]
      tmp_tab$AA_pos<-as.numeric(str_sub(tmp_tab$mutation,2, -2))
      
      ## read the sample tsv.txt file 
      if(infos$origin[j]=="bichat" & infos$sample[j]%in%samples12$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_1-2/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="bichat" & infos$sample[j]%in%samples23$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_2-3/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="bichat" & infos$sample[j]%in%samples4$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_4/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="bichat" & infos$sample[j]%in%samples5$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_5/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="bichat" & infos$sample[j]%in%samples6$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_6/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="bichat" & infos$sample[j]%in%samples7$V1){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_7/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      } else if (infos$origin[j]=="pitie"){
        dataset<-readr::read_tsv(paste(rep,"/Accurate/Seq_Pitie/",genes$GeneId[i],"/",infos$sample[j],"_",genes$GeneId[i],".tsv.txt",sep=""), col_names = T, na=character())
      }
       
      # For each mutation of the panel     
       for( k in 1:nrow(mut_tab)){
         # if this mutation is prensent in the sample --> write the infos from tmp_tab 
        if(mut_tab$mutation[k]%in%tmp_tab$mutation){
          mut_tab[k,]<-tmp_tab[which(tmp_tab$mutation==mut_tab$mutation[k]),]
          # else read in the tsv.txt file the depth and the reason it is not a mutation ("discarded", "no_mutation")
        } else {
          mut<-dataset[which(dataset$pos_spike_AA==mut_tab$AA_pos[k]),]
          mut_tab$depth[k]<-mut$depth_total[3]
          mut_tab$obs[k]<-mut$mutation[3]
          mut_tab$freq[k]<-0 # set the frequency at 0 (the cases of too low coverage will be converted into NA at the end of the loop) 
        }
       }
      if (j==1 & i==1){
        mut_tab_full<-mut_tab
      } else {
        mut_tab_full<-rbind(mut_tab_full,mut_tab)
        }
     }
  }
  ## set the frequency of low coverage mutations as NA
  mut_tab_full$freq[which(mut_tab_full$obs=="discarded" | mut_tab_full$obs=="discarded_minor")]<-NA
  ## write the table
  setwd(paste(rep,"/Accurate/",sep=""))
  write.table(mut_tab_full,file=paste("Mutations_no_indel_",pop,"_full_seq_1-2-3-4-5-6-7_avec_pitie2.txt",sep=""),sep="\t",row.names = F,quote = F)
}

mut_tab_with_no_cov("J0")
mut_tab_with_no_cov("Suivi")
mut_tab_with_no_cov("J0_Suivi")







