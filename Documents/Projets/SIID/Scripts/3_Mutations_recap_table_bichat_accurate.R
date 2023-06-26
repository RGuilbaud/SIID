#### Recap mutations
library(dplyr)
library(readxl)

make_recap<-function(input_name, output_name, pop){
  Start_time<-Sys.time()
  
  
  ############# CHOOSE THE SEQUENCING RUN ##############
  SeqRun<-"Seq_7"
  ######################################################
  
  
  # Working directory
  rep<-paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/",sep="")
  setwd("~/Projets/SIID/data")
  
  
  # Files to import
  genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
  infos<-read_xlsx("emergen_romane_lineage.xlsx")
  
  setwd(rep)
  
  if(pop == "J0"){
    infos<-infos[which(infos$individual=="patient_J0" | infos$individual=="personnel"),]
  } else if (pop == "Suivi"){
    infos<-infos[which(infos$individual=="patient_suivi"),]
  }
  
  # ##### ONLY TEST
  # rep<-"C:/Users/romain.coppee/Documents/test_romane/test"
  # setwd(rep)
  # #####
  
  for (i in 1:nrow(genes)){
    print(genes$GeneId[i])
    # gene name
    gene_name<-genes$GeneId[i]
    # for each gene read the "results" file if it exists
    if (file.exists(paste(rep,"/",gene_name,"/",input_name,".txt",sep=""))){
      dataset<-read.table(file = paste(gene_name,"/",input_name,".txt",sep=""),sep="\t",header = T)
      #dataset<-dataset[which(dataset$Id_sample%in%infos$sample),]
      #filenames = dir(gene_name, pattern="*.tsv.txt")
      ## list of the samples tsv.txt files and names
      filenames<-paste(dataset$Id_sample,"_",gene_name,".tsv.txt",sep="")
      nbsamples<-length(filenames)
      sample_name = gsub(paste0("\\_",gene_name,".*"), "", filenames)
      
      mut_tab<-data.frame(Id_sample=sample_name,nb_major=NA,nb_minor=NA,nb_total=NA,gene_cov=NA)
      
      # for each sample, if the sample has mutations for this gene, count the major and minor mutations
      # if none write a line of zeros
      for (ind in mut_tab$Id_sample){
        if (ind %in% dataset$Id_sample){
          mut_tab$nb_major[mut_tab$Id_sample==ind]<-length(which(dataset[dataset$Id_sample==ind,]=="major"))
          mut_tab$nb_minor[mut_tab$Id_sample==ind]<-length(which(dataset[dataset$Id_sample==ind,]=="minor"))
          mut_tab$nb_total[mut_tab$Id_sample==ind]<- mut_tab$nb_major[mut_tab$Id_sample==ind] + mut_tab$nb_minor[mut_tab$Id_sample==ind]
        } else {
          mut_tab$nb_major[mut_tab$Id_sample==ind]<-0
          mut_tab$nb_minor[mut_tab$Id_sample==ind]<-0
          mut_tab$nb_total[mut_tab$Id_sample==ind]<-0
        }
        
      }
      ## if there is no "results" file for the gene, write lines for the samples with 0 mutations 
    } else {
      setwd(rep)
      filenames = dir(gene_name, pattern="*.tsv.txt")
      nbsamples<-length(filenames)
      sample_name = gsub(paste0("\\_",gene_name,".*"), "", filenames)
      mut_tab<-data.frame(Id_sample=sample_name,nb_major=0,nb_minor=0,nb_total=0)
    }
    
    ## For ecah sample read the tsv.txt file 
    for(my_file in filenames){
      tmp_data<-readr::read_tsv(paste(rep,"/",gene_name,"/",my_file,sep=""), col_names = T, na=character())
      my_file_name<-gsub(paste0("\\_",gene_name,".*"), "", my_file)
      # calc gene coverage = nb codons with enough coverage / total nb of codons
      nb_codons<- nrow(tmp_data)/3
      nb_too_low<-length(which(tmp_data$obs=="too_low"))
      pct_cov<-((nb_codons - nb_too_low) * 100) / nb_codons
      mut_tab$gene_cov[which(mut_tab$Id_sample==my_file_name)]<-pct_cov
      ## write table with mutants from each sample
      mutants<-tmp_data[which(tmp_data$mutation%in%colnames(dataset)),]
      if(nrow(mutants!=0)){
        mutants<-select(mutants,Id_sample,depth,mutation,obs,freq)
        mutants$gene<-gene_name
        
        if (my_file==filenames[1] & gene_name==genes$GeneId[1]){
          mutants_tab<-mutants
        } else {
          mutants_tab<-rbind(mutants_tab,mutants)
        }
      }
    }
    

    
    mut_tab$Gene<-gene_name
    if (i == 1){
      recap<-mut_tab
    } else {
      recap<-rbind(recap,mut_tab)
    }
   
  }
  
  
  recap<-recap[which(recap$Id_sample%in%infos$sample),]
  mutants_tab<-mutants_tab[which(mutants_tab$Id_sample%in%infos$sample),]
  
  write.table(recap,file=paste("Recap_mutations",output_name,".txt"),sep="\t",row.names = F,quote=F)
  write.table(mutants_tab,file = paste("Mutations",output_name,".txt" ),sep="\t",row.names = F,quote=F)
  
  print(Sys.time() - Start_time)
}

#make_recap("results","all")
#make_recap("results_no_indel","no_indel","all")
make_recap("results_no_indel","no_indel_J0","J0")
make_recap("results_no_indel","no_indel_Suivi","Suivi")
make_recap("results","J0","J0")



#### Merge outputs from different seq + merge with pitie 
### JO

## import the bichat data from the different runs
setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_1-2/")
tab_B12<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B12<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_2-3/")
tab_B23<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B23<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_4/")
tab_B4<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B4<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_5/")
tab_B5<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B5<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_6/")
tab_B6<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B6<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_7/")
tab_B7<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_B7<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)

## import the list of pitié samples
setwd("~/Projets/SIID/data")
infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos<-infos[which(infos$origin=="pitie"),]

## import the data from pitié + bichat_fast and extract only the pitié samples 
setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/")
tab_P<-read.table(file="Mutations no_indel_J0 .txt",sep="\t",header = T)
tab_P<-tab_P[which(tab_P$Id_sample%in%infos$sample),]
rec_tab_P<-read.table(file="Recap_mutations no_indel_J0 .txt",sep="\t",header = T)
rec_tab_P<-rec_tab_P[which(rec_tab_P$Id_sample%in%infos$sample),]

## Mutations no_indel_J0 table
# Elim the doubles of the samples reseq +  Merge the bichat runs
tab_B4<-tab_B4[which((tab_B4$Id_sample%in%tab_B6$Id_sample)==F),]
tab_B6<-tab_B6[which((tab_B6$Id_sample%in%tab_B23$Id_sample)==F),]
tab_Btmp<-rbind(tab_B12,rbind(tab_B23, tab_B4))
tab_Btmp<-tab_Btmp[which((tab_Btmp$Id_sample%in%tab_B5$Id_sample)==F),]
tab_B<-rbind(tab_Btmp, rbind(tab_B5,rbind(tab_B6,tab_B7)))

# merge pitié and bichat
tab<-rbind(tab_B,tab_P)


## Recap_mutations no_indel_J0 table
# Elim the doubles of the samples reseq +  Merge the bichat runs
rec_tab_B4<-rec_tab_B4[which((rec_tab_B4$Id_sample%in%rec_tab_B6$Id_sample)==F),]
rec_tab_B6<-rec_tab_B6[which((rec_tab_B6$Id_sample%in%rec_tab_B23$Id_sample)==F),]
rec_tab_Btmp<-rbind(rec_tab_B12,rbind(rec_tab_B23, rec_tab_B4))
rec_tab_Btmp<-rec_tab_Btmp[which((rec_tab_Btmp$Id_sample%in%rec_tab_B5$Id_sample)==F),]
rec_tab_B<-rbind(rec_tab_Btmp, rbind(rec_tab_B5, rbind(rec_tab_B6,rec_tab_B7)))

# merge pitié and bichat
rec_tab<-rbind(rec_tab_B,rec_tab_P)


## write the new tables
# bichat all seq
setwd("~/Projets/SIID/2_Filter_trim/Accurate")
write.table(rec_tab_B,file=paste("Recap_mutations no_indel_J0_seq1-2-3-4-5-6-7 .txt"),sep="\t",row.names = F,quote=F)
write.table(tab_B,file = paste("Mutations no_indel_J0_seq1-2-3-4-5-6-7 .txt" ),sep="\t",row.names = F,quote=F)

# bichat + pitié
write.table(rec_tab,file=paste("Recap_mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt"),sep="\t",row.names = F,quote=F)
write.table(tab,file = paste("Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt" ),sep="\t",row.names = F,quote=F)




### Suivis (Bichat)
## import the bichat data from the different runs
setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_1-2/")
tab_B12<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B12<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_2-3/")
tab_B23<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B23<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_4/")
tab_B4<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B4<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_5/")
tab_B5<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B5<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_6/")
tab_B6<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B6<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)

setwd("~/Projets/SIID/2_Filter_trim/Accurate/Seq_7/")
tab_B7<-read.table(file="Mutations no_indel_Suivi .txt",sep="\t",header = T)
rec_tab_B7<-read.table(file="Recap_mutations no_indel_Suivi .txt",sep="\t",header = T)


setwd("~/Projets/SIID/data")
infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos<-infos[which(infos$origin=="pitie"),]

## Mutations no_indel_J0 table
# Elim the doubles of the samples reseq +  Merge the bichat runs
tab_B4<-tab_B4[which((tab_B4$Id_sample%in%tab_B6$Id_sample)==F),]
tab_B6<-tab_B6[which((tab_B6$Id_sample%in%tab_B23$Id_sample)==F),]
tab_Btmp<-rbind(tab_B12,rbind(tab_B23, tab_B4))
tab_Btmp<-tab_Btmp[which((tab_Btmp$Id_sample%in%tab_B5$Id_sample)==F),]
tab_B<-rbind(tab_Btmp, rbind(tab_B5,rbind(tab_B6,tab_B7)))
##


## Recap_mutations no_indel_J0 table
# Elim the doubles of the samples reseq +  Merge the bichat runs
rec_tab_B4<-rec_tab_B4[which((rec_tab_B4$Id_sample%in%rec_tab_B6$Id_sample)==F),]
rec_tab_B6<-rec_tab_B6[which((rec_tab_B6$Id_sample%in%rec_tab_B23$Id_sample)==F),]
rec_tab_Btmp<-rbind(rec_tab_B12,rbind(rec_tab_B23, rec_tab_B4))
rec_tab_Btmp<-rec_tab_Btmp[which((rec_tab_Btmp$Id_sample%in%rec_tab_B5$Id_sample)==F),]
rec_tab_B<-rbind(rec_tab_Btmp, rbind(rec_tab_B5, rbind(rec_tab_B6,rec_tab_B7)))
##

## write the new tables
# bichat all seq
setwd("~/Projets/SIID/2_Filter_trim/Accurate")
write.table(rec_tab_B,file=paste("Recap_mutations no_indel_Suivi_seq1-2-3-4-5-6-7 .txt"),sep="\t",row.names = F,quote=F)
write.table(tab_B,file = paste("Mutations no_indel_Suivi_seq1-2-3-4-5-6-7 .txt" ),sep="\t",row.names = F,quote=F)


