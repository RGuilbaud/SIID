Start_time<-Sys.time()
library(readr)
library(data.table)
library(reshape2)
library(Biostrings)
library(stringr)
library(dplyr)

setwd("~/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
setwd("~/Projets/SIID/Scripts/")
source(file = "Fonction_code_genetique.R")


############# CHOOSE THE SEQUENCING RUN ##############
SeqRun<-"Seq_7"
######################################################



for (gene in 1:nrow(genes)) {
  gene_name<-genes$GeneId[gene]
  print(gene_name)
  # gene start position
  start<-genes$Start[gene]
  # gene end position 
  end<-genes$End[gene]
  
  # work repertory
  setwd(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/",gene_name,sep=""))
  
  #PARAMETERS
  min_coverage_major = 20
  min_coverage_minor = 50
  min_freq_alt = 5
  max_freq_alt = 49.99
  min_freq_indel = 30
  start_spike = start
  end_spike = end
  
  # files to use
  ref_seq <- read.table(paste("~/Projets/SIID/Seq_ref/",genes$Seq_name[gene],sep=""), header = TRUE, sep = "\t")
  filenames = dir(pattern="*.tsv")
  
  # # OPTINAL! if a subset is needed 
  # samplelist<-read.csv(file = "~/Projets/SIID/2_Filter_trim/Accurate/samplelist_seq5_redo.csv",sep=";",header = F)
  # samplelist$filename<-paste(samplelist$V1,"_",gene_name,".tsv",sep="")
  # filenames<-filenames[which(filenames%in%samplelist$filename)]
  
  # Fonction n'est pas inclus/pr?sent dans 
  '%!in%'<- function(x,y)!('%in%'(x,y))
  
  
  # Importation et formatage fichier d'entr?e
  for (my_file in filenames){
    print(my_file)
    #add lacking rows
    file_name = as.numeric(sub("\\_.*", "", my_file))
    my_sample <- read_tsv(my_file, col_names = TRUE)
    for(i in start_spike:end_spike){
      if(i %!in% my_sample$Position){
        my_sample[nrow(my_sample) + 1,] = list(file_name, i, 0, 0, 0, 0, 0, "N", "N", NA)
      }
    }
    my_sample = my_sample[order(my_sample$Position),]
    
    #fixing zero problem
    for(i in 1:nrow(my_sample)){
      if(my_sample$n_A[i] == 0 && my_sample$n_T[i] == 0 && my_sample$n_G[i] == 0 && my_sample$n_C[i] == 0){
        my_sample$n_A[i] = 1
      }
    }
    
    
    my_sample[my_sample=="0"]<-1
    my_sample$complement <- NULL
    my_sample$qual <- NULL
    
    #table combination seq ref et sample
    full_data <- cbind(my_sample, ref_seq)
    full_data$depth_total <- rowSums(full_data[3:7])
    full_data$major_depth <- apply(full_data[3:7], 1, FUN=max)
    
    # calc des % ? chaque pos
    for (i in 1:nrow(full_data)){
      #print(i)
      full_data$freq_A[i] <- full_data$n_A[i]/full_data$depth_total[i]*100
      full_data$freq_T[i] <- full_data$n_T[i]/full_data$depth_total[i]*100
      full_data$freq_G[i] <- full_data$n_G[i]/full_data$depth_total[i]*100
      full_data$freq_C[i] <- full_data$n_C[i]/full_data$depth_total[i]*100
      full_data$freq_indel[i] <- full_data$n_indel[i]/full_data$depth_total[i]*100
      full_data$freq_alt[i] <- max(full_data[i,14:18][full_data[i,14:18] != max(full_data[i,14:18])])
      full_data$freq_maj[i] <- max(full_data[i,14:18])
      
      # def all?le major
      if(full_data$depth_total[i] >= min_coverage_major){
        if(full_data$freq_A[i] == full_data$freq_maj[i]){
          full_data$major_seq[i] = "A"
        } else if(full_data$freq_T[i] == full_data$freq_maj[i]){
          full_data$major_seq[i] = "T"
        } else if(full_data$freq_G[i] == full_data$freq_maj[i]){
          full_data$major_seq[i] = "G"
        } else if(full_data$freq_C[i] == full_data$freq_maj[i]){
          full_data$major_seq[i] = "C"
        } else {
          full_data$major_seq[i] = "del"
        }
      } else {
        full_data$major_seq[i] = "N"
      }
      
      # cas allele alt est indel
      if(full_data$freq_alt[i] == full_data$freq_indel[i]){
        # Si freq indel sous seuil
        if(full_data$freq_alt[i] < min_freq_indel && full_data$depth_total[i] >= min_coverage_major){
          # Si on a une freq alt autre que l'indel au dessus du seuil
          if (max(full_data[i,14:17][full_data[i,14:17] != full_data$freq_maj[i]])>= min_freq_alt && full_data$depth_total[i] >= min_coverage_minor){
            if(full_data$freq_A[i] >= min_freq_alt && full_data$freq_A[i] < max_freq_alt){
              full_data$minor_seq[i] = "A"
              full_data$freq_alt[i]<-full_data$freq_A[i]
            } else if (full_data$freq_T[i] >= min_freq_alt && full_data$freq_T[i] < max_freq_alt){
              full_data$minor_seq[i] = "T"
              full_data$freq_alt[i]<-full_data$freq_T[i]
            } else if (full_data$freq_G[i] >= min_freq_alt && full_data$freq_G[i] < max_freq_alt){
              full_data$minor_seq[i] = "G"
              full_data$freq_alt[i]<-full_data$freq_G[i]
            } else if (full_data$freq_C[i] >= min_freq_alt && full_data$freq_C[i] < max_freq_alt){
              full_data$minor_seq[i] = "C"
              full_data$freq_alt[i]<-full_data$freq_C[i]
            } 
          } else {
            # Si indel et autre allèle tous sous le seuil
            full_data$minor_seq[i] = full_data$major_seq[i]
          }
          # Si freq allèle au dessus du seuil
        } else if (full_data$freq_indel[i] >= min_freq_indel && full_data$freq_indel[i] < max_freq_alt){
          full_data$minor_seq[i] = "del"
        } else {
          full_data$minor_seq[i] = "N"
        }
      } else { # cas allele alt non indel 
        if(full_data$freq_alt[i] < min_freq_alt && full_data$depth_total[i] >= min_coverage_major){
          full_data$minor_seq[i] = full_data$major_seq[i]
        } else if (full_data$depth_total[i] >= min_coverage_minor && full_data$freq_alt[i] >= min_freq_alt){
          if(full_data$freq_A[i] >= min_freq_alt && full_data$freq_A[i] < max_freq_alt){
            full_data$minor_seq[i] = "A"
          } else if (full_data$freq_T[i] >= min_freq_alt && full_data$freq_T[i] < max_freq_alt){
            full_data$minor_seq[i] = "T"
          } else if (full_data$freq_G[i] >= min_freq_alt && full_data$freq_G[i] < max_freq_alt){
            full_data$minor_seq[i] = "G"
          } else if (full_data$freq_C[i] >= min_freq_alt && full_data$freq_C[i] < max_freq_alt){
            full_data$minor_seq[i] = "C"
          } 
        } else {
          full_data$minor_seq[i] = "N"
        }
      }
      
      
      ### AUTOM  
      codon_check<-function(it,data_seq){
        if(it%%3 == 0){
          if(data_seq[it-2] == "N" || data_seq[it-1] == "N" || data_seq[it] == "N"){
            output = "too_low"
            return(output)
          } else if (data_seq[it-2] == "del" || data_seq[it-1] == "del" || data_seq[it] == "del"){
            output = "del"
            return(output)
          } else {
            output = "pass"
            return(output)
          }
        } else {
          output = NA
          return(output)
        }
      }  
      
      full_data$N_major[i]<-codon_check(i,full_data$major_seq)
      full_data$N_minor[i]<-codon_check(i,full_data$minor_seq)
      ###
      
      
      #CODE GENETIQUE
      
      
      full_data<-code_genet(full_data,i)
      
      # ### AUTOM
      # if(i%%3 == 0){
      #   if (full_data$N_minor[i] == "del"){
      #     full_data$minor[i] = "del"
      #   } else if (full_data$N_minor[i] == "too_low"){
      #     full_data$minor[i] = "too_low"
      #   } else {
      #     s_minor <- paste(full_data$minor_seq[(i-2):i], collapse="")
      #     s_minor = DNAString(s_minor)
      #     full_data$minor[i] = as.character(translate(s_minor, if.fuzzy.codon="X", no.init.codon = T))
      #   }
      #   
      #   if (full_data$N_major[i] == "del"){
      #     full_data$major[i] = "del"
      #   } else if (full_data$N_major[i] == "too_low"){
      #     full_data$major[i] = "too_low"
      #   } else {
      #     s_major <- paste(full_data$major_seq[(i-2):i], collapse="")
      #     s_major = DNAString(s_major)
      #     full_data$major[i] = as.character(translate(s_major, if.fuzzy.codon="X", no.init.codon = T))
      #   }
      # }
      # ##### 
      
      
      if(i%%3 == 0){
        if(full_data$N_minor[i] == "pass"){
          full_data$final_minor[i] = full_data$minor[i]
        } else if (full_data$N_minor[i] == "del") {
          full_data$final_minor[i] = "del"
        } else {
          full_data$final_minor[i] = "too_low"
        }
        
        if(full_data$N_major[i] == "pass"){
          full_data$final_major[i] = full_data$major[i]
        } else if (full_data$N_major[i] == "del") {
          full_data$final_major[i] = "del"
        } else {
          full_data$final_major[i] = "too_low"
        }
      }
      

      
      
      # Ecrit de quel type de mutation ou non mutation il s'agit
      if(i%%3 == 0){
        if(full_data$final_major[i] == "too_low"){
          full_data$mutation[i] = "discarded"
        } else if(full_data$final_major[i] != full_data$AA_ref[i]) {
          full_data$mutation[i] = gsub(" ", "", paste(full_data$AA_ref[i],full_data$pos_spike_AA[i],full_data$final_major[i]))
        } else if(full_data$final_minor[i] == "too_low"){
          full_data$mutation[i] = "discarded_minor"
        } else if(full_data$final_minor[i] != full_data$AA_ref[i]) {
          full_data$mutation[i] = gsub(" ", "", paste(full_data$AA_ref[i],full_data$pos_spike_AA[i],full_data$final_minor[i]))
        } else {
          full_data$mutation[i] = "no_mutation"
        }
      }
      
      
      
      if(i%%3 == 0){
        if(full_data$mutation[i] != "no_mutation"){
          if(full_data$mutation[i] == "discarded"){
            full_data$obs[i] = "too_low"
          } else if(full_data$mutation[i] == "discarded_minor") {
            full_data$obs[i] = "too_low_minor"
          #} else if (full_data$major_seq[i] != full_data$nucl_ref[i] || full_data$major_seq[i-1] != full_data$nucl_ref[i-1] || full_data$major_seq[i-2] != full_data$nucl_ref[i-2]){
          } else if (full_data$final_major[i] != full_data$AA_ref[i]){
            full_data$obs[i] = "major"
          } else {
            full_data$obs[i] = "minor"
          }
        } else {
          full_data$obs[i] = "-"
        }
      }
      
      # ?crit les freq major et minor
      if(i%%3 == 0){
        if(full_data$obs[i] == "minor"){
          full_data$freq[i] = full_data[(i-2):i,] %>% filter(minor_seq != nucl_ref) %>% select(freq_alt) %>% max()
          #full_data$freq[i] = max(full_data$freq_alt[i-2], full_data$freq_alt[i-1], full_data$freq_alt[i])
        } else if(full_data$obs[i] == "major") {
          full_data$freq[i] = full_data[(i-2):i,] %>% filter(major_seq != nucl_ref) %>% select(freq_maj) %>% max()
          #full_data$freq[i] = max(full_data$freq_maj[i-2], full_data$freq_maj[i-1], full_data$freq_maj[i])
        } else if(full_data$obs[i] == "too_low_minor"){
          full_data$freq[i] = "too_low_minor"
        } else if(full_data$obs[i] == "too_low"){
          full_data$freq[i] = "too_low"
        } else {
          full_data$freq[i] = NA
        }
      }
      
      # on ?crit la couverture du codon 
      if(i%%3 == 0){
        if(full_data$obs[i] == "minor" || full_data$obs[i] == "major"){
          if(full_data$major_seq[i-2] != full_data$nucl_ref[i-2] || full_data$minor_seq[i-2] != full_data$nucl_ref[i-2]){
            full_data$depth[i] = full_data$depth_total[i-2]
          } else if (full_data$major_seq[i-1] != full_data$nucl_ref[i-1] || full_data$minor_seq[i-1] != full_data$nucl_ref[i-1]){
            full_data$depth[i] = full_data$depth_total[i-1]
          } else if (full_data$major_seq[i] != full_data$nucl_ref[i] || full_data$minor_seq[i] != full_data$nucl_ref[i]){
            full_data$depth[i] = full_data$depth_total[i]
          }
          
        } else {
          full_data$depth[i] = NA
        }
      } else {
        full_data$depth[i] = NA
      }
      
      
    }
    write.table(x = full_data, file = paste(my_file, ".txt", sep=""), sep="\t", row.names=F,quote = F)
    
  }
  
  
  # Cr?e matrice des mutations
  temp = list.files(pattern="*.tsv.txt")
  myfiles = lapply(temp, read.delim)
  full_data_combine = rbindlist(myfiles, fill=FALSE, idcol=NULL)
  
  AA_data = full_data_combine[complete.cases(full_data_combine), ]
  if (nrow(AA_data)!=0){
    AA_data_clean = AA_data[,c(1,30,31,32,27)]
    matrix_mutation = dcast(AA_data_clean, Id_sample ~ mutation, value.var="obs")
    
    write.table(matrix_mutation, file="results.txt", col.names = TRUE, sep="\t", row.names = FALSE)
    
    # matrix mutations sans les indels
    AA_data_clean_NOindel<-AA_data_clean[which(str_detect(AA_data_clean$mutation,"del")==F),]
    AA_data_clean_NOindel<-AA_data_clean_NOindel[which(str_detect(AA_data_clean_NOindel$mutation,"ins")==F),]
    
    if (nrow(AA_data_clean_NOindel)!=0){
      matrix_mutation_NOindel = dcast(AA_data_clean_NOindel, Id_sample ~ mutation, value.var="obs")
      write.table(matrix_mutation_NOindel, file="results_no_indel.txt", col.names = TRUE, sep="\t", row.names = FALSE)
    }
  }


}

print(Sys.time() - Start_time)
