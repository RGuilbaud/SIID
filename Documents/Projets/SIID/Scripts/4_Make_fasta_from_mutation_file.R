############## Make fasta from mutation file (tsv.txt) ##################

library(readxl)

#### import gene table
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
genes<-genes[order(genes$Start),]

#### import the sample list from each seq
setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
infos1_2<-read.csv(file = "samplelist_seq1-2.csv",sep=";",header = F)
infos2_3<-read.csv(file = "samplelist_seq2-3.csv",sep=";",header = F)
infos4<-read.csv(file = "samplelist_seq4.csv",sep=";",header = F)
infos5<-read.csv(file = "samplelist_seq5.csv",sep=";",header = F)
infos6<-read.csv(file = "samplelist_seq6.csv",sep=";",header = F)
infos7<-read.csv(file = "samplelist_seq7.csv",sep=";",header = F)

#### Filter samples redone: remove the first seq that did not work (or the best if both worked)
infos1_2<-infos1_2[which((infos1_2$V1%in%infos5$V1)==F),]
infos2_3<-infos2_3[which((infos2_3$V1%in%infos5$V1)==F),]
infos4<-infos4[which((infos4$V1%in%infos5$V1)==F),]

infos4<-infos4[which((infos4$V1%in%infos6$V1)==F),]
infos6<-infos6[which((infos6$V1%in%infos2_3$V1)==F),]


#### combine all lists into one table and add the path to their directory
infos<-as.data.frame(c(infos1_2$V1,infos2_3$V1,infos4$V1, infos5$V1, infos6$V1, infos7$V1))
colnames(infos)<-"V1"
infos$V2<-NA
infos$V2[which(infos$V1%in%infos1_2$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_1-2/"
infos$V2[which(infos$V1%in%infos2_3$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_2-3/"
infos$V2[which(infos$V1%in%infos4$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_4/"
infos$V2[which(infos$V1%in%infos5$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_5/"
infos$V2[which(infos$V1%in%infos6$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_6/"
infos$V2[which(infos$V1%in%infos7$V1)]<-"C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Seq_7/"


#### make multi fasta 

### write the fasta sample per sample 
for(i in 1:nrow(infos)){
  
  sample<-infos$V1[i]
  data_rep<-infos$V2[i]
  
  make_multifasta<-function(sample,data_rep, i){
    ## Full genome sequence
    # write a full sequence of N
    seq_tab<-data.frame(pos=seq(genes$Start[1],genes$End[12],by=1),seq="N")
    
    # for each gene replace the N by the major sequence 
    for (j in 1:nrow(genes)){
      gene<-genes$GeneId[j]
      setwd(paste(data_rep,gene,sep=""))
      full_tab<-read.table(file= paste(sample,"_",gene,".tsv.txt",sep=""),sep="\t",header = T)
      seq_major<-as.data.frame(full_tab[,c(2,21)])
      
      seq_major$major_seq<-gsub("del","-",seq_major$major_seq)
      seq_tab$seq[which(seq_tab$pos>=genes$Start[j] & seq_tab$pos<=genes$End[j])]<-seq_major$major_seq[match(seq_tab$pos[which(seq_tab$pos>=genes$Start[j] & seq_tab$pos<=genes$End[j])],seq_major$Position)]
    }
    
    ## store the info that will go in the two lines of the fasta 
    entete<-paste(">",sample,sep="")
    sequence_major<-t(as.data.frame(seq_tab$seq))
    
    
    setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate")
    ## if it is the first line : create the file 
    if (i ==1 ){
      write.table(entete, file=paste("Consensus_major.fasta"),col.names = F,row.names = F, quote = F, sep="")
      write.table(sequence_major, file=paste("Consensus_major.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
    } else { ## else : add the line at the end of the file
      write.table(entete, file=paste("Consensus_major.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
      write.table(sequence_major, file=paste("Consensus_major.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
    }
  }
  
  make_multifasta(sample,data_rep,i)
}



