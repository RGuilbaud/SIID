############## Make fasta from mutation file (tsv.txt) ##################

library(readxl)

#### import gene table
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")
genes<-genes[order(genes$Start),]

#### import pitié sample list 
infos<-read_xlsx("emergen_romane_lineage.xlsx")
infos<-infos[which(infos$origin=="pitie"),]

#### work rep
data_rep<-("~/Projets/SIID/2_Filter_trim/Accurate/Seq_Pitie/")
setwd(data_rep)


#### multi fasta 
### write the fasta sample per sample 
for(i in 1:nrow(infos)){
  sample<-infos$sample[i]
  
  make_multifasta<-function(sample, i){
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
      write.table(entete, file=paste("Consensus_major_pitie.fasta"),col.names = F,row.names = F, quote = F, sep="")
      write.table(sequence_major, file=paste("Consensus_major_pitie.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
    } else { ## else : add the line at the end of the file
      write.table(entete, file=paste("Consensus_major_pitie.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
      write.table(sequence_major, file=paste("Consensus_major_pitie.fasta"),append = T,col.names = F,row.names = F, quote = F, sep="")
    }
  }
  
  make_multifasta(sample,i)
}



