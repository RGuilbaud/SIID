#### Tsv files for each ORF ###

# Library import
library(readr)


# Files to import
setwd("~/Projets/SIID/data")
genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")


############# CHOOSE THE SEQUENCING RUN ##############
SeqRun<-"Seq_7"
######################################################



setwd(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/tsv_files/",sep=""))
tsv_files<-dir(pattern="*.tsv")

# # OPTINAL! if a subset is needed 
# samplelist<-read.csv(file = "~/Projets/SIID/2_Filter_trim/Accurate/samplelist_seq5_redo.csv",sep=";",header = F)
# samplelist$filename<-paste(samplelist$V1,".tsv",sep="")
# tsv_files<-tsv_files[which(tsv_files%in%samplelist$filename)]

### remove the samples without a median coverage of 50
if(file.size(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/Samples_med_cov_inf50.txt",sep="")) > 0){
  setwd(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/",sep=""))
  filtered_samples<-read.table(file ="Samples_med_cov_inf50.txt",sep="\t",header = F)
  
  tsv_files<-tsv_files[which((tsv_files%in%filtered_samples[,1])==F)]
}


# for each sample, we write a new file containing the postions inside the gene (= between the start and the end positions)
for (sample_file in tsv_files){
  setwd(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/tsv_files/",sep=""))
  data_name<-gsub(".tsv","",sample_file)
  data_table<-readr::read_tsv(sample_file, col_names = T)
  
  for (i in 1:nrow(genes)){
    print(genes$GeneId[i])
    # gene name
    gene_name<-genes$GeneId[i]
    # gene start position
    start<-genes$Start[i]
    # gene end position 
    end<-genes$End[i]
    
    # create new dir for the gene (if it doesn't already exist)
    setwd(paste("~/Projets/SIID/2_Filter_trim/Accurate/",SeqRun,"/",sep=""))
    dir.create(gene_name)

      
    subdata_table<-data_table [which(data_table$Position>=start & data_table$Position<=end),]
    write.table(subdata_table,file = paste(gene_name,"/",data_name,"_",gene_name,".tsv",sep=""), row.names = F, quote = F,sep="\t")
  }

}