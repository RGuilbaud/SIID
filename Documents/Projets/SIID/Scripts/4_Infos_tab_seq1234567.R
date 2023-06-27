############# Create info/metadata tables  
library(readxl)

### import existing metadata 
setwd("~/Projets/SIID/data")
infos<-read_xlsx("emergen_romane_lineage.xlsx")
metadata<-read_excel("Metadata_VF.xlsx",sheet="Metadata_copie")
metadata$DDN2<-as.Date(metadata$DDN, format =  "%d/%m/%Y")


############################# Bichat resequencing summary #####################################

### import the sample lists from the bichat runs 
setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
samples12<-read.csv(file="samplelist_seq1-2.csv",sep=";",header = F)
samples12$seq<-"Seq_1-2"
samples23<-read.csv(file="samplelist_seq2-3.csv",sep=";",header = F)
samples23$seq<-"Seq_2-3"
samples4<-read.csv(file="samplelist_seq4.csv",sep=";",header = F)
samples4$seq<-"Seq_4"
samples5<-read.csv(file="samplelist_seq5.csv",sep=";",header = F)
samples5$seq<-"Seq_5"
samples6<-read.csv(file="samplelist_seq6.csv",sep=";",header = F)
samples6$seq<-"Seq_6"
samples7<-read.csv(file="samplelist_seq7.csv",sep=";",header = F)
samples7$seq<-"Seq_7"

### import the variants obtained in our analysis 
nextclade<-read.csv(file = "nextclade_1-2-3-4-5-6-7.csv",sep = ";", header = T)

### create a table with all the samples and add infos (if reseq, variant, date, type, etc....)
samples<-rbind(samples12,rbind(samples23,rbind(samples4,rbind(samples5, rbind(samples6,samples7)))))
colnames(samples)<-c("Nom","barcode","Ct","num_barcode","seq")
samples$reseq<-"NON"
samples$reseq[which(duplicated(samples$Nom)==T | duplicated(samples$Nom, fromLast=TRUE)==T)]<-"OUI"
samples$variant_Nextstrain<-nextclade$clade[match(samples$Nom,nextclade$seqName)]
samples$variant_pango<-nextclade$Nextclade_pango[match(samples$Nom,nextclade$seqName)]

samples$date_prelev<-metadata$`DATE_PRELE`[match(samples$Nom,metadata$STAT_DOSSIER)]
samples$type<-infos$individual[match(samples$Nom,infos$sample)]

### gives the reason why I filtered these samples (the samples reseq because their Ct was a little high do not have any tag)
setwd("~/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
samples$Filtre<-"NON"
NA20<-read.table(file ="samples_more_20pctNA.txt",sep="\t" )
for (i in 1:nrow(samples)){
  
  if(samples$Nom[i]%in%NA20$V1){
    samples$Filtre[i]<-"Données manquantes sup 20%"
  }
  
  ## if there is names in the "Samples_med_cov_inf50.txt" file
  if (file.size(paste("~/Projets/SIID/2_Filter_trim/Accurate/",samples$seq[i],"/Samples_med_cov_inf50.txt",sep="")) > 0){
    med_inf50<-read.table(file = paste("~/Projets/SIID/2_Filter_trim/Accurate/",samples$seq[i],"/Samples_med_cov_inf50.txt",sep=""),sep="\t")
    options(scipen=999)
    if(samples$Nom[i]%in%med_inf50$V1){
      samples$Filtre[i]<-"Mediane inf 50X"
    }
  }

}

## samples done twice with good quality both times --> keep the best 
samples$Filtre[which(samples$Nom=="112206101397" & samples$seq=="Seq_6")]<-"Doublon non gardé"
samples$Filtre[which(samples$Nom=="112206103360" & samples$seq=="Seq_4")]<-"Doublon non gardé"

write.csv(samples,file="Infos_samples_seq1-2-3-4-5-6-7.csv",row.names = F)



#### keep only samples not reseq
infos_b<-infos[which(infos$origin=="bichat"),]
failed<-samples[which(samples$Filtre!="NON"),]
failed_nreseq<-failed[which((failed$Nom%in%samples$Nom[which(samples$Filtre=="NON")]==F)),]
infos_b_nreseq<-infos_b[which((infos_b$sample%in%samples$Nom)==F ),]



########################### Bichat Tableau envoi pitié ################################################
#### keep only the samples actually used in the analysis : remove the reseq samples
samples_analysis<-samples[which(samples$Filtre=="NON"),]
samples_reseq<-samples_analysis[which((samples_analysis$reseq=="OUI" & samples_analysis$seq=="Seq_5") | samples_analysis$Nom=="112206103360" | samples_analysis$Nom=="112206101397" ),]

samples_analysis<-samples_analysis[which(samples_analysis$reseq=="NON"),]
samples_analysis<-rbind(samples_analysis,samples_reseq)

write.csv(samples_analysis,"Infos_samples_seq1-2-3-4-5-6-7_analysis.csv",row.names = F)

#### add the metadata 
samples_analysis2<-samples_analysis[,c(1,3,7,8,9,10)]
samples_analysis2$sexe<-metadata$Sexe[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$Cardio<-metadata$Cardio[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$Greffe_rein<-metadata$`Greffé rein`[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$IS_med_int<-metadata$`IS Med Int`[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$CT_med_int<-metadata$`CT Med Int`[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$transplante_pneumo<-metadata$`Transplantés Pneumo`[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]
samples_analysis2$VIH<-metadata$VIH[match(samples_analysis2$Nom,metadata$STAT_DOSSIER)]

samples_analysis2$Cardio<-gsub("Néphro", "Nephro",samples_analysis2$Cardio)
samples_analysis2$Greffe_rein<-gsub("Néphro", "Nephro",samples_analysis2$Greffe_rein)

write.csv(samples_analysis2,file="Infos_samples_reseq_bichat.csv")



################################### Bichat + pitié: Create a common metadata table from the metadata table from pitié ##############
#### import metadata and nexclade output for pitié data
setwd("~/Projets/SIID/data")
metadata_p<-read_xlsx("Metadata_PSL_Ct rectifié.xlsx")
metadata_p$Dossier<-paste("fastq_total_",metadata_p$Dossier,sep="")
setwd("~/Projets/SIID/2_Filter_trim/Accurate/")
nextclade_p<-read.csv("nextclade_pitie.csv",sep=";")
metadata_p$Variant<-nextclade_p$clade_nextstrain[match(metadata_p$Dossier,nextclade_p$seqName)]
metadata_p<-metadata_p[,1:8]

pitie_abs<-nextclade_p[which((nextclade_p$seqName%in%metadata_p$Dossier)==F),]

#### format the bichat metadata to match pitié metadata format 
metadata_b<-metadata[which(metadata$STAT_DOSSIER%in%samples_analysis$Nom),]
metadata_b_format<-data.frame(Dossier=metadata_b$STAT_DOSSIER,DDP=metadata_b$DATE_PRELE,Categorie=NA,Sex=metadata_b$Sexe,DDN=metadata_b$DDN2,Categorie_ID=metadata_b$Type,Variant=NA,Ct=NA)
metadata_b_format$Categorie<-samples_analysis$type[match(metadata_b_format$Dossier,samples_analysis$Nom)]
metadata_b_format$Variant<-samples_analysis$variant_Nextstrain[match(metadata_b_format$Dossier,samples_analysis$Nom)]
metadata_b_format$Ct<-samples_analysis$Ct[match(metadata_b_format$Dossier,samples_analysis$Nom)]
metadata_b_format$Categorie_ID<-gsub("Néphro", "Nephro",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("CARDIO", "TRANSPL. CARD.",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("CT Med Int", "MED. INT.",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("IS Med Int", "MED. INT.",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("Nephro dialyse", "TRANSPL. REIN",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("Nephro", "TRANSPL. REIN",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("Ritux", "RITUX",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("VIH1", "VIH",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("VIH2", "VIH",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("\\+", "",metadata_b_format$Categorie_ID)
metadata_b_format$Categorie_ID<-gsub("VIH  VIH", "VIH",metadata_b_format$Categorie_ID)

#### add a column "origin" and merge the two tables 
metadata_b_format$Origin<-"Bichat"
metadata_p$Origin<-"Pitie"

metadata_new<-rbind(metadata_b_format,metadata_p)

##### New column for immunodepression 
metadata_new$Immunodepression<-NA
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="HEMATO")]<-"Hemato_Oncology"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="HEPATO-GREFFE")]<-"Liver_transplantation"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="MED. INT.")]<-"Autoimmune_or_inflamatory_diseases"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="ONCO")]<-"Oncology"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="Pneumo")]<-"Pneumology"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="RHUMATO")]<-"Autoimmune_or_inflamatory_diseases"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="RITUX")]<-"RITUXIMAB"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="TRANSPL. CARD.")]<-"Heart_transplantation"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="TRANSPL. CARD.VIH")]<-"Heart_transplantation"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="TRANSPL. REIN")]<-"Kidney_transplantation"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="TRANSPL. REIN / VIH")]<-"Kidney_transplantation"
metadata_new$Immunodepression[which(metadata_new$Categorie_ID=="VIH")]<-"HIV_infection"
metadata_new$Immunodepression[which(is.na(metadata_new$Categorie_ID) & metadata_new$Categorie=="patient_J0")]<-"Autoimmune_or_inflamatory_diseases"


#### write the full metadata table
setwd("~/Projets/SIID/data")
write.csv(metadata_new, file="Metadata_complet.csv",row.names = F)
