####### Patients suivis

library(readxl)
library(dplyr)
library(ggplot2)
library(paletteer) 
library(gridExtra)
library(ggforce)
library(ggrepel)
library(forcats)

## Ggplot2 theme importation
setwd("~/Projets/SIID/Scripts/")
source(file = "Function_theme_Publication.R")

## infos on sars cov 2 genes
setwd("C:/Users/Romane/Documents/Projets/SIID/data")

genes<-read.table(file = "Genes_SARS-CoV2.txt",header=T,sep="\t")

## samples file (bichat)
suivi<-read_xlsx("Metadata_VF.xlsx",sheet="Main_update")
colnames(suivi)<- c(colnames(suivi[,1:3]),"Date_heure",colnames(suivi[,5:ncol(suivi)]))

suivi$jour<-as.numeric(gsub("J","",suivi$J_x))

suivi$nom_prenom<-paste(suivi$Nom,suivi$Prénom,sep="_")
list_noms<-unique(suivi$nom_prenom)

suivi$num_patient<-NA

## creation of a code for each patient
for (i in 1:length(list_noms)){
  num_patient<-paste("Patient",i,sep = "")
  suivi$num_patient[which(suivi$nom_prenom==list_noms[i])]<-num_patient
}

tab_suivi<-select(suivi, num_dossier, num_patient, Date_heure, jour,Clade, Pango)
rm(suivi)
rm(list_noms)


## metadata 
metadata<-read_xlsx("Metadata_VF.xlsx",sheet="Metadata_copie")
#metadata<-metadata[which(metadata$STAT_DOSSIER%in%infos$dossier),]


## infos on samples (from the analysis)
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/")
infos<-read.csv(file = "Infos_samples_seq1-2-3-4-5-6-7.csv")
# elim doublons
infos$supr<-"NON"
infos$supr[which(infos$Filtre!="NON" | (infos$reseq=="OUI" & (infos$seq=="Seq_1-2" | infos$seq=="Seq_2-3" | infos$seq=="Seq_4")))]<-"OUI"
infos<-infos[which(infos$supr=="NON"),]

metadata<-metadata[which(metadata$STAT_DOSSIER%in%infos$Nom),]

tab_suivi<-tab_suivi[which(tab_suivi$num_dossier%in%infos$Nom),]
tab_suivi$type_prelev<-infos$type[match(tab_suivi$num_dossier,infos$Nom)]
tab_suivi$Clade<-infos$variant_Nextstrain[match(tab_suivi$num_dossier,infos$Nom)]
tab_suivi$Pango<-infos$variant_pango[match(tab_suivi$num_dossier,infos$Nom)]

list_patients<-unique(tab_suivi$num_patient)

## keep only suivi
for (i in 1:length(list_patients)){
  tmp<-tab_suivi[which(tab_suivi$num_patient==list_patients[i]),]
  if(nrow(tmp)>1){
    if(length(which(tmp$type_prelev=="patient_suivi"))!=0){
      if (exists("tab_suivi2")==F){
        tab_suivi2<-tmp
      } else {
        tab_suivi2<-rbind(tab_suivi2,tmp)
      }
    }
  }
}


## nbe jours suivi
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Suivi/")

png(filename = "Nombre_jours_suivi.png",width = 800,height = 600)
ggplot(tab_suivi2,aes(x=jour))+
  geom_bar(stat = "bin", binwidth = 1)+
  scale_x_continuous(breaks = seq(0,max(tab_suivi$jour) + 5, by = 10))+
  xlab("Jours de prelevement") + ylab("nombre de patients")
dev.off()

median(tab_suivi2$jour[which(tab_suivi2$jour!=0)])
summary(tab_suivi2$jour[which(tab_suivi2$jour!=0)])


### nbe prelevements (nombre de fois ou un num patient est present)
nb_prelev<-tab_suivi2 %>% dplyr :: count(num_patient)

png(filename = "Nombre_prelevement_par_patient.png",width = 800,height = 600)
ggplot(nb_prelev,aes(x=n))+
  geom_bar(stat = "bin", binwidth = 1)+
  scale_x_continuous(breaks = seq(0,max(nb_prelev$n) + 1, by = 1))+
  xlab("nombre de prelevements") + ylab("nombre de patients")
dev.off()  



## mutations table
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
mut_J0_suivi<-read.table(file="Mutations_no_indel_J0_Suivi_full_seq_1-2-3-4-5-6-7_min5-20pct.txt",sep = "\t", header = T)
#mut_suivi<-read.table(file="Mutations_no_indel_suivi_full.txt",sep = "\t", header = T)
mut_J0_suivi<-mut_J0_suivi[which(mut_J0_suivi$Id_sample%in%tab_suivi2$num_dossier),]
#mut_suivi<-mut_suivi[which(mut_suivi$Id_sample%in%tab_suivi2$num_dossier),]

list_patients2<-unique(tab_suivi2$num_patient)

mut_J0_suivi$variant<-infos$variant_Nextstrain[match(mut_J0_suivi$Id_sample,infos$Nom)]


## import mutations specific of each variant (extracted from covariants)
setwd("C:/Users/Romane/Documents/Projets/SIID/data")
mut_20I<-read.table(file="Mut_20I_Alpha_V1.txt",sep=":",header = F)
mut_21J<-read.table(file="Mut_21J_Delta.txt",sep=":",header = F)
mut_21K<-read.table(file="Mut_21K_Omicron.txt",sep=":",header = F)
mut_21L<-read.table(file="Mut_21L_Omicron.txt",sep=":",header = F)
mut_22B<-read.table(file="Mut_22B_Omicron.txt",sep=":",header = F)
mut_22E<-read.table(file="Mut_22E_Omicron.txt",sep=":",header = F)
mut_22F<-read.table(file="Mut_22F_Omicron.txt",sep=":",header = F)



### identify the mutations characteristic of variants in our data
mut_J0_suivi$mut_spe<-NA
for (n in 1:nrow(mut_J0_suivi)){
  if(mut_J0_suivi$mutation[n] %in% mut_20I$V2[which(mut_20I$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"20I_Alpha_V1"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "20I_Alpha_V1", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_21J$V2[which(mut_21J$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"21J_Delta"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "21J_Delta", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_21K$V2[which(mut_21K$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"21K_Omicron"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "21K_Omicron", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_21L$V2[which(mut_21L$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"21L_Omicron"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "21L_Omicron", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_22B$V2[which(mut_22B$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"22B_Omicron"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "22B_Omicron", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_22E$V2[which(mut_22E$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"22E_Omicron"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "22E_Omicron", sep=";")
    }
  }
  if(mut_J0_suivi$mutation[n] %in% mut_22F$V2[which(mut_22F$V1==mut_J0_suivi$gene[n])]){
    if(is.na(mut_J0_suivi$mut_spe[n])){
      mut_J0_suivi$mut_spe[n]<-"22F_Omicron"
    } else {
      mut_J0_suivi$mut_spe[n]<-paste(mut_J0_suivi$mut_spe[n], "22F_Omicron", sep=";")
    }
  }
}


#################################################################################################
###### Table/matrix of mutation frequency on all suivis from a patient ######


### Filter the patients : keep only patients with at least two days "suivi"
pls_suivis<-nb_prelev[which(nb_prelev$n>2),]
tab_suivi_filter<-tab_suivi2[which(tab_suivi2$num_patient%in%pls_suivis$num_patient),]

list_patients_filter<-unique(tab_suivi_filter$num_patient)


for (patient in list_patients_filter){
  tmp_tab_filter<-tab_suivi_filter[which(tab_suivi_filter$num_patient==patient),]

  # extract J0 (or the lowest suivi if no J0)
  #tmp_tab_J0_filter<-tmp_tab_filter[which(tmp_tab_filter$type_prelev=="patient_J0"),]
  tmp_tab_J0_filter<-tmp_tab_filter[which(tmp_tab_filter$jour==min(tmp_tab_filter$jour)),]
  J0_filter<-mut_J0_suivi[which(mut_J0_suivi$Id_sample == tmp_tab_J0_filter$num_dossier),]
  J0_filter$mutation_gene<-paste(J0_filter$mutation,J0_filter$gene,sep="_")
  
  #J0_filter$jour<-0
  J0_filter$jour<-min(tmp_tab_filter$jour)
  
  mut_tab<-J0_filter
  
  # extract suivi
  tmp_tab_suivi_filter<-tmp_tab_filter[which(tmp_tab_filter$type_prelev=="patient_suivi" & tmp_tab_filter$jour!=min(tmp_tab_filter$jour)),]
  tmp_tab_suivi_filter<-tmp_tab_suivi_filter[order(tmp_tab_suivi_filter$jour),]
  for(j in 1:nrow(tmp_tab_suivi_filter)){
    #check variant covid --> if the varaint is not the same than J0 the suivi is skipped
    if(tmp_tab_suivi_filter$Clade[j]==tmp_tab_J0_filter$Clade[1]){
      
      suivi_filter<-mut_J0_suivi[which(mut_J0_suivi$Id_sample == tmp_tab_suivi_filter$num_dossier[j]),]
      suivi_filter$Id_sample<-as.character(suivi_filter$Id_sample)
      suivi_filter$mutation_gene<-paste(suivi_filter$mutation,suivi_filter$gene,sep="_")
      
      if((tmp_tab_suivi_filter$jour[j] %in%mut_tab$jour)==F){
        suivi_filter$jour<-tmp_tab_suivi_filter$jour[j]
      } else {
        # if two samples from the same partient the same day, the second is called day.2
        suivi_filter$jour<-paste(tmp_tab_suivi_filter$jour[j],".2",sep="")
      }
      
      mut_tab<-rbind(mut_tab,suivi_filter)
    }
  }
  mut_tab$jour<-as.numeric(mut_tab$jour)
  mut_patient<-unique(mut_tab$mutation_gene[which(mut_tab$obs=="major" | mut_tab$obs=="minor")])
  mut_tab<-mut_tab[which(mut_tab$mutation_gene%in%mut_patient),]
  mut_tab_full<-mut_tab
  mut_tab<-mut_tab[which(mut_tab$freq!=0 | is.na(mut_tab$freq)==T),]
  
  ### mut that appear in suivi (only if J0 available)
  J0_mut<-mut_tab$mutation_gene[which(mut_tab$jour==0)]
  if(length(J0_mut>0)){
    suivi_mut<-mut_tab$mutation_gene[which(mut_tab$jour!=0)]
    new_mut<-suivi_mut[which((suivi_mut%in%J0_mut)==F)]
    if(patient==list_patients_filter[1]){
      new_mut_list<-new_mut
    } else {
      new_mut_list<-c(new_mut_list, new_mut)
    }
  }

  ## type of sample info : Ecouvillon nasal or 	L. broncho-alvéolaire
  tmp_tab_filter$prelevement<-metadata$VARIABLE[match(tmp_tab_filter$num_dossier,metadata$STAT_DOSSIER)]
  tmp_tab_filter<-tmp_tab_filter[order(tmp_tab_filter$jour),]
  prelev <- ifelse(tmp_tab_filter$prelevement == "Ecouvillon nasal", "red", "blue")
  
  list_mut<-unique(mut_tab$mutation_gene)

  if(length(unique(mut_tab$jour))>1){

    ## graphs of freq by day for each sample
    setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Suivi/")
    
    mut_tab_graph<-mut_tab
    mut_tab_graph$gene_start<-genes$Start[match(mut_tab_graph$gene,genes$GeneId)]
    mut_tab_graph$gene_order<-mut_tab_graph$gene_start * 10000
    mut_tab_graph$order<-mut_tab_graph$gene_order + mut_tab_graph$AA_pos
    
   # png(filename = paste(patient,"_heatmap_freq_mut_suivi_publi.png"), width=17, height=7,units="in",res = 200)
    pdf(file = paste(patient,"_heatmap_freq_mut_suivi_publi.pdf"), width=17, height=7)
    p1<-ggplot(data = mut_tab_graph, aes(x=fct_reorder(mutation, order), y=as.factor(jour), fill=freq)) +
      geom_tile() +
      ylab("Day")+
      scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray")+
      xlab("Mutation")+
      theme_Publication()+
      labs(fill='Frequency')+
      facet_grid(. ~ reorder(gsub("ORF","",gene),gene_start),scales="free",space = "free")+
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90,size = 10),axis.text.y = element_text(colour = prelev),legend.key.size= unit(0.5, "cm"))
    print(p1)
    dev.off()

    png(filename = paste(patient,"_lineplot_freq_mut_suivi.png"), width=10, height=10,units="in",res = 200)
    p2<-ggplot(mut_tab_full, aes(x=as.factor(jour), y=freq, group=mutation_gene)) +
      geom_line(aes( color=mutation_gene),size=0.6)+
      geom_point(aes(color=mutation_gene))+
      theme(legend.position='none', panel.background = element_rect(fill = "darkgray"),axis.text.x = element_text(colour = prelev))+
      xlab("jour")+
      scale_color_paletteer_d("palettesForR::Named")
    print(p2)
    dev.off()

    # Zoom mutants freq<50
    mut_tab_low_freq<-mut_tab[which(mut_tab$freq<50),]

    png(filename = paste(patient,"_lineplot_LOW_freq_mut_suivi_NOzero.png"), width=10, height=10,units="in",res = 200)
    p3<-ggplot(mut_tab_low_freq, aes(x=as.factor(jour), y=freq, group=mutation_gene)) +
      geom_line(aes( color=mutation_gene),size=0.6)+
      geom_point(aes(color=mutation_gene))+
      theme(legend.position='none', panel.background = element_rect(fill = "darkgray"),axis.text.x = element_text(colour = prelev))+
      xlab("jour")+
      scale_color_paletteer_d("palettesForR::Named")
    print(p3)
    dev.off()


    # Zoom mutants freq<50
    mut_tab_low_freq_full<-mut_tab_full[which(mut_tab_full$freq<50),]

    png(filename = paste(patient,"_lineplot_LOW_freq_mut_suivi.png"), width=10, height=10,units="in",res = 200)
    p4<-ggplot(mut_tab_low_freq_full, aes(x=as.factor(jour), y=freq, group=mutation_gene)) +
      geom_line(aes( color=mutation_gene),size=0.6)+
      geom_point(aes(color=mutation_gene))+
      theme(legend.position='none', panel.background = element_rect(fill = "darkgray"),axis.text.x = element_text(colour = prelev))+
      xlab("jour")+
      scale_color_paletteer_d("palettesForR::Named")
    print(p4)
    dev.off()


    ## recup cas J0= minor ou NA et on a un major dans les suivis
    for (l in 1:length(list_mut)){
        mutation<-mut_tab_full[which(mut_tab_full$mutation_gene%in%list_mut[l]),]
        mutation<-mutation[order(mutation$jour),]
        #if(mutation$obs[which(mutation$jour==0)]=="minor" | is.na(mutation$freq[which(mutation$jour==0)])==T | mutation$freq[which(mutation$jour==0)]==0){
        if(mutation$obs[which(mutation$jour==min(tmp_tab_filter$jour))]=="minor" | is.na(mutation$freq[which(mutation$jour==min(tmp_tab_filter$jour))])==T | mutation$freq[which(mutation$jour==min(tmp_tab_filter$jour))]==0){
          if(length(which(mutation$obs=="major"))>=1){
            if(exists("min_to_maj")==F){
              min_to_maj<-mutation
            } else {
              min_to_maj<-rbind(min_to_maj,mutation )
            }

          }
        }
      }

    if(exists("min_to_maj")==T){
      png(filename = paste(patient,"_lineplot_freq_mut_min_to_maj_suivi.png"), width=10, height=10,units="in",res = 200)
      p5<-ggplot(min_to_maj, aes(x=as.factor(jour), y=freq, group=mutation_gene)) +
        geom_line(aes( color=mutation_gene),size=0.6)+
        geom_point(aes(color=mutation_gene))+
        theme(legend.position='none', panel.background = element_rect(fill = "darkgray"),axis.text.x = element_text(colour = prelev))+
        xlab("jour")+
        scale_color_paletteer_d("ggsci::default_igv")
      print(p5)
      dev.off()

      min_to_maj$patient<-patient
      if(exists("min_to_maj_tab")==F){
        min_to_maj_tab<-min_to_maj
      } else {
        min_to_maj_tab<-rbind(min_to_maj_tab,min_to_maj)
      }

      rm(min_to_maj)
    }

    
    #### New mut in suivi
    if(length(J0_mut>0)){
    mut_tab_full_new<-mut_tab_full[which(mut_tab_full$mutation_gene%in%new_mut),]
    color_set<-c("#edd229", "#cb279e", "#3f942c", "#dc8382", "#e73c1e", "#80a626", "#792611", "#8848e6", "#377b9b", "#afc3e7", "#3babc4", "#52d8b4", "#ec87b0", "#beef7d", "#649ad5", "#76262a", "#8dc35b", "#629f53", "#d8e6c3", "#ccb875", "#eae1ad", "#809fb3", "#9ab59b", "#54563d", "#506591", "#536bbb", "#578761", "#462419", "#8de990", "#b1c434", "#bf9a4b", "#e36a55", "#7deb36", "#a53c9c", "#b8e4e3", "#704213", "#415418", "#a2a94e", "#9e5eec", "#e438d8", "#9993e1", "#c64d96", "#9a7475", "#5aeb75", "#724b3c", "#aa5344", "#e6ea5d", "#ddbe37", "#3f1a39", "#cab49f", "#75768e", "#df3891", "#67497f", "#f0b233", "#e05e93", "#6b23b2", "#4691eb", "#6a2859", "#22433d", "#a7491a", "#8fe05b", "#e2ccdb", "#eb3077", "#af6595", "#a72628", "#e06329", "#e23a43", "#5c172d", "#96296e", "#b3a025", "#6bbd82", "#937226", "#e45ae2", "#9a71a0", "#1d233b", "#35732b", "#d29c7a", "#df75cc", "#336aaa", "#8b72e2", "#c66ee0", "#355f61", "#d8a0aa", "#dae588", "#3d4daf", "#a62eb1", "#263118", "#c1b5ea", "#439ba2", "#ef4bc2", "#59e547", "#da3d60", "#f128a4", "#8cbbb8", "#6f7a21", "#2f2059", "#351b70", "#9b2edb", "#30417c", "#a62754", "#51827c", "#5ee7f2", "#a1e231", "#c944e6", "#a6b782", "#cdf060", "#87d0f0", "#e9d76d", "#808451", "#f0b66e", "#7285ec", "#c99ec0", "#c17225", "#efcd9a", "#3c1a8f", "#39a289", "#0e2420", "#67521f", "#8885b9", "#eb8d1e", "#8d5dad", "#0d353f", "#55f4df", "#d39231", "#6a1c6c", "#433cd4", "#61efaf", "#4a4ac9", "#53b0e0", "#e38c5f", "#e1ec2c", "#342d18", "#dd9ee3", "#877c64", "#3fc068", "#a9edca", "#384c68", "#4fcdcf", "#9e5166", "#372c37", "#14391a", "#536ff5", "#9f6943", "#c0ea9c", "#2f5d3a", "#674a5b", "#72c3a7", "#d45f6d", "#47b52f", "#73349b")
    
    png(filename = paste(patient,"_lineplot_freq_New_mut_suivi.png"), width=10, height=10,units="in",res = 200)
    p6<-ggplot(mut_tab_full_new, aes(x=as.factor(jour), y=freq, group=mutation_gene)) +
      geom_line(aes( color=mutation_gene),size=0.6)+
      geom_point(aes(color=mutation_gene))+
      theme(legend.position='none',axis.text.x = element_text(colour = prelev))+
      xlab("jour")+
      #scale_color_paletteer_d("palettesForR::Named")
      scale_color_manual(values=color_set)
    print(p6)
    dev.off()
    }
  }
  
  
  ### distance matrix between the suivis
  mut_tab_full$jour_sample<-paste("J",mut_tab_full$jour,"_",mut_tab_full$Id_sample,sep="")
  
  matrix_mut<-as.data.frame(matrix(0,length(unique(mut_tab_full$jour_sample)),length(unique(mut_tab_full$mutation))))
  colnames(matrix_mut)<-unique(mut_tab_full$mutation)
  rownames(matrix_mut)<-unique(mut_tab_full$jour_sample)
  
  
  for(i in 1:ncol(matrix_mut)){
    tmp<-mut_tab_full[which(mut_tab_full$mutation==colnames(matrix_mut[i])),]
    tmp_maj<-tmp[which(tmp$obs=="major"),]
    tmp_min<-tmp[which(tmp$obs=="minor"),]
    tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
    matrix_mut[which(rownames(matrix_mut)%in%tmp_maj$jour_sample),i]<-1
    matrix_mut[which(rownames(matrix_mut)%in%tmp_min$jour_sample),i]<-0.5
    matrix_mut[which(rownames(matrix_mut)%in%tmp_nocov$jour_sample),i]<- NA
    
  }
  
  rm(tmp)
  rm(tmp_maj)
  rm(tmp_min)
  rm(tmp_nocov)
  
  # library(adegenet)
  # library(dartR)
  library(proxy)
  
  # matrix_mut2<-df2genind(matrix_mut,ncode=1) %>% gi2gl()
  # 
  # 
  # dist_mat<-as.matrix(dist(matrix_mut,method = "minkowski"))
  # dist_mat<-as.matrix(gl.dist.ind(matrix_mut2,method="Absolute"))
  
  simil_mat<-as.matrix(simil(matrix_mut))
  write.csv(simil_mat,file=paste(patient,"taux_simil.csv",sep="_"))
}


###############################################################################################################################
### New mut in suivis
new_mut_list<-unique(new_mut_list)
new_mut_tab<-data.frame(mut=new_mut_list,nb_J0=NA,nb_suivi=NA,nb_patients_J0=NA,nb_patients_suivi=NA)

# add infos to mut_J0_suivi
mut_J0_suivi<-mut_J0_suivi[which(mut_J0_suivi$Id_sample%in%tab_suivi_filter$num_dossier),]
mut_J0_suivi$patient<-tab_suivi_filter$num_patient[match(mut_J0_suivi$Id_sample,tab_suivi_filter$num_dossier)]
mut_J0_suivi$jour<-tab_suivi_filter$jour[match(mut_J0_suivi$Id_sample,tab_suivi_filter$num_dossier)]
mut_J0_suivi$mut_gene<-paste(mut_J0_suivi$mutation,mut_J0_suivi$gene,sep="_")

#extract major and minor (+ the others in an other table)
mut_J0_suivi_maj_min<-mut_J0_suivi[which(mut_J0_suivi$obs=="major" | mut_J0_suivi$obs=="minor"),]
mut_J0_suivi_not_maj_min<-mut_J0_suivi[which(mut_J0_suivi$obs!="major" & mut_J0_suivi$obs!="minor"),]

## count for each mutation the number of samples J0 and suivi that have it + nb of patients they represent 
for (m in 1:nrow(new_mut_tab)){
  new_mut_tab$nb_J0[m]<-length(unique(mut_J0_suivi_maj_min$Id_sample[which(mut_J0_suivi_maj_min$jour==0 & mut_J0_suivi_maj_min$mut_gene==new_mut_tab$mut[m])]))
  new_mut_tab$nb_suivi[m]<-length(unique(mut_J0_suivi_maj_min$Id_sample[which(mut_J0_suivi_maj_min$jour!=0 & mut_J0_suivi_maj_min$mut_gene==new_mut_tab$mut[m])]))
  new_mut_tab$nb_patients_J0[m]<-length(unique(mut_J0_suivi_maj_min$patient[which(mut_J0_suivi_maj_min$jour==0 & mut_J0_suivi_maj_min$mut_gene==new_mut_tab$mut[m])]))
  new_mut_tab$nb_patients_suivi[m]<-length(unique(mut_J0_suivi_maj_min$patient[which(mut_J0_suivi_maj_min$jour!=0 & mut_J0_suivi_maj_min$mut_gene==new_mut_tab$mut[m])]))
  }

## extract the mutations carried only by at least 2 suivis 
mut_J0_suivi_new<-mut_J0_suivi_maj_min[which(mut_J0_suivi_maj_min$mut_gene%in%new_mut_tab$mut[which(new_mut_tab$nb_J0==0 & new_mut_tab$nb_patients_suivi>1)]),]

## add the non major or minor lines for those mutations and samples
mut_J0_suivi_not_maj_min<-mut_J0_suivi_not_maj_min[which(mut_J0_suivi_not_maj_min$mut_gene%in%mut_J0_suivi_new$mut_gene & mut_J0_suivi_not_maj_min$Id_sample%in%mut_J0_suivi_new$Id_sample),]
mut_J0_suivi_new<-rbind(mut_J0_suivi_new,mut_J0_suivi_not_maj_min)

## graph
mut_J0_suivi_new$pos_gene<-genes$Start[match(mut_J0_suivi_new$gene, genes$GeneId)]
mut_J0_suivi_new$gene_order<-mut_J0_suivi_new$pos_gene * 10000
mut_J0_suivi_new$order<-mut_J0_suivi_new$gene_order + mut_J0_suivi_new$AA_pos

png(filename = "mutations_nouvelles_dans_suivis.png", width=15, height=10,units="in",res = 200)
ggplot(mut_J0_suivi_new, aes(x=as.factor(reorder(mut_gene,order)), y=freq, group=patient)) +
  geom_line(aes(color=patient))+
  geom_point(aes(color=patient))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))+
  scale_color_manual(values=c("#9cd245", "#665cd7", "#65d863", "#d053cd", "#d2c844", "#8737a7", "#62dc9c", "#db4795", "#50923a", "#cb83cf", "#808134", "#6f8eda", "#d57231", "#5d5597", "#d29e46", "#82bde1", "#cb3f31", "#7addd3", "#c53e5e", "#5d9d7c", "#934376", "#bbd18f", "#4f6280", "#d57f78", "#498f9c", "#82542f", "#cea6ca", "#45623c", "#86535c", "#d0b292"))
dev.off()


mut_J0_suivi_new$patient_jour<-paste(mut_J0_suivi_new$patient,mut_J0_suivi_new$jour, sep = "_j")
mut_J0_suivi_new$num_patient<-as.numeric(gsub("Patient","",mut_J0_suivi_new$patient))
mut_J0_suivi_new$num_patient<-mut_J0_suivi_new$num_patient * 100
mut_J0_suivi_new$sample_order<-mut_J0_suivi_new$num_patient + mut_J0_suivi_new$jour
list_patients_mut_new<-as.data.frame(unique(mut_J0_suivi_new$patient))

ggplot(data = mut_J0_suivi_new, aes(x=as.factor(reorder(mut_gene,order)), y=reorder(patient_jour,sample_order), fill=freq)) + 
  geom_tile(color = "black") +
  ylab("jour")+
  scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray",limits=c(5,100))+
  theme(axis.text.x = element_text(angle=90,size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


suivis_new_mut<-unique(mut_J0_suivi_new$patient)
for (p in 1:length(suivis_new_mut)){
  tmp_mut_J0_suivi_new<-mut_J0_suivi_new[which(mut_J0_suivi_new$patient==suivis_new_mut[p]),]
  png(filename = paste(suivis_new_mut[p],"_lineplot_freq_New_mut_not_J0_suivi.png"), width=10, height=10,units="in",res = 200)
  p7<-ggplot(tmp_mut_J0_suivi_new, aes(x=mut_gene, y=freq, group=as.factor(jour))) +
    geom_line(aes( color=as.factor(jour)),size=0.6)+
    geom_point(aes(color=as.factor(jour)))+
    theme(axis.text.x = element_text(angle=90))
  print(p7)
  dev.off()
}


mut_J0_suivi_new_mut_spe<-mut_J0_suivi_new[which(is.na(mut_J0_suivi_new$mut_spe)==F),]
#mut_J0_suivi_new_mut_spe$patient_jour<-paste(mut_J0_suivi_new_mut_spe$patient,mut_J0_suivi_new_mut_spe$jour, sep = "_j")

png(filename = "mutations_nouvelles_spe_variants_dans_suivis.png", width=15, height=10,units="in",res = 200)
ggplot(mut_J0_suivi_new_mut_spe, aes(x=as.factor(reorder(mut_gene,order)), y=freq, group=patient_jour)) +
  geom_line(aes(color=patient_jour))+
  geom_point(aes(color=patient_jour))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))+
  scale_color_manual(values=c("#e26bb0", "#f1cf39", "#e0451e", "#312a55", "#747c27", "#98d05d", "#4352d0", "#5a5282", "#6ba7d3", "#e9c28d", "#41973b", "#2d5777", "#a3e53e", "#b538be", "#cea850", "#5aeb75", "#8b6f2f", "#e39084", "#7a2daf", "#67e842", "#e48824", "#81a431", "#b48e69", "#e93280", "#2d1b88", "#c4a82a", "#393a86", "#825169", "#e0a7bc", "#2c6c67", "#626fd2", "#984f93", "#9a5f55", "#4e74a0", "#9f3344", "#e047e2", "#ac7ecd", "#594422", "#e39155", "#d0d22c", "#813490", "#592562", "#8761eb", "#ea7255", "#3f742b", "#56eaa0", "#4f86d6", "#59233c", "#92ddd5", "#7e245c", "#e9a1e3", "#c67386", "#5332cc", "#64c886", "#bbd940", "#b573a0", "#accc7e", "#cbce9f", "#582386", "#c580e6", "#d74460", "#ec3eb6", "#ae356c", "#49bc4e", "#8bc49c", "#d83a3e", "#958ac5", "#bfbbeb", "#ba2d95", "#ab7223", "#418760", "#7b3c19", "#96ec83", "#788654", "#e9688f", "#314921", "#1e3753", "#4ea3a1", "#d7d56b", "#55e2c7", "#d56bd5", "#622423", "#082742", "#912419", "#a242e3", "#b2562b", "#982092", "#58be2c", "#9c71e2", "#63cde8"))
dev.off()


# ### focus on the new mut seen in several patients (>4)
# focus<-new_mut_tab$mut[which(new_mut_tab$nb_patients_suivi>=4 & new_mut_tab$nb_J0==0)]
# mut_focus<-mut_J0_suivi_new[which(mut_J0_suivi_new$mut_gene%in%focus),]
# list_mut_focus<-unique(mut_focus$mut_gene)
# 
# for (f in 1:length(focus)){
#   focus_tab<-mut_focus[which(mut_focus$mut_gene==focus[f]),]
#   png(filename =paste(focus[f], "_dans_suivis.png",sep = ""), width=15, height=10,units="in",res = 200)
#   p8<-ggplot(focus_tab, aes(x=jour, y=freq, group=patient)) +
#     geom_line(aes(color=patient))+
#     geom_point(aes(color=patient))+
#     ggtitle(focus[f])+
#     theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))+
#     scale_color_manual(values=c("#9cd245", "#665cd7", "#65d863", "#d053cd", "#d2c844", "#8737a7", "#62dc9c", "#db4795", "#50923a", "#cb83cf", "#808134", "#6f8eda", "#d57231", "#5d5597", "#d29e46", "#82bde1", "#cb3f31", "#7addd3", "#c53e5e", "#5d9d7c", "#934376", "#bbd18f", "#4f6280", "#d57f78", "#498f9c", "#82542f", "#cea6ca", "#45623c", "#86535c", "#d0b292"))
#   print(p8)
#   dev.off()
# }
##############################################################################################################################################

### are variants that start min and become maj seen outside of the suivi pannel?
setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/")

mut_J0_non_suivi<-read.table(file="Mutations no_indel_J0_seq1-2-3-4-5-6-7_pitie .txt",sep = "\t", header = T)
mut_J0_non_suivi<-mut_J0_non_suivi[which((mut_J0_non_suivi$Id_sample%in%tab_suivi_filter$num_dossier)==F),]
min_maj_mut<-unique(min_to_maj_tab$mutation_gene)
mut_J0_non_suivi$mutation_gene<-paste(mut_J0_non_suivi$mutation,mut_J0_non_suivi$gene,sep="_")
mut_J0_min_maj<-mut_J0_non_suivi[which(mut_J0_non_suivi$mutation_gene%in%min_maj_mut),]
length(unique(mut_J0_min_maj$Id_sample))
length(unique(mut_J0_min_maj$mutation_gene))



##############################################################################################################################
## Heatmaps mutations min to maj

setwd("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/Suivi/")
list_min_maj<-unique(min_to_maj_tab$patient)
min_to_maj_tab<-min_to_maj_tab[which(min_to_maj_tab$freq!=0 | is.na(min_to_maj_tab$freq)==T),]

min_to_maj_tab$gene_start<-genes$Start[match(min_to_maj_tab$gene,genes$GeneId)]
min_to_maj_tab$gene_order<-min_to_maj_tab$gene_start * 10000
min_to_maj_tab$order<-min_to_maj_tab$gene_order + min_to_maj_tab$AA_pos

#min_to_maj_tab$freq[which(is.na(min_to_maj_tab$gene))]<-NA

for (m in 1:length(list_min_maj)){
  # makes the graphs 4 by 4 and print them in one figure
  if((m +3) %%4 == 0){
    # writes the data into 4 tables
    tab1<-min_to_maj_tab[which(min_to_maj_tab$patient==list_min_maj[m]),]
    rec1<-tab1 %>% group_by(jour,Id_sample) %>% summarise
    rec1$prelev<-metadata$VARIABLE[match(rec1$Id_sample, metadata$STAT_DOSSIER)]
    prelev1 <- ifelse(rec1$prelev == "Ecouvillon nasal", "red", "blue")
    tab2<-min_to_maj_tab[which(min_to_maj_tab$patient==list_min_maj[m+1]),]
    rec2<-tab2 %>% group_by(jour,Id_sample) %>% summarise
    rec2$prelev<-metadata$VARIABLE[match(rec2$Id_sample, metadata$STAT_DOSSIER)]
    prelev2 <- ifelse(rec2$prelev == "Ecouvillon nasal", "red", "blue")
    tab3<-min_to_maj_tab[which(min_to_maj_tab$patient==list_min_maj[m+2]),]
    rec3<-tab3 %>% group_by(jour,Id_sample) %>% summarise
    rec3$prelev<-metadata$VARIABLE[match(rec3$Id_sample, metadata$STAT_DOSSIER)]
    prelev3 <- ifelse(rec3$prelev == "Ecouvillon nasal", "red", "blue")
    tab4<-min_to_maj_tab[which(min_to_maj_tab$patient==list_min_maj[m+3]),]
    rec4<-tab4 %>% group_by(jour,Id_sample) %>% summarise
    rec4$prelev<-metadata$VARIABLE[match(rec4$Id_sample, metadata$STAT_DOSSIER)]
    prelev4 <- ifelse(rec4$prelev == "Ecouvillon nasal", "red", "blue")
    
    # create the 4 graphs 
    plot1<-ggplot(data = tab1, aes(x=fct_reorder(mutation_gene, order), y=as.factor(jour), fill=freq)) + 
      geom_tile() +
      ylab("jour")+
      scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray",limits=c(5,100))+
      theme(axis.text.x = element_text(angle=90,size = 7),
            axis.title.x = element_blank(),
            legend.position='none',
            axis.text.y = element_text(colour = prelev1)) + 
      ggtitle(tab1$patient[1])
    plot2<-ggplot(data = tab2, aes(x=fct_reorder(mutation_gene, order), y=as.factor(jour), fill=freq)) + 
      geom_tile() +
      ylab("jour")+
      scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray",limits=c(5,100))+
      theme(axis.text.x = element_text(angle=90,size = 7),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour = prelev2)) + 
      ggtitle(tab2$patient[1])
    plot3<-ggplot(data = tab3, aes(x=fct_reorder(mutation_gene, order), y=as.factor(jour), fill=freq)) + 
      geom_tile() +
      ylab("jour")+
      scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray",limits=c(5,100))+
      theme(axis.text.x = element_text(angle=90,size = 7),
            legend.position='none',
            axis.text.y = element_text(colour = prelev3)) +
      xlab("Mutation")+
      ggtitle(tab3$patient[1])
    plot4<-ggplot(data = tab4, aes(x=fct_reorder(mutation_gene, order), y=as.factor(jour), fill=freq)) + 
      geom_tile() +
      ylab("jour")+
      scale_fill_gradientn(colours = c("blue","cyan","green","yellow","red"), na.value = "darkgray",limits=c(5,100))+
      theme(axis.text.x = element_text(angle=90,size = 7),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour = prelev4)) + 
      xlab("Mutation")+
      ggtitle(tab4$patient[1])
    
    name<-paste("Patients",gsub("Patient","",tab1$patient[1]),gsub("Patient","",tab2$patient[1]),gsub("Patient","",tab3$patient[1]),gsub("Patient","",tab4$patient[1]),sep="_")
    
    # write the figure containg the 4 graphs
    png(filename = paste("Heatmap_freq_mut_min_to_maj_suivi_",name,".png",sep=""), width=10, height=10,units="in",res = 200)
    grid.arrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2) 
    dev.off()
  }
}

###################################### Comp cases where two different samples (nasal & bronch-pahryngé) the same day ###############################
#### count the number of differences between each days of a patient
comp_suivis<-function(patient){
  mut_patient<-mut_J0_suivi[which(mut_J0_suivi$patient==patient),]
  mut_patient$prelev<-metadata$VARIABLE[match(mut_patient$Id_sample, metadata$STAT_DOSSIER)]
  mut_patient$jour_prelev<-paste(mut_patient$jour,mut_patient$prelev,sep="_")
  jours_suivis<-unique(mut_patient$jour_prelev)
  mut_patient$mut<-"No_mut"
  mut_patient$mut[which(mut_patient$obs=="major" | mut_patient$obs=="minor")]<-"Mut"
  print(length(unique(mut_patient$mut_gene[which(mut_patient$mut=="Mut")])))
  comp_mat<-as.data.frame(matrix(NA,length(jours_suivis),length(jours_suivis)))
  colnames(comp_mat)<-jours_suivis
  rownames(comp_mat)<-jours_suivis
  for(a in 1:length(jours_suivis)){
    for(b in 1:length(jours_suivis)){
      comp_tab<-data.frame(Ja = mut_patient[which(mut_patient$jour_prelev==jours_suivis[a]),c(12,11,15,4)], Jb = mut_patient[which(mut_patient$jour_prelev==jours_suivis[b]),c(12,11,15,4)])
      colnames(comp_tab)<-gsub("Ja", paste("J",jours_suivis[a],sep=""),colnames(comp_tab))
      colnames(comp_tab)<-gsub("Jb", paste("J",jours_suivis[b],sep=""),colnames(comp_tab))
      comp_tab$diff<-"NO"
      comp_tab$diff[which(comp_tab[,3]!=comp_tab[,7])]<-"Diff"
      comp_mat[a,b]<-length(comp_tab$diff[which(comp_tab$diff=="Diff")])
    }
  }
  write.csv(comp_mat,file=paste("nb_mut_diff_",patient,".csv",sep=""))
}

comp_suivis("Patient119")
comp_suivis("Patient103")
comp_suivis("Patient77")


##### OLDER: Look case by case 
mut_J0_suivi_771<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112207096440"),]
mut_J0_suivi_771$mutation_gene<-paste(mut_J0_suivi_771$mutation,mut_J0_suivi_771$gene,sep="_")
length(which(mut_J0_suivi_771$obs=="major" | mut_J0_suivi_771$obs=="minor"))
mut_J0_suivi_772<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112207097950"),]
mut_J0_suivi_772$mutation_gene<-paste(mut_J0_suivi_772$mutation,mut_J0_suivi_772$gene,sep="_")
length(which(mut_J0_suivi_772$obs=="major" | mut_J0_suivi_772$obs=="minor"))

mut_J0_suivi_77<-full_join(mut_J0_suivi_771,mut_J0_suivi_772,by="mutation_gene",suffix = c(".1", ".2"))
mut_J0_suivi_77<-mut_J0_suivi_77[which(mut_J0_suivi_77$obs.1=="major" | mut_J0_suivi_77$obs.1=="minor" | mut_J0_suivi_77$obs.2=="major" | mut_J0_suivi_77$obs.2=="minor"),]

mut_J0_suivi_77$comp<-NA
for (c in 1:nrow(mut_J0_suivi_77)){
  if(mut_J0_suivi_77$obs.1[c]==mut_J0_suivi_77$obs.2[c]){
    mut_J0_suivi_77$comp[c]<-"Id"
  } else {
    mut_J0_suivi_77$comp[c]<-"Diff"
  }
}

mut_J0_suivi_77diff<-mut_J0_suivi_77[which(mut_J0_suivi_77$comp=="Diff"),]



mut_J0_suivi_1031<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112201055247"),]
mut_J0_suivi_1031$mutation_gene<-paste(mut_J0_suivi_1031$mutation,mut_J0_suivi_1031$gene,sep="_")
length(which(mut_J0_suivi_1031$obs=="major" | mut_J0_suivi_1031$obs=="minor"))
mut_J0_suivi_1032<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112201056353"),]
mut_J0_suivi_1032$mutation_gene<-paste(mut_J0_suivi_1032$mutation,mut_J0_suivi_1032$gene,sep="_")
length(which(mut_J0_suivi_1032$obs=="major" | mut_J0_suivi_1032$obs=="minor"))

mut_J0_suivi_103<-full_join(mut_J0_suivi_1031,mut_J0_suivi_1032,by="mutation_gene",suffix = c(".1", ".2"))
mut_J0_suivi_103<-mut_J0_suivi_103[which(mut_J0_suivi_103$obs.1=="major" | mut_J0_suivi_103$obs.1=="minor" | mut_J0_suivi_103$obs.2=="major" | mut_J0_suivi_103$obs.2=="minor"),]

mut_J0_suivi_103$comp<-NA
for (c in 1:nrow(mut_J0_suivi_103)){
  if(mut_J0_suivi_103$obs.1[c]==mut_J0_suivi_103$obs.2[c]){
    mut_J0_suivi_103$comp[c]<-"Id"
  } else {
    mut_J0_suivi_103$comp[c]<-"Diff"
  }
}

mut_J0_suivi_103diff<-mut_J0_suivi_103[which(mut_J0_suivi_103$comp=="Diff"),]



mut_J0_suivi_1191<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112206103360"),]
mut_J0_suivi_1191$mutation_gene<-paste(mut_J0_suivi_1191$mutation,mut_J0_suivi_1191$gene,sep="_")
length(which(mut_J0_suivi_1191$obs=="major" | mut_J0_suivi_1191$obs=="minor"))
mut_J0_suivi_1192<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112207100753"),]
mut_J0_suivi_1192$mutation_gene<-paste(mut_J0_suivi_1192$mutation,mut_J0_suivi_1192$gene,sep="_")
length(which(mut_J0_suivi_1192$obs=="major" | mut_J0_suivi_1192$obs=="minor"))
mut_J0_suivi_1193<-mut_J0_suivi[which(mut_J0_suivi$Id_sample=="112207102446"),]
mut_J0_suivi_1193$mutation_gene<-paste(mut_J0_suivi_1193$mutation,mut_J0_suivi_1193$gene,sep="_")
length(which(mut_J0_suivi_1193$obs=="major" | mut_J0_suivi_1193$obs=="minor"))

#comp 2 samples broncho alveolaire J149 vs J182
mut_J0_suivi_119ba<-full_join(mut_J0_suivi_1191,mut_J0_suivi_1193,by="mutation_gene",suffix = c(".1", ".3"))
mut_J0_suivi_119ba<-mut_J0_suivi_119ba[which(mut_J0_suivi_119ba$obs.1=="major" | mut_J0_suivi_119ba$obs.1=="minor" | mut_J0_suivi_119ba$obs.3=="major" | mut_J0_suivi_119ba$obs.3=="minor"),]

mut_J0_suivi_119ba$comp<-NA
for (c in 1:nrow(mut_J0_suivi_119ba)){
  if(mut_J0_suivi_119ba$obs.1[c]==mut_J0_suivi_119ba$obs.3[c]){
    mut_J0_suivi_119ba$comp[c]<-"Id"
  } else {
    mut_J0_suivi_119ba$comp[c]<-"Diff"
  }
}

summary(as.factor(mut_J0_suivi_119ba$comp))
length(which(mut_J0_suivi_119ba$comp=="Diff" & (mut_J0_suivi_119ba$obs.1=="major" | mut_J0_suivi_119ba$obs.1=="minor")))
length(which(mut_J0_suivi_119ba$comp=="Diff" & (mut_J0_suivi_119ba$obs.3=="major" | mut_J0_suivi_119ba$obs.3=="minor")))
length(which(mut_J0_suivi_119ba$comp=="Diff" & ((mut_J0_suivi_119ba$obs.3=="major" & mut_J0_suivi_119ba$obs.1=="minor")| (mut_J0_suivi_119ba$obs.1=="major" & mut_J0_suivi_119ba$obs.3=="minor"))))

#comp 2 samples consecutive days J181 vs J182
mut_J0_suivi_119c<-full_join(mut_J0_suivi_1192,mut_J0_suivi_1193,by="mutation_gene",suffix = c(".2", ".3"))
mut_J0_suivi_119c<-mut_J0_suivi_119c[which(mut_J0_suivi_119c$obs.2=="major" | mut_J0_suivi_119c$obs.2=="minor" | mut_J0_suivi_119c$obs.3=="major" | mut_J0_suivi_119c$obs.3=="minor"),]

mut_J0_suivi_119c$comp<-NA
for (c in 1:nrow(mut_J0_suivi_119c)){
  if(mut_J0_suivi_119c$obs.2[c]==mut_J0_suivi_119c$obs.3[c]){
    mut_J0_suivi_119c$comp[c]<-"Id"
  } else {
    mut_J0_suivi_119c$comp[c]<-"Diff"
  }
}

summary(as.factor(mut_J0_suivi_119c$comp))
length(which(mut_J0_suivi_119c$comp=="Diff" & (mut_J0_suivi_119c$obs.2=="major" | mut_J0_suivi_119c$obs.2=="minor")))
length(which(mut_J0_suivi_119c$comp=="Diff" & (mut_J0_suivi_119c$obs.3=="major" | mut_J0_suivi_119c$obs.3=="minor")))
length(which(mut_J0_suivi_119c$comp=="Diff" & ((mut_J0_suivi_119c$obs.3=="major" & mut_J0_suivi_119c$obs.2=="minor")| (mut_J0_suivi_119c$obs.2=="major" & mut_J0_suivi_119c$obs.3=="minor"))))


#comp the other 2 samples J149 vs J181
mut_J0_suivi_119d<-full_join(mut_J0_suivi_1192,mut_J0_suivi_1193,by="mutation_gene",suffix = c(".1", ".2"))
mut_J0_suivi_119d<-mut_J0_suivi_119d[which(mut_J0_suivi_119d$obs.2=="major" | mut_J0_suivi_119d$obs.2=="minor" | mut_J0_suivi_119d$obs.1=="major" | mut_J0_suivi_119d$obs.1=="minor"),]

mut_J0_suivi_119d$comp<-NA
for (c in 1:nrow(mut_J0_suivi_119d)){
  if(mut_J0_suivi_119d$obs.1[c]==mut_J0_suivi_119d$obs.2[c]){
    mut_J0_suivi_119d$comp[c]<-"Id"
  } else {
    mut_J0_suivi_119d$comp[c]<-"Diff"
  }
}

summary(as.factor(mut_J0_suivi_119d$comp))
length(which(mut_J0_suivi_119d$comp=="Diff" & (mut_J0_suivi_119d$obs.1=="major" | mut_J0_suivi_119d$obs.1=="minor")))
length(which(mut_J0_suivi_119d$comp=="Diff" & (mut_J0_suivi_119d$obs.2=="major" | mut_J0_suivi_119d$obs.2=="minor")))
length(which(mut_J0_suivi_119d$comp=="Diff" & ((mut_J0_suivi_119d$obs.2=="major" & mut_J0_suivi_119d$obs.1=="minor")| (mut_J0_suivi_119d$obs.1=="major" & mut_J0_suivi_119d$obs.2=="minor"))))


############################################ PCA suivi #######################################################################

## create the table for the pca 

#mut_JO_suivi<-rbind(mut_J0,mut_suivi)
#use mutation_gene to avoid issues if two genes have the same mutation
mut_J0_suivi$mutation_gene<-paste(mut_J0_suivi$mutation,mut_J0_suivi$gene,sep="_")
mut_J0_suivi<-mut_J0_suivi[which(mut_J0_suivi$Id_sample%in%tab_suivi_filter$num_dossier),]

tab_pca<-as.data.frame(matrix(0,length(unique(mut_J0_suivi$Id_sample)),length(unique(mut_J0_suivi$mutation_gene))))
colnames(tab_pca)<-unique(mut_J0_suivi$mutation_gene)
rownames(tab_pca)<-unique(mut_J0_suivi$Id_sample)


for(i in 1:ncol(tab_pca)){
  tmp<-mut_J0_suivi[which(mut_J0_suivi$mutation_gene==colnames(tab_pca[i])),]
  tmp_maj<-tmp[which(tmp$obs=="major"),]
  tmp_min<-tmp[which(tmp$obs=="minor"),]
  tmp_nocov<-tmp[which(tmp$obs=="discarded" | tmp$obs=="discarded_minor" ),]
  tab_pca[which(rownames(tab_pca)%in%tmp_maj$Id_sample),i]<-1
  tab_pca[which(rownames(tab_pca)%in%tmp_min$Id_sample),i]<-0.5
  tab_pca[which(rownames(tab_pca)%in%tmp_nocov$Id_sample),i]<- NA
  
  # tab_pca_mino[which(rownames(tab_pca_mino)%in%tmp_min$Id_sample),i]<-1
}

rm(tmp)
rm(tmp_maj)
rm(tmp_min)
rm(tmp_nocov)

# remove mutations with no variance in the pannel
tab_pca<-tab_pca[ , which(apply(tab_pca, 2, var,na.rm=T) != 0)]  


## pca
library(nipals)

pca<-nipals(tab_pca,ncomp=3)

barplot(pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))


li<-as.data.frame(pca$scores)
li$Patient<-tab_suivi2$num_patient[match(rownames(li),tab_suivi2$num_dossier)]
li$jour<-tab_suivi2$jour[match(rownames(li),tab_suivi2$num_dossier)]

# color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
#                       "#FC4E2A", "#E7298A", "#BD0026", "#800026",
#                       "#88B2BEFF", "#2F8FA8FF", "#4A879CFF", "#316283FF", "#345A78FF" ,"#0F3B6CFF" ,"#081F39FF" ,
#                       "#984EA3")
                      

# png(filename = "PCA_all_genes_Axes12_suivis.png", width=8, height=8,units="in",res = 400)
# ggplot(li, aes(x=Axis1 , y=Axis2, color=Patient , label= rownames(li)))+
#   geom_point(size=2)+
#   scale_color_paletteer_d("pals::polychrome")+
#   theme(panel.background = element_rect(fill="gray70"))
# dev.off()


png(filename = "PCA_all_genes_Axes12_suivis_better_zoom.png", width=16, height=12,units="in",res = 400)
ggplot(li, aes(x=PC1 , y=PC2, color=Patient , label= rownames(li)))+
  geom_point(size=2)+
  scale_color_paletteer_d("pals::polychrome")+
  theme(panel.background = element_rect(fill="gray70"))+
  facet_zoom(xlim = c(-0.03, 0.005))+
  geom_text_repel(aes(label = li$jour),size = 3,max.overlaps=50) 
dev.off()

# # no outliers
# li2<-li[which(li$Axis1>-10),]
# png(filename = "PCA_all_genes_Axes12_suivis_zoom.png", width=8, height=8,units="in",res = 400)
# ggplot(li2, aes(x=Axis1 , y=Axis2, color=Patient , label= rownames(li2)))+
#   geom_point(size=2)+
#   scale_color_paletteer_d("pals::polychrome")+
#   theme(panel.background = element_rect(fill="gray70"))
# dev.off()

### per patient

for (p in 1:length(list_patients_filter)){
  tab_pca_p<-tab_pca[which(rownames(tab_pca)%in%tab_suivi_filter$num_dossier[which(tab_suivi_filter$num_patient==list_patients_filter[p])]),]
  
  # remove mutations with no variance in the pannel
  tab_pca_p<-tab_pca_p[ , which(apply(tab_pca_p, 2, var,na.rm=T) != 0)]  
  
  
  ## pca
  library(nipals)
  
  pca_p<-nipals(tab_pca_p,ncomp=3)
  
  barplot(pca_p$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
  
  
  li_p<-as.data.frame(pca_p$scores)
  li_p$Patient<-tab_suivi2$num_patient[match(rownames(li_p),tab_suivi2$num_dossier)]
  li_p$jour<-tab_suivi2$jour[match(rownames(li_p),tab_suivi2$num_dossier)]
  li_p$prelev<-metadata$VARIABLE[match(rownames(li_p),metadata$STAT_DOSSIER)]
  li_p$prelev<-gsub("Aspiration bronchique", "L. broncho-alvéolaire", li_p$prelev)
  
  png(filename = paste("PCA_all_genes_Axes12_suivis_",list_patients_filter[p],".png",sep=""), width=8, height=8,units="in",res = 200)
  plotpca<-ggplot(li_p, aes(x=PC1 , y=PC2, color=prelev , label= rownames(li_p)))+
    geom_point(size=4)+
    #scale_color_paletteer_d("ggthemes::Classic_Cyclic")+
    scale_color_manual(values=c("red","blue"))+
    ggtitle(list_patients_filter[p])+
    #theme(panel.background = element_rect(fill="gray70"))+
    geom_text_repel(aes(label = li_p$jour),size = 3,max.overlaps=50) 
  print(plotpca)
  dev.off()

}

