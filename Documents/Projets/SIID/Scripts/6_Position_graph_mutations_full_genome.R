############################## Graph position mutations on the genome ############################

# 1) Packages input 
library(ggplot2)
library(gridExtra)
library(stringr)
library(readxl)
library(forcats)
library(dplyr)
library(BioCircos)

# 2) input files 

#### Repertories ####
INREP<-("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
INFOREP<-("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/")
DATAREP<-("C:/Users/Romane/Documents/Projets/SIID/")
OUTREP<-("C:/Users/Romane/Documents/Projets/SIID/2_Filter_trim/Accurate/minor_5_20pct/")
setwd(INREP)


#### mutations ####
input1<-read.table(file = "Mutations_no_indel_J0_full_seq_1-2-3-4-5-6-7_avec_pitie2_min5-20pct.txt",sep="\t",header=T,stringsAsFactors = F)
#input1<-input1[which(input1$obs=="minor" | input1$obs=="major"),]


#### variants ####
infos_p<-read.csv(paste(INFOREP,"nextclade_pitie.csv",sep=""),sep=";")
infos_b<-read.csv(paste(INFOREP,"nextclade_1-2-3-4-5-6-7.csv",sep=""),sep=";")
infos<-rbind(infos_b,infos_p)
infos$clade_legacy<- gsub(" ", "_",
                   gsub("\\(","",
                        gsub("\\)","",
                             gsub(",","", infos$clade_legacy))))


input1$variant<-infos$clade_legacy[match(input1$Id_sample,infos$seqName)]
rm(infos_b)
rm(infos_p)


#### Genes ####
setwd(paste(DATAREP,"/data",sep=""))
info_ref<-read.table(file="Genes_SARS-CoV2.txt",header=T)
info_ref$aa_length<-(info_ref$End - info_ref$Start +1) / 3

input1$gene_start<-info_ref$Start[match(input1$gene,info_ref$GeneId)]
input1$gene_aa<-paste(input1$gene,input1$AA_pos,sep="_")


# 3) count the number of samples with a mutation + nb of covered samples for each variant 

mut<-unique(input1$gene_aa)

for (i in 1:length(mut)){
  tmp<-input1[which(input1$gene_aa==mut[i]),]
  variants<-unique(tmp$variant)
  for (j in 1:length(variants)){
    graph_tab_tmp<-data.frame(gene=tmp$gene[1],aa_pos=tmp$AA_pos[1],variant=variants[j],nb_mut=NA,nb_cov=NA, gene_start=tmp$gene_start[1])
    graph_tab_tmp$nb_mut<-length(unique(tmp$Id_sample[which(tmp$obs[which(tmp$variant==variants[j])]=="major" | tmp$obs[which(tmp$variant==variants[j])]=="minor")]))
    graph_tab_tmp$nb_cov<-length(unique(tmp$Id_sample[which(tmp$obs[which(tmp$variant==variants[j])]=="major" | tmp$obs[which(tmp$variant==variants[j])]=="minor" | tmp$obs[which(tmp$variant==variants[j])]=="no_mutation")]))
    
    if(j==1){
      graph_tab<-graph_tab_tmp
    } else {
      graph_tab<-rbind(graph_tab,graph_tab_tmp)
    }
  }
  if(i==1){
    graph_tab_full<-graph_tab
  } else {
    graph_tab_full<-rbind(graph_tab_full, graph_tab)
  }
  
}

graph_tab_full$nb_mut_cov<-graph_tab_full$nb_mut/graph_tab_full$nb_cov


# 4) Graph nb samples with mut / nb samples covered for each variant at each position

### all variants
color_set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                      "#FC4E2A", "#E7298A", "#BD0026", "#800026",
                      "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5" ,"#08519C" ,"#08306B" ,
                      "#984EA3","grey")
setwd(OUTREP)             
  
png(filename = "Mutations_positions.png", width=25, height=8,units="in",res = 400)
p1<-ggplot(data=graph_tab_full, aes(x=aa_pos, y=nb_mut_cov, fill=variant)) +
  geom_bar(stat="identity", width=5)+
  facet_grid(. ~ fct_reorder(gene, gene_start), scales = "free_x", space = "free_x")+
  scale_fill_manual(values = color_set)
print(p1)
dev.off()


### alpha, delta & omicron 
color_set2<-c("#B3B3B3",  "#800026",
                      "#C6DBEF", "#9ECAE1",  "#2171B5" )

graph_tab_full2<-graph_tab_full[which(graph_tab_full$variant%in%c("20I_Alpha_V1","21J_Delta","21K_Omicron","21L_Omicron","22B_Omicron")),]

png(filename = "Mutations_positions_sel_variants.png", width=25, height=8,units="in",res = 400)
p2<-ggplot(data=graph_tab_full2, aes(x=aa_pos, y=nb_mut_cov, fill=variant)) +
  geom_bar(stat="identity", width=5)+
  facet_grid(. ~ fct_reorder(gene, gene_start), scales = "free_x", space = "free_x")+
  scale_fill_manual(values = color_set2)
print(p2)
dev.off()


### 5) circos

### Interactive circos

info_ref<-info_ref[order(info_ref$Start),]
tmp_genome<-as.data.frame(matrix(info_ref$aa_length,ncol = nrow(info_ref)))
colnames(tmp_genome)<-info_ref$GeneId
MyGenome<-as.list(tmp_genome)
graph_tab_full2<-arrange(graph_tab_full2 ,gene_start, aa_pos)


tracklist<-BioCircosSNPTrack("22B_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="22B_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="22B_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="22B_Omicron")],
                                         colors = "#08306B", labels = "", size = 2, shape = "circle",
                                         opacities = 1, maxRadius = 0.95, minRadius = 0.85, range = 0)

tracklist<-tracklist + BioCircosSNPTrack("21L_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="21L_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21L_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21L_Omicron")],
                                         colors = "#2171B5", labels = "", size = 2, shape = "circle",
                                         opacities = 1, maxRadius = 0.80, minRadius = 0.70, range = 0)

tracklist<-tracklist + BioCircosSNPTrack("21K_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="21K_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21K_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21K_Omicron")],
                                         colors = "#6BAED6", labels = "", size = 2, shape = "circle",
                                         opacities = 1, maxRadius = 0.65, minRadius = 0.55, range = 0)

tracklist<-tracklist + BioCircosSNPTrack("Delta", graph_tab_full2$gene[which(graph_tab_full2$variant=="21J_Delta")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21J_Delta")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21J_Delta")],
                                         colors = "#800026", labels = "", size = 2, shape = "circle",
                                         opacities = 1, maxRadius = 0.50, minRadius = 0.4, range = 0)

tracklist<-tracklist + BioCircosSNPTrack("Alpha", graph_tab_full2$gene[which(graph_tab_full2$variant=="20I_Alpha_V1")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="20I_Alpha_V1")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="20I_Alpha_V1")],
                  colors = "#B3B3B3", labels = "", size = 2, shape = "circle",
                  opacities = 1, maxRadius = 0.35, minRadius = 0.15, range = 0)


#png(filename = "Ciros_mutations_positions_sel_variants.png", width=8, height=8,units="in",res = 400)
p3<-BioCircos(tracklist,genome = MyGenome, genomeFillColor = c("#cc4a6d", "#c95733", "#c9a443", "#a4db52", "#5d8d37", "#61d89d", "#5c8ed9", "#605ec0", "#7f47d2", "#d082d0", "#923a85", "#d746c2"),
          genomeTicksScale = 100,genomeLabelOrientation = 90, genomeLabelDy = 30, genomeLabelTextSize = "11pt")
print(p3)
#dev.off()
## Save with Export -> save as web page


## lines
tracklist2<-BioCircosLineTrack("22B_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="22B_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="22B_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="22B_Omicron")],
                                         color = "#08306B", opacities = 1, maxRadius = 1, minRadius = 0.8, range = 0)

tracklist2<-tracklist2 + BioCircosLineTrack("21L_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="21L_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21L_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21L_Omicron")],
                                         color = "#2171B5", opacities = 1, maxRadius = 0.8, minRadius = 0.6, range = 0)

tracklist2<-tracklist2 + BioCircosLineTrack("21K_Omicron", graph_tab_full2$gene[which(graph_tab_full2$variant=="21K_Omicron")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21K_Omicron")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21K_Omicron")],
                                         color = "#6BAED6",opacities = 1, maxRadius = 0.6, minRadius = 0.4, range = 0)

tracklist2<-tracklist2 + BioCircosLineTrack("Delta", graph_tab_full2$gene[which(graph_tab_full2$variant=="21J_Delta")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="21J_Delta")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="21J_Delta")],
                                         color = "#800026", opacities = 1, maxRadius = 0.40, minRadius = 0.20, range = 0)

tracklist2<-tracklist2 + BioCircosLineTrack("Alpha", graph_tab_full2$gene[which(graph_tab_full2$variant=="20I_Alpha_V1")], graph_tab_full2$aa_pos[which(graph_tab_full2$variant=="20I_Alpha_V1")], graph_tab_full2$nb_mut[which(graph_tab_full2$variant=="20I_Alpha_V1")],
                                         color = "#B3B3B3", opacities = 1, maxRadius = 0.20, minRadius = 0, range = 0)


#png(filename = "Ciros_lines_mutations_positions_sel_variants.png", width=8, height=8,units="in",res = 400)
p4<-BioCircos(tracklist2,genome = MyGenome, genomeFillColor = c("#cc4a6d", "#c95733", "#c9a443", "#a4db52", "#5d8d37", "#61d89d", "#5c8ed9", "#605ec0", "#7f47d2", "#d082d0", "#923a85", "#d746c2"),
              genomeTicksScale = 100,genomeLabelOrientation = 90,genomeTicksTextSize = 12, genomeLabelDy = 40, genomeLabelTextSize = "12pt")
print(p4)
#dev.off()
## Save with Export -> save as web page