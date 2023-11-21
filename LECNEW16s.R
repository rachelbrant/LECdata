#### 16S slice and dice LEC data 
## step 1. Import data created in qiime2 to R as a phyloseq object. This is done with the command qza_to_phyloseq from the package qiime2R
## The clean and good files are as follows: table = NEW_tableAUG29.qza, tree=rooted_NEW_UNITE.qza, taxonomic assignment is NEWESTtaxonomyUNITE.qza and metadata is FULLmetadata.txt
install.packages("qiime2R")
library(qiime2R)
library(MicrobeDS)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
BiocManager::install("cowplot")
#install.packages("labdsv")
library(labdsv)
BiocManager::install("ALDEx2")
library(ALDEx2)
library(qiime2R)
install.packages("DESeq2")
library(DESeq2)
library(qiimer)
library(phyloseq)
library(vegan)
library(microbiome)
library(ggplot2)
library(RColorBrewer) #pretty colour pallettes
library(gplots)   #venn diagrams
library(cowplot) #saving graphs
library(pheatmap)

physeqLEC16S2<-qza_to_phyloseq(features="/Users/rachelbrant/Desktop/LEC_data/table_LEC16S.qza",tree = "/Users/rachelbrant/Desktop/LEC_data/rooted-tree_LEC16S.qza", taxonomy="/Users/rachelbrant/Desktop/LEC_data/taxonomy_LEC16S.qza", metadata = "/Users/rachelbrant/Desktop/LEC_data/LEC-metadata4.txt")
physeqLEC16S2
physeqLEC16S22<-prune_taxa(as.logical(tax_table(physeqLEC16S2)[,1]=="d__Bacteria"),physeqLEC16S2)
physeqLEC16S22
table1<-tax_table(physeqLEC16S2)
table2<-otu_table(physeqLEC16S2)
write.csv(table1,"/Users/rachelbrant/Desktop/taxonomybac.csv")
sample_names(physeqLEC16S2)


##slice and dice 

##innoculate soil no plants 
ControlSoil<-prune_samples(sample_data(physeqLEC16S22)$Species=="Control",physeqLEC16S22)
ControlSoil1<-prune_samples(sample_data(ControlSoil)$Soil_type!="Control",ControlSoil)
ControlSoil
ControlSoil1
Soil<-prune_samples(sample_data(physeqLEC16S22)$Species!="Control",physeqLEC16S22)
Soil2<-prune_samples(sample_data(Soil)$Species!="GreenHouse",Soil)
Soil13<-prune_samples(sample_data(Soil2)$Soil_Age!="Control",Soil2)
Soil33<-prune_taxa(as.logical(tax_table(Soil2)[,5]!="Mitochondria"),Soil2)
Soil4<-prune_taxa(as.logical(tax_table(Soil33)[,3]!="NB1-j"),Soil33)
sample_names(ControlSoil)
range(sample_sums(ControlSoil))
rareSoil3<-rarefy_even_depth(Soil33,900,rngseed = 356789)
rareSoil3

rareSoil<-rarefy_even_depth(physeqLEC16S2,800,rngseed = 351590)
rareSoil
##Young soil all plants
YoungSoil<-prune_samples(sample_data(Soil2)$Soil_Age=="Young",Soil2)
YoungSoil<-prune_samples(sample_data(YoungSoil)$Species!="Control",YoungSoil)
YoungSoil
range(sample_sums(YoungSoil))
rareYoungSoil<-rarefy_even_depth(YoungSoil,1000,rngseed = 552597)
rareYoungSoil

##Middle soil all plants
MidSoil<-prune_samples(sample_data(Soil2)$Soil_Age=="Middle",Soil2)
MidSoil<-prune_samples(sample_data(MidSoil)$Species!="Control",MidSoil)
MidSoil
range(sample_sums(MiddleSoil))
rareMiddleSoil<-rarefy_even_depth(MiddleSoil,1000,rngseed = 351590)
rareMiddleSoil

##Old soil all plants
OldSoil<-prune_samples(sample_data(Soil2)$Soil_Age=="Old",Soil2)
OldSoil<-prune_samples(sample_data(OldSoil)$Species!="Control",OldSoil)
OldSoil
range(sample_sums(OldSoil))
rareOldSoil<-rarefy_even_depth(OldSoil,1000,rngseed = 351590)
rareOldSoil
###### plant species now#####
##CP species 
CPphylo<-prune_samples(sample_data(Soil2)$Species=="CP",Soil2)
CPphylo1<-prune_samples(sample_data(CPphylo)$Soil_Age!="Control",CPphylo)
range(sample_sums(CPphylo))
rareCPphylo<-rarefy_even_depth(CPphylo,1000,rngseed = 351590)
rareCPphylo
##GV species
GVphylo<-prune_samples(sample_data(Soil2)$Species=="GV",Soil2)
GVphylo1<-prune_samples(sample_data(GVphylo)$Soil_Age!="Control",GVphylo)
GVphylo
range(sample_sums(GVphylo))
rareGVphylo<-rarefy_even_depth(GVphylo,1000,rngseed = 353490)
rareGVphylo
###SA Species
SAphylo<-prune_samples(sample_data(Soil2)$Species=="SA",Soil2)
SAphylo1<-prune_samples(sample_data(SAphylo)$Soil_Age!="Control",SAphylo)
range(sample_sums(SAphylo))
rareSAphylo<-rarefy_even_depth(SAphylo,1000,rngseed = 549290)
rareSAphylo
##SC species
SCphylo<-prune_samples(sample_data(Soil2)$Species=="SC",Soil2)
SCphylo1<-prune_samples(sample_data(SCphylo)$Soil_Age!="Control",SCphylo)
range(sample_sums(SCphylo))
rareSCphylo<-rarefy_even_depth(SCphylo,1000,rngseed = 341390)
rareSCphylo

#####venn diagrams of shared ASVS
###################   SHARED AND UNIQUE SEQUENCES BY SITE
library(gplots) #required for Venn diagram. Online Venn diagram tools like Venny are also available. 
rareControlSoil
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(rareControlSoil)$Soil_Age))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]

#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
#Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(rareControlSoil)$Soil_Age %in% poplist[[k]],rareControlSoil)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(rareControlSoil))
##prune down to only those ASVS
phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)
install.packages("VennDiagram")



##Cp ven
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
library(phyloseq)
pops<-as.character(unique(sample_data(CPphylo)$Soil_Age))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]

#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(CPphylo)$Soil_Age %in% poplist[[k]],CPphylo)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(CPphylo))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)
install.packages("VennDiagram")
library(VennDiagram) 
venn.diagram(venn.list)



##gv ven
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
library(phyloseq)
pops<-as.character(unique(sample_data(GVphylo)$Soil_Age))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]


#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(GVphylo)$Soil_Age %in% poplist[[k]],GVphylo)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(GVphylo))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)



##SA ven 
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
library(phyloseq)
pops<-as.character(unique(sample_data(SAphylo)$Soil_Age))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]


#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(SAphylo)$Soil_Age %in% poplist[[k]],SAphylo)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SAphylo))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)


#####SC ven
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
library(phyloseq)
pops<-as.character(unique(sample_data(SCphylo)$Soil_Age))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]


#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(SCphylo)$Soil_Age %in% poplist[[k]],SCphylo)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SCphylo))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)

################richness
ControlSoil1$Soil_Age <- factor(ControlSoil1$Soil_Age,
                       levels = c('Young','Middle',"Old"),ordered = TRUE)
ControlRichness<-estimate_richness(ControlSoil1,measures=c("Observed","Shannon","InvSimpson"))
head(ControlRichness)
pl<-plot_richness(ControlSoil1,x="Soil_Age",measures=c("Shannon"))+ geom_boxplot(aes(fill=Soil_Age))+ theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8))
plot3<-pl+scale_fill_manual(values =c("#6699CC","#009E73","#999999"),breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
plot4<-plot3+theme(legend.position = "none")+scale_x_discrete(limits=c("Young","Middle","Old"))
plot4
plot_richness(rareControlSoil,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=8))

YoungRichness<-estimate_richness(rareYoungSoil,measures=c("Observed","Shannon","InvSimpson"))
head(YoungRichness)


MiddleRichness<-estimate_richness(rareMiddleSoil,measures=c("Observed","Shannon","InvSimpson"))
head(MiddleRichness)

OldRichness<-estimate_richness(rareOldSoil,measures=c("Observed","Shannon","InvSimpson"))
head(OldRichness)

########### species richness
CPRichness<-estimate_richness(rareCPphylo,measures=c("Observed","Shannon","InvSimpson"))
head(CPRichness)
plot_richness(rareCPphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plot_richness(rareCPphylo,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))


GVRichness<-estimate_richness(rareGVphylo,measures=c("Observed","Shannon","InvSimpson"))
head(GVRichness)
plot_richness(rareGVphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plot_richness(rareGVphylo,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))


SARichness<-estimate_richness(rareSAphylo,measures=c("Observed","Shannon","InvSimpson"))
head(SARichness)
plot_richness(rareSAphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plot_richness(rareSAphylo,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))


SCRichness<-estimate_richness(rareSCphylo,measures=c("Observed","Shannon","InvSimpson"))
head(SCRichness)
plot_richness(rareSCphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plot_richness(rareSCphylo,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))



##################beta and ordination
library(vegan)
ord.nmds.brayControl <- ordinate(ControlSoil1, method="NMDS",k=3, distance="bray",trymax=100)
ord1<-plot_ordination(ControlSoil1, ord.nmds.brayControl, color="Soil_Age", title="Bray NMDS")
ord1+stat_ellipse(geom="polygon",aes(fill=Soil_Age),type="norm",alpha=0.4) + theme_bw()
ord1
ord.nmds.bray <- ordinate(ControlSoil1, method="NMDS",k=2, distance="bray",trymax=50)
ordinationLECITS <- ordinate(ControlSoil1, "PCoA", "bray")

beta_div_bray_cleanLEC16S <- plot_ordination(ControlSoil1, ord.nmds.bray, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLEC16S <- beta_div_bray_cleanLEC16S + stat_ellipse() + ggtitle("A")  + theme_classic() + scale_color_brewer("Age", palette = "Set2")
braybac<-bray_clean_presenceLEC16S+scale_color_manual(values =c("#6699CC","darkorange1","#999999"),breaks=c('Young', 'Middle', 'Old'),name = "Soil Age",labels=c('Young'="Young", 'Middle'="Intermediate", 'Old'="Old")) + theme(axis.text=element_text(size=8),axis.title = element_text(size=8))
braybac2<-braybac+xlim(-1,1.5)+ylim(-0.6,0.6)
braybac3<-braybac2+theme(legend.position = "none")
braybac4<-braybac3+annotate("text", x=1, y=-0.55, label= "stress = 0.0867",size=5)
braybac5<-braybac4+theme(axis.text.x = element_text(size=13),
                axis.text.y = element_text(size=13))+theme(axis.title=element_text(size=14))
braybac5
##bray curtis species 
###Bray curtis Cp 
ord.nmds.brayCP <- ordinate(CPphylo, method="NMDS",k=3, distance="bray",trymax=500)
ord1<-plot_ordination(CPphylo, ord.nmds.brayCP, color="Soil_Age", title="Bray NMDS")
ord1+stat_ellipse(geom="polygon",aes(fill=Soil_Age),type="norm",alpha=0.4) + theme_bw()

ord.nmds.brayCP16 <- ordinate(CPphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECCP16 <- ordinate(CPphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSCP16 <- plot_ordination(CPphylo, ordinationLECCP16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECCP16<-beta_div_bray_cleanLECITSCP16 + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYCP<-bray_clean_presenceLECCP16+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYCP+ylim(-1,1)+xlim(-1,1)

ord.nmds.brayGV16 <- ordinate(GVphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECGV16 <- ordinate(GVphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSGV16 <- plot_ordination(GVphylo, ordinationLECGV16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECGV16<-beta_div_bray_cleanLECITSGV16 + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYGV<-bray_clean_presenceLECGV16+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYGV+ylim(-1,1)+xlim(-1,1)

ord.nmds.braySA16 <- ordinate(SAphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECSA16 <- ordinate(SAphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSSA16 <- plot_ordination(SAphylo, ordinationLECSA16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECSA16<-beta_div_bray_cleanLECITSSA16 + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYSA<-bray_clean_presenceLECSA16+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYSA+ylim(-1,1)+xlim(-1,1)

ord.nmds.braySC16 <- ordinate(SCphylo, method="NMDS",k=6, distance="bray",trymax=500)
ordinationLECSC16 <- ordinate(SCphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSSC16 <- plot_ordination(SCphylo, ordinationLECSC16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECSC16<-beta_div_bray_cleanLECITSSC16 + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYSC<-bray_clean_presenceLECSC16+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYSC+ylim(-1.2,1)+xlim(-1.2,1)

##bray curtis GV
ord.nmds.brayGV16 <- ordinate(GVphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECGV16 <- ordinate(GVphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSGV16 <- plot_ordination(GVphylo, ordinationLECGV16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECGV16<-beta_div_bray_cleanLECITSGV16 + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECGV16+scale_color_brewer(palette = "Set1", direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))






##bray curtis SA
ord.nmds.braySA16 <- ordinate(SAphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECSA16 <- ordinate(SAphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSSA16 <- plot_ordination(SAphylo, ordinationLECSA16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECSA16<-beta_div_bray_cleanLECITSSA16 + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECSA16+scale_color_brewer(palette = "Set1", direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))




##bray curtis SC
ord.nmds.braySC16 <- ordinate(SCphylo, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECSC16 <- ordinate(SCphylo, "PCoA", "bray")
beta_div_bray_cleanLECITSSC16 <- plot_ordination(SCphylo, ordinationLECSC16, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECSC16<-beta_div_bray_cleanLECITSSC16 + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECSC16+scale_color_brewer(palette = "Set1", direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))




library(dplyr)
soil_meta<-data.frame(sample_data(Soil2))
soil_meta$sampleid<-rownames(soil_meta)
soil_meta
soil_richnessControl<-estimate_richness(Soil2,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessControl$sampleid<-rownames(soil_richnessControl)
soil_richness132<-left_join(soil_richnessControl,soil_meta,"sampleid")
soil_richness132
plotrichness1<-plot_richness(ControlSoil,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
library(lme4)
m1<-lm(Shannon ~ Soil_Age*Species, data=soil_richness132)
summary(m1)


RichnessFungi<-estimate_richness(rareSoilFungi2,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessFungi)
plotrichness1<-plot_richness(SoilFungi2,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Species)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1+scale_fill_brewer(palette="RdBu",direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
plot_richness(rareControlSoilFungi2,x="Soil_Age",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))

library(dplyr)
soil_metaGRU<-data.frame(sample_data(rareSoil3))
soil_metaGRU$sampleid<-rownames(soil_metaGRU)
soil_metaGRU
soil_richnessGRU<-estimate_richness(rareSoil3,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessGRU$sampleid<-rownames(soil_richnessGRU)
soil_richnessGRU1<-left_join(soil_richnessGRU,soil_metaGRU,"sampleid")
soil_richnessGRU1
library(lme4)
bacteriaSoilSpecies<-lm(Shannon~ Species*Soil_Age, data=soil_richnessGRU1)
summary(bacteriaSoilSpecies)


library(MuMIn)
m1_ml<-update(m1,REML=F)
summary(m1_ml)
options(na.action = "na.fail")
library(MuMIn)
m1_aic_rank<-dredge(m1_ml)
m1_aic_rank



##per species richness differences 
library(dplyr)
library(ggplot2)
soil_metaCP<-data.frame(sample_data(rareCPphylo))
soil_metaCP$sampleid<-rownames(soil_metaCP)
soil_metaCP
soil_richnessCP<-estimate_richness(rareCPphylo,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessCP$sampleid<-rownames(soil_richnessCP)
soil_richnessCP1<-left_join(soil_richnessCP,soil_metaCP,"sampleid")
soil_richnessCP1
plotrichness1<-plot_richness(rareCPphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
library(lme4)
m1<-aov(Shannon ~ Soil_Age, data=soil_richnessCP1)
summary(m1)


library(dplyr)
library(ggplot2)
soil_metaGV<-data.frame(sample_data(rareGVphylo))
soil_metaGV$sampleid<-rownames(soil_metaGV)
soil_metaGV
soil_richnessGV<-estimate_richness(rareGVphylo,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessGV$sampleid<-rownames(soil_richnessGV)
soil_richnessGV<-left_join(soil_richnessGV,soil_metaGV,"sampleid")
soil_richnessGV
plotrichness1GV<-plot_richness(rareGVphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1GV+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
library(lme4)
m1<-aov(Shannon ~ Soil_Age, data=soil_richnessGV)
summary(m1)

library(dplyr)
library(ggplot2)
soil_metaSA<-data.frame(sample_data(rareSAphylo))
soil_metaSA$sampleid<-rownames(soil_metaSA)
soil_metaSA
soil_richnessSA<-estimate_richness(rareSAphylo,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessSA$sampleid<-rownames(soil_richnessSA)
soil_richnessSA1<-left_join(soil_richnessSA,soil_metaSA,"sampleid")
soil_richnessSA1
plotrichness1SA<-plot_richness(rareSAphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1SA+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
library(lme4)
m1<-aov(Shannon ~ Soil_Age, data=soil_richnessSA1)
summary(m1)


library(dplyr)
library(ggplot2)
soil_metaSC<-data.frame(sample_data(rareSCphylo))
soil_metaSC$sampleid<-rownames(soil_metaSC)
soil_metaSC
soil_richnessSC<-estimate_richness(rareSCphylo,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessSC$sampleid<-rownames(soil_richnessSC)
soil_richnessSC<-left_join(soil_richnessSC,soil_metaSC,"sampleid")
soil_richnessSC
plotrichness1SC<-plot_richness(rareSCphylo,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichness1SC+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
library(lme4)
m1<-aov(Shannon ~ Soil_Age, data=soil_richnessSC)
summary(m1)
