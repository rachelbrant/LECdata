###LEC ITS phyloseq downstream analyses
## sept 30 2022

library(qiime2R)
library(MicrobeDS)
library(labdsv)
library(ALDEx2)
library(qiime2R)
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

physeqLECFungi12<-qza_to_phyloseq(features="/Users/rachelbrant/Desktop/LEC_data/mergedLECFungi.qza",tree = "/Users/rachelbrant/Desktop/LEC_data/merged_rooted-treeLECFungi.qza", taxonomy="/Users/rachelbrant/Desktop/LEC_data/merged-taxaFungiNew.qza", metadata = "/Users/rachelbrant/Desktop/LEC_data/LEC-metadata4.txt")
physeqLECFungi12
table<-tax_table(physeqLECFungi12)
write.csv(table,"/Users/rachelbrant/Desktop/taxtable.csv")
sample_names(physeqLECFungi12)

##innoculate soil no plants 
ControlSoilFungi3<-prune_samples(sample_data(physeqLECFungi12)$Species=="Control",physeqLECFungi12)
ControlSoilFungi2<-prune_samples(sample_data(ControlSoilFungi3)$Soil_Age!="Control",ControlSoilFungi3)
ControlSoilFungi4<-prune_taxa(as.logical(tax_table(ControlSoilFungi2)[,1]=="Fungi"),ControlSoilFungi2)
sample_names(ControlSoilFungi2)
range(sample_sums(ControlSoilFungi4))
rareControlSoilFungi24<-rarefy_even_depth(ControlSoilFungi2,1000,rngseed = 354307)
rareControlSoilFungi24
ControlSoilFungi5<-prune_taxa(as.logical(tax_table(ControlSoilFungi4)[,6]!="unidentified"),ControlSoilFungi4)
ControlSoilFungi7<-prune_taxa(as.logical(tax_table(ControlSoilFungi5)[,6]!="NA"),ControlSoilFungi5)
ControlSoilFungi8<-prune_taxa(as.logical(tax_table(ControlSoilFungi4)[,5]!="unidentified"),ControlSoilFungi4)
ControlSoilFungi9<-prune_taxa(as.logical(tax_table(ControlSoilFungi8)[,5]!="NA"),ControlSoilFungi8)
ControlSoilFungi9
ControlSoilFungi10<-prune_taxa(as.logical(tax_table(physeqLECFungi12)[,2]=="Glomeromycota"),physeqLECFungi12)
ControlSoilFungi10
otu_table(ControlSoilFungi10) <- otu_table(ControlSoilFungi10) + 1
##prune down to only those ASVS
phy.sub.subFun<-prune_taxa(phy.sub.keepFun,ControlSoilFungi10)
###samples with plants 
SoilFungi<-prune_samples(sample_data(physeqLECFungi12)$Species!="Control",physeqLECFungi12)
SoilFungi2<-prune_samples(sample_data(SoilFungi)$Soil_type!="Control",SoilFungi)
SoilFungi13<-prune_samples(sample_data(SoilFungi)$Species!="GreenHouse",SoilFungi)
SoilFungi4<-prune_taxa(as.logical(tax_table(SoilFungi13)[,1]=="Fungi"),SoilFungi13)
SoilFungi5<-prune_taxa(as.logical(tax_table(SoilFungi4)[,2]=="Ascomycota"),SoilFungi4)
SoilFungi6<-prune_taxa(as.logical(tax_table(SoilFungi4)[,5]!="NA"),SoilFungi4)
SoilFungi4
SoilFungi7<-prune_taxa(as.logical(tax_table(SoilFungi6)[,5]!="unidentified"),SoilFungi6)
SoilFungi7

###This is For the new analysis for just classified fungi to Genera. This is January 2023 from MA instructions 
SoilFungi8<-prune_taxa(as.logical(tax_table(SoilFungi7)[,6]!="unidentified"),SoilFungi7)
SoilFungi8
otu1<-otu_table(ControlSoilFungi10)
meta<-sample_data(ControlSoilFungi4)
na.omit(meta)
write.csv(otu1,"/Users/rachelbrant/Desktop/otusoil2.csv")
write.csv(meta,"/Users/rachelbrant/Desktop/otucontrolmeta.csv")
sample_names(SoilFungi5)
range(sample_sums(SoilFungi2))
rareSoilFungiNew<-rarefy_even_depth(SoilFungi13,800,rngseed = 995368)
rareSoilFungiNew

#####################richness for just inocula 
RichnessFungiCP<-estimate_richness(rareSoilCP1,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessFungiCP)
write.csv(RichnessFungiCP,"/Users/rachelbrant/Desktop/richnessCPfungi.csv")

RichnessFungiGV<-estimate_richness(rareSoilGV1,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessFungiGV)
write.csv(RichnessFungiGV,"/Users/rachelbrant/Desktop/richnessGVfungi.csv")

RichnessFungiSA<-estimate_richness(rareSoilSA1,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessFungiSA)
write.csv(RichnessFungiSA,"/Users/rachelbrant/Desktop/richnessSAfungi.csv")

RichnessFungiSC<-estimate_richness(rareSoilSC1,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessFungiSC)
write.csv(RichnessFungiSC,"/Users/rachelbrant/Desktop/richnessSCfungi.csv")


RichnessBacCP<-estimate_richness(rareCPphylo,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessBacCP)
write.csv(RichnessBacCP,"/Users/rachelbrant/Desktop/richnessCPbac.csv")

RichnessBacGV<-estimate_richness(rareGVphylo,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessBacGV)
write.csv(RichnessBacGV,"/Users/rachelbrant/Desktop/richnessGVbac.csv")

RichnessBacSA<-estimate_richness(rareSAphylo,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessBacSA)
write.csv(RichnessBacSA,"/Users/rachelbrant/Desktop/richnessSAbac.csv")


RichnessBacSC<-estimate_richness(rareSCphylo,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessBacSC)
write.csv(RichnessBacSC,"/Users/rachelbrant/Desktop/richnessSCbac.csv")

plotrichness<-plot_richness(ControlSoilFungi2,x="Species",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age),breaks=c('Young', 'Middle', 'Old',name = "Soil Age")) + theme_bw() + labs(x="Species",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plot_richness<-plot_richness(ControlSoilFungi2,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age),breaks=c('Young', 'Middle', 'Old',name = "Soil Age")) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8))
plotfungi<-plotrichness+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',name = "Soil Age"))
plot_richness(rareSoil3,x="Species",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
SoilFungiCP1$Soil_Age <- factor(SoilFungiCP1$Species,                                    # Change ordering manually
                  levels = c( "Young", "Middle", "Old"))

ControlRichnessCP<-estimate_richness(ControlRichness,measures=c("Observed","Shannon","InvSimpson"))
head(ControlRichness16s)
plotrichness<-plot_richness(ControlSoilFungi2,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8))
plotrichness2<-plotrichness+scale_fill_manual(values =c("#6699CC","#009E73","#999999"),breaks=c("Young","Middle","Old"))
plotrichness3<-plotrichness2+theme(legend.position = "none")
plotrichness3<-plotrichness3+scale_x_discrete(limits=c("Young","Middle","Old"))
plotrichness3
newSTorder = c("Young","Middle","Old")
plotrichness2$ControlSoilFungi2$Soil_Age <- as.character(plotrichness2$ControlSoilFungi2$Soil_Age)
plotrichness2$ControlSoilFungi2$Soil_Age <- factor(plotrichness2$ControlSoilFungi2$Soil_Age, levels=newSTorder)



library(dplyr)
soil_metaFun3<-data.frame(sample_data(ControlSoilFungi2))
soil_metaFun3$sampleid<-rownames(soil_metaFun3)
soil_metaFun3
soil_richnessControlFungi3<-estimate_richness(ControlSoilFungi2,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessControlFungi3$sampleid<-rownames(soil_richnessControlFungi3)
soil_richnessControl4<-left_join(soil_richnessControlFungi3,soil_metaFun3,"sampleid")
soil_richnessControl4
library(lme4)
m11<-aov(InvSimpson~ Soil_Age, data=soil_richnessControl4)
summary(m11)
TukeyHSD(m11)
library(effectsize)
options(es.use_symbols = TRUE) # get nice symbols when printing! (On Windows, requires R >= 4.2.0)
eta_squared(m11, partial = FALSE)
sample_data(ControlSoilFungi2)


################richness for plant x soil age interaction ### 
library(ggplot2)
library(ggpubr)
RichnessFun<-estimate_richness(SoilFungi3,measures=c("Observed","Shannon","InvSimpson"))
head(RichnessBac)
plotrichnessFun<-plot_richness(SoilFungi3,x="Species",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age),position=position_dodge(1)) + theme_bw() + labs(x="Species",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichnessFun+scale_fill_brewer(palette="GnBu",direction=-1)+ylim(0,5)
plot_richness(rareControlSoilFungi2,x="Species",measures=c("Shannon")) + geom_point(size=5,pch=21,aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))

library(dplyr)
soil_metaFunBLE<-data.frame(sample_data(SoilFungi3))
soil_metaFunBLE$sampleid<-rownames(soil_metaFunBLE)
soil_metaFunBLE
soil_richnessFungiBLE<-estimate_richness(SoilFungi3,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessFungiBLE$sampleid<-rownames(soil_richnessFungiBLE)
soil_richnessFungiBLE1<-left_join(soil_richnessFungiBLE,soil_metaFunBLE,"sampleid")
soil_richnessFungiBLE1
library(lme4)
fungiSoilPlantBLEB<-lm(Shannon~ Species*Soil_Age, data=soil_richnessFungiBLE1)
summary(fungiSoilPlantBLEB)



####################  SHARED AND UNIQUE SEQUENCES BY SITE
library(gplots) #required for Venn diagram. Online Venn diagram tools like Venny are also available. 
ControlSoilFungi2
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(ControlSoilFungi2)$Soil_Age))
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
  phy.sub<-prune_samples(sample_data(rareControlSoilFungi2)$Soil_Age %in% poplist[[k]],rareControlSoilFungi2)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(rareControlSoilFungi2))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)

bob<-otu_table(ControlSoilFungi10)
library(wesanderson)
####### Beta diversity for innoculate soil fungi 
bray_dist <- vegdist(bob, distance="bray",trymax=500)
beta_nmds_bray = metaMDS(bray_dist, distance="bray", k=2)
ordinationLECITS <- ordinate(ControlSoilFungi10, "PCoA", "bray")
beta_div_bray_cleanLECITS <- plot_ordination(ControlSoilFungi10, ord.nmds.bray) 
bray_clean_presenceLECITS <- beta_div_bray_cleanLECITS + geom_point(aes(color=Soil_Age))+stat_ellipse(aes(color= Soil_Age)) + theme_classic() 
bray22<-bray_clean_presenceLECITS+scale_color_manual(values=c("#6699CC","darkorange1","#999999"),breaks=c('Young', 'Middle', 'Old'),name = "Soil Age",labels=c('Young'="Young", 'Middle'="Intermediate", 'Old'="Old"))+ theme(axis.text=element_text(size=8),axis.title = element_text(size=8))
bray23<-bray22+ylim(-0.6,0.62)+xlim(-1,1.5)
bray24<-bray23+theme(legend.title=element_text(size=14),legend.text = element_text(size =12))
bray24
bray50<-bray24+ylim(-0.6,0.62)+xlim(-1,1.5)
bray52<-bray50+annotate("text", x=0.75, y=-0.5, label= "stress = 0.04519",size=5)
bray53<-bray52+theme(axis.text.x = element_text(size=13),
              axis.text.y = element_text(size=13))+theme(axis.title=element_text(size=14))
bray53
####beta diversity with all plants fungi 
ord.nmds.brayplants <- ordinate(rareSoilFungiNew, method="NMDS",k=6, distance="bray",trymax=500)
ord.nmds.brayplants
ordinationLECITSplants <- ordinate(SoilFungi13, "PCoA", "bray")
beta_div_bray_cleanLECITSplants <- plot_ordination(rareSoilFungiNew, ord.nmds.brayplants, type= "samples", color= "Soil") + geom_point(size=3)
bray_clean_presenceLECITSplants <- beta_div_bray_cleanLECITSplants + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECITSplants+scale_color_brewer(palette = "RdBu", direction=-1,breaks=c('CP', 'GV', 'SA',"SC",name = "Species"))
bray_clean_presenceLECITSplants

SoilFungi<-prune_samples(sample_data(physeqLECFungi12)$Species!="Control",physeqLECFungi12)
SoilFungi2<-prune_samples(sample_data(SoilFungi)$Soil_type!="Control",SoilFungi)
SoilFungi2<-prune_samples(sample_data(SoilFungi2)$Species!="GreenHouse",SoilFungi2)
SoilFungiCP11<-prune_samples(sample_data(SoilFungi4)$Species=="CP",SoilFungi4)
SoilFungiGV11<-prune_samples(sample_data(SoilFungi4)$Species=="GV",SoilFungi4)
SoilFungiSA11<-prune_samples(sample_data(SoilFungi4)$Species=="SA",SoilFungi4)
SoilFungiSC11<-prune_samples(sample_data(SoilFungi4)$Species=="SC",SoilFungi4)
SoilFungiOld<-prune_samples(sample_data(SoilFungi4)$Soil_Age=="Old",SoilFungi4)
SoilFungiMiddle<-prune_samples(sample_data(SoilFungi4)$Soil_Age=="Middle",SoilFungi4)
SoilFungiYoung<-prune_samples(sample_data(SoilFungi4)$Soil_Age=="Young",SoilFungi4)
library(ggplot2)

rareSoilCP1<-rarefy_even_depth(SoilFungiCP11,1000,rngseed = 354300)
rareSoilGV1<-rarefy_even_depth(SoilFungiGV11,1000,rngseed = 354350)
rareSoilSA1<-rarefy_even_depth(SoilFungiSA11,700,rngseed = 354360)
rareSoilSC1<-rarefy_even_depth(SoilFungiSC11,1000,rngseed = 354377)

####beta diversity with CP across sites fungi 
ord.nmds.brayCP <- ordinate(SoilFungiCP1, method="NMDS",k=7, distance="bray",trymax=500)
ordinationLECITSCP <- ordinate(rareSoilCP1, "PCoA", "bray")
beta_div_bray_cleanLECITSCP <- plot_ordination(rareSoilCP1, ordinationLECITSCP, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECITSCP <- beta_div_bray_cleanLECITSCP  + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type"))+ theme_classic()
BrayCPfun<-bray_clean_presenceLECITSCP+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BrayCPfun+ylim(-1,1)+xlim(-1,1)

####beta diversity with GV across sites fungi 
ord.nmds.brayGV<-ordinate(rareSoilGV1, method="NMDS",k=6, distance="bray",trymax=500)
ordinationLECITSGV <- ordinate(rareSoilGV1, "PCoA", "bray")
beta_div_bray_cleanLECITSGV <- plot_ordination(rareSoilGV1, ordinationLECITSGV, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECITSGV <- beta_div_bray_cleanLECITSGV + stat_ellipse() +ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYGVfun<-bray_clean_presenceLECITSGV+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYGVfun+ylim(-1,1)+xlim(-1,1)

####beta diversity with SA across sites fungi 
ord.nmds.braySA <- ordinate(rareSoilSA1, method="NMDS",k=8, distance="bray",trymax=500)
ordinationLECITSSA <- ordinate(rareSoilSA1, "PCoA", "bray")
beta_div_bray_cleanLECITSSA <- plot_ordination(rareSoilSA1, ordinationLECITSSA, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECITSSA <- beta_div_bray_cleanLECITSSA + stat_ellipse() + ggtitle("Bray Curtis") + guides(color = guide_legend(title = "Soil Type")) + theme_classic()
BRAYSAfun<-bray_clean_presenceLECITSSA+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
BRAYSAfun+ylim(-1,1)+xlim(-1,1)
####beta diversity with SC across sites fungi 
ord.nmds.braySC <- ordinate(rareSoilSC1, method="NMDS",k=10, distance="bray",trymax=500)
ordinationLECITSSC <- ordinate(rareSoilSC1, "PCoA", "bray")
beta_div_bray_cleanLECITSSC <- plot_ordination(rareSoilSC1, ordinationLECITSSC, type= "samples", color= "Soil_Age") + geom_point(size=3)
bray_clean_presenceLECITSSC <- beta_div_bray_cleanLECITSSC  + stat_ellipse() + ggtitle("Bray Curtis")  + guides(color = guide_legend(title = "Soil Type"))  + theme_classic()
BRAYSCfun<-bray_clean_presenceLECITSSC+scale_color_manual(values=c("black","steelblue","darkgrey","purple"),breaks=c('Young', 'Middle', 'Old',"Control", name = "Soil Age"))
BRAYSCfun+ylim(-1,1)+xlim(-1,1)
###beta diversity with old across species of plant 
ord.nmds.brayOld <- ordinate(SoilFungiOld, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECITSOld <- ordinate(SoilFungiOld, "PCoA", "bray")
beta_div_bray_cleanLECITSOld <- plot_ordination(SoilFungiOld, ordinationLECITSOld, type= "samples", color= "Species") + geom_point(size=3)
bray_clean_presenceLECITSOld <- beta_div_bray_cleanLECITSOld  + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECITSOld+scale_color_brewer(palette = "Set1", direction=-1,breaks=c('CP', 'GV', 'SA',"SC",name = "Plants"))

ord.nmds.brayYoung <- ordinate(SoilFungiYoung, method="NMDS",k=6, distance="bray",trymax=50)
ordinationLECITSYoung <- ordinate(SoilFungiYoung, "PCoA", "bray")
beta_div_bray_cleanLECITSYoung<- plot_ordination(SoilFungiYoung, ordinationLECITSYoung, type= "samples", color= "Species") + geom_point(size=3)
bray_clean_presenceLECITSYoung <- beta_div_bray_cleanLECITSYoung  + ggtitle("Bray Curtis")  + theme_classic()
bray_clean_presenceLECITSYoung+scale_color_brewer(palette = "Set1", direction=-1,breaks=c('CP', 'GV', 'SA',"SC",name = "Plants"))


###species venns for fungi CP
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(SoilFungiCP)$Soil_Age))
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
  phy.sub<-prune_samples(sample_data(SoilFungiCP)$Soil_Age %in% poplist[[k]],SoilFungiCP)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SoilFungiCP))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)

###species venns for fungi GV
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(SoilFungiGV)$Soil_Age))
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
  phy.sub<-prune_samples(sample_data(SoilFungiGV)$Soil_Age %in% poplist[[k]],SoilFungiGV)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SoilFungiGV))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)


###Species venns for fungi SA
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(SoilFungiSA)$Soil_Age))
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
  phy.sub<-prune_samples(sample_data(SoilFungiSA)$Soil_Age %in% poplist[[k]],SoilFungiSA)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SoilFungiSA))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)


###Species venns for fungi SC
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(SoilFungiSC)$Soil_Age))
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
  phy.sub<-prune_samples(sample_data(SoilFungiSC)$Soil_Age %in% poplist[[k]],SoilFungiSC)
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SoilFungiSC))
  ##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
plot<-venn(venn.list)



####Fixed plot
soil_richnessControl2$Soil_Age <- factor(soil_richnessControl2$Soil_Age , levels=c("Young", "Middle","Old"))
plot1<-ggplot(soil_richnessControl2) +
  aes(x = Soil_Age, y = Shannon,fill=Soil_Age) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette="RdBu",direction=1) +theme_minimal()+xlab("Soil Age")+ylab("Prokaryotic Shannon Diversity")
plot1+ theme(legend.position = "right")






#PERMANOVA
library(vegan)
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


OTU<-otu_table(SoilFungiCP11)
OTU1<-as.table(OTU)
class(OTU1)
dim(OTU1)
#Convert Sample Data to     
soil.sr22S<-as(sample_data(SoilFungiCP11),"matrix")
soil.sr22S
###### NMDS Ordination
soil.nmds11<-metaMDS(soil.vr22S,distance="bray",trymax = 600,halfchange = TRUE)   
soil.adonisr1322S<-adonis2(soil.vr22S ~ factor(Soil_Age)*factor(Species) ,data=soil.sr22S,permutations=999,method="bray")
soil.adonisr1322S
soil.nmds11
betadisper(soil.nmds11)


##Permdisp
otu=read.table(,header=TRUE)
OTU2<-otu_table(rareControlSoilFungi24)
head(OTU)
rdp2=OTU2[,ncol(OTU2)]
names(rdp2)=row.names(OTU2)
otu=OTU2[,-ncol(OTU2)]
dim(otu)
dim(OTU2)


class(soil.sr22S)
map<-as.table(soil.sr22S)
dim(map)
View(map)
braycurtis.d2=vegdist(t(OTU2), method="bray")
braycurtis.d2
sorenson.d=vegdist(t(otu), method="bray", binary=TRUE)
head(sorenson.d) 
braycurtis.mds <- metaMDS(braycurtis.d)
plot(braycurtis.mds, type="t")
b2=betadisper(braycurtis.d2, group=map[,"Soil_Age"], type="median")
b2
b.perm2=permutest(b2, group=map[,"Soil_Age"], type="median", permutations=999, pairwise=TRUE)
b.perm2



##For control soil too
braycurtis.d=vegdist(t(OTU1), method="bray")
braycurtis.d
OTU<-otu_table(SoilFungiCP11)
head(OTU)
rdp=OTU1[,ncol(OTU1)]
names(rdp)=row.names(OTU1)
otu=OTU1[,-ncol(OTU1)]
dim(otu)
dim(OTU1)
braycurtis.d3=vegdist(t(OTU1), method="bray")
braycurtis.d3
sorenson.d=vegdist(t(otu), method="bray", binary=TRUE)
head(sorenson.d) 
braycurtis.mds3 <- metaMDS(braycurtis.d3)
braycurtis.mds3
plot(braycurtis.mds3, type="t")
b3=betadisper(braycurtis.d3, group=map[,"Soil_Age"], type="median")
b3
b.perm2=permutest(b3, group=map[,"Soil_Age"], type="median", permutations=999, pairwise=TRUE)
b.perm2



taxo<-function(resultsobject,physeqobject,alpha){
  sigtab<-resultsobject[which(resultsobject$padj<alpha),]
  sigtab<- cbind(as(sigtab, "data.frame"), as(tax_table(physeqobject)[rownames(sigtab), ], "matrix"))
  colnames(sigtab)[7:12]<-c("Kingdom","Phylum","Class","Order","Family","Genus")
  return(sigtab)
}          


###### Function to parse significant data from DESeq2 results
deseqplot_data<-function(sigtab){
# Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
# Copy Across Genus Labels and Fill in Any Unassigned
  sigtab$Genus.long<-as.character(sigtab$Genus)
  sigtab$Genus.long[grep("unclassified",sigtab$Genus)] <-paste0("[",as.character(sigtab$Family[grep("unclassified",sigtab$Genus)]),"]")
  
# Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus.long, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus.long = factor(as.character(sigtab$Genus.long), levels=names(x))
  return(sigtab)
}
###Aldex2
##Diff abundance
library(ALDEx2)
library(vegan)
library(dplyr)
install.packages("DESeq2")
library(DESeq2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2",force=TRUE)
gp_soilocean3<-subset_samples(ControlSoilFungi4,Soil_Age %in% c("Young","Middle"))
gp_soilocean3
otu_table(gp_soilocean3) <- otu_table(gp_soilocean3) + 1
sample_sums(gp_soilocean3)
gpmod<-phyloseq_to_deseq2(gp_soilocean3, ~ Soil_Age)
gp.deseq = DESeq(gpmod, test="Wald", fitType="mean")
resultsNames(gp.deseq)
gp_results<-results(gp.deseq,contrast=c("Soil_Age","Young","Middle"))
gp_results_taxa<-taxo(gp_results,gp_soilocean3,0.05) #set threshold signfiicance 
gp_plot_data<-deseqplot_data(gp_results_taxa)
head(gp_plot_data)

gp_plot<-ggplot(gp_plot_data,aes(x=Genus.long, y=log2FoldChange)) + geom_point(shape=21,size=6,aes(fill=Phylum)) + theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=12)) +  scale_fill_brewer(palette="Set2")+geom_hline(yintercept = 0,linetype="dashed") + labs(x="Genus") + scale_shape_manual(values=c(21,22,23,24,25))
gp_plot



gp_soilocean3<-subset_samples(SoilFungi4,Soil_Age %in% c("Young","Middle"))
gp_soilocean3
otu_table(gp_soilocean3) <- otu_table(gp_soilocean3) + 1
gp_soilocean3<-subset_taxa(gp_soilocean3, Kingdom=="Fungi")
sample_sums(gp_soilocean3)
gpmod1<-phyloseq_to_deseq2(gp_soilocean3, ~ Soil_Age)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

gp_soilocean_v<-t(vegan_otu(gp_soilocean3))
gp_soilocean_sample<-sample_data(gp_soilocean3)$Soil_Age
gp_soilocean_aldexmodel<- aldex(gp_soilocean_v, gp_soilocean_sample, mc.samples=128, test="t", effect=TRUE,
                                include.sample.summary=TRUE, denom="all", verbose=FALSE)
write.csv(gp_soilocean_aldexmodel,"/Users/rachelbrant/Desktop/AldexFungi.csv")
aldex.plot(gp_soilocean_aldexmodel, type="MW", test="welch", xlab="Dispersion",ylab="Difference")
conds <- c(rep("Young", 6), rep("Old",7))
library(DESeq2)
gp.deseq = DESeq(gpmod1, test="Wald", fitType="mean")
resultsNames(gp.deseq)
gp_results<-results(gp.deseq,contrast=c("Soil_Age","Young","Middle"))
gp_results_taxa<-taxo(gp_results,gp_soilocean3,0.05) #set threshold signficance 
gp_plot_data<-deseqplot_data(gp_results_taxa)
gp_plot_data
head(gp_plot_data)
library(wesanderson)
write.csv(gp_plot_data,"/Users/rachelbrant/Desktop/deseqoutputYOnew.csv")
gp_plotYO<-ggplot(gp_plot_data,aes(x=Genus, y=log2FoldChange)) + geom_point(shape=21,size=6,aes(fill=Family)) + theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=12)) +  scale_fill_manual(values = c("#E5F5F9", "#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C" ,"#A6D854",
                                                                                                               "#D4B9DA", "#6A51A3" ,"#7F0000", "#D9D9D9", "#FFF7BC" ,"#000000", "#F0F0F0"))+geom_hline(yintercept = 0,linetype="dashed") + labs(x="Family") + scale_shape_manual(values=c(21,22,23,24,25))
gp_plotYO

gp_plotYM<-ggplot(gp_plot_data,aes(x=Genus, y=log2FoldChange)) + geom_point(shape=21,size=6,aes(fill=Family)) + theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=12)) +  scale_fill_manual(values=c("#E5F5F9", "#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C" ,"#A6D854",
                                                                                                             "#D4B9DA", "#6A51A3" ,"#7F0000", "#D9D9D9", "#FFF7BC" ,"#000000", "#F0F0F0", "#C7EAE5", "#003C30" ,"#F16913",
                                                                                                             "#FFF7FB", "#8C6BB1" ,"#C7E9B4"))+geom_hline(yintercept = 0,linetype="dashed") + labs(x="Family") + scale_shape_manual(values=c(1,26))
gp_plotYM


library(RColorBrewer)
n <- 30
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
col <- sample(col_vec, n)
area <- rep(1,n)
pie(area, col = col)

####Redo of venn by soil type
###Species venns for fungi SA
#Empty lists - we have three groups so we need three empty lists
venn.list<-rep(list(NA),5)
poplist<-rep(list(NA),5)
pops<-as.character(unique(sample_data(SoilFungiOld)$Species))
pops
#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]
poplist[[4]]<-pops[4]
poplist[[5]]<-pops[5]
#Name the Groups in the Vennlist
names(venn.list)<-pops



#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:5){
  
#Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(SoilFungiOld)$Species %in% poplist[[k]],SoilFungiOld)
#Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  phy.sub.keep
  rownames(otu_table(SoilFungiOld))
##prune down to only those ASVS
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
#Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
library(gplots)
plot<-venn(venn.list)
##This plot is really complicated. Should we include it? Not sure. Ask 


library(ggplot2) # to create data visualizations 
install.packages("RColorBrewer")
library(ggthemes) # use themes to clean up the data visualizations 
library(RColorBrewer)
library(cowplot) #load in library to create plot grids

#multi-panel side by side graphs
library(cowplot)
multi_plot<- ggarrange(plot2+bray22+ braybac, #plots that are going to be included in this multipanel figure
                       labels = c("A", "B", "C","D"), #labels given each panel 
                       ncol = 2, nrow = 2, #adjust plot space 
                       common.legend = T)
multi <- (plot2+ braybac+bray22) + 
  plot_layout(widths = c(3,1))+
  plot_annotation(tag_levels = 'I') #add figure labels
multi #view multi-panel figure 
library(patchwork)
install.packages("patchwork")
library(patchwork)
design <- "
  122
  153
  443
"
((braybac2|bray222)/(plot22|plotFUN)+plot_layout(nrow=2,byrow=TRUE,heights =c(10,10,6,6)))
plot5<-plot4+ggtitle("A")+theme(axis.title=element_text(size=10,face="bold"))
plot5
plotrichness4<-plotrichness3+ggtitle("A")+theme(plot.title=element_text(size=14,face="bold"))
plotrichness4
braybac2<-braybac5+ggtitle("A")+theme(plot.title=element_text(size=10,face="bold"))
braybac2
bray222<-bray53+ggtitle("B")+theme(plot.title=element_text(size=10,face="bold"))
bray222
plot22<-plotSimper333+ggtitle("C")+theme(plot.title=element_text(size=10,face="bold"))
plot22
plotFUN<-plotSimper5+ggtitle("D")+theme(plot.title=element_text(size=10,face="bold"))
plotFUN
layout <- c(
  area(1,1,1,3),
  area(2,1),
  area(2,3),
  area(3,1),
  area(3,3))
plot(layout)
bray222


### Alpha per Sp
library(dplyr)
library(ggplot2)
soil_metafun<-data.frame(sample_data(rareSoilFungiNew))
soil_metafun$sampleid<-rownames(soil_metafun)
soil_metafun
soil_richnessfun<-estimate_richness(rareSoilFungiNew,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessfun$sampleid<-rownames(soil_richnessfun)
soil_richnessfun<-left_join(soil_richnessfun,soil_metafun,"sampleid")
soil_richnessfun
plotrichnessCpfun<-plot_richness(rareSoilCP1,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichnessCpfun+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
library(lme4)
m11good<-lm(Shannon ~ Soil_Age*Species, data=soil_richnessfun)
summary(m11good)

soil_metaGVfun<-data.frame(sample_data(SoilFungiGV11))
soil_metaCPfun$sampleid<-rownames(soil_metaGVfun)
soil_metaGVfun
soil_richnessGVfun<-estimate_richness(SoilFungiGV11,measures=c("Observed","Shannon","InvSimpson"))

#Add Richness Onto Our Metadata
soil_richnessGVfun$sampleid<-rownames(soil_richnessGVfun)
soil_richnessGVfun<-left_join(soil_richnessGVfun,soil_metaGVfun,"sampleid")
soil_richnessGVfun
plotrichnessGVfun<-plot_richness(rareSoilGV1,x="Soil_Age",measures=c("Shannon")) + geom_boxplot(aes(fill=Soil_Age)) + theme_bw() + labs(x="Soil_Age",y="Alpha Diversity (Shannon)") + theme(axis.text=element_text(size=12),axis.title = element_text(size=12))
plotrichnessGVfun+scale_fill_brewer(palette="GnBu",direction=-1,breaks=c('Young', 'Middle', 'Old',"Control",name = "Soil Age"))
library(lme4)
m111<-aov(Shannon ~ Soil_Age, data=soil_richnessGVfun)
summary(m111)
