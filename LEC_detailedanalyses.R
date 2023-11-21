rareSoilFungiCP<-rarefy_even_depth(SoilFungiCP,2000,rngseed = 354390)

soil_Fun<-data.frame(sample_data(rareSoilFungiCP))
soil_Fun$sampleid<-rownames(soil_Fun)
soil_Fun
soil_richnessCP<-estimate_richness(rarSoilFungiCP,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessCP$sampleid<-rownames(soil_richnessCP)
soil_richnessCP1<-left_join(soil_richnessCP,soil_Fun,"sampleid")
soil_richnessCP1
library(lme4)
m11<-lm(Shannon~Soil_Age, data=soil_richnessCP1)
summary(m11)


rareSoilFungiGV<-rarefy_even_depth(SoilFungiGV,2000,rngseed = 354390)

soil_Fun<-data.frame(sample_data(rareSoilFungiGV))
soil_Fun$sampleid<-rownames(soil_Fun)
soil_Fun
soil_richnessGV<-estimate_richness(rareSoilFungiGV,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessGV$sampleid<-rownames(soil_richnessGV)
soil_richnessGV1<-left_join(soil_richnessGV,soil_Fun,"sampleid")
soil_richnessGV1
library(lme4)
m11<-lm(Shannon~Soil_Age, data=soil_richnessGV1)
summary(m11)


plot<-ggplot(soil_richnessControl4) +
  aes(x = Soil_Age, y = Shannon, fill = Soil_Age) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal()
plot


rareSoilFungiSA<-rarefy_even_depth(SoilFungiSA,700,rngseed = 354390)

soil_FunSA<-data.frame(sample_data(rareSoilFungiSA))
soil_FunSA$sampleid<-rownames(soil_FunSA)
soil_FunSA
soil_richnessSA<-estimate_richness(rareSoilFungiSA,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessSA$sampleid<-rownames(soil_richnessSA)
soil_richnessSA1<-left_join(soil_richnessSA,soil_FunSA,"sampleid")
soil_richnessSA1
library(lme4)
m11<-lm(Shannon~Soil_Age, data=soil_richnessSA1)
summary(m11)


rareSoilFungiSC<-rarefy_even_depth(SoilFungiSC,2000,rngseed = 354390)

soil_FunSC<-data.frame(sample_data(rareSoilFungiSC))
soil_FunSC$sampleid<-rownames(soil_FunSC)
soil_FunSC
soil_richnessSC<-estimate_richness(rareSoilFungiSC,measures=c("Observed","Shannon","InvSimpson"))
#Add Richness Onto Our Metadata
soil_richnessSC$sampleid<-rownames(soil_richnessSC)
soil_richnessSC1<-left_join(soil_richnessSC,soil_FunSC,"sampleid")
soil_richnessSC1
library(lme4)
m11<-lm(Shannon~Soil_Age, data=soil_richnessSC1)
summary(m11)



############# indicator species
library(labdsv)
library(vegan)
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

ControlSoilFungiYoung<-prune_samples(sample_data(ControlSoilFungi4)$Soil_Age=="Young",ControlSoilFungi4)
ControlSoilFungiMid<-prune_samples(sample_data(ControlSoilFungi4)$Soil_Age=="Middle",ControlSoilFungi4)
ControlSoilFungiOld<-prune_samples(sample_data(ControlSoilFungi4)$Soil_Age=="Old",ControlSoilFungi4)

otuZ<-vegan_otu(ControlSoilFungi4)
SoilFungiCP11<-prune_samples(sample_data(SoilFungi13)$Species=="CP",SoilFungi13)
SoilFungiGV11<-prune_samples(sample_data(SoilFungi13)$Species=="GV",SoilFungi3)
SoilFungiSA11<-prune_samples(sample_data(SoilFungi13)$Species=="SA",SoilFungi3)
SoilFungiSC11<-prune_samples(sample_data(SoilFungi13)$Species=="SC",SoilFungi3)
otuCP<-vegan_otu(SoilFungiCP11)
metaCP<-sample_data(SoilFungiCP11)
otuGV<-vegan_otu(SoilFungiGV11)
metaGV<-sample_data(SoilFungiGV11)
otuSA<-vegan_otu(SoilFungiSA11)
metaSA<-sample_data(SoilFungiSA11)
otuSC<-vegan_otu(SoilFungiSC11)
metaSC<-sample_data(SoilFungiSC11)



indicator_matrixControlY<- vegan_otu(ControlSoilFungiYoung)
indicator_matrixControlY
write.csv(otuZ,"/Users/rachelbrant/Desktop/otuZ.csv")

indicator_matrixControlM<- vegan_otu(ControlSoilFungiMid)
indicator_matrixControlM
write.csv(indicator_matrixControlM,"/Users/rachelbrant/Desktop/NewmatrixMid.csv")

indicator_matrixControlO<- vegan_otu(ControlSoilFungiOld)
indicator_matrixControlO
write.csv(indicator_matrixControlO,"/Users/rachelbrant/Desktop/NewmatrixOld.csv")





indicatormatrixGV1<-vegan_otu(SoilFungiGV1)
write.csv(indicatormatrixGV1,"/Users/rachelbrant/Desktop/otutableGVNEW.csv")

indicatormatrixSA1<-vegan_otu(SoilFungiSA1)
write.csv(indicatormatrixSA1,"/Users/rachelbrant/Desktop/otutableSANEW.csv")

indicatormatrixSC1<-vegan_otu(SoilFungiSC1)
write.csv(indicatormatrixSC1,"/Users/rachelbrant/Desktop/otutableSCNEW.csv")


indicatormatrixCP1<-vegan_otu(SoilFungiCP1)
write.csv(indicatormatrixCP1,"/Users/rachelbrant/Desktop/otutableCPNEW.csv")

indiactormatrixCPbac1<-vegan_otu(CPphylo)
write.csv(indiactormatrixCPbac1,"/Users/rachelbrant/Desktop/otutableCPbacCC.csv")

indiactormatrixGVbac<-vegan_otu(GVphylo)
write.csv(indiactormatrixGVbac,"/Users/rachelbrant/Desktop/otutableGVbacCC.csv")

indiactormatrixSAbac<-vegan_otu(SAphylo)
write.csv(indiactormatrixSAbac,"/Users/rachelbrant/Desktop/otutableSAbacCC.csv")


indiactormatrixSCbac<-vegan_otu(SCphylo)
write.csv(indiactormatrixSCbac,"/Users/rachelbrant/Desktop/otutableSCbacCC.csv")
library(phyloseq)
indicatormatrixcontrolbac<-vegan_otu(ControlSoil1)
write.csv(indicatormatrixcontrolbac,"/Users/rachelbrant/Desktop/outtableControlbacNEW.csv")


veg_numeric12<-with(sample_data(ControlSoilFungi4),ifelse(Soil_Age=="Young",1,2))
veg_numeric12
indicator_matrix12
veg_indic13<-indval(indicator_matrix12,veg_numeric12)
summary(veg_indic13)
veg_indic_output1<-data.frame(indval=veg_indic13$indcls,pval=veg_indic13$pval,cluster=veg_indic13$maxcls)
veg_indic_output1$taxon<-gsub("^X","",rownames(veg_indic_output1))
veg_indic_significant1<-subset(veg_indic_output1,indval>=0.5 & pval <=0.05)
head(veg_indic_significant1)
ps_rare_tax<-data.frame(tax_table(rareSoilFungiCP))
ps_rare_tax$taxon<-rownames(ps_rare_tax)
veg_indic_significant<-left_join(veg_indic_significant,ps_rare_tax,"taxon")
head(veg_indic_significant)
write.csv(otuGeneraonly,"/Users/rachelbrant/Desktop/OTUgeneraonly.csv")



detach(otutableUGH)
detach(Newmatrix2)
detach(otutableCPNEW)
attach(otutableGVNEW)
detach(otutableSA1)
detach(otutableSANEW)
detach(otutableCPbacCC)
detach(otutableSCbacCC)
detach(outtableControlbacNEW)
detach(Numberleaf_data)
###control indicator species analyses to get core taxa 
invO = multipatt(abundO,timeO, func = "IndVal", duleg = TRUE,control = how(nperm=999))
abundO = otutableSANEW[,3:ncol(otutableSANEW)]
timeO = otutableSANEW$Age
abundO
timeO
fresh<-na.omit(otutableSAbacCC)
summary(invO,alpha=0.05)

invN$sign



install.packages("indicspecies")
library(indicspecies)

write.csv(wetland,"/Users/rachelbrant/Desktop/wetland.csv")




##dominant taxa
install.packages("compare") 
library(compare)
library(ggpubr)
top4taxaFUNGICONTROL<-top_taxa(rareControlSoilFungi2, 4)
top4taxaFUNGIALL<-top_taxa(rareSoilFungi2, 4)
top4taxaBACCONTROL<-top_taxa(rareControlSoil, 4)
top4taxaBACALL<-top_taxa(rareSoil3, 4)
top4taxaBACCONTRO + stat_compare_means(
  method = "wilcox.test")
clean_rare_all_compare_taxa <- aggregate_taxa(rareControlSoil,"Phylum")
clean_rare_all_compare_taxa
##dominant taxa for whole data clean_rare
dom.tax.allITS <- dominant_taxa(rareControlSoil,level = "Family", group="Soil_Age")
head(dom.tax.allITS$dominant_overview)
print(dom.tax.allITS$dominant_overview)

##subsetsamplesmore
SoilFungiCPY<-prune_samples(sample_data(SoilFungiCP)$Soil_Age=="Young",SoilFungiCP)
SoilFungiCPM<-prune_samples(sample_data(SoilFungiCP)$Soil_Age=="Middle",SoilFungiCP)
SoilFungiCPO<-prune_samples(sample_data(SoilFungiCP)$Soil_Age=="Old",SoilFungiCP)














plotabundance<-prune_samples(sample_sums(SoilFungi5)>0,SoilFungi5)
plotabundance <- plotabundance %>%
  aggregate_taxa(level = "Family") %>%  
  microbiome::transform(transform = "compositional")
plotabundance
physeq_sample_plot2 <- plotabundance %>%
  plot_composition(average_by = "Soil_Age")
library(wesanderson)
physeq_sample_plot2+scale_fill_brewer(palette="GnBu")+theme_bw()






library(ALDEx2)
library(vegan)
#Strip Out Matrix using functiion above 
gp_soilocean_v<-t(vegan_otu(ControlSoilFungi4))
gp_soilocean_v
type(gp_soilocean_v)
#Strip Out a vector of covariates from the phyloseq object (needs to be two levels)
gp_soilocean_sample<-sample_data(gp_soilocean)$SampleType

#Run ALDEX model
gp_soilocean_aldexmodel<- aldex(gp_soilocean_v, gp_soilocean_sample, mc.samples=128, test="t", effect=TRUE)
                                
### indicator species analysis
library(indicspecies)
wetkm = kmeans(pc, centers=3)
groupskm = pc$Age
groupskm
attach(pc)
class(pc)
pc2<-as.matrix(pc)
head(pc2)
abund = pc[,3:ncol(pc)]
time = pc$Age
inv = multipatt(abund, groupskm, control = how(nperm=9999))
summary(inv)
meta
attach(LEC.metadata22)

####SIMPER NOV 22

vibes<-simper(otuZ, meta$Soil_Age, permutations=100)
summary(vibes,ordered=TRUE)
summary(vibes,ordered=TRUE)$Young_Old

bacControl<-vegan_otu(ControlSoil)
metabac<-sample_data(ControlSoil)

bacOld<-simper(bacControl,metabac$Soil_Age, permutations=100)
summary(bacOld)
summary(bacOld,ordered=TRUE)$Middle_Young

bacMid<-vegan_otu(MidSoil)
metabacMid<-sample_data(MidSoil)

bacYoung<-vegan_otu(YoungSoil)
metabacYoung<-sample_data(YoungSoil)

bacOld11<-vegan_otu(OldSoil)
metabacOld1<-sample_data(OldSoil)

vibesbacYoung<-simper(bacYoung,metabacYoung$Species, permutations=100)
summary(vibesbacYoung)
summary(vibesbacYoung,ordered=TRUE)$CP_GV

vibesbacMid<-simper(bacMid,metabacMid$Species, permutations=100)
summary(vibesbacMid)
summary(vibesbacMid,ordered=TRUE)$SA_SC

vibesbacOld1<-simper(bacOld11,metabacOld1$Species, permutations=100)
summary(vibesbacOld1)
summary(vibesbacOld1,ordered=TRUE)$SA_SC

##by plant
otuCP<-vegan_otu(SoilFungiCP11)
metaCP<-sample_data(SoilFungiCP11)
otuGV<-vegan_otu(SoilFungiGV11)
metaGV<-sample_data(SoilFungiGV11)
otuSA<-vegan_otu(SoilFungiSA11)
metaSA<-sample_data(SoilFungiSA11)
otuSC<-vegan_otu(SoilFungiSC11)
metaSC<-sample_data(SoilFungiSC11)
metaCP

otuOld<-vegan_otu(SoilFungiOld)
metaOld<-sample_data(SoilFungiOld)

otuMid<-vegan_otu(SoilFungiMiddle)
metaMid<-sample_data(SoilFungiMiddle)


otuGeneraonly<-vegan_otu(ControlSoilFungi6)
metaGeneraonly
otuGeneraonly

install.packages("remotes")
remotes::install_github("vmikk/metagMisc",force=TRUE)
library(metagMisc)
otu_table(ControlSoilFungi6) <- otu_table(ControlSoilFungi6) + 1
function1<-phyloseq_average(ControlSoilFungi7, avg_type = "aldex", acomp_zero_impute = NULL,
                 aldex_samples = 128, aldex_denom = "all", group = "Soil_Age",
                 drop_group_zero = FALSE)

library(ALDEx2)
gp_soilocean<-subset_samples(ControlSoilFungi7,Soil_Age %in% c("Young","Old"))
gp_soilocean
#Strip Out Matrix using functiion above 
gp_Genera2<-t(vegan_otu(gp_soilocean))
gp_Genera2
#Gap rip Out a vector of covariates from the phyloseq object (needs to be two levels)
gp_Genera_sample<-sample_data(gp_soilocean)

#Run ALDEX model
gp_Genera_aldexmodel<-aldex(gp_Genera2, gp_Genera_sample, mc.samples=128, test="t", effect=TRUE,
                                include.sample.summary=TRUE, denom="all", verbose=FALSE)

##inocula genera only simper 
metaGeneraonly<-sample_data(ControlSoilFungi9)
otuGenera1<-(vegan_otu(ControlSoilFungi9))
vibesGenera<-simper(otuGenera1,metaGeneraonly$Soil_Age,permutations=100)
summary(vibesGenera,ordered=TRUE)$Young_Old
summary(vibesGenera,ordered=TRUE)$Middle_Young
otu_table(ControlSoilFungi9) <- otu_table(ControlSoilFungi9) + 1
otu_table(ControlSoilFungi9)


#inocula plant genera only simper 
metaGeneraonly7<-sample_data(SoilFungi7)
metaGeneraonly7
otuGenera7<-(vegan_otu(SoilFungi7))
otuGenera7
vibesGenera7<-simper(otuGenera7,metaGeneraonly7$Soil_Age,permutations=100)
summary(vibesGenera7,ordered=TRUE)$Young_Old
summary(vibesGenera7,ordered=TRUE)$Middle_Young
summary(vibesGenera7,ordered=TRUE)$Middle_Old
otu_table(SoilFungi7) <- otu_table(SoilFungi7) + 1
otu_table(SoilFungi7)


otu_table(ControlSoilFungi9) <- otu_table(ControlSoilFungi9) + 1
otu_table(ControlSoilFungi9)
write.csv(otuGenera1,"/Users/rachelbrant/Desktop/otugeneraonly1.csv")
attach(otugeneraonly1)
class(otugeneraonly1)
class(otuGenera1)
otu<-as.matrix(otugeneraonly1)
class(otu)
class(metaGeneraonly)
otuSoil<-vegan_otu(SoilFungi4)
metaSoil<-sample_data(SoilFungi4)

vibesSoilNew<-simper(otuOld,metaOld$Species, permutations=100)
summary(vibesSoilNew)
summary(vibesSoilNew,ordered=TRUE)$SA_SC

otuSoilbac<-vegan_otu(Soil4)
metaSoilbac<-sample_data(Soil4)
vibesSoilBac<-simper(otuSoilbac,metaSoilbac$Soil_Age,permutations=100)
summary(vibesSoilBac)
summary(vibesSoilBac,ordered=TRUE)$Middle_Old
taxtableSoilBac<-tax_table(Soil4)
write.csv(taxtableSoilBac,"/Users/rachelbrant/Desktop/taxtablebacAll.csv")


taxtableSoilFung<-tax_table(SoilFungi4)
write.csv(taxtableSoilFung,"/Users/rachelbrant/Desktop/taxtableFungiAllsoil.csv")
vibesMid<-simper(otuMid,metaMid$Species, permutations=100)
summary(vibesMid)
summary(vibesMid,ordered=TRUE)$GV_SC

vibesYoung<-simper(otuYoung,metaYoung$Species, permutations=100)
summary(vibesYoung)
summary(vibesYoung,ordered=TRUE)$SA_SC
citation("vegan")
vibesCP<-simper(otuCP, metaCP$Soil_Age, permutations=100)
summary(vibesCP)
summary(vibesCP,ordered=TRUE)$Control_Old
summary(vibesCP,ordered=TRUE)

vibesGV<-simper(otuGV, metaGV$Soil_Age, permutations=100)
summary(vibesGV)
summary(vibesGV,ordered=TRUE)$Control_Old
summary(vibesGV,ordered=TRUE)

vibesSA<-simper(otuSA, metaSA$Soil_Age, permutations=100)
summary(vibesSA)
summary(vibesSA,ordered=TRUE)$Control_Middle
summary(vibesSA,ordered=TRUE)

vibesSC<-simper(otuSC,metaSC$Soil_Age, permutations=100)
summary(vibesSC)
summary(vibesSC,ordered=TRUE)$Control_Old
summary(vibesSC,ordered=TRUE)
library(ggplot2)
detach(simpergraph)
attach(simpergraph)
simpergraph$Comparison = factor(simpergraph$Comparison, levels=c('Young_Old','Middle_Young','Middle_Old'),labels=c('Young_Old'="Young vs Old",'Middle_Young'="Young vs Intermediate",'Middle_Old'="Intermediate vs Old"))
simpergraph$Direction=factor(simpergraph$Direction , levels=c('Old', 'Intermediate', 'Young'))
ggplot1<-ggplot(simpergraph,aes(x = reorder(Phylum, +Contribution), y = Contribution, fill = Direction)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplot11<-ggplot1+coord_flip()
ggplot11
ggplot2<-ggplot11+scale_y_continuous(limits=c(0,1))+theme(axis.text.x = element_text(angle = 45,hjust=0.8,size=8))
ggplot2
plotSiMPER1<-ggplot2+scale_fill_manual(values=c("#6699CC","darkorange1","#999999"),breaks=c("Young","Intermediate","Old"))+guides(fill=guide_legend(title="Soil Type"))+ylab("Contribution Proportion")+xlab("Phylum")
plotSiMPER1
plotSimper22<-plotSiMPER1+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))+guides(color = guide_legend(override.aes = list(size = 0.4)))
plotSimper333<-plotSimper22+theme(axis.text=element_text(size=12),axis.text.x = element_text(size=13),
                                axis.title=element_text(size=14))
plotsim33<-plotSimper333+xlab("Bacterial Phylum")


detach(simpergraph)
detach(SimpergraphFung)
simpergraph$Comparison = factor(simpergraph$Comparison, levels=c('Young_Old','Middle_Young','Middle_Old'),labels=c('Young_Old'="Young vs Old",'Middle_Young'="Young vs Intermediate",'Middle_Old'="Intermediate vs Old"))
simpergraph$Direction=factor(simpergraph$Direction , levels=c('Old', 'Intermediate', 'Young'))
SimpergraphFung$Direction <- factor(SimpergraphFung$Direction , levels=c('Old', 'Intermediate', 'Young'))
ggplot122<-ggplot(SimpergraphFung,aes(x = reorder(Phylum, +Contribution), y = Contribution, fill = Direction)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplot222<-ggplot122+coord_flip()
ggplot222
ggplot223<-ggplot222+scale_y_continuous(limits=c(0,1))+theme(axis.text.x = element_text(angle = 45,hjust=0.8,size=8))
ggplot223
plotSiMPER2<-ggplot223+scale_fill_manual(values=c("#6699CC","darkorange2","#999999"),breaks=c("Young","Intermediate","Old"))+guides(fill=guide_legend(title="Soil Type"))+ylab("Contribution Proportion")+xlab("Phylum")
plotSiMPER2
plotSimper4<-plotSiMPER2+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))+guides(color = guide_legend(override.aes = list(size = 0.4)))
plotSimper5<-plotSimper4+theme(axis.text=element_text(size=12),axis.text.x = element_text(size=13),
                               axis.title=element_text(size=14))
plotsim6<-plotSimper5+xlab("Fungal Phylum")





##########simper per species and for fungi guilds ##########################
library(ggplot2)
cbp1 <- c(  "#009E73","#E69F00", "#56B4E9","#999999",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                   p + scale_color_manual(values = cbp1)

attach(Young_GUILD)
ggplotY<-ggplot(Young_GUILD, aes(x=reorder(Guild, +Contribution), y = Contribution, fill = Species)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplotY
ggplot21Y<-ggplotY+coord_flip()
ggplot212Y<-ggplot21Y+theme(axis.text.x = element_text(angle = 45,hjust=0.8,size=8))
plotSiMPERY<-ggplot212Y+scale_fill_manual(values=cbp1,breaks=c("SC","SA","GV","CP"))+guides(fill=guide_legend(title="Species"))+ylab("Contribution Proportion")+xlab("Guild")
plotSiMPERY
plotSimper2Y<-plotSiMPERY+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))
plotSimper2Y
plotSimper4Y<-plotSimper2Y+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),
                                axis.title=element_text(size=12))
plotSimper4Y+labs(title="Young Guild Contribution")

detach(Young_GUILD)
attach(Mid_GUILD)
ggplotM<-ggplot(Mid_GUILD, aes(x=reorder(Guild, +Contribution), y = Contribution, fill = Species)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplotM
ggplot21M<-ggplotM+coord_flip()
ggplot212M<-ggplot21M+theme(axis.text.x = element_text(angle = 45,hjust=0.8,size=8))
plotSiMPERM<-ggplot212M+scale_fill_manual(values=cbp1,breaks=c("SC","SA","GV","CP"))+guides(fill=guide_legend(title="Species"))+ylab("Contribution Proportion")+xlab("Guild")
plotSiMPERM
plotSimper2M<-plotSiMPERM+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))
plotSimper2M
plotSimper4M<-plotSimper2M+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),
                                 axis.title=element_text(size=12))
plotSimper4M+labs(title="Intermediate Guild Contribution")

detach(Mid_GUILD)
attach(Old_GUILD)
ggplotO<-ggplot(Old_GUILD, aes(x=reorder(Guild, +Contribution), y = Contribution, fill = Species)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplotO
ggplot21O<-ggplotO+coord_flip()
ggplot212O<-ggplot21O+theme(axis.text.x = element_text(angle = 45,hjust=0.8,size=8))
plotSiMPERO<-ggplot212O+scale_fill_manual(values=cbp1,breaks=c("SC","SA","GV","CP"))+guides(fill=guide_legend(title="Species"))+ylab("Contribution Proportion")+xlab("Guild")
plotSiMPERO
plotSimper2O<-plotSiMPERO+theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))
plotSimper2O
plotSimper4O<-plotSimper2O+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),
                                 axis.title=element_text(size=12))
plotSimper4O+labs(title="Old Guild Contribution")




detach(Guild1)
attach(Guild1)
summary(Guild1)
Guild1$Direction <- factor(Guild1$Direction , levels=c('Old', 'Intermediate', 'Young'))
ggplotG<-ggplot(Guild1,aes(x = factor(Guild,level=c("EMF","Undefined","AMF","Saprotroph","Pathogen")), y = Contribution, fill = Direction)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplotG
ggplot21G<-ggplotG+coord_flip()
ggplot21G
plotSiMPERG<-ggplot21G+scale_fill_manual(values=c("#6699CC","darkorange2","#999999"),breaks=c("Young","Intermediate","Old"))+guides(fill=guide_legend(title="Soil Type"))+ylab("Contribution Proportion")
plotSiMPERG
plotSimper2G<-plotSiMPERG+theme(legend.title = element_text(size = 14),legend.text = element_text(size = 14))+guides(color = guide_legend(override.aes = list(size = 0.8)))
plotSimper4G<-plotSimper2G+theme(axis.text=element_text(size=14),
                               axis.title=element_text(size=12))
plotSimper4G
plot55<-plotSimper4G+xlab("Functional Guild")
plot56<-plot55+theme(axis.text=element_text(size=14),
             axis.title=element_text(size=14))
plot56

detach(GuildInocula)
attach(GuildInocula)
GuildInocula$Guild<-level=c('value1', 'value2', 'value3')
GuildInocula$Direction <- factor(GuildInocula$Direction , levels=c('Old', 'Intermediate', 'Young'))
ggplotInoc<-ggplot(GuildInocula,aes(x = factor(Guild,level=c('EMF','Undefined', 'AMF', 'Saprotroph','Pathogen')), y = Contribution, fill = Direction)) +
  geom_col(position = "stack") +
  scale_fill_hue(direction = 1) +
  theme_minimal()
ggplotInoc
ggplot21Inoc<-ggplotInoc+coord_flip()
ggplot21Inoc
plotSiMPERInoc<-ggplot21Inoc+scale_fill_manual(values=c("#6699CC","darkorange2","#999999"),breaks=c("Young","Intermediate","Old"))+guides(fill=guide_legend(title="Soil Type"))+ylab("Contribution Proportion")
plotSiMPERInoc
plotSimper2Inoc<-plotSiMPERInoc+theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12))+guides(color = guide_legend(override.aes = list(size = 0.8)))
plotSimper4Inoc<-plotSimper2Inoc+theme(axis.text=element_text(size=14),
                                 axis.title=element_text(size=12))
plotInoc55<-plotSimper4Inoc+ylim(0,1)+xlab("Functional Guild")
plotinoc66<-plotInoc55+theme(legend.position="none")         
plotinoc67<-plotinoc66+theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=14))
plotinoc67
library(patchwork)
((plotinoc67|plot56)+plot_layout(nrow=1,byrow=TRUE,heights =c(10,10)))





####quick Basal Area Anova test to see if vary AM associating trees BA varies
attach(AMtrees)
Amonly<-subset(AMtrees,AMtrees$MF_Type=="AM")
Amonly
anovaBA<-aov(BA~Unit,data=Amonly)
summary(anovaBA)
TukeyHSD(anovaBA)




  

SubsettedFungiLEC<-subset_taxa(SoilFungi44, Phylum=="Glomeromycota" | Phylum=="Ascomycota" | Phylum=="Basidiomycota"| Phylum=="Chytridiomycota"|Phylum=="Mucoromycota"|Phylum=="unidentified")
SubsettedFungiLECGV<-prune_samples(sample_data(SubsettedFungiLEC)$Species=="GV",SubsettedFungiLEC)
SubsettedFungiLECSA<-prune_samples(sample_data(SubsettedFungiLEC)$Species=="SA",SubsettedFungiLEC)
SubsettedFingiLECSC<-prune_samples(sample_data(SubsettedFungiLEC)$Species=="SC",SubsettedFungiLEC)
library(dplyr)
b1<-phyloseq::psmelt(SubsettedFungiLECGV) %>%
  ggplot(data = ., aes(x = Soil_Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Soil Age", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
b1

attach(MID_SpSimpGuild)
library(ggplot2)
ggplot(MID_SpSimpGuild) +
  aes(x = Guild, y = Contribution, fill = Species) +
  geom_col() +
  scale_fill_manual(values=cbp1,breaks=c("SC","SA","GV","CP"))+guides(fill=guide_legend(title="Species"))+ylab("Contribution Proportion")+xlab("Guild")+
  theme_minimal()



######### abundance plots updated
GV.updated <- tax_glom(SoilFungiGV11, taxrank = "Phylum")
# Top N taxa
N <- 20
top <- names(sort(taxa_sums(GV.updated), decreasing = TRUE))[1:N]

# Calculate relative abundance
GP.genus.prop <- transform_sample_counts(GV.updated, function(x) x / sum(x) )

# Subset object to top N taxa
GP.genus.prop.top <- prune_taxa(top, GV.updated)


sample_data(SoilFungi46)$Soil_Age <- factor(sample_data(SoilFungi46)$Soil_Age, levels =c("Young","Middle","Old"))
levels(sample_data(SoilFungi46)$Soil_Age)
table(phyloseq::tax_table(SoilFungiGV11)[, "Phylum"])
ps_rel_abund = phyloseq::transform_sample_counts(GP.genus.prop.top , function(x){x / sum(x)})
phyloseq::otu_table(GP.genus.prop.top)[1:5, 1:5]
plotGVphy<-phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotGVphy1<-plotGVphy+theme(legend.position = "none")
plotGVphy1
table(phyloseq::tax_table(SoilFungiSA11)[, "Phylum"])
ps_rel_abundSA = phyloseq::transform_sample_counts(SoilFungiSA11, function(x){x / sum(x)})
phyloseq::otu_table(SoilFungiSA11)[1:5, 1:5]
plotSAphy<-phyloseq::plot_bar(ps_rel_abundSA, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotSAphy1<-plotSAphy+theme(legend.position = "none")
plotfunsa<-plotSAphy1+theme(axis.title.y = element_blank())
plotSAphy
table(phyloseq::tax_table(SoilFungiSC11)[, "Phylum"])
ps_rel_abundSC = phyloseq::transform_sample_counts(SoilFungiSC11, function(x){x / sum(x)})
phyloseq::otu_table(SoilFungiSC11)[1:5, 1:5]
plotSCphy<-phyloseq::plot_bar(ps_rel_abundSC, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotSCphy1<-plotSCphy+theme(legend.title = element_text(size=13),
                               legend.text = element_text(size=12))
plotfunsc<-plotSCphy1+theme(axis.title.y = element_blank())

##bacteria now 
sample_data(GVphylo1)$Soil_Age <- factor(sample_data(GVphylo1)$Soil_Age, levels =c("Young","Middle","Old"))
######### abundance plots updated
GV.updatedbac <- tax_glom(GVphylo1, taxrank = "Phylum")
# Top N taxa
N <- 20
top <- names(sort(taxa_sums(GV.updatedbac), decreasing = TRUE))[1:N]

# Calculate relative abundance
GP.genus.propbac <- transform_sample_counts(GP.genus.prop.topbac, function(x) x / sum(x) )

# Subset object to top N taxa
GP.genus.prop.topbac <- prune_taxa(top, GV.updatedbac)
plotGVphybac<-phyloseq::plot_bar(GP.genus.propbac, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotGVphybac1<-plotGVphybac+theme(legend.position = "none")
plotGVphybac1
######### abundance plots updated SA
sample_data(SAphylo1)$Soil_Age <- factor(sample_data(SAphylo1)$Soil_Age, levels =c("Young","Middle","Old"))
SA.updatedbac <- tax_glom(SAphylo1, taxrank = "Phylum")
# Top N taxa
N <- 20
top <- names(sort(taxa_sums(SA.updatedbac), decreasing = TRUE))[1:N]

# Calculate relative abundance
SA.genus.propbac <- transform_sample_counts(SA.genus.prop.topbac, function(x) x / sum(x) )
SA.genus.prop.topbac <- prune_taxa(top, SA.updatedbac)
# Subset object to top N taxa
SA.genus.prop.topbac <- prune_taxa(top, SA.updatedbac)
plotSAphybac<-phyloseq::plot_bar(SA.genus.propbac, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")+scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotSAphybac1<-plotSAphybac+theme(legend.position = "none")
plotSAAA<-plotSAphybac1+theme(axis.title.y = element_blank())

######### abundance plots updated
sample_data(SCphylo1)$Soil_Age <- factor(sample_data(SCphylo1)$Soil_Age, levels =c("Young","Middle","Old"))
SC.updatedbac <- tax_glom(SCphylo1, taxrank = "Phylum")
# Top N taxa
N <- 20
top <- names(sort(taxa_sums(SC.updatedbac), decreasing = TRUE))[1:N]
library(viridis)
# Calculate relative abundance
SC.genus.propbac <- transform_sample_counts(SC.genus.prop.topbac, function(x) x / sum(x) )

# Subset object to top N taxa
SC.genus.prop.topbac <- prune_taxa(top, SC.updatedbac)
plotSCphybac12<-phyloseq::plot_bar(SC.genus.propbac,fill="Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")+scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  labs(x = "Soil Age", y = "Relative Abundance\n") +
  facet_wrap(~ Soil_Age, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
plotSCphybac12
plotSCphybac11<-plotSCphybac12+theme(legend.title = element_text(size=12),
legend.text = element_text(size=11))
plotSCCC<-plotSCphybac11+theme(axis.title.y = element_blank())

library(patchwork)
((plotGVphy1|plotfunsa|plotfunsc)+plot_layout(nrow=1,byrow=TRUE,heights =c(10,10)))

sample_data(SoilFungi44)
ps_phylum12 <- phyloseq::tax_glom(SoilFungi44, "Phylum")
phyloseq::taxa_names(ps_phylum12) <- phyloseq::tax_table(ps_phylum12)[, "Phylum"]
phyloseq::otu_table(ps_phylum12)[1:3, 1:3]
b1912<-phyloseq::psmelt(ps_phylum12) %>%
  ggplot(data = ., aes(x = Soil, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(labels = c("Fusarium","AMF","PGPB"),values = c(Ascomycota = "#56B4E9",
                                                                    Glomeromycota = "#E69F00",
                                                                    Proteobacteria = "#999999"))+
  labs(x = "Soil", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
b111<-b1912+theme_bw()
b121<-b111+theme(axis.text=element_text(size=15),axis.title = element_text(size=15))
b131<-b121+theme(legend.text = element_text(size =15))+theme(axis.text.y = element_text(size=15))+theme(axis.text.x = element_text(size=15))
b141<-b131+theme( strip.background = element_blank() )
b141