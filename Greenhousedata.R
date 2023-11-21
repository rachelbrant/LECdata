###attempting to analyze greenhouse plant growth data
##file is called Greenhouse_dataall, it includes all plants in three soil types
attach(Greenhouse_dataall)
growthLast<-subset(Greenhouse_dataall,week=="5")
shapiro.test(growthLast$length)
growthLast$loglength<-log(growthLast$length)
View(growthLast)
##growthLast is the main dataset used because it is solely week 5, no remnant
growthLast<-subset(Greenhouse_dataall,week=="5")
growthLast<-subset(growthLast,inoculant!="Remnant")
## subsetting growthLast into Species if they end up being needed 
growthLastCP<-subset(growthLast,species=="corpub")
growthLastGV<-subset(growthLast,species=="geuvir")
growthLastSA<-subset(growthLast,species=="solarg")
growthLastSC<-subset(growthLast,species=="solcae")
growthLastCP$logCP=log(growthLastCP$leafarea)
growthLastGV$logGV=log(growthLastGV$leafarea)
growthLastSA$logSA=log(growthLastSA$leafarea)
growthLastSC$logSC=log(growthLastSC$leafarea)
summary(growthLastGV)

##ANOVAs of each species in each inoculant,
#UPDATAD we ended up doing  lsmeans instead, see below  
trythisAOVCP<-aov(leafarea~inoculant,data=growthLastCP)
summary(trythisAOVCP)
TukeyHSD(trythisAOVCP)

trythisAOVGV<-aov(leafarea~inoculant,data=growthLastGV)
summary(trythisAOVGV)
TukeyHSD(trythisAOVGV)

trythisAOVSA<-aov(growthLastSA$logSA~inoculant,data=growthLastSA)
summary(trythisAOVSA)
TukeyHSD(trythisAOVSA)

trythisAOVSC<-aov(leafarea~inoculant,data=growthLastSC)
summary(trythisAOVSC)
TukeyHSD(trythisAOVSC)


########Lsmeans instead of individual ANOVA 3-22-23 ##########
library(lsmeans)
# adjusted means and comparisons, treating machine C as control
attach(Greenhouse_dataall)
growthLast<-subset(Greenhouse_dataall,week=="5")
growthLast
leaf.lm3CP5 <- lmer(loglength ~ inoculant+ (1|block), data = growthLastCP)
anova(leaf.lm3CP5)
print(lsmeans (leaf.lm3CP5,  list(pairwise ~inoculant )),omit=10)

leaf.lm3GV5 <- lmer(loglength ~ inoculant+ (1|block), data = growthLastGV)
anova(leaf.lm3GV5)
print(lsmeans (leaf.lm3GV5,  list(pairwise ~inoculant )),omit=10)

leaf.lm3SA5 <- lmer(loglength ~ inoculant+ (1|block), data = growthLastSA)
anova(leaf.lm3SA5)
print(lsmeans (leaf.lm3SA5,  list(pairwise ~inoculant )),omit=10)

leaf.lm3SC5 <- lmer(loglength ~ inoculant+ (1|block), data = growthLastSC)
anova(leaf.lm3SC5)
print(lsmeans (leaf.lm3SC5,  list(pairwise ~inoculant )),omit=10)

##  How length could be a proxy for overall growth is another importnat avenue of research
##  Add a paragraph about how putting things in at the early stages and transplanting, but need 
##  To check on other morphometric data because we couldnt destructively harvest the plants
##  We don't know how they are doing I the field but there is an ongoing study to assess this 
##  Length is highly correlated to overall area, conducted proxies for growth, all highly correlated so we just present this one
### One of the weaknesses of this study is that we couldnt perform any destructive analyes, so it would be interesting to see how they performed in the field, which is an ongoing study. 

qqnorm(residuals(leaf.lm3))
##checking residuals with DHARMa package 
#install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = leaf.lm3, plot = T)


###Summary Statistic SE for graphs 
library(plyr)           
library(ggplot2)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
attach(Greenhouse_dataall)

Greenhoues<-na.rm(Greenhouse_dataall)
library(dplyr)
df.summary2 <- week12Area %>%
  group_by(week12Area$Unit,week12Area$Species) %>%
  summarise(
    sd = sd(week12Area$leafLength),
    len = mean(week12Area$leafLength)
  )
df.summary2


tgc <- summarySE(week12Area, measurevar="leafLength", groupvars=c("Unit","Species"))
tgc
write.csv(tgc,"/Users/rachelbrant/Desktop/Selength12.csv")
View(growthLast)
df <- na.omit(Greenhouse_dataall)
gorwthlast2<-na.omit(growthLast)
View(df)
data_col_NA <- growthLast[!is.na(growthLast$length), ]    # Drop NAs in only one column
data_col_NA                               # Print data frame without NA in one column
View(data_col_NA)
df.summary3 <- data_col_NA %>%
  group_by(data_col_NA$inoculant,data_col_NA$species) %>%
  summarise(
    sd = sd(data_col_NA$length),
    len = mean(data_col_NA$length)
  )
df.summary3
tgc2 <- summarySE(data_col_NA, measurevar="length", groupvars=c("inoculant","species"))
tgc2
write.csv(tgc2,"/Users/rachelbrant/Desktop/Selength5.csv")

##added the tgc (SE values) with the means to create boxes for the graph, data is 
##called newgreenhousePlant
######## full script for graph used in MS to display growth data##### 
attach(newgreenhousePlant)
attach(newse)
newse$Soil2 <- factor($Soil2, levels = c("Young", "Intermediate", "Old"))
p<-ggplot(newse,aes(x = Soil2, y=leafArea,fill=Soil2)) + 
  geom_col()+scale_fill_manual(values=c("#6699CC","#009E73","#999999"))+
  geom_errorbar(aes(x=Soil2, ymin=leafArea-se, ymax=leafArea+se), width=0.2, colour="black", alpha=0.4)
p
p1<-p+facet_grid(~Species)
pL<-p1+theme_bw()
pL
p2<-pL+xlab("Soil Age")+ylab("Leaf Area")
p2
p3<-p2
p3
p4<-p3+scale_x_discrete(labels=c("Young" = "Young", "Middle" = "Intermediate",
                                     "Old" = "Old"))
p5<-p4+ theme(axis.text.x = element_text(
                                     hjust=1,size = 13, angle = 45))

p6<-p5+theme(strip.text.x = element_text(size = 15))+theme(axis.title.y = element_text(size = 13))
p6 + guides(fill=guide_legend(title="Soil Age"))+theme(axis.title.x=element_text(size=13))
dev.off()

View(Numberlef)
Numberlef<-subset(Numberleaf_data,Numberleaf_data$Species!="CP")

nLeafthree$Soil2 <- factor(nLeafthree$Soil2, levels = c("Young", "Intermediate", "Old"))
nLeafSA$Soil2 <- factor(nLeafSA$Soil2, levels = c("Young", "Intermediate", "Old"))
nLeafSC$Soil2 <- factor(nLeafSC$Soil2, levels = c("Young", "Intermediate", "Old"))

pNum<-ggplot(nLeafthree,aes(x = Soil2, y=nLeaf,fill=Soil2)) + 
  geom_boxplot()+scale_fill_manual(values=c("darkorange2","#999999","#6699CC"))+scale_color_manual(values=c("darkorange2","#999999","#6699CC"))
pNum
p1Num<-pNum+facet_grid(~Species)
pLNum<-p1Num+theme_bw()
pLNum
p2Num<-pLNum+xlab("Soil Age")+ylab("Number of Leaves")
p2Num
p3Num<-p2Num
p3Num
p4Num<-p3Num+scale_x_discrete(labels=c("Young" = "Young", "Middle" = "Intermediate",
                                 "Old" = "Old"))
p5Num<-p4Num+ theme(axis.text.x = element_text(
  hjust=1,size = 13, angle = 45))

p6Num<-p5Num+theme(strip.text.x = element_text(size = 15))+theme(axis.title.y = element_text(size = 13))
p6Num + guides(fill=guide_legend(title="Soil Age"))+theme(axis.title.x=element_text(size=13))











attach(Selength5)
Selength5$inoculant <- factor(Selength5$inoculant, levels = c("Young", "Intermediate", "Old"))
p<-ggplot(Selength5,aes(x = inoculant, y=length,fill=inoculant)) + 
  geom_col()+scale_fill_manual(values=c("#6699CC","darkorange2","#999999"))+
  geom_errorbar(aes(x=inoculant, ymin=length-se, ymax=length+se), width=0.2, colour="black", alpha=0.4)
p
p1<-p+facet_grid(~species)
pL<-p1+theme_bw()
pL
p2<-pL+xlab("Soil Age")+ylab("Leaf length (cm)")
p2
p3<-p2
p3
p4<-p3+scale_x_discrete(labels=c("Young" = "Young", "Middle" = "Intermediate",
                                 "Old" = "Old"))
p5<-p4+ theme(axis.text.x = element_text(
  hjust=1,size = 13, angle = 45))

p6<-p5+theme(strip.text.x = element_text(size = 15))+theme(axis.title.y = element_text(size = 13))
p7<-p6 + guides(fill=guide_legend(title="Soil Age"))+theme(axis.title.x=element_text(size=13))
p8<-p7+theme(axis.text.x = element_blank()) + theme(axis.title.x = element_blank())
p9<-p8+ylim(0,15)
p10<-p9+theme(legend.title=element_text(size=14),legend.text = element_text(size =12))
p10


attach(Selength12)
Selength12$inoculant <- factor(Selength12$inoculant, levels = c("Young", "Intermediate", "Old"))
p8<-ggplot(Selength12,aes(x = inoculant, y=length,fill=inoculant)) + 
  geom_col()+scale_fill_manual(values=c("#6699CC","darkorange2","#999999"))+
  geom_errorbar(aes(x=inoculant, ymin=length-se, ymax=length+se), width=0.2, colour="black", alpha=0.4)
p8
p9<-p8+facet_grid(~species)
p9L<-p9+theme_bw()
p9L
p92<-p9L+xlab("Soil Age")+ylab("Leaf length (cm)")
p92
p93<-p92
p93
p94<-p93+scale_x_discrete(labels=c("Young" = "Young", "Middle" = "Intermediate",
                                 "Old" = "Old"))
p95<-p94+ theme(axis.text.x = element_text(
  hjust=1,size = 13, angle = 45))

p96<-p95+theme(strip.text.x = element_text(size = 15))+theme(axis.title.y = element_text(size = 13))
p97<-p96 + guides(fill=guide_legend(title="Soil Age"))+theme(axis.title.x=element_text(size=13))
p98<-p97+ylim(0,25)
p99<-p98+theme(legend.title=element_text(size=14),legend.text = element_text(size =12))
library(patchwork)
((p10)/(p99)+plot_layout(nrow=2,byrow=TRUE,heights =c(15,15)))

### adjusted means and comparisons  ##### ##copy of above script line 42 ### 
library(lsmeans)
attach(Greenhouse_dataall)
growthLast<-subset(Greenhouse_dataall,week=="5")
growthLast
leaf.lm <- lm(leafarea ~ inoculant*species, data = growthLast)
print(lsmeans (leaf.lm,  list(pairwise ~ inoculant |species)),omit=10)
plot(leaf.lm)
##check for package to plot residuals of LMER ## 

#####AM tree density comparison ########
attach(AMtrees)
summary(AMtrees)
AMonly<-subset(AMtrees, MF_Type=="AM",data=AMtrees)
amtreetest<-aov(BA~Unit,data=AMonly)
EMonly<-subset(AMtrees, MF_Type=="EM",data=AMtrees)
Emtreetest<-aov(BA~Unit,data=EMonly)
summary(Emtreetest)
TukeyHSD(Emtreetest)


##A quick tukey HSD test on the Soil Vars to see which direction P and K are sig####
attach(SoilVariables)
library(dplyr)
install.packages("corrr")
library(corrr)
newsoilvas<-SoilVariables%>%
  dplyr::select(pH,OM,Bray_IP,CA_lb,Mg_lb,K_lb,CEC_lb,Silt)
summary(newsoilvas)
cor.vas<-correlate(newsoilvas,method="spearman",diagonal=1)
cor.vas
install.packages("ggcorrplot")
library(ggcorrplot)
corr <- cor(newsoilvas)
head(corr[, 1:6])

ggcorrplot(corr, hc.order = TRUE)

###run pca on soil variables
###PCA ONE WAY TO DO IT 
pca_soilvars<-prcomp(newsoilvas,scale=TRUE,center=TRUE)
summary(pca_soilvars)
newsoilvas1<-as.matrix(newsoilvas)
 install.packages("devtools")
devtools::install_github("arleyc/PCAtest")
PCAtest(newsoilvas, varcorr=FALSE, counter=FALSE, plot=TRUE)
library(vegan)
library(PCAtest)
library(ggfortify)
levels(SoilVariables$Site)
plotsoilvas1<-autoplot(pca_soilvars, data =SoilVariables,colour="Site",size = 3, loadings = TRUE,loadings.label.repel=T,loadings.label = TRUE, loadings.label.size = 6,loadings.colour="black",loadings.label.colour = "black",frame.type="norm")
pl2<-plotsoilvas1+theme_bw()+scale_fill_manual(values=c("darkorange2","#999999","#6699CC"))+scale_color_manual(values=c("darkorange2","#999999","#6699CC"))
plotsoilvas1
pl3<-pl2+theme(axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16))+theme(axis.text=element_text(size=16),
                                                   axis.title=element_text(size=16))+theme(legend.title = element_text(size=16), #change legend title font size
                                                                                           legend.text = element_text(size=14)) #change legend text font size)

pl3
library(patchwork)


anovaSoil<-aov(K_lb~Site,data=SoilVariables)
summary(anovaSoil)
TukeyHSD(anovaSoil)
(K_lb~Site)
##A quick anova on soil C and N for new table 3-26#####
attach(Cnsoil)
anovaCN<-aov(N....~Site,data=Cnsoil)
summary(anovaCN)



###################### Number of leaves##################################

###number of leaves data, we dont end up using this because it is too highly variable
#attach(Numberleaf_data)
#library(moments)
#sapply(Numberleaf_data, mode)
#shapiro.test(Numberleaf_data$nLeaf)
#shapiro.test(completerecords$variable_log)
#skewness(Numberleaf_data$nLeaf, na.rm = TRUE)
#library(dplyr)
#Numberleaf_data %>%
#mutate(log_var = log(nLeaf))
#head(Numberleaf_data)
#Numberleaf_data$variable_log <- log(Numberleaf_data$nLeaf)
#completerecords <- na.omit(Numberleaf_data) 
#head(completerecords)
#summary(completerecords)

##More plots, but didnt use these 
#library(ggpubr)
#ggboxplot(growthLast, x = "inoculant", y = "leafarea",
#   color = "inoculant", palette = "jco")+
#  stat_compare_means()

##Plots for number leaves, but dont use 
#pnumber <- ggboxplot(completerecords, x = "inoculant", y = "nLeaf",
#            color="inoculant",palette = "jco",
#           add = "jitter", short.panel.labs = FALSE)
#pnumber1<-pnumber+stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(430,420,410))+ylim(0,800)
#pnumber1+ theme(legend.position="none")
#my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))



#library(MASS)
##attempt at binomial dis for leaf number, need to subset all the data first
attach(Numberleaf_data)
summary(Numberleaf_data)
Numberleaf_data$Soil2<-as.factor(Numberleaf_data$Soil2)
library(MASS)
nLeafthree<-subset(Numberleaf_data,Numberleaf_data$Species!="Cp")
nLeafGV<-subset(Numberleaf_data,Numberleaf_data$Species=="Gv")
nLeafSA<-subset(Numberleaf_data,Numberleaf_data$Species=="Sa")
nLeafSC<-subset(Numberleaf_data,Numberleaf_data$Species=="Sc")
shapiro.test($logLeaf)

#Numberleaf_data$logLeaf=log(Numberleaf_data$nLeaf)
#redologdata<-subset(Numberleaf_data,Numberleaf_data$logLeaf!="NA")
#write.csv(redologdata,"/Users/rachelbrant/Desktop/logleaves.csv")
nleavesCP <- glm.nb(Height ~Soil2, data = nLeafCP)
summary(nleavesCP)
nleavesGV <- glm.nb(nLeaf ~Soil2, data = nLeafGV)
summary(nleavesGV)
nleavesSA <- glm.nb(nLeaf ~Soil2, data = nLeafSA)
summary(nleavesSA)
nleavesSC <- aov(nLeaf ~Soil2, data = nLeafSC)
anova(nleavesSC)
summary(nleavesSC)
TukeyHSD(nleavesSC)
(nleavesSA)
library(dplyr)
library(corrr)
newleaves<-nLeafCP%>%
  dplyr::select(nLeaf,leafWidth,leafLength)
summary(newleaves)
cor.leaves<-correlate(newleaves,method="spearman",diagonal=1)
cor.leaves

newleaves<-nLeafGV%>%
  dplyr::select(nLeaf,leafWidth,leafLength)
cor.leaves<-correlate(newleaves,method="spearman",diagonal=1)
cor.leaves

newleaves<-nLeafSA%>%
  dplyr::select(nLeaf,leafWidth,leafLength)
cor.leaves<-correlate(newleaves,method="spearman",diagonal=1)
cor.leaves

newleaves<-nLeafSC%>%
  dplyr::select(nLeaf,leafWidth,leafLength)
cor.leaves<-correlate(newleaves,method="spearman",diagonal=1)
cor.leaves

###glmmtb model fit w random effect of block 
##use block
summary(nleavesSC)
attach(logleaves)
summary(logleaves)
leavesGlm.nb<-glm.nb(logLeaf ~ Unit*Species, data = logleaves)
summary(leavesGlm.nb)
#mSA<-glm.nb(nLeaf ~ inoculant, data = nLeafSA)
#summary(mSA)
#mSC<-glm.nb(nLeaf ~ inoculant, data = nLeafSC)
#summary(mSC)
#poisCP <- glm(nLeaf~ inoculant, "poisson", data = nLeafCP)
#summary(poisCP)
#poisGV <- glm(nLeaf~ inoculant, "poisson", data = nLeafGV)
#summary(poisGV)
#poisSA<-glm(nLeaf~ inoculant, "poisson", data = nLeafSA)
#summary(poisSA)
#poisSC<-glm(nLeaf~ inoculant, "poisson", data = nLeafSC)
#summary(poisSC)

#ggboxplot(nLeafCP, x = "inoculant", y = "variable_logn",
#          color = "inoculant", palette = "jco")+
#  stat_compare_means()

#pCP <- ggboxplot(nLeafCP, x = "inoculant", y = "nLeaf",
#        color="inoculant",palette = "jco",
#          add = "jitter")
#pCP1<-pCP+stat_compare_means(comparisons = my_comparisonsCP, label = "p.signif",label.y = c(500,480))+ylim(0,600)
#pCP1+ theme(legend.position="none")
#my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))
#my_comparisonsCP <- list(c("Young", "Old"))

#GVp <- ggboxplot(nLeafGV, x = "inoculant", y = "nLeaf",
#        color="inoculant",palette = "jco",
#        add = "jitter")
#GVp
#GVp1<-GVp+stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y=c(360,300,340))+ylim(0,400)
#GVp1
#GVp1+ theme(legend.position="none")
#my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))
#my_comparisonsCP <- list(c("Young", "Old"))

#SAp <- ggboxplot(nLeafSA, x = "inoculant", y = "nLeaf",
#         color="inoculant",palette = "jco",
#         add = "jitter")
#SAp1<-SAp+stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(360,300,340))+ylim(0,400)
#SAp1+ theme(legend.position="none")
#my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))
#my_comparisonsCP <- list(c("Young", "Old"))

#SCp <- ggboxplot(nLeafSC, x = "inoculant", y = "nLeaf",
#  color="inoculant",palette = "jco",
#   add = "jitter")
#SCp1<-SCp+stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(520,500,480))+ylim(0,600)
#SCp1+ theme(legend.position="none")
#my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))
#my_comparisonsCP <- list(c("Young", "Old"))


###fix box plot to facet wrap the species and have x axis be for all 
library(esquisse)
library(moments)
na.omit(growthLast$leafarea)
datafixed<-log10(growthLast$leafarea)
my_comparisons <- list(c("Young", "Old"), c("Old", "Middle"), c("Young","Middle"))
library(lme4)
library(nlme)
install.packages("ggpubr")
library(ggpubr)
install.packages("lmerTest")
library(lmerTest)


##MA suggested try a mixed model to control for block, ended up not using  
#model1<-lm(logweight~inoculant+species,data=growthLast)
#model2<-lm(logweight~inoculant*species,data=growthLast)
#anova(lmer(logweight~inoculant+species+(1|block),data=growthLast),lmer(logweight~inoculant*species+(1|block),data=growthLast))



attach(week12Area)
week12Area<-na.omit(week12Area)
View(week12Area)
week12Area$leafLength<-log(week12Area$leafLength)
View(week12Area)
shapiro.test(week12Area$leafLength)    

########Lsmeans week 12 ############################
library(lsmeans)
# adjusted means and comparisons
library(lmerTest)
lengthCP<-subset(week12Area,week12Area$Species=="Cp")
leaf.lm3CP <- lmer(leafLength ~ Unit + (1|Block), data = lengthCP)
anova(leaf.lm3CP)
print(lsmeans (leaf.lm3CP,  list(pairwise ~ Unit)),omit=10)

lengthGV<-subset(week12Area,week12Area$Species=="Gv")
leaf.lm3GV <- lmer(leafLength ~ Unit + (1|Block), data = lengthGV)
anova(leaf.lm3GV)
print(lsmeans (leaf.lm3GV,  list(pairwise ~ Unit)),omit=10)

lengthSA<-subset(week12Area,week12Area$Species=="Sa")
leaf.lm3SA <- lmer(leafLength ~ Unit + (1|Block), data = lengthSA)
anova(leaf.lm3SA)
print(lsmeans (leaf.lm3SA,  list(pairwise ~ Unit)),omit=10)

lengthSC<-subset(week12Area,week12Area$Species=="Sc")
leaf.lm3SC <- lmer(leafLength ~ Unit + (1|Block), data = lengthSC)
anova(leaf.lm3SC)
print(lsmeans (leaf.lm3SC,  list(pairwise ~ Unit)),omit=10)

anova(leaf.lm3)
print(lsmeans (leaf.lm3,  list(pairwise ~ Unit |Species)),omit=10)
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = leaf.lm3SC, plot = T)




###### new box plots for the tree densities #####
attach(Tree_boxplots)
Tree_boxplots$Location <- factor(Tree_boxplots$Location, levels = c("Young", "Intermediate", "Old"))
Canopy<-subset(Tree_boxplots,Tree_boxplots$Category=="Canopy")
Trees<-subset(Tree_boxplots,Tree_boxplots$Association!="Unknown")
plotTree<-ggplot(Trees,aes(x = Association, y = Density,fill=Location)) +
  geom_boxplot(position = position_dodge(width = 1)) +scale_fill_manual(values=c("#6699CC","#009E73","#999999"))+
theme_minimal()
PlotTree1<-plotTree+ylab("Average Density")+ theme(axis.text.x = element_text(
  hjust=1,size = 10, angle = 45))
plotTree2<-PlotTree1+facet_grid(~Category)
plotTree45<-PlotTree1+ggtitle('A')
class(Canopy$se)
Canopy$se<-as.numeric(Canopy$se)


library(dplyr)
AMtrees$Location <- factor(AMtrees$Location , levels = c("Young", "Intermediate", "Old"))
Plot54<-ggplot(AMtrees, aes(x = MF_Type, y = BA, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 1))+scale_fill_manual(values=c("#6699CC","#009E73","#999999"))
               
Plot55<-Plot54+theme_minimal()
Plot56<-Plot55+xlab("Association")+ylab("Basal Area (cm)")
Plot57<-Plot56+ggtitle("B")
install.packages("patchwork")
library(patchwork)
((plotTree45)/(Plot57)+plot_layout(nrow=2,byrow=TRUE,heights =c(10,10,6,6)))

Amsub<-subset(Tree_boxplots,Tree_boxplots$Association=="AM")
EMsum<-subset(Tree_boxplots,Tree_boxplots$Association=="EM")
AMBA<-subset(AMtrees,AMtrees$MF_Type=="AM")
EMBA<-subset(AMtrees,AMtrees$MF_Type=="EM")
amanova<-aov(BA~Location,data=AMBA)
summary(amanova)
TukeyHSD(amanova)




########################## canopy and midstory NMDS for June LEC updates ##########
detach(Tree_NMDSinfo)
detach(AMtrees)
dttach(DataTrees)

attach(envitrees2)
as.matrix(envitrees2)
as.matrix(DataTrees)
#data(dune)
#data("dune.env")
View(dune.env)
View(dune)
library(vegan) # for vegetation and community analysis
library(dplyr)
type(tree.mds)
 set.seed(1019)
tree.mds <- metaMDS(DataTrees, distance = "bray", autotransform = TRUE)

plot(tree.mds)
data(dune)
data("dune.env")
View(dune.env)


dune.mds <- metaMDS(DataTrees, distance = "bray", autotransform = FALSE)
plot(dune.mds,type="n")
ordiplot(dune.mds,type="n")
orditorp(dune.mds,display="species",col="red",air=0.01)
orditorp(dune.mds,display="sites",cex=1.25,labels = T, pch = c(16, 8, 17, 6) [(envitrees2$Plot)], col = c( "darkorange2","#6699CC","#999999") [(envitrees2$Unit)])



library(ggplot2)
install.packages("BiodiversityR")
library(BiodiversityR)
library(readxl)
library(ggsci)
library(ggrepel)
install.packages("ggforce")
library(ggforce)
library(vegan)
as.list(envitrees2)
class(dune.env$Use)
envitrees2$Age<-as.factor(envitrees2$Age)
summary(envitrees2)
with(envitrees2, levels(Age))
ordiplot(tree.mds,type="none",scaling=2,ylim = c(-2,2),air=0.01,xlim = c(-1,1),cex.lab=1.4,cex.axis=1.3)

orditorp(tree.mds, display = "species", labels = T, pch = c(16, 9, 17) [(envitrees2$Age)], col = c( "darkorange2","#999999","#6699CC") [(envitrees2$Age)], cex = 2)
plot(dune.spp.fit, p.max = 0.06, col = "black", cex = 1.2)
ordiellipse(tree.mds, groups = envitrees2$Age, lty = 1, col = c("darkorange2","#999999","#6699CC"))
legend(x="bottomleft", legend=levels(envitrees2$Age), pch = c(16, 9, 17), col = c( "darkorange2","#999999","#6699CC"))
dune.spp.fit <- envfit(tree.mds, DataTrees, permutations = 999)
head(dune.spp.fit)
plot(dune.spp.fit, p.max = 0.05, col = "black", cex = 0.7)

library(dplyr)
View(DataTrees)
soil.nmds11<-metaMDS(DataTrees,distance="bray",trymax = 600,halfchange = TRUE)   
soil.adonisr1322S<-adonis2(DataTrees ~ factor(Unit)+(1|Plot) ,data=envitrees2,permutations=999,method="bray")
(soil.adonisr1322S)






##density 
attach(CanopyDensity)
densitytrees<-lm(AvgDensity~Genus*Location,data=CanopyDensity)
summary(densitytrees)




#### growthrate instead using leaf length over time 
attach(cleaned_growthdata)
summary(cleaned_growthdata)
class(cleaned_growthdata$length)
cleanedgeu<-subset(cleaned_growthdata,cleaned_growthdata$species=="geuvir")
cleanedsa<-subset(cleaned_growthdata,cleaned_growthdata$species=="solarg")
cleanedsc<-subset(cleaned_growthdata,cleaned_growthdata$species=="solcae")
library(gplots)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


tgcgeu <- summarySE(cleanedgeu, measurevar="length", groupvars=c("week","inoculant"))
tgcgeu
tgcsa <- summarySE(cleanedsa, measurevar="length", groupvars=c("week","inoculant"))
tgcsa
tgcsc <- summarySE(cleanedsc, measurevar="length", groupvars=c("week","inoculant"))
tgcsc
tgcgeu$inoculant <- factor(tgcgeu$inoculant , levels = c("Young", "Intermediate", "Old"))
levels(tgcgeu$inoculant)
cleanedgeu$week<-as.factor(cleanedgeu$week)
plotgeu<-ggplot(tgcgeu, aes(x=week, y=length, colour=inoculant)) + 
  geom_errorbar(aes(ymin=length-se, ymax=length+se), width=.1) +
  geom_line()+scale_color_manual(values=c("#6699CC","#009E73","#999999"))+
  geom_point()+theme_bw()
pgeu<-plotgeu+ scale_x_continuous(breaks=c(1,2, 3, 4,5,12))
pgeu1<-pgeu+theme(legend.position = "none")+ylab("Leaf Length (cm)")+ylim(5,25)
pgeu1
tgcsa$inoculant <- factor(tgcsa$inoculant , levels = c("Young", "Intermediate", "Old"))
levels(tgcsa$inoculant)
cleanedsa$week<-as.factor(cleanedsa$week)
plotsa<-ggplot(tgcsa, aes(x=week, y=length, colour=inoculant)) + 
  geom_errorbar(aes(ymin=length-se, ymax=length+se), width=.1) +
  geom_line()+scale_color_manual(values=c("#6699CC","#009E73","#999999"))+
  geom_point()+theme_bw()
psa<-plotsa+ scale_x_continuous(breaks=c(1,2, 3, 4,5,12))
psa1<-psa+theme(legend.position = "none")+theme(axis.title.y=element_blank())+ylim(5,25)
psa1
tgcsc$inoculant <- factor(tgcsc$inoculant , levels = c("Young", "Intermediate", "Old"))
levels(tgcsc$inoculant)
cleanedsc$week<-as.factor(cleanedsc$week)
plotsc<-ggplot(tgcsc, aes(x=week, y=length, colour=inoculant)) + 
  geom_errorbar(aes(ymin=length-se, ymax=length+se), width=.1) +
  geom_line()+scale_color_manual(values=c("#6699CC","#009E73","#999999"))+
  geom_point()+theme_bw()
psc<-plotsc+ scale_x_continuous(breaks=c(1,2, 3, 4,5,12))
psc1<-psc+guides(color=guide_legend(title="Restoration Age"))+theme(axis.title.y=element_blank())+ylim(5,25)

library(patchwork)
((pgeu1|psa1|psc1)+plot_layout(nrow=1,byrow=TRUE,heights =c(10,10)))





### repeated measures anova 
########################Repeated Measures Anova############################
if(!require(psych)){install.packages("psych")}
if(!require(nlme)){install.packages("nlme")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(rcompanion)){install.packages("rcompanion")}
library(nlme)

###########More appropriate way############
qqnorm(cleaned_growthdata$loglength, pch = 1, frame = FALSE)
qqline(cleaned_growthdata$loglength, col = "steelblue", lwd = 2)
library("car")
qqPlot(cleaned_growthdata$loglength)
repeated.aov<-with(cleanedsa,aov(length ~ factor(inoculant)*factor(week))) 
summary(repeated.aov)
cleaned_growthdata$length <-(cleaned_growthdata$length)
shapiro.test(loglength)
cleaned_growthdata$length 
loglength = log(cleaned_growthdata$length)
summary(loglength) 

