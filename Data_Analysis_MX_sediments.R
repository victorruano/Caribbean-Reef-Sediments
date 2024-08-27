
################# Parrotfish size and abundance analysis ##################
library(ggplot2)
library(viridis)
library(hrbrthemes)

setwd("~/Documents/Research/Sediments")
pfsize<-read.csv("Parrotfish_size.csv")
pf.abund<-read.csv("Parrotfish_total_abundance.csv")

dplot.pf.site<-ggplot(data=pfsize, 
                      aes(x=Size_class, 
                          group=Site, 
                          fill=Site))+
  geom_density(adjust=1.5, 
               alpha=.4)+
  theme_classic()+ 
  labs(y = "Density", 
       x = "Size class")+
  theme(legend.position = "top",
        axis.ticks.length=unit(.25, "cm"),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
dplot.pf.site

dplot.pf.local<-ggplot(data=pfsize, 
                       aes(x=Size_class, 
                           group=Locality, 
                           fill=Locality))+
  geom_density(adjust=1.5, 
               alpha=0.5)+
  scale_fill_manual(values=c("firebrick3","goldenrod1", "steelblue3"))+
  theme_classic()+ 
  labs(y = "Density", 
       x = "Parrotfish size class")+
  theme(legend.position = "top",
        axis.ticks.length=unit(.25, "cm"),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
dplot.pf.local

ggsave("PFish_Size.tiff", plot = dplot.pf.local, 
       width = 6, height = 5, units = "in", dpi = 600)

# KS test
Akumal<-subset(pfsize, pfsize$Locality=="Akumal")
PAllen<-subset(pfsize, pfsize$Locality=="Punta Allen")
PMaroma<-subset(pfsize, pfsize$Locality=="Punta Maroma")

ks.test(Akumal$Size_class, PAllen$Size_class)
ks.test(PAllen$Size_class, PMaroma$Size_class)

# Nested ANOVA for abundance
shapiro.test(pf.abund$Abundance)

sq.pfabund<-sqrt(pf.abund$Abundance)
shapiro.test(sq.pfabund)

bartlett.test(sq.pfabund~pf.abund$Site)
bartlett.test(sq.pfabund~pf.abund$Locality)

# Nested ANOVA
summary(aov(sq.pfabund~pf.abund$Locality+Error(pf.abund$Site)))
summary(aov(sq.pfabund~pf.abund$Locality+pf.abund$Site))

################## Seidments multivariate analysis ######################
library(vegan)
library(BiodiversityR)
library(dplyr)
library(ggrepel)

mx.seds<-read.csv("Sediments_MX_allsites.csv")

mx.seds$Site<-factor(mx.seds$Site, 
                     levels = c("Dicks",
                                "Langosta",
                                "Yal Ku",
                                "Punta Allen N", 
                                "Punta Allen C", 
                                "San Antonio", 
                                "Mar F5"))

# Different magnitudes and relatively large ranges
sapply(mx.seds[,7:15], FUN = range)
sq.mx<-sqrt(mx.seds[,7:15])

# Same magnitudes and relatively short ranges
sapply(sq.mx, FUN = range)

# PERMANOVA
adonis2(sq.mx~mx.seds$Locality/mx.seds$Site, method = "bray", permutations = 9999)

# Simper

dissimilarity <- vegdist(sq.mx, method = "bray")
summary(simper(sq.mx,mx.seds$Site, permutations = 9999))
simper_results <- simper(dissimilarity,mx.seds$Site, permutations = 9999)

# MDS and stress
sedmx.mds=metaMDS(sq.mx, distance = "bray")
min(sedmx.mds$stress)

# Vectors for MDS
mxsed.fit <- envfit(sedmx.mds, mx.seds, permutations = 999)

site.scrs <- as.data.frame(scores(sedmx.mds, display = "sites")) 
#save NMDS results into dataframe

site.scrs <- cbind(site.scrs, Site = mx.seds$Site) 
#add grouping variable "Site" to dataframe

head(site.scrs)

spp.scrs <- as.data.frame(scores(mxsed.fit, display = "vectors")) 
#save species intrinsic values into dataframe

spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) 
#add species names to dataframe

spp.scrs <- cbind(spp.scrs, pval = mxsed.fit$vectors$pvals) 
#add pvalues to dataframe to select species which are significant

sig.spp.scrs <- subset(spp.scrs, pval<=0.05) 
#subset data to show species significant at 0.05

head(spp.scrs)
head(sig.spp.scrs)

nmds.plot.mxsed <- ggplot(site.scrs, 
                          aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(NMDS1, NMDS2, 
                 colour = factor(site.scrs$Site), 
                 shape = factor(site.scrs$Site)), 
             size = 3, 
             alpha = 0.6)+
  scale_shape_manual(values = c(16,16,16,17,17,17,15))+
  scale_color_manual(values=c("brown3","coral", "goldenrod", 
                              "olivedrab", "chartreuse", "springgreen",
                              "steelblue4"))+
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 1, 
                                        linetype = "solid"))+
  labs(colour = "Site", 
       shape = "Site")+
  theme(legend.title = element_blank(), 
        legend.position = "top", 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 10))

nmds.plot.mxsed

# Add significant species vectors to ordination plot
nmds.plot.mxsed +
  geom_segment(data = sig.spp.scrs, 
               aes(x = 0, xend=NMDS1, 
                   y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "grey10", 
               lwd=0.3) + 
  #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, 
                               label = Species), 
                           cex = 3, 
                           direction = "both", 
                           segment.size = 0.25) 
#add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

################## Linear mixed effects models ###################
sed.eros<-read.csv("Sediment_erosion_lm.csv")
sed.eros.binomal<-read.csv("Sediment_erosion_binomial_V2.csv")

# Bioerosion
sed.lme.be<-lme(Bioerosion~Coral_grain, 
                data=sed.eros, random = ~Locality|Site)
summary(sed.lme.be)

# Coral cover
sed.lme.cc<-lme(Coral_cover~Coral_grain, 
                data=sed.eros, random = ~Locality|Site)
summary(sed.lme.cc)

# Carbonate production
sed.lme.cc<-lme(C_production~Coral_grain, 
                data=sed.eros, random = ~Locality|Site)
summary(sed.lme.cc)

# Halimeda
hal.lme<-lme(Halimeda_cover~Halimeda_grain, data=sed.eros, random = list(~1|Site))
summary(hal.lme)

# Logistic regression
sed.eros.glme2<-glmer(Bioerosion~Coral_grain+(Locality|Site),
                      family=binomial,data=sed.eros.binomal)

###################### Nested ANOVAs ########################

# Bioerosion
summary(aov(log(sed.eros$Bioerosion) ~ sed.eros$Locality + 
              Error(sed.eros$Site)))
summary(aov(log(sed.eros$Bioerosion) ~ sed.eros$Locality + 
              sed.eros$Site))

TukeyHSD(aov(sed.eros$Bioerosion ~ sed.eros$Locality))

# Coral grains

shapiro.test(sed.eros$Coral_grain)
bartlett.test(sed.eros$Coral_grain~sed.eros$Locality)
bartlett.test(sed.eros$Coral_grain~sed.eros$Site)

summary(aov(sed.eros$Coral_grain ~ sed.eros$Locality + 
              Error(sed.eros$Site)))
summary(aov(sed.eros$Coral_grain ~ sed.eros$Locality + 
              sed.eros$Site))

TukeyHSD(aov(sed.eros$Coral_grain ~ sed.eros$Locality))

# Carbonate production
shapiro.test(log(sed.eros$C_production))
bartlett.test(log(sed.eros$C_production)~sed.eros$Locality)
bartlett.test(log(sed.eros$C_production)~sed.eros$Site)

summary(aov(log(sed.eros$C_production) ~ sed.eros$Locality + 
              Error(sed.eros$Site)))
summary(aov(log(sed.eros$C_production) ~ sed.eros$Locality + 
              sed.eros$Site))

# Net c production
shapiro.test(log(sed.eros$Net_c))
bartlett.test(log(sed.eros$Net_c)~sed.eros$Locality)
bartlett.test(log(sed.eros$Net_c)~sed.eros$Site)

summary(aov(log(sed.eros$Net_c) ~ sed.eros$Locality + 
              Error(sed.eros$Site)))
summary(aov(log(sed.eros$Net_c) ~ sed.eros$Locality + 
              sed.eros$Site))

# Coral cover
shapiro.test(sed.eros$Coral_cover)
bartlett.test(sed.eros$Coral_cover~sed.eros$Locality)
bartlett.test(sed.eros$Coral_cover~sed.eros$Site)

summary(aov(sed.eros$Coral_cover ~ sed.eros$Locality + 
              Error(sed.eros$Site)))
summary(aov(sed.eros$Coral_cover ~ sed.eros$Locality + 
              sed.eros$Site))



