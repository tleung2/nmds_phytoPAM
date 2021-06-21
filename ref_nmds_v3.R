
  ## LOAD PACKAGES
library(vegan) 
library(viridis)
library(NbClust) #for finding optimal clusters
library(cluster) #for computing cluster silhouettes
library(tidyverse)
library(factoextra) #computing heirarchical cluster
library(vegan3d) # For plotting 3D images, autoload vegan package
library(plotly)
library(plot3D)

  ## Use 2018-2020 fluorescence data with references all in one dataset. 
  ## Then plot the ordination to see how the fluorescence samples  
  ## compare with the references (default and lab measured)

  ## Read data----from Fluorescence Ref Folder--------------------
data.all=read.csv("all_data_2018-2020_v3.csv",header=T)
data.all$fit_error=as.numeric(data.all$fit_error)


#######################################################################
  ########################  TIDY DATASET   ##########################

  ## dataset without non-default references
data.all2 <- data.all %>%
  filter(source !='non_Default')

########################################################################
  ########################  NORMALIZE DATA  ########################
  
  ## Normalize samples realtive to max (max F values scaled to 1)
  ## First subset samples then scale to max F 
data.norm <- as.data.frame(t(apply(data.all2[2:6], 1, 
                                   (function(x) round((x/(max(x))),3)))))

  ## Normalize checked columns
check.norm <- as.data.frame(t(apply(data.all2[14:18], 1, 
                                   function(x) round((x/(max(x)))*1000,2))))
  ## Add references back to normalized data
data.norm2<-bind_cols(data.norm,data.all2[1])
data.norm3<-bind_cols(data.norm2,data.all2[8])
data.norm4<-bind_cols(data.norm3,data.all2[c(9,10,13,19:24)])

check.norm2<-bind_cols(check.norm,data.all[1])
check.norm3<-bind_cols(check.norm2,data.all[8])
check.norm4<-bind_cols(check.norm3,data.all[9,10,19:24])


  ## Remove samples with Fit error = 32000 (user error)
data.norm4 <- data.norm4 %>%
  filter(fit_error !='32000')

check.norm4 <- check.norm4 %>%
  filter(fit_error !='32000')

  ## Subset 2019-2020 datapoints
data.subset<-subset(data.norm4, source == "2019" | source == "2020" | source == "Default" | source == "2018")
data.subset2<-subset(check.norm4, source == "2019" | source == "2020" | source == "2018" | source == "Default")

  
#######################################################################
  ##########################  HISTOGRAMS  #########################

  ## Histrogram of fit errors for each year
hist.1<-ggplot(data.subset2, aes(x = fit_error, fill = source)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Range of Fit Errors") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 18),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text = element_text(size = 18), legend.key=element_rect(fill='white'),
        legend.title = element_blank()) +
  facet_wrap(.~source, scales = "free_x")
hist.1

########################################################################
  #########################   BOXPLOTS
  ### Boxplot only lakes with consisitent sampling (not random testing) 
  ### and ignoring null values

  ### Convert z-off background fluorescence to numeric
data.subset$zoff_F =as.numeric(data.subset$zoff_F)
  ### filter our random lakes and lakes with 2 sample sites
  ### Use pipe (%>%) and (!) to do this and then graph boxplot
  ### this eliminates making another dataframe

data.subset %>% filter(!site %in% c("Pleasant Creek", "Honey Creek Resort",
                                    "Lacey-Keosauqua", "Lake Wapello",
                                    "Marble Beach", "Red Haw", "Crandall's",
                                    "Union Grove 10x Dilution", "Default",
                                    "Denison", "McIntosh Woods",
                                    "North Twin East")) %>%
  #ggplot(aes(x = fct_rev(as_factor(site)), y = zoff_F, fill = site)) +
  ggplot(aes(x = site, y = zoff_F, fill = site)) +
  geom_boxplot(color = "black", fill = "gray", alpha = 0.7) +
  scale_x_discrete(limits = rev) +  ## reverse order of flipped x axis
  labs(y = 'Background Fluoresence (%)') +
  coord_flip() +
  #scale_fill_manual(values = c("#66CCCC", "#FFFF00", "#99CC00", "#FF9900"),
                    #labels = c("'Blue' group", "'Brown' group",
                               #"'Green' group","'Red' group")) +
  #facet_wrap(.~order, scale = "free", ncol = 2) +
  #scale_fill_viridis_d(option = "turbo") +
  theme(panel.background = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(color = "white"))

########################################################################
  #########################  RUN nMDS  ############################
  
  ## 1) start mds: 
  ## Is there a difference in fluorescence of default fluorescence signatures
  ## vs lab calibrated fluorescence vs field samples?
  ### note:run more than once if first attempt is no convergence
set.seed(123)
  ## mds of default ref and cyano cultures
  ## Import PAM-ref
  ## mds.ref2 => k = 3 whereas mds.ref => k = 2
mds.ref2<-metaMDS(PAM_ref[, c(2:6)], distance = "bray", k = 3, 
                  maxit = 999, trace =2)
  ## mds of samples and default ref
mds.subset3<-metaMDS(data.subset[, c(1:5)], distance = "bray", k = 2,
                    maxit = 999, trace =2)
  ## mds of deconvolution F (check F) of samples and default ref
mds.subset2<-metaMDS(data.subset2[, c(1:5)], distance = "bray", k = 3,
                    maxit = 999, trace =2)
stressplot(mds.data)
mds.data
mds.subset
mds.subset3

  ## 2) Saving nMDS output
save(mds.subset3b, file = "mds_subset3b.rda")
save(mds.ref2, file = "mds_ref2.rda")

  ## 3) Plot nMDS output using base R
plot(mds.subset4)  ## Base R will automatically recognize ordination
plot(mds.data)

########################################################################
  #################  PREPARING OUTPUT FOR GGPLOT  ##################
  
  ## 1) Extract nMDS output into a dataframe 
  ## Use the score () to extraxt site scores and convert to data frame
mds.scores3<-as.data.frame(scores(mds.subset3))
mds.scores<-as.data.frame(scores(mds.data))

  ## 2) create solumn of site names from row names of meta.scores
mds.scores3$site<-rownames(mds.subset3)
mds.scores$site<-rownames(mds.scores)

  ## 3) add details to mds.scores dataframe
grp.fit<-round(data.subset$fit_error, digits = 0)  ## Round Fit error to 1 decimal place
mds.scores3$fit<-data.subset$fit_error
mds.scores3$source<-data.subset$source
mds.scores3$sample<-data.subset$sample
mds.scores3$location<-data.subset$site
mds.scores3$tot_chla<-data.subset$tot_chla
mds.scores3$cyano_chla<-data.subset$cyano_chla
mds.scores3$green_chla<-data.subset$green_chla
mds.scores3$brown_chla<-data.subset$brown_chla
mds.scores3$pe_chla<-data.subset$pe_chla
mds.scores3$class<-data.all2$fit_class

  
  ## 4) Extract Species scores into dataframe 
  ## Use score () to extract species score from mds output 
  ## and convert to data frame
species.score<-as.data.frame(scores(mds.data, "species"))

  ## 5) create columns of species from the species score dataframe
species.score$species<-rownames(species.score)
  ##check species dataframe
head(species.score)

  ## 6) Create polygons for default (factory) reference spectra
  ## Renaming default references
mds.scores3$sample[which(mds.scores3$sample == "syleo")] <- "'Blue' group"  #Synechococcus leopoliensis.
mds.scores3$sample[which(mds.scores3$sample == "chlorella")] <- "'Green' group" #Chorella vulgaris"
mds.scores3$sample[which(mds.scores3$sample == "phaeo")] <- "'Brown' group" #Phaeodactylum tricornutum"
mds.scores3$sample[which(mds.scores3$sample == "crypto")] <- "'Red' group" #Cryptomonas ovata"
  ## hull values for grp Default
grp.default3<-mds.scores3[mds.scores3$source=="Default",][chull(mds.scores3[mds.scores3$source=="Default", 
                                             c("NMDS1", "NMDS2")]),]
 ## hull values for grp Customized
grp.cust<-mds.scores3[mds.scores3$grp.source == "non_Default",][chull(mds.scores3[mds.scores3$grp.source == "non_Default", 
                                             c("NMDS1", "NMDS2")]),]
 ## hull values for mixed cultures
grp.mix<-mds.scores[mds.scores$grp.source == "mixed",][chull(mds.scores[mds.scores$grp.source == "mixed", 
                                            c("NMDS1", "NMDS2")]),]

 ## combine grp default and grp cust
hull.data<-rbind(grp.default,grp.cust)
hull.data
  
  ## Do same for ref nmds
  ## Add PAM-ref details to mds.scores nmds output
mds.scores$sample<-PAM_ref$sample
mds.scores$source<-PAM_ref$source

mds.scores$sample[which(mds.scores$sample == "syleo")] <- "'Blue' group"  #Synechococcus leopoliensis.
mds.scores$sample[which(mds.scores$sample == "chlorella")] <- "'Green' group" #Chorella vulgaris"
mds.scores$sample[which(mds.scores$sample == "phaeo")] <- "'Brown' group" #Phaeodactylum tricornutum"
mds.scores$sample[which(mds.scores$sample == "crypto")] <- "'Red' group" #Cryptomonas ovata"

  ## Save mds scores dataframe
save(mds.scores3, file = "mds_scores3.rda")

  ## Do same for ref and cyano mds output
  ## 1) Extract nMDS output into a dataframe 
ref2.scores<-as.data.frame(scores(mds.ref2))
  ## 2) create solumn of site names from row names of meta.scores
ref2.scores$site<-rownames(mds.ref)
  ## 3) add details to mds.scores dataframe
ref2.scores$source<-PAM_ref$source
ref2.scores$sample<-PAM_ref$sample
  ## 4) rename refs
ref2.scores$sample[which(ref.scores$sample == "syleo")] <- "'Blue' group"  #Synechococcus leopoliensis.
ref2.scores$sample[which(ref.scores$sample == "chlorella")] <- "'Green' group" #Chorella vulgaris"
ref2.scores$sample[which(ref.scores$sample == "phaeo")] <- "'Brown' group" #Phaeodactylum tricornutum"
ref2.scores$sample[which(ref.scores$sample == "crypto")] <- "'Red' group" #Cryptomonas ovata"

  ## Save ref mds scores dataframe
save(ref2.scores, file = "ref_scores3.rda")


########################################################################
  ###################  Plot nMDS using ggplot  #####################

  ## -------------- Plot nMDS of references and cyanos  -----------------
  ## hulling only default refs, alpha = transparency w/ 0 being transparent
  ## Use grp.ref for references 
  ## subset for cyanos as grp.cust (group customized)
grp.ref<-ref.scores[ref.scores$source=="Default",][chull(ref.scores[ref.scores$source=="Default",
                                                                    c("NMDS1", "NMDS2")]),]
grp.cust<-subset(ref.scores, source == "non-Default")
  ## plot nmds with ggplot
p1<-ggplot() +
  # this adds default refs scores
  geom_polygon(data = grp.ref, aes(x = NMDS1, y = NMDS2, group = source), 
               fill = NA, color = "gray", size = 1) +
  geom_point(data = grp.ref, aes(x=NMDS1, y=NMDS2), size = 2) +
  geom_text(data = grp.ref, aes(x = NMDS1, y = NMDS2, label = sample), 
            size = 3.5, nudge_x = 0.2) +
  # this adds cyanos scores
  geom_point(data = grp.cust, aes(x=NMDS1, y=NMDS2, shape = sample), 
             size = 3) +
  #geom_text(data = grp.cust, aes(x = NMDS1, y = NMDS2, label = sample), 
            #size = 3.5, nudge_x = 0.2) +
  scale_shape_manual(values = c(0:9)) + ## assign multiple shapes
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.key=element_rect(fill='white'),
        legend.title = element_blank(),
        legend.position = "na")
p1

  ## 3D plot for nmds (because k =3)

plot_ly(ref2.scores, x=ref2.scores$NMDS1, 
        y=ref2.scores$NMDS2, 
        z=ref2.scores$NMDS3, 
        type="scatter3d", mode="markers", color=source)
ggplot(ref2.scores, aes(x=NMDS1, y=NMDS2, z=NMDS3, color=source)) + 
  theme_void() +
  axes_3D() +
  stat_3D()
 ## --------------  Plot references and samples  -------------------
 ## Subset 2018 data from the mds output
mds.2019<-subset(mds.scores3, grp.source == "2019")
mds.2020<-subset(mds.scores3, grp.source == "2020")
mds.2018<-subset(mds.scores3, grp.source == "2018")
mds.sample<- mds.scores3 %>%
  filter(source !='Default')
mds.zerogreen<-subset(mds.scores3, green_chla == "0")
mcintosh.2018<-subset(mds.2018, location == "McIntosh Woods")
beeds.2018<-subset(mds.2018, location == "Beed's Lake")

  ## subset by fit error
fit0<-subset(mds.scores3, grp.fit == 0)
fit1<-subset(mds.scores3, grp.fit == 1)
fit2<- subset(mds.scores3, grp.fit == 2)
fit3<-subset(mds.scores3, grp.fit == 3)
fit4<-subset(mds.scores3, grp.fit > 3)

p3<-ggplot() +
  geom_polygon(data = grp.default3, 
               aes(x = NMDS1, y = NMDS2, group = grp.source), 
               fill = NA, color = "gray", size = 1) +
  geom_point(data = grp.default3, aes(x=NMDS1, y=NMDS2), size =2) +
  #geom_point(data = mds.sample, aes(x=NMDS1, y=NMDS2, color = class), size =2) +
  #geom_point(data = mds.2020, aes(x=NMDS1, y=NMDS2), size = 6, color = "#9999FF") +
  geom_point(data = beeds.2018, aes(x=NMDS1, y=NMDS2, color = grp.fit), size = 6) +
  #geom_point(data = mds.2018, aes(x=NMDS1, y=NMDS2), color = "#33CC99", size = 6) +
  #geom_point(data = mds.2019, aes(x=NMDS1, y=NMDS2), size = 6, color = "#FF9900") +
  #geom_text(data = grp.cust, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 5, vjust = 0, nudge_x = 0.07) +
  geom_text(data = grp.default3, aes(x = NMDS1, y = NMDS2, label = sample), size = 4, vjust = 0, nudge_x = 0.07) +
  #geom_text(data = mds.2020, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = mds.2018, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = mds.2019, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  geom_text(data = beeds.2018, aes(x = NMDS1, y = NMDS2, label = sample), size = 3, vjust = 0.1, nudge_x = 0.05) +
  scale_colour_viridis(option = "plasma") + 
  #scale_color_manual(values = c("#de4968", "#929596"), #de4968
                     #label = c("> 1", "0")) +
  #facet_wrap(.~ location) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.key=element_rect(fill='white'),
        legend.position = c(0.9,0.9),
        legend.title = element_blank())
p3

  ## Turn on interactive graph
library(plotly)
ggplotly(p3)
  ## Subset and plot green valley 
mds.gval<-subset(mds.scores3, location == "Green Valley")
mds.twin<-subset(mds.scores3, location == "North Twin East" | location == "North Twin West")
p5<-ggplot() +
  geom_polygon(data = grp.default3, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  geom_point(data = grp.default3, aes(x=NMDS1, y=NMDS2), color = "black", size =7) +
  geom_point(data = mds.twin, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 7) +
  geom_text(data = mds.twin, aes(x = NMDS1, y = NMDS2, label = sample), size = 5, vjust = 0, nudge_x = 0.01) +
  geom_text(data = grp.default3, aes(x = NMDS1, y = NMDS2, label = sample), size = 7, vjust = 0, nudge_x = 0.06) +
  scale_color_viridis_d(option = "plasma", direction = -1, breaks = c("2018", "2019", "2020")) + 
  #scale_color_viridis(option = "plasma", direction = -1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p5

ggplotly(p5)

## Plot references and 2019 samples
## Subset 2018 data from the mds output
mds.2019<-subset(mds.scores, grp.source == "2019")

p5<-ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  #geom_point(data = hull.data, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_point(data = mds.2019, aes(x=NMDS1, y=NMDS2, color = week), size = 3) +
  geom_text(data = mds.2019, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 4, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = hull.data, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  coord_equal() +
  scale_color_viridis_d() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())

p6<- p5 + facet_wrap(.~location, ncol = 3)
p6


## Plot Black Hawk Lake 2018 + 2019 samples
## Subset black Hawk data from the mds output
mds.chase<-subset(mds.scores, location %in% c('Brushy Creek','North Twin East','North Twin West',
                                                  'Black Hawk','Denison'))

p7<-ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  #geom_point(data = hull.data, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_point(data = mds.chase, aes(x=NMDS1, y=NMDS2, color = week), size = 3) +
  geom_text(data = mds.chase, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 4, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = hull.data, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  coord_equal() +
  scale_color_viridis_d() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p8 <- p7 + facet_wrap(.~location, ncol = 2)
p8


  ##2018-2019 samaples
mds.all<-subset(mds.scores, grp.source %in% c('2018','2019'))
p9<-ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  #geom_point(data = hull.data, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_point(data = mds.all, aes(x=NMDS1, y=NMDS2, color = month), size = 4) +
  geom_text(data = mds.all, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = hull.data, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  coord_equal() +
  scale_color_viridis_d(breaks=c("May", "June", "July", "August")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.key=element_rect(fill='white'),
        legend.title = element_blank(), strip.text = element_text(size = 16)) 
p9

p10<- p9 + facet_wrap(.~location, ncol = 3)
p10

 ##Subset for green valley
mds.gv<-subset(mds.scores, location == "Green Valley")
p11<-ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  #geom_point(data = hull.data, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_point(data = mds.gv, aes(x=NMDS1, y=NMDS2, color = month), size = 7) +
  #geom_text(data = mds.gv, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = hull.data, aes(x = NMDS1, y = NMDS2, label = site), size = 4, vjust = 0, nudge_x = 0.01) +
  #coord_equal() +
  scale_color_viridis_d(breaks=c("May", "June", "July", "August")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 22),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=22), axis.title = element_text(size = 22),
        legend.text = element_text(size = 22), legend.key=element_rect(fill='white'),
        legend.title = element_blank(), strip.text = element_text(size = 22)) 
p11

p10<- p9 + facet_wrap(.~location, ncol = 3)
p10

  ## Clustering: Determine the optimal number of clusters
  ## the NbClust package used to do this
  ## use the scaled or normalized data 
nb<-NbClust(data.norm[,c(1:5)], distance = "euclidean", min.nc = 2, max.nc = 10,
            method = "average", index = "all")

  ## K-means Clustering Analysis
km.res<-eclust(data.norm[,c(1:5)], "kmeans", k = 2,
               nstart = 25, graph = FALSE)
  ## visualize k-means cluster
  ## Notes: ellipse.alpha = ellipse color transparency, show.clus.cent = cluster center
p9 <- fviz_cluster(km.res, geom = "point", ellipse.type = "norm",
                   pointsize = 4, ellipse.alpha = 0, ellipse.width = 2,
                   show.clust.cent = FALSE) +
  scale_color_viridis_d() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 16), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p9  

  ## Dendrogram of cluster analysis
  ## similar samples will be closer together
  ## First, create new dataframe with normalized/scaled values and sample info
data.hc<-cbind(data.norm,data.all[,c(8:12)])
data.hc<-cbind(data.hc,data.all$sample)
  ## Start hierarchical clustering analysis
res.hc <- eclust(data.hc[,c(1:5)], "hclust", k = 2, method = "average", graph = FALSE)
plot(res.hc) #view hclust results
  ## Make dendrogram from hclust results
dend<-as.dendrogram(res.hc)
plot(dend)
  ## this is very difficult to see because there are so many samples
par(mfrow = c(1,1)) #2,1 = two trees laying on top of each other
plot(dend, main = "Original dendrogram")
dend_list <- get_subdendrograms(dend, 2) #creates list of 2 clusters
sapply(dend_list, plot) #applies list of clusters to a plot

## cut the dendrogram into different groups
grp<-as.data.frame(cutree(dend, k = 2, order_clusters_as_data = FALSE))
grp #check order of samples and grouped clusters
write.csv(grp,'grp_cluster.csv') #save this as a .csv file
order_clust=read.csv("grp_cluster.csv",header=T) #bring it back in

IDs<- data.all %>% mutate(order = row_number())
order_clust <-order_clust %>% mutate(row_name = row_number())
final<- merge(order_clust, IDs, by = "row_name") #Joins IDs and order_clust by "row_name" 

test.trans<-as.data.frame(t(data.all))
test<- test.trans %>% mutate(order = rownames(test.trans))


############################################################################################
  ##################    MAKING 3D NMDS WITH GGPLOT    #################

   ### This uses plot3D package so turn it on if you haven't
   ### this package is follows the same coding as ggplot
library(rgl)  ## loads package for 3D plotting

   ### Set the x,y,z coordinates
   ### For this plot, will use the nmds scores since 
   ### stress was lowest at 3 dimensions (because 2D did not converge)

   ### Create df of ref scores
default.3d<-mds.scores3 %>%  
  subset(source == "Default")
   ### create df without references
mds.3d<-mds.scores3 %>%  
  subset(!source %in% "Default")
   ### Assign dimensions to x,y,z
x <- mds.3d$NMDS1
y <- mds.3d$NMDS2
z <- mds.3d$NMDS3
   

   ### Basic scatter plot
rgl.open()
rgl.bg(color = "white") # Setup the background color
rgl.points(x, y, z, size = 5,
           color = get_colors(mds.3d$fit)) #color = "blue"

   ### Color sample points by group
   ### Create function
#' Get colors for the different levels of 
#' a factor variable
#' 
#' @param groups a factor variable containing the groups
#'  of observations
#' @param colors a vector containing the names of 
#   the default colors to be used
get_colors <- function(groups, group.col = palette()){
  #groups <- as.factor(groups)
  #ngrps <- length(levels(groups))
  #if(ngrps > length(group.col)) 
    #group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
# Add x, y, and z Axes
rgl.lines(c(0, 1), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(0,0.5), c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), c(0,1), color = "green")

  ### Add plane (use references)
# define function to calculate axes limit
lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1}
# xlab, ylab, zlab: axis labels
# show.plane : add axis planes
# show.bbox : add the bounding box decoration
# bbox.col: the bounding box colors. The first color is the
# the background color; the second color is the color of tick marks
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add plane
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = c(-0.3526627, 0.3076674,0.4029682,0.1305049), 
             y = c(-0.05494934, -0.1255089, 0.03109923, 0.15533561),
             z = c(0.02211692, -0.01647502, 0.03821252, 0.0314686))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}

rgl_add_axes(x, y, z)




   ### clab is used to change legend title
   ### The points are automatically colored according to variable Z
colVar <- sapply(mds.3d$fit,function(a){ifelse(a==0,'gray','red')})
colVar <- factor(colVar,levels=c('gray','red'))
   ### viewing angle: theta = azimuthal direction (rotate up/down?)
   ### viewing angle: phi = co-latitude (horizontal??)
   ### Use extreme for each and then 45 for both to test it out

scatter3D(x=x,y=y,z=z, cex = 1.3,
          phi = 45,theta = 45, ticktype = "detailed",
          xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
          clab = c("Least", "Square Fit"),
          colvar=as.integer(colVar),
          colkey=list(at=c(0,1),side=4),
          col=as.character(levels(colVar)),
          pch=19)

   ### Add ref scores to 3D plot
   ### Cyano ref
scatter3D(x = -0.3526627, y = -0.05494934, z = 0.02211692, add = TRUE, colkey = FALSE, 
          pch = 18, cex = 3, col = "black")

  #############   BARPLOT: CHLA COMMUNITY COMPOSITON   ##############
  
   ### select only chla data from PhytoPAM for year 2018
chla.comm<-data.all2[c(5:151),c(1,9,10,11,20:23)]

   ### Pivot longer to make column for taxa and chla values
chla.comm2<-pivot_longer(chla.comm,5:8,names_to = "taxa", values_to = "chla")


   ### Plot stacked barplot
ggplot(chla.comm2, aes(fill = taxa, y = chla, x = week)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(.~site)
