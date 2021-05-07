
library(vegan) 
library(viridis)
library(NbClust) #for finding optimal clusters
library(cluster) #for computing cluster silhouettes
library(tidyverse)
library(factoextra) #computing heirarchical cluster

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
                                   (function(x) round((x/(max(x)))*1000,3)))))

  ## Normalize checked columns
check.norm <- as.data.frame(t(apply(data.all[13:17], 1, 
                                   function(x) round((x/(max(x)))*1000,2))))
  ## Add references back to normalized data
data.norm2<-bind_cols(data.norm,data.all2[1])
data.norm3<-bind_cols(data.norm2,data.all2[8])
data.norm4<-bind_cols(data.norm3,data.all2[c(9,10,19)])

check.norm2<-bind_cols(check.norm,data.all[1])
check.norm3<-bind_cols(check.norm2,data.all[8])
check.norm4<-bind_cols(check.norm3,data.all[9,10])


  ## Remove samples with Fit error = 32000 (user error)
data.norm4 <- data.norm4 %>%
  filter(fit_error !='32000')

check.norm4 <- check.norm4 %>%
  filter(fit_error !='32000')

  ## Subset 2019-2020 datapoints
data.subset<-subset(data.norm4, source == "2019" | source == "2020" | source == "Default" | source == "2018")
data.subset2<-subset(check.norm4, source == "2019" | source == "2020" | source == "2018" | source == "Default")

  ## Remove 
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
  #########################  RUN nMDS  ############################
  
  ## 1) start mds: 
  ## Is there a difference in fluorescence of default fluorescence signatures
  ## vs lab calibrated fluorescence vs field samples?
  ### note:run more than once if first attempt is no convergence
set.seed(123)
  ## mds of samples, default ref, cyano cultures
mds.data<-metaMDS(PAM_ref[, c(2:6)], distance = "bray", k = 2, 
                  maxit = 999, trace =2)
  ## mds of samples and default ref
mds.subset3<-metaMDS(data.subset[, c(1:5)], distance = "bray", k = 3,
                    maxit = 999, trace =2)
  ## mds of deconvolution F (check F) of samples and default ref
mds.subset2<-metaMDS(data.subset2[, c(1:5)], distance = "bray", k = 3,
                    maxit = 999, trace =2)
stressplot(mds.subset)
mds.data
mds.subset
mds.subset2

  ## 2) Saving nMDS output
save(mds.subset3, file = "mds_subset3.rda")

  ## 3) Plot nMDS output using base R
plot(mds.subset4)  ## Base R will automatically recognize ordination
rm(mds.subset)

########################################################################
  #################  PREPARING OUTPUT FOR GGPLOT  ##################
  
  ## 1) Extract nMDS output into a dataframe 
  ## Use the score () to extraxt site scores and convert to data frame
mds.scores3<-as.data.frame(scores(mds.subset3))

  ## 2) create solumn of site names from row names of meta.scores
mds.scores3$site<-rownames(mds.subset3)

  ## 3) add details to mds.scores dataframe
grp.fit<-round(data.subset$fit_error, digits = 0)  ## Round Fit error to 1 decimal place
mds.scores3$grp.fit<-data.subset$fit_error
mds.scores3$grp.source<-data.subset$source
mds.scores3$sample<-data.subset$sample
mds.scores3$location<-data.subset$site
mds.scores3$tot_chla<-data.subset$tot_chla

#mds.scores$week<-as.factor(week)
#mds.scores$month<-month
  
 ##look at the data
head(mds.scores)

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
grp.default3<-mds.scores3[mds.scores3$grp.source=="Default",][chull(mds.scores3[mds.scores3$grp.source=="Default", 
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
  
  

########################################################################
  ###################  Plot nMDS using ggplot  #####################

  ## -------------- Plot nMDS of references only  -----------------
  ## hulling only default refs, alpha = transparency w/ 0 being transparent
p1<-ggplot() +
  geom_polygon(data = grp.default, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = NA, alpha = 0.5) +
  geom_point(data = grp.default, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 5) +
  geom_text(data = grp.default, aes(x = NMDS1, y = NMDS2, label = sample), size = 7, nudge_x = 0.1) +
  scale_color_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p1

  ## subset of non-default 
cust.scores<-subset(mds.scores, grp.source == "non-Default")
  ## Begin plot using ggplot--hulling only non-default refs, but there is error in the chull
  ## chull for non-default is only hulling 7 species and not all of them.
p2<-ggplot() +
  geom_polygon(data = grp.cust, aes(x = NMDS1, y = NMDS2, group = grp.source,), fill = c("gray"), alpha = 0.5) +
  geom_point(data = cust.scores, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_text(data = cust.scores, aes(x = NMDS1, y = NMDS2, label = sample), size = 5, vjust = 0, nudge_x = 0.1) +
  coord_equal() +
  scale_color_viridis_d() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p2

 ## --------------- Plot references and samples  -------------------
 ## Subset 2018 data from the mds output
mds.2019<-subset(mds.scores3, grp.source == "2019")
mds.2020<-subset(mds.scores3, grp.source == "2020")
mds.2018<-subset(mds.scores3, grp.source == "2018")
mds.sample<- mds.scores3 %>%
  filter(grp.source !='Default')

  ## subset by fit error
fit0<-subset(mds.scores3, grp.fit == "0")
fit1<-subset(mds.scores3, grp.fit == 1)
fit2<- subset(mds.scores3, grp.fit == 2)
fit3<-subset(mds.scores3, grp.fit == 3)
fit4<-subset(mds.scores3, grp.fit > 3)

p3<-ggplot() +
  geom_polygon(data = grp.default3, 
               aes(x = NMDS1, y = NMDS2, group = grp.source), 
               fill = "gray", alpha =0.3, linetype = 2) +
  geom_point(data = grp.default3, aes(x=NMDS1, y=NMDS2), size =9) +
  #geom_point(data = mds.2020, aes(x=NMDS1, y=NMDS2), size = 6, color = "#9999FF") +
  geom_point(data = mds.sample, aes(x=NMDS1, y=NMDS2, color = tot_chla), size = 6) +
  #geom_point(data = mds.2018, aes(x=NMDS1, y=NMDS2), color = "#33CC99", size = 6) +
  #geom_point(data = mds.2019, aes(x=NMDS1, y=NMDS2), size = 6, color = "#FF9900") +
  geom_text(data = grp.default3, aes(x = NMDS1, y = NMDS2, label = sample), size = 7, vjust = 0, nudge_x = 0.07) +
  #geom_text(data = mds.2020, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = mds.2018, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  #geom_text(data = mds.2019, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
  geom_text(data = mds.sample, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0.1, nudge_x = 0.01) +
  scale_colour_viridis_c(option = "turbo") + 
  scale_shape_manual(values = c(8:14)) + ## assign multiple shapes
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.key=element_rect(fill='white'),
        legend.position = "right",
        legend.title = element_blank())
p3
p4<- p3 + facet_wrap(.~location)
p4

  ## Subset and plot green valley 
mds.gval<-subset(mds.scores3, location == "Green Valley")
mds.twin<-subset(mds.scores3, location == "North Twin East" | location == "North Twin West")
p5<-ggplot() +
  geom_polygon(data = grp.default3, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  geom_point(data = grp.default3, aes(x=NMDS1, y=NMDS2), color = "black", size =7) +
  geom_point(data = mds.twin, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 7) +
  geom_text(data = mds.twin, aes(x = NMDS1, y = NMDS2, label = grp.fit), size = 7, vjust = 0, nudge_x = 0.01) +
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
  ###############   TESTING FOR SIG DIFF BTWN GROUPS   ################


