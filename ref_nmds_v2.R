library(ggplot2)
library(GGally)
library(dplyr)
library(vegan)
library(viridis)
install.packages("ggiraph")
library(ggiraph)


  ## Use 2018 fluorescence data with references all in one dataset. 
  ## Then plot the ordination to see how the fluorescence samples 2018 
  ## compare with the references (default and lab measured)

  ##Read data----from Fluorescence Ref Folder
data.ref=read.csv("PAM references.csv",header=T)
data.all=read.csv("all_data_2018-2019.csv",header=T)
data.2019=read.csv("all_data_2019.csv",header=T)

  ## compare leds
ggpairs(data.ref, columns = c("led440","led480","led540","led590","led625"),
        diag = list(continuous = "density"), axisLabels = "show")
  ## create group variable (assigns phytoplankton type to each row)
grp.ref<-data.ref$sample
grp.source<-data.ref$source
 
  ## start mds: Is there a difference in fluorescence of default fluorescence signatures
  ## vs lab calibrated fluorescence signatures?
mds.ref<-metaMDS(data.ref[, c(2:6)], distance = "bray", k = 2, maxit = 999)
stressplot(mds.ref)
mds.ref

  ##Use the score () to extraxt site scores and convert to data frame------
mds.scores<-as.data.frame(scores(mds.ref))
  ##create solumn of site names from rowmanes of meta.scores
mds.scores$site<-rownames(mds.scores)
  ## add group variables created
mds.scores$grp.ref<-grp.ref
mds.scores$grp.source<-grp.source

  #look at the data
head(mds.scores)
  ## Use score () to extract species score and convert to data frame
species.score<-as.data.frame(scores(mds.ref, "species"))
  ## create columns of species from the species score dataframe
species.score$species<-rownames(species.score)
  ##check species dataframe
head(species.score)
  ## hull values for grp Default
grp.default<-mds.scores[mds.scores$grp.source=="Default",][chull(mds.scores[mds.scores$grp.source=="Default", 
                                             c("NMDS1", "NMDS2")]),]
 ##hull values for grp Customized
grp.cust<-mds.scores[mds.scores$grp.source == "non-Default",][chull(mds.scores[mds.scores$grp.source == "non-Default", 
                                             c("NMDS1", "NMDS2")]),]
##hull values for mixed cultures
grp.mix<-mds.scores[mds.scores$grp.source == "mixed",][chull(mds.scores[mds.scores$grp.source == "mixed", 
                                            c("NMDS1", "NMDS2")]),]

 ## combine grp default and grp cust
hull.data<-rbind(grp.default,grp.cust,grp.mix)
hull.data
  
## Begin plot using ggplot--hulling only default refs, alpha = transparency w/ 0 being transparent
p1<-ggplot() +
  geom_polygon(data = grp.default, aes(x = NMDS1, y = NMDS2, group = grp.source,), fill = c("gray"), alpha = 0.5) +
  geom_point(data = grp.default, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_text(data = grp.default, aes(x = NMDS1, y = NMDS2, label = grp.ref), size = 5, nudge_x = 0.1) +
  coord_equal() +
  scale_color_viridis_d() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p1

  ##subset of non-default 
cust.scores<-subset(mds.scores, grp.source == "non-Default")
  ## Begin plot using ggplot--hulling only non-default refs, but there is error in the chull
  ## chull for non-default is only hulling 7 species and not all of them.
p2<-ggplot() +
  geom_polygon(data = grp.cust, aes(x = NMDS1, y = NMDS2, group = grp.source,), fill = c("gray"), alpha = 0.5) +
  geom_point(data = cust.scores, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_text(data = cust.scores, aes(x = NMDS1, y = NMDS2, label = grp.ref), size = 5, vjust = 0, nudge_x = 0.1) +
  coord_equal() +
  scale_color_viridis_d() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p2

## Begin plot using ggplotplot both default and customized references
p3<-ggplot() +
  geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.source), fill = c("gray"), alpha = 0.5) +
  geom_point(data = mds.scores, aes(x=NMDS1, y=NMDS2, color = grp.source), size = 3) +
  geom_text(data = mds.scores, aes(x = NMDS1, y = NMDS2, label = grp.ref), size = 3, vjust = 0, nudge_x = 0.0) +
  coord_equal() +
  scale_color_viridis_d() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p3
p4 <- p3 + geom_point_interactive(aes(tooltip = model, data_id = model), size = 2)  
x<-girafe(code = print(p4))

identify(x = hull.data$NMDS1, y = hull.data$NMDS2)

 ##
plot(mds.ref, display = "sites")
text(nmds1, nmds2, row.names(mds.scores), cex=0.6, pos=4, col="red")
ordihull(mds.ref, grp.source,lty = 2, col = "red", label = T)
## which distance method is best for the 2018 dataset------
data<-na.omit(data)
data<- data[data$`440nm` >= 0, ]
data.2018f<-data[c(1:162),c(6:10)]
rank.totus <- rankindex(as.matrix(ref2), data.2018f, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.totus, decreasing = TRUE)[1]), "method."))


## create group variable (assigns fit error to each row)
grp.error<-data$fit_error
grp.type<-data$grp

  ## standardize data from 0 to 1 since there are negative values
data.norm<- as.data.frame(apply(data[, 6:10], 2, 
                                 function(x) (x - min(x))/(max(x)-min(x))))
data.norm<-cbind(data[,c(1:5,11,12)],data.norm[,c(1:5)])
## start mds 
data<-data[c(11,1:10,12)]
mds.2018<-metaMDS(data[c(1:151), c(7:11)], distance = "euclidean", k = 2, maxit = 999, wascores = TRUE)
mds.2018
stressplot(mds.2018)
summary(mds.2018)

##Use the score () to extraxt site scores and convert to data frame------
mds.2018scores<-as.data.frame(scores(mds.2018))
##create solumn of site names from rowmanes of meta.scores
mds.2018scores$site<-rownames(mds.2018scores)
## add group variables created
mds.2018scores$grp.error<-grp.error
mds.2018scores$grp.type<-grp.type
#look at the data
head(mds.2018scores)
## Use score () to extract species score and convert to data frame
species.2018score<-as.data.frame(scores(mds.2018, "species"))
## create columns of species from the species score dataframe
species.2018score$species<-rownames(species.2018score)
##check species dataframe
head(species.score)
## hull values for grp Default
grp.default<-mds.2018scores[mds.2018scores$grp.type=="Default",
                        ][chull(mds.2018scores[mds.2018scores$grp.type=="Default", 
                                           c("NMDS1", "NMDS2")]),]
##hull values for grp non-Default
grp.nondefault<-mds.2018scores[mds.2018scores$grp.type=="non-Default",
                     ][chull(mds.2018scores[mds.2018scores$grp.type=="non-Default", 
                                        c("NMDS1", "NMDS2")]),]
##hull values for grp Sample
grp.sample2018<-mds.2018scores[mds.2018scores$grp.type=="Sample_2018",
                               ][chull(mds.2018scores[mds.2018scores$grp.type=="Sample_2018", 
                                                      c("NMDS1", "NMDS2")]),]

## combine grp default and grp cust
hull.data<-rbind(grp.default,grp.nondefault)
hull.data
hull.data<-rbind(hull.data,grp.sample2018)
hull.data

## Begin plot using ggplot
p1<-ggplot() +
  geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = grp.type,), fill = c("gray"), alpha = 0.5) +
  geom_point(data = mds.2018scores, aes(x=NMDS1, y=NMDS2, color = grp.type), size = 3) +
  #geom_text(data = mds.2018scores, aes(x = NMDS1, y = NMDS2, label()), size = 6, vjust = 0, nudge_x = 0.2) +
  coord_equal() +
  scale_color_viridis_d() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.key=element_rect(fill='white'),
        legend.title = element_blank())
p1


plot(mds.ref, display = "sites")