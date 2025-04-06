# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Script: Landscape habitat amount x Patch size ####
# # # # # # # # # by Ricardo Spillari # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

### Library required packages ###
library(sp)
library(sf)
library(raster)
library(vegan)
library(betapart)
library(bipartite)
library(ggplot2)
library(ggpubr)
library(TITAN2)
library(ggrepel)
library(reshape2)
library(glmmTMB)
library(MuMIn)

# # # # # # # # # # # #
#### Data Fitting #####
# # # # # # # # # # # #

#### Response Variables ####

# Ant data Vasconcelos & Delabie 2000 # 
# Updated by Ricardo Spillari 2025 #

## Richness ##
ants<-read.csv("Cephalotes_data.csv",h=T, sep = ",")
head(ants)
# Species Abundance/Frequency Table
abund_table<-table(ants$Fragment,ants$Species)
write.csv(abund_table,"Abundance_table.csv",row.names = F) # Save abundance table #
# Species Occurrence Table
occ <- decostand(abund_table, "pa")
write.csv(occ,"Occurence_table.csv",row.names = F) # Save occurence table #
rowSums(occ)->Richness # Species Richness - Total number of species #
# Habitat-use guilds richness - Forest specialists
subset(ants,ants$Habitat_guild=="Forest")->fores
table(fores$Fragment,fores$Species)->abund_fores
occ_fores <- decostand(abund_fores, "pa")
rowSums(occ_fores)->Richness_forest # Species richness of forest specialist #
# Habitat-use guilds richness - Generalists of Habitat
subset(ants, ants$Habitat_guild=="Generalist")->gen
table(gen$Fragment,gen$Species)->abund_gen
occ_gen <- decostand(abund_gen, "pa")
rowSums(occ_gen)->Richness_generalist # Species richness of generalist of habitat #
# Habitat-use guilds richness - Open habitat specialists
subset(ants,ants$Habitat_guild=="Open habitat")->op
table(op$Fragment,op$Species)->abund_op
occ_op <- decostand(abund_op, "pa")
rowSums(occ_op)->Richness_opens
Richness_open<-c(1,3,3,2,1,1,1,0,1) # Species richness of open habitat specialist #

## Fragment ID ##
ID<-rownames(abund_table)

## Composition - PCoA ##
vegdist(occ,"jaccard")->dist_jac # Using occurrence data 
betadisper((dist_jac),ID)->pcoa
pcoa$vectors[,1]->eixo1
pcoa$vectors[,2]->eixo2 
# Species Scores PCoA 
wascores(eixo1, occ)->species_scores
wascores(eixo2, occ)->species_scores2
data.frame(PCoA1=species_scores,PCoA2=species_scores2)->spp_scores
ants[match(rownames(spp_scores),ants$Species),]$Habitat_guild->guildas
guildas[is.na(guildas)]<-"Morphospecies"
data.frame(spp_scores,guildas=guildas)->spp_scores # Position of species in the PCoA plot

## Response Variable Table ##
respvars_table<-cbind(ID,Richness,Richness_forest,Richness_generalist,Richness_open,eixo1,eixo2)
resp_table<-as.data.frame(respvars_table,row.names = F)
colnames(resp_table)<-c("Fragment","Richness","Richness_forest","Richness_generalist","Richness_open","PCoA1","PCoA2")
head(resp_table) # Response variables table #
write.csv(resp_table,"Response_table.csv",row.names = F) # Save #
as.numeric(resp_table$Richness)->resp_table$Richness
as.numeric(resp_table$Richness_forest)->resp_table$Richness_forest
as.numeric(resp_table$Richness_generalist)->resp_table$Richness_generalist
as.numeric(resp_table$Richness_open)->resp_table$Richness_open
as.numeric(resp_table$PCoA1)->resp_table$PCoA1
as.numeric(resp_table$PCoA2)->resp_table$PCoA2

#### Predictor variables ####

## Landscape habitat amount ##
# Sample points coordinates 
coordenadas<-read.csv("Centroids_reserves.csv",h=T, sep = ",")
# Spatial data frame 
SpatialPointsDataFrame(coords = coordenadas[,2:3],data = coordenadas,proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))->coords_spatial
st_as_sf(coords_spatial)->coords_spatial
# Transform coordinates
st_transform(coords_spatial,"+proj=utm +zone=21 +south +ellps=GRS80 +units=m +no_defs")->coords_spatial_t
# Buffer max extent
buff<-st_buffer((coords_spatial),5000)
# Raster of land use e cover
cover<-raster("brasil_coverage_1993.tif")
# Cut raster in the max extent 
rastercortado<-crop(cover,extent(buff))
# Transform raster
projectRaster(rastercortado,crs = "+proj=utm +zone=21 +south +ellps=GRS80 +units=m +no_defs")->rastercortado_t
# Diff size buffers
buffer_distance_meters <- c(500,1000,1500,2000,2500,3000,3500,4000)
# Landscape covers for each buffer size
porcentage<-{}
for (x in 1:9){
  for (y in 1:8){
    buff_shp <- st_buffer(coords_spatial_t[x,], buffer_distance_meters[y])
    crop(rastercortado_t,buff_shp)->crop2
    mask1<-mask(crop2,buff_shp)
    plot(crop2)
    plot(mask1)
    values(mask1)->val
    true<-is.na(val)
    val<-val[!true]
    floresta<-length(val[val==3])
    pastagem<-length(val[val==15])
    total<-length(val)
    forestcove<-floresta/total
    pastagemcover<-pastagem/total
    porcentage[[length(porcentage)+1]]<-data.frame(frag=coordenadas$frag[x],forest =forestcove,pastagem=pastagemcover, scale =buffer_distance_meters[y])
  }
}
landscape_data<-do.call(rbind.data.frame,porcentage)
write.csv(landscape_data,"Landscape_coves.csv", row.names = F)
# Correlation: Forest cover x Pasture cover
{}->value_cor
for (x in buffer_distance_meters) {
  landscape<-landscape_data[landscape_data$scale==x,]
  cor(landscape$forest,landscape$pastagem)->a
  a->value_cor[[length(value_cor)+1]]
}
value_cor # Correlation for each buff size #
write.csv(value_cor,"cor_forest_pasture.csv",row.names = F) # Save #
# Use only forest cover data
habitat<-landscape_data[,c(1,2,4)] # Landscape habitat amount for each sample plot in each buffer size #

## Patch size ##
patch<-rep(c(rep(1,8),rep(10,8),rep(100,8)),3)# Size of each fragment #

## Farm ##
farm<-c(rep("Esteio",24),rep("Dimona",24),rep("Porto Alegre",24)) # Farm where each fragment is located #
as.data.frame(farm)->farm1
as.factor(farm)->farm2

## Predictor Variables Table ##
pred_table<-cbind(habitat$frag,farm1,habitat$scale,habitat$forest,patch)
colnames(pred_table)<-c("Fragment","Farm","Scale","Forest","Size")
head(pred_table) # Predictor variables table #
write.csv(pred_table,"Predictor_table.csv",row.names = F) # Save #
as.data.frame(pred_table)->pred_table
pred_table$Farm<-farm2
# Predictor variables padronization
pred_pad<-scale(data.frame(pred_table$Forest,pred_table$Size))
pred_table_pad<-cbind(pred_table[1:3],pred_pad)
colnames(pred_table_pad)<-c("Fragment","Farm","Scale","Forest","Size")
head(pred_table_pad) # Predictor variables table patronized #
# Correlation: Landscape habitat amount x Patch size #
{}->value_cor2
for (x in buffer_distance_meters) {
  landscape<-pred_table_pad[pred_table_pad$Scale==x,]
  cor(landscape$Forest,landscape$Size)->a
  a->value_cor2[[length(value_cor2)+1]]
}
value_cor2 # Correlation for each buff size #
write.csv(value_cor2,"cor_forest_patch.csv",row.names = F) # Save #

# # # # # # # # # # # #
#### Data Analysis ####
# # # # # # # # # # # # 

# Effects of landscape habitat amount and patch size over total species richness and forest specialist richness #

# Richness ####

#### Scale of effect - Species Richness ####
values<-{}
for (h in buffer_distance_meters ) {
  paisagem<-pred_table_pad[pred_table_pad$Scale==h,]
  paisagem$Richness<-resp_table$Richness
  glmmTMB(Richness~Forest+Size+ (1|Farm), family = "gaussian",data=paisagem)->result
  AIC(result)->akai
  summary(result)[[6]]->coefs
  data.frame(coefs$cond,rsquare=r.squaredGLMM(result),aic=akai,scale=h)->values[[length(values)+1]]
}
do.call(rbind,values)->valuesR
valuesR$Pr...z..<0.05->valuesR$p ;valuesR 
write.csv(valuesR,"values_Richness.csv",row.names = T) # Save #
paisagem_rich<-pred_table_pad[pred_table_pad$Scale=="2000",] # Scale of effect = 2000

#### Final Model - Species Richness ####
cbind(paisagem_rich,Richness=resp_table$Richness)->table_richall
model_rich<-glmmTMB(Richness ~ Forest + Size +(1|Farm), family = "gaussian", data=table_richall) 
summary(model_rich)
r.squaredGLMM(model_rich)
# Correlations
# Forest amount x Farm
aov(paisagem_rich$Forest~paisagem_rich$Farm)->anova
summary(anova) 
TukeyHSD(anova)
# Richness x Farm
aov(table_richall$Richness~paisagem_rich$Farm)->resu
summary(resu)
TukeyHSD(resu)

# Composition ####

#### Scale of effect - Species Composition ####
values3<-{}
for (h in buffer_distance_meters ) {
  paisagem<-pred_table_pad[pred_table_pad$Scale==h,]
  paisagem$Composition<-resp_table$PCoA1
  glmmTMB(Composition~Forest+Size+(1|Farm), family = "gaussian",data=paisagem)->result
  AIC(result)->akai
  summary(result)[[6]]->coefs
  data.frame(coefs$cond,rsquare=r.squaredGLMM(result),aic=akai,scale=h)->values3[[length(values3)+1]]
}
do.call(rbind,values3)->valuesCOM 
valuesCOM$Pr...z..<0.05->valuesCOM$p; valuesCOM
write.csv(valuesCOM,"values_Composition.csv",row.names = T) # Save 
paisagem_com<-pred_table_pad[pred_table_pad$Scale=="1500",] # Scale of effect = 1500  

#### Final Models - Species Composition ####
cbind(paisagem_com,Composition=resp_table$PCoA1)->table_comp
model_comp<-glmmTMB(Composition ~ Forest + Size + (1|Farm) , family = "gaussian", data=table_comp) 
summary(model_comp)
r.squaredGLMM(model_comp)
# Correlations
# Forest amount x Farm
aov(paisagem_com$Forest~paisagem_com$Farm)->anova
summary(anova) 
TukeyHSD(anova)
# Composition x Farm
aov(table_comp$Composition~paisagem_com$Farm)->resu
summary(resu)
TukeyHSD(resu)

#### Nestedness x Turnover ####
occ # Occurrence table #
paisagem_com # Predictor table #

## NODF ##
ant.for <- occ[order(paisagem_com$Forest, decreasing = TRUE), ] # Rows ordered based on the landscape Forest amount #
ant.pat <- occ[order(paisagem_com$Size, decreasing = TRUE), ] # Rows ordered based on the patch Size #
null_for <- oecosimu(as.matrix(ant.for), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_for
null_pat <- oecosimu(as.matrix(ant.pat), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_pat

## Baselga ##
overal   <- beta.multi(occ, index.family="jaccard") # Multisite #
pair.ant <- beta.pair(occ, index.family="jaccard")  # Pairwise #
# Effect of habitat amount over turnover and nestedness, while controling by patch size 
mantel.partial(pair.ant$beta.jtu, dist(paisagem_com$Forest), dist(paisagem_com$Size)) # Turnover #
mantel.partial(pair.ant$beta.jne, dist(paisagem_com$Forest), dist(paisagem_com$Size)) # Nestedness #

# # # # # # # # #
#### Figures ####
# # # # # # # # #

## Figure scale of effects ####
rsquare<-c(unique(valuesR$adj.r.square),unique(valuesCOM$adj.r.square))
scale<-rep(buffer_distance_meters,2)
cbind(rsquare,scale)->datar
as.data.frame(datar)->datarsquare
datarsquare$cores<-c("black","black","black","red","black","black","black","black","black","black","red","black","black","black","black","black")
datarsquare$plot<-c(rep("(a) Species Richness",8),rep("(b) Species Composition",8))
figscale<-ggplot(data=datarsquare,mapping=aes(x=scale,y=rsquare))+
  geom_line(linetype="dashed")+
  geom_point(size=3,aes(colour = cores))+
  scale_color_manual(values = c("black","#ff3131"))+
  guides(color="none")+
  labs(x="Landscape Size", y="Adj. R²")+
  scale_x_continuous(n.breaks = 8,)+
  scale_y_continuous(n.breaks = 6)+
  facet_wrap(.~plot,scales="free_y",nrow=2,ncol=1)+
  theme_minimal()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "right",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "top"); figscale
# Save
tiff("Scales_of_effect_all.tif",width = 4,height = 5, units = "in",res = 600)
figscale
dev.off()

## Figure Richness x Habitat amount X Patch size ####
# Plot Species Richness
pred_table[pred_table$Scale==2000,]->pred_data_rich
pred_data_rich$Richness<-resp_table$Richness
richplot<-ggplot(data=pred_data_rich,mapping = aes(x = Forest, y = Richness,size=Size))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  geom_point(data= pred_data_rich,aes(shape = Farm, size = Size))+ 
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
  labs(x="", y="Species\nRichness", shape= "Farm", size= "Patch size")+
  theme_classic()+
  scale_y_continuous("Species\nRichness",n.breaks = 8)+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "bottom" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "center",legend.box = "vertical"); richplot

# Plot Composition
pred_table[pred_table$Scale==1500,]->pred_data_comp
pred_data_comp$Composition<-resp_table$PCoA1
pred_data_comp$Composition2<-resp_table$PCoA2
compplot<-ggplot(data=pred_data_comp,mapping = aes(x = Forest, y = pred_data_comp$Composition*-1,size = Size))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  geom_point(data= pred_data_rich,aes(shape = Farm, size = Size))+ 
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
  labs(x="Landscape habitat amount", y="Composition",shape= "Farm", size= "Patch size")+
  theme_classic()+
  scale_y_continuous("Species \nComposition",n.breaks = 8)+
  scale_x_continuous("Landscape Habitat Amount", n.breaks = 6)+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "center",legend.position = "bottom" ,legend.box = "vertical"); compplot
# Plot Final
figure<-ggarrange(richplot,compplot, 
          common.legend = T, legend="bottom",
          labels= c("(a)","(b)"),
          ncol= 1, nrow=2,
          align = "hv", 
          font.label = list(size = 14, color = "black", face = "bold", family = "sans", position = "top")); figure
# Save #
tiff("Habitat amount.tiff",width = 5,height = 7, units = "in",res = 600)
figure
dev.off()

## Figure PCoA ##
pcoaplot<-ggplot()+
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0,linetype="dashed",alpha=0.5)+
  geom_point(data=pred_data_comp, aes(x=Composition,y=Composition2, colour=Forest, size=Size,shape=Farm))+
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3, 9)) +
  labs(x="PCoA 1",y="PCoA 2")+
  labs(size= "Patch size",shape="Farm",color = "Habitat amount")+
  scale_colour_gradient(low="#b0c900", high="#0A7029",)+
  geom_point(data=spp_scores,aes(x=PCoA1,y=PCoA2),shape=21)+
  geom_text_repel(data=spp_scores[spp_scores$guildas==c("Forest","Generalist"),],aes(x=PCoA1,y=PCoA2,label = gsub("_", " ",rownames(spp_scores[spp_scores$guildas==c("Forest","Generalist"),]))),family="sans",fontface="italic",size=2,max.overlaps = 10,force= 10,box.padding = 0.5,point.padding = 0.5 )+
  theme_classic()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold") ,legend.justification = "top" );pcoaplot
# Save
tiff("PCoA.tiff",width = 8,height = 6,units = "in",res = 600)
pcoaplot
dev.off()
# # # # # # # # # # # # # # # # # #