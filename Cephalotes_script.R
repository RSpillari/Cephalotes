# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Script: Landscape habitat amount x Patch size ####
# # # # # # # # # # by: Ricardo Ruaro # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

### Load required packages ###
library(betareg)
library(ggpubr)
library(cowplot) 
library(landscapemetrics)
library(pbapply)
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
library(iNEXT)
library(purrr)
library(stars)
library(units)
library(betaC)
library(ggnewscale)
library(ggspatial)
#### Response Variables ----

## Load Data - Presence table
read.csv("Cephalotes_data.csv",h=T, sep = ",") -> ants
head(ants)
#
length(unique(ants$Species)) # Total number of collected species
strsplit(unique(ants$Species), "_" )->striped
matrix(unlist(striped),ncol=2,byrow=T)->table_strip
length(unique(table_strip[,1])) # Total number of genera collected
#
(table_strip[,2])->species_name
length(species_name[!grepl("^sp\\d+$", species_name)]) # Number of species identified at species level
length(species_name[grepl("^sp\\d+$", species_name)]) # Number of species identified at morphospecies level
#
rowSums(table(table_strip[,1],table_strip[,2]))->nsp_gen # Number of species per genera
write.csv(nsp_gen, "nsp_gen.csv")
#
table(paste(ants$Fragment,ants$Sample,sep = "_"),ants$Species)->abun_point 
occ_point <- decostand(abun_point, "pa") # Occurrence table per sampling point
colSums(occ_point)-> ocpoint # Number of sampling point occurrence for each species
#
## Abundance table - Count of occurrence points for each species in each fragment/plot
frag_names <- sapply(strsplit(rownames(occ_point), "_"), `[`, 1)
occ_fragment_count <- aggregate(occ_point, by = list(Fragment = frag_names), FUN = sum)
rownames(occ_fragment_count)<-occ_fragment_count$Fragment
occ_fragment_count$Fragment<-NULL
abund_plot<-occ_fragment_count 
abund_plot # Abundance per sampling plot
#
decostand(abund_plot, "pa")-> occ_sp #Occurrence table per sampling point
colSums(occ_sp)->ocplot
cbind(ocpoint,ocplot)->records
as.data.frame(records)->records
records[order(records$ocpoint,decreasing = T),]->recordsfinal
write.csv(recordsfinal, "records.csv")

# Singletons, doubletons, ...
table(colSums(occ_sp)) # Species records per sampling plot
table(colSums(abund_plot)) # Species records for sampling points
#
esteio<-(occ_sp[1:3,]) 
dimona<-(occ_sp[4:6,])
poa<-(occ_sp[7:9,])
#
in_esteio <- colSums(esteio) > 0
in_dimona <- colSums(dimona) > 0
in_poa <- colSums(poa) > 0
#
sum(in_esteio) # Number of species registered in the Esteio Farm
sum(in_dimona) # Number of species registered in the Dimona Farm
sum(in_poa)    # Number of species registered in the Poa Farm
#
sum(in_esteio & !in_dimona & !in_poa) # Number of species unique to the Esteio Farm
sum(!in_esteio & in_dimona & !in_poa) # Number of species unique to Dimona Farm
sum(!in_esteio & !in_dimona & in_poa) # Number of species unique to Poa Farm

## Total number of species per fragment ##
rowSums(occ_sp)->sp_obs # Richness observed 

## ENS_fun - Individual based Rarefaction ##
ENS_fun <- function(Sn, n) {
  equation <- function(S, Sn, n) {
    abs(S * (1 - (1 - (1 / S)) ^ n) - Sn)
  }
  out <- optimize(
    equation,
    interval = c(Sn, 10000000),
    Sn = Sn,
    n = n
  )
  out = out$minimum
  if (n == 1)
    out = NA
  return(out)
}
#
full_data<-rarefy(abund_plot,sample= min(rowSums(abund_plot)))
map_dbl(full_data, ENS_fun, n=min(rowSums(abund_plot)))-> SAD_effect # Richness Rarefaction

## Habitat-use guilds - Forest specialists richness ##
subset(ants,ants$Habitat_guild=="Forest")->fores
table(paste(fores$Fragment,fores$Sample,sep="_"),fores$Species)->abund_fores_point
occ_fores_point<- decostand(abund_fores_point, "pa") # Occurence per sampling point
frag_names <- sapply(strsplit(rownames(occ_fores_point), "_"), `[`, 1)
occ_fores_count <- aggregate(occ_fores_point, by = list(Fragment = frag_names), FUN = sum) # Abundance table per sampling plot
rownames(occ_fores_count)<-occ_fores_count$Fragment
occ_fores_count$Fragment<-NULL
forest_abun<-occ_fores_count # Forest specialist abundance table
length(colnames(forest_abun)) # Number of species identified as forest specialists

## Total number of forest species observed
occ_sp_fores<-decostand(occ_fores_count,"pa")
rowSums(occ_sp_fores)->sp_fores_total

## Individual based  Rarefaction 
full_data2<-vegan::rarefy(forest_abun,sample= min(rowSums(forest_abun)))
map_dbl(full_data2, ENS_fun, n=min(rowSums(forest_abun)))-> SAD_effect_forest 

## Habitat-use guilds - Habitat generalists richness ##
subset(ants, ants$Habitat_guild=="Generalist")->gen
table(paste(gen$Fragment,gen$Sample,sep="_"),gen$Species)->abund_gen_point 
occ_gen_point<- decostand(abund_gen_point, "pa") # Occurencce table per sampling point
frag_names <- sapply(strsplit(rownames(occ_gen_point), "_"), `[`, 1)
occ_gen_count <- aggregate(occ_gen_point, by = list(Fragment = frag_names), FUN = sum) # Abundance table per sampling plot
rownames(occ_gen_count)<-occ_gen_count$Fragment
occ_gen_count$Fragment<-NULL
gen_abun<-occ_gen_count # Forest specialist abundance
length(colnames(gen_abun)) # Number of species identified as generalists

## Total number of habitat generalist species observed
occ_sp_gen<-decostand(gen_abun,"pa")
rowSums(occ_sp_gen)->sp_gen_total

## Habitat-use guilds - Open habitat specialist richness ##
subset(ants, ants$Habitat_guild=="Open habitat")->op
table(paste(op$Fragment,op$Sample,sep="_"),op$Species)->abund_op_point
occ_op_point<- decostand(abund_op_point, "pa") # Occurencce per sampling point
frag_names <- sapply(strsplit(rownames(occ_op_point), "_"), `[`, 1)
occ_op_count <- aggregate(occ_op_point, by = list(Fragment = frag_names), FUN = sum) # Abundance per sampling
rownames(occ_op_count)<-unique(frag_names)
occ_op_count$Fragment<-NULL
op_abun<-occ_op_count # Open habitat specialist abundance 
length(colnames(op_abun)) # Species identified as open-habitat specialists

## Total number of open habitat specialist species observed
occ_sp_op<-decostand(op_abun,"pa")
rowSums(occ_sp_op)->sp_op_total
# Fragment 3209 with no open habitat specialist species

## Beta diversity between sample points of each fragment - Heterogeneity ##
#  Cover based Rarefaction 
unique(substr(rownames(occ_point),1,4))->fragz
data.frame(frag=fragz,C=NA)->C_data
for(x in fragz){
  C_data[C_data$frag==x,2] = min(beta_stand(occ_point[substr(rownames(occ_point),1,4)==x,],setsize = 30,list("C_target"),summarise = F, resamples = 9999))  
}
round(min(C_data$C),2)->tagC
allBetaC<-lapply(fragz,function(x){
  betaC_frag<-beta_stand(occ_point[substr(rownames(occ_point),1,4)==x,],setsize = 30,list("beta_C"),args = list(C = tagC),summarise = F, resamples = 9999)
  mean(betaC_frag)->betaC
  return(data.frame(betaC,frag=x))
})
do.call("rbind",allBetaC)->allBetaC ## Beta diversity within each fragment while considering the sample effort

## Fragment ID ##
ID<-rownames(abund_plot)

## Response Variable Table ##
respvars_table<-cbind(ID,SAD_effect,SAD_effect_forest,allBetaC$betaC)
resp_table<-as.data.frame(respvars_table,row.names = F)
colnames(resp_table)<-c("Fragment","Richness","Richness_forest","Within_fragment_β_diversity")
head(resp_table) # Response variables table 
#write.csv(resp_table,"Response_table.csv",row.names = F) # Save #
as.numeric(resp_table$Richness)->resp_table$Richness
as.numeric(resp_table$Richness_forest)->resp_table$Richness_forest
as.numeric(resp_table$Within_fragment_β_diversity)->resp_table$Within_fragment_β_diversity
rownames(resp_table)<-resp_table$Fragment
resp_table<-resp_table[,2:4]

#### Predictor Variables ----
read.csv("pdbff_coord.csv",h=T,sep=",")->pdbff
st_as_sf(pdbff,coords=c("long","lat"),crs=4326)->pdbff
#
read_sf("arie_pdbff_kml.kml")->frags
st_cast(frags, "POLYGON")->frags
st_intersects(pdbff,frags)->int
frags[unlist(int),]->frags
#
read_stars("brasil_coverage_1993.tif") -> ant_rast
st_crop(ant_rast,st_bbox(st_buffer((pdbff),5000))) -> ant_rast_crop
st_warp(ant_rast_crop,crs=4326,use_gdal = T)->ant_rast_wgs
st_transform(ant_rast_wgs, crs = st_crs(32621))->ant_rast_utm
#
st_transform(pdbff,32621)->pdbff
st_transform(frags,32621)->frags

## Patch size ##
st_intersection(frags,pdbff)->pdbff_temp
frags$frag<-pdbff_temp$frag
frags[order(frags$frag),]->frags
st_area(frags)->frag_area
frag_area_ha <- set_units(frag_area, "ha")
as.numeric(frag_area_ha)->patch_size
log(patch_size)->logPS

##Landscape habitat amount ##
buffer_distance_meters <- c(250,500,750,1000,1500,2000,2500,3000)
all_esc<-pblapply(buffer_distance_meters,function(y){
  st_buffer(pdbff,dist = y)->ant_buffer
  ant_buffer$buffer_distance_meters<-y
  return(ant_buffer)
})
do.call("rbind",all_esc)->all_esc
#
st_warp(ant_rast_utm,crs = st_crs((ant_rast_utm))) -> ant_rast_utm
SoE_ant<-pblapply(buffer_distance_meters,function(y){
  st_buffer(pdbff,dist = y)->ant_buffer
  data.frame(site=ant_buffer$frag,HA=NA,PA=NA)->ant_ha
  for(x in 1:nrow(ant_ha)){
    st_crop(ant_rast_utm,ant_buffer[x,],mask=T)->buf_temp
    plot(buf_temp)
    lsm_c_ca(buf_temp)->cover
    cover[cover$class==3,]$value/sum(cover$value)->HA
    HA[1]->ant_ha$HA[x]
    cover[cover$class==15,]$value/sum(cover$value)->PA
    PA[1]->ant_ha$PA[x]
  }
  ant_ha$HA[is.na(ant_ha$HA)] <- 0
  ant_ha$PA[is.na(ant_ha$PA)] <- 0
  return(ant_ha[,c("HA", "PA")])
})
do.call("cbind",SoE_ant)->SoE_ant
colnames(SoE_ant)<-paste0(colnames(SoE_ant),sprintf("%04d",c(0250,0250,0500,0500,0750,0750,1000,1000,1500,1500,2000,2000,2500,2500,3000,3000)))
rownames(SoE_ant)<-ID
# Correlation: Forest cover x Pasture cover
cor(SoE_ant$HA0250,SoE_ant$PA0250)
cor(SoE_ant$HA0500,SoE_ant$PA0500)
cor(SoE_ant$HA0750,SoE_ant$PA0750)
cor(SoE_ant$HA1000,SoE_ant$PA1000)
cor(SoE_ant$HA1500,SoE_ant$PA1500)
cor(SoE_ant$HA2000,SoE_ant$PA2000)
cor(SoE_ant$HA2500,SoE_ant$PA2500)
cor(SoE_ant$HA3000,SoE_ant$PA3000)
# Landscape habitat amount for each sample plot in each buffer size
SoE_ant<-SoE_ant[,c(1,3,5,7,9,11,13,15)]

## Farm ##
farm<-c(rep("Esteio",3),rep("Dimona",3),rep("Porto Alegre",3)) # Farm where each fragment is located #
as.factor(farm)->farm

## Predictor Variables Table ##
pred_table<-cbind(ID,farm,patch_size,logPS,SoE_ant)
colnames(pred_table)<-c("Fragment","Farm","Size","Size_log",colnames(SoE_ant))
pred_table # Predictor variables table # For figure construction

## Corelation between patch size and forest cover ##
cor(logPS,SoE_ant$HA0250)
cor(logPS,SoE_ant$HA0500)
cor(logPS,SoE_ant$HA0750)
cor(logPS,SoE_ant$HA1000)
cor(logPS,SoE_ant$HA1500)
cor(logPS,SoE_ant$HA2000)
cor(logPS,SoE_ant$HA2500)
cor(logPS,SoE_ant$HA3000)

# # # # # # # # # # # #
#### Data Analysis ####
# # # # # # # # # # # #
SoE_ant2<-SoE_ant[,4:8]

isolation<-c(150,700,0,600,225,150,300,900,450)
as.data.frame(isolation)->is
#
final_resu<-{}
for (y in 1:ncol(resp_table)) {
  lms_ant<-lapply(1:ncol(SoE_ant2),function(x){
  data.frame(S=resp_table[,y],Frag=ID,Farm=farm,PS=logPS,isolation=is,HA=SoE_ant2[,x])->temp_data
  scale(temp_data[,c("HA","PS","isolation")])->temp_data[,c("HA","PS","isolation")]
  glmmTMB(S~HA+PS+(1|Farm),family = "gaussian", data=temp_data)->lm_temp
  data.frame(summary(lm_temp)[[6]]$cond,
             var=rownames(summary(lm_temp)[[6]]$cond),
             escala=colnames(SoE_ant2)[x],
             AIC=AIC(lm_temp),
             R2m=r.squaredGLMM(lm_temp)[1],
             R2c=r.squaredGLMM(lm_temp)[2],
             COR=cor(temp_data[,c("HA","PS")])[1,2]) -> coef_temp
  colnames(coef_temp)=c("Estimate","Error","t","P","var","escala","AIC","R2m","R2c","COR")
  return(coef_temp)
})
do.call("rbind",lms_ant)->lms_ant
as.numeric(gsub("HA","",lms_ant$escala))->lms_ant$escala
lms_ant[lms_ant$var!="(Intercept)",]->ant_coefs
#ant_coefs[order(abs(ant_coefs$R2m),decreasing = T),]
final_resu[[colnames(resp_table)[y]]] <- ant_coefs  
}; final_resu

final_table <- do.call(rbind, final_resu)
final_table$resp_var <- rep(names(final_resu), sapply(final_resu, nrow))

write.csv(final_table,"model_coefs.csv", row.names = F)

#### Partial regression ####
### All ants - Richness 2500m 
SoE_ant2[,4]->Forest
as.data.frame(cbind(Forest,logPS))->tab
scale(tab)->table_pad
cbind(table_pad,resp_table[,1])->table_pad
colnames(table_pad)<-c("Forest","Size","S")
as.data.frame(table_pad)->table_partial
# Regress Forest without size effect
model_aux <- lm(S ~ Size, data = table_partial)
resid_richness <- residuals(model_aux)
#
model <- lm(resid_richness ~ Forest, data = table_partial)
summary(model)  
#
# Regress Patch without forest effect
model_aux1 <- lm(S ~ Forest, data = table_partial)
resid_richness2 <- residuals(model_aux1)
# 
model2 <- lm(resid_richness2 ~ Size, data = table_partial)
summary(model2) 

# Forest ants - Richness 3000m
SoE_ant2[,5]->Forest
as.data.frame(cbind(Forest,logPS))->tab
scale(tab)->table_pad
cbind(table_pad,resp_table[,2])->table_pad
colnames(table_pad)<-c("Forest","Size","Sfor")
as.data.frame(table_pad)->table_partial
#
#Regress Forest without size effect
model_aux2 <- lm(Sfor ~ Size, data = table_partial)
resid_richness3 <- residuals(model_aux2)
#
model3 <- lm(resid_richness3 ~ Forest, data = table_partial)
summary(model3)
#
# Regress Patch without forest effect
model_aux3 <- lm(Sfor ~ Forest, data = table_partial)
resid_richness4 <- residuals(model_aux3)
#
model3 <- lm(resid_richness4 ~ Size, data = table_partial)
summary(model3)  

####  Beta diversity between fragments ####
## All ants
lms_ant<-lapply(1:ncol(SoE_ant2),function(x){
    data.frame(Frag=ID,Farm=farm,PS=logPS,isolation=is,HA=SoE_ant2[,x])->temp_data
    scale(temp_data[,c("HA","PS")])->temp_data[,c("HA","PS")]
    vegdist(temp_data[,5],"euclidean")->distforest
    vegdist(temp_data[,3],"euclidean")->distpatch
    vegdist(as.numeric(farm),"euclidean")->diffarm
    vegdist(abund_plot,"bray")->dist_bray
    diffs<-data.frame(unlist(melt(dist_bray)),unlist(melt(distforest)),unlist(melt(distpatch)),as.factor(unlist(melt(diffarm))))
    colnames(diffs)<-c("diffbeta","difforest","diffpatch","diffarm")
    glmmTMB(diffbeta~difforest+diffpatch+(1|diffarm), family = beta_family(link = "logit"), data=diffs)->lm_temp
    data.frame(summary(lm_temp)[[6]]$cond,
               var=rownames(summary(lm_temp)[[6]]$cond),
               escala=colnames(SoE_ant2)[x],
               AIC=AIC(lm_temp),
               R2m=r.squaredGLMM(lm_temp)[1],
               R2c=r.squaredGLMM(lm_temp)[2],
               COR=cor(temp_data[,c("HA","PS")])[1,2]) -> coef_temp
    colnames(coef_temp)=c("Estimate","Error","t","P","var","escala","AIC","R2m","R2c","COR")
    return(coef_temp)
    final_resu[[length(lms_ant)+1]] <-coef_temp 
  })
############
#
vegdist(abund_plot,"bray")->dist_bray
#
finalbet<-{}
for (w in 1:ncol(SoE_ant2)) {
  cbind(SoE_ant2[,w],logPS)->dat_temp
  colnames(dat_temp)<-c("HA","PS")
  as.data.frame(dat_temp)->data_temp
  scale(data_temp)->data_tem_pad
  vegdist(data_tem_pad[,1],"euclidean")->distforest
  vegdist(data_tem_pad[,2],"euclidean")->distpatch
  diffs<-data.frame(unlist(melt(dist_bray)),unlist(melt(distforest)),unlist(melt(distpatch)))
  colnames(diffs)<-c("diffbeta","difforest","diffpatch")
  modelbeta<-betareg(diffbeta ~ difforest+diffpatch,data=diffs)
  sum_model <- summary(modelbeta)
  coefs <- sum_model$coefficients$mean  # Coeficientes e p-valores
  pseudo_r2 <- sum_model$pseudo.r.squared  # Pseudo R²
  akai<-AIC(modelbeta)
  scale<-w
  finalbet[[length(finalbet)+1]]<-data.frame(coefs,pseudo_r2,akai,scale)
}
do.call("rbind",finalbet)->finalbeta
rep(c("Intercept","HA","PS"),5)->finalbeta$vars
finalbeta[finalbeta$vars!="Intercept",]->tabeta
rep(c(1000,1500,2000,2500,3000),each=2)->tabeta$scale
rep("Between-fragments β-diversity",5)->tabeta$resp_var
#
write.csv(tabeta,"model_coefs_beta.csv",row.names = F)
#### Nestedness x Turnover ####
## All Ants
occ_all <- decostand(abund_plot, "pa")
logPS
SoE_ant$HA3000->ha_data

## NODF ##
ant.for1 <- occ_all[order(ha_data, decreasing = TRUE), ] # Rows ordered based on the landscape Forest amount #
ant.pat1<- occ_all[order(logPS, decreasing = TRUE), ] # Rows ordered based on the patch Size #
null_for1 <- oecosimu(as.matrix(ant.for1), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_for1
null_pat1<- oecosimu(as.matrix(ant.pat1), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_pat1

## Baselga ##
overal_all  <- beta.multi(occ_all, index.family="jaccard") # Multisite #
pair.ant_all <- beta.pair(occ_all, index.family="jaccard")  # Pairwise #
# Effect of habitat amount over turnover and nestedness, while controling by patch size 
mantel.partial(pair.ant_all$beta.jtu, dist(ha_data), dist(logPS)) # Turnover #
mantel.partial(pair.ant_all$beta.jne, dist(ha_data), dist(logPS)) # Nestedness #

## Forest Specialists 
occ_fores <- decostand(forest_abun, "pa")

## NODF ##
ant.for2 <- occ_fores[order(ha_data, decreasing = TRUE), ] # Rows ordered based on the landscape Forest amount #
ant.pat2 <- occ_fores[order(logPS, decreasing = TRUE), ] # Rows ordered based on the patch Size #
null_for2 <- oecosimu(as.matrix(ant.for2), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_for2
null_pat2 <- oecosimu(as.matrix(ant.pat2), nestednodf, method = "swap", thin = 10, burnin = 1000, nsimul = 9999, order = FALSE); null_pat2

## Baselga ##
overal_fores   <- beta.multi(occ_fores, index.family="jaccard") # Multisite #
pair.ant_fores <- beta.pair(occ_fores, index.family="jaccard")  # Pairwise #
# Effect of habitat amount over turnover and nestedness, while controling by patch size 
mantel.partial(pair.ant_fores$beta.jtu, dist(ha_data), dist(logPS)) # Turnover #
mantel.partial(pair.ant_fores$beta.jne, dist(ha_data), dist(logPS)) # Nestedness #

# # # # # # # # #
#### Figures ####
# # # # # # # # #

## Figure scale of effects ####
rsquare<- c(unique(final_table$R2m),unique(tabeta$pseudo_r2))
aic<-c(unique(final_table$AIC),unique(tabeta$akai))
estimateHA<-c(final_table$Estimate[final_table$var=="HA"],tabeta$Estimate[tabeta$vars=="HA"])
estimatePS<-c(final_table$Estimate[final_table$var=="PS"],tabeta$Estimate[tabeta$vars=="PS"])
pvalHA<-c(final_table$P[final_table$var=="HA"],tabeta$Pr...z..[tabeta$vars=="HA"])
pvalPS<-c(final_table$P[final_table$var=="PS"],tabeta$Pr...z..[tabeta$vars=="PS"])
scale<-rep(unique(final_table$escala),4)
cbind(rsquare,scale,aic,estimateHA,estimatePS,pvalHA,pvalPS)->datar
as.data.frame(datar)->datarsquare
datarsquare$plot<- rep(c("Ant richness","Forest specialist richness","Within-fragment β-diversity", "Between-fragments β-diversity"), each = 5)
figscale<-ggplot(data=datarsquare,mapping=aes(x=scale,y=rsquare ))+
  geom_line(linewidth=0.8,aes(linetype = plot, colour = plot))+
  geom_point(size=3,aes(pch=plot,color=plot))+
  scale_y_continuous(n.breaks = 6)+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  scale_shape_manual(values = c(16,17,15,18))+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f","#f8961e"))+
  labs(x = "Landscape size",y = "R²", color = "Response Variable",shape = "Response Variable",linetype = "Response Variable") +
  theme_minimal()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "none",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.justification = "center"); figscale
#
tiff("Scales_r².tif",width = 10,height = 6, units = "in",res = 600)
figscale
dev.off()
#
figaic<-ggplot(data=datarsquare,mapping=aes(x=scale,y=aic))+
  geom_line(linewidth=0.8,aes(linetype = plot, colour = plot))+
  geom_point(size=3,aes(pch=plot,color=plot))+
  scale_y_continuous(n.breaks = 8)+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  scale_shape_manual(values = c(16,17,15,18))+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f","#f8961e"))+
  labs(x = "Landscape size",y = "AIC", color = "Response Variable",shape = "Response Variable",linetype = "Response Variable") +
  theme_minimal()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "none",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "center"); figaic

figestimateHA<-ggplot(data=datarsquare,mapping=aes(x=scale,y=estimateHA))+
  geom_line(linewidth=0.8,aes(linetype = plot, colour = plot))+
  geom_point(size=3,aes(pch=plot,color=plot))+
  scale_y_continuous(n.breaks = 6)+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  scale_shape_manual(values = c(16,17,15,18))+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f","#f8961e"))+
  labs(x = "Landscape size",y = "Effect size HA", color = "Response Variable",shape = "Response Variable",linetype = "Response Variable") +
  theme_minimal()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "none",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "center"); figestimateHA

figpHA<-ggplot(data=datarsquare,mapping=aes(x=scale,y=pvalHA))+
  geom_line(linewidth=0.8,aes(linetype = plot, colour = plot))+
  geom_point(size=3,aes(pch=plot,color=plot))+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  scale_shape_manual(values = c(16,17,15,18))+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f","#f8961e"))+
  labs(x = "Landscape size",y = "p-value HA", color = "Response Variable",shape = "Response Variable",linetype = "Response Variable") +
  scale_y_continuous(n.breaks = 6)+
  theme_minimal()+
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5)+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "right",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "center"); figpHA

figpPS<-ggplot(data=datarsquare,mapping=aes(x=scale,y=pvalPS))+
  geom_line(linewidth=0.8,aes(linetype = plot, colour = plot))+
  geom_point(size=3,aes(pch=plot,color=plot))+
  scale_linetype_manual(values=c("solid","dashed","dotted","dotdash"))+
  scale_shape_manual(values = c(16,17,15,18))+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f","#f8961e"))+
  labs(x = "Landscape size",y = "p-value FS", color = "Response Variable",shape = "Response Variable",linetype = "Response Variable") +
  scale_y_continuous(n.breaks = 6)+
  theme_minimal()+
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5)+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "right",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "center"); figpPS

get_legend(figpPS)->legendar

plots<-ggarrange(figscale,figaic,figestimateHA,figpHA,figpPS,
                  common.legend = F,legend = "none",labels = c("a)","b)","c)","d)","e)"),
                    ncol= 3, nrow=2,
                    align = "hv",
                    font.label = list(size = 14, color = "black", face = "bold", family = "sans",
                                      position = "top"))
plots+annotation_custom(legendar,xmin = 0.65,ymax=0.55)->figure_coefs_models

tiff("coefs_final.tif",width = 10,height = 6, units = "in",res = 600)
figure_coefs_models
dev.off()

## Figure Richness x Habitat amount X Patch size ####
data_fig<-cbind(resp_table,pred_table)
#
# Plot Species Richness x Habitat amount
richplot<-ggplot(data=data_fig,mapping = aes(x = HA2500*100, y = Richness))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  geom_point(data= data_fig,aes(color = Farm,size = Size))+ 
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
  labs(x="", y="", color= "Farm", size= "Fragment size (ha)")+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic()+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f"))+
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous("Species\nrichness",n.breaks = 8)+
  theme(strip.text=element_text(family="sans",size=14,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "right" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "top",legend.box = "vertical"); richplot
#
# Plot Forest Specilist Richnness x Habitat amount
rich_forest<-ggplot(data=data_fig,mapping = aes(x = HA3000*100, y = Richness_forest))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  geom_point(data= data_fig,aes(color = Farm,size=Size))+ 
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
  labs(x="Habitat amount (%)", y="", color= "Farm", size= "Fragment size (ha)")+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f"))+
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic()+
  scale_y_continuous("Forest specialist\nrichness",n.breaks = 8)+
  theme(strip.text=element_text(family="sans",size=14,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "right" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "top",legend.box = "vertical"); rich_forest
# 
# ## Plot Species Richness x Patch size
#  S_patch<-ggplot(data=data_fig,mapping = aes(x = logPS, y = Richness))+
#   #geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
#    geom_point(data= data_fig,aes(shape = Farm,size=Size))+
#    scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
#    labs(x=" ", y=" ", shape= "Farm",)+
#    theme_classic()+
#    scale_y_continuous("",n.breaks = 8)+
#    theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "none" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "top",legend.box = "vertical"); S_patch
# #
# ## Plot Forest Specilist Richnness X Patch size
# S_patch_forest<-ggplot(data=data_fig,mapping = aes(x = logPS, y = Richness_forest))+
#    geom_point(data= data_fig,aes(shape = Farm,,size=Size))+ 
#   scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
#    guides(size="none")+
#    labs(x="Patch size", y=" ", shape= "Farm",size="")+
#    theme_classic()+
#    scale_y_continuous(" ",n.breaks = 8)+
#    theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "none" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "top",legend.box = "vertical"); S_patch_forest
# #
#
painel_final<-ggarrange(richplot,rich_forest,
                  common.legend = T,legend = "bottom",
                  ncol= 1, nrow=2,
                  labels = c("a)","b)"),
                  align = "hv", 
                  font.label = list(size = 14, color = "black", face = "bold", family = "sans"));
#                
tiff("richness.tiff",width = 5,height = 8.5 , units = "in",res = 600)
painel_final
dev.off()

#### Figure Heterogeneity -Beta diversity within each fragment ####
hetero<-ggplot(data=data_fig,mapping = aes(x = HA3000*100, y = Within_fragment_β_diversity,size=Size))+
  #geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  geom_point(data= data_fig,aes(color = Farm,size = Size))+ 
 scale_size_continuous(breaks = c(1, 10, 100), range = c(3,6.5)) +
  labs(x="Habitat amount (%)", y="",color= "Farm", size= "Fragment size (ha)")+
  theme_classic()+
  scale_y_continuous("Within-fragment β-diversity",n.breaks = 8)+
  scale_colour_manual(values=c("#1f8d49","#90be6d","#f9c74f"))+
  scale_x_continuous("Habitat amount (%)", n.breaks = 6)+
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "top",legend.position = "bottom" ,legend.box = "vertical"); hetero
tiff("whithinfragment.tif",width = 4.5,height =5 ,units = "in",res = 600)
hetero
dev.off()

#### Figure Beta diversity ####
cbind(SoE_ant2[,5],logPS)->dat
colnames(dat)<-c("HA","PS")
as.data.frame(dat)->data_
scale(data_)->data_pad
vegdist(data_tem_pad[,1],"euclidean")->distforest
vegdist(data_tem_pad[,2],"euclidean")->distpatch
diffs<-data.frame(unlist(melt(dist_bray)),unlist(melt(distforest)),unlist(melt(distpatch)))
colnames(diffs)<-c("diffbeta","distforest","diffpatch")

forestdif<-ggplot(data = diffs,mapping = aes(x = distforest, y =diffbeta))+
  geom_point(size=2.5)+
  geom_smooth(method=glm, color='black',size=1,linewidth=0.6, method.args=list(family="Gamma"),se=F)+
  labs(x="Diferences in habitat amount", y="")+
  theme_classic()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",))

patchdif<-ggplot(data = diffs,mapping = aes(x = diffpatch, y =diffbeta))+
  geom_jitter(size = 2.5, width = 0.1, height = 0) +
  labs(x="Diferences in fragment size", y="Between-fragments β-diversity")+
  theme_classic()+
  scale_y_continuous()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",))

beta_figure<-ggarrange(patchdif,forestdif,
          labels= c("a)","b)"),
          ncol= 2, nrow=1,
          align = "hv", 
          font.label = list(size = 14, color = "black", face = "bold", family = "sans", position = "top"))

tiff("beta_fig.tif",width = 8,height =4,units = "in",res = 600)
beta_figure
dev.off()

#### Figure PCoA ####
## Composition - PCoA ##
vegdist(abund_plot,"bray")->dist_bray 
betadisper((dist_bray),ID)->pcoa
pcoa$vectors[,1]->eixo1
pcoa$vectors[,2]->eixo2 
# Species Scores PCoA 
wascores(eixo1, abund_plot)->species_scores
wascores(eixo2, abund_plot)->species_scores2
data.frame(PCoA1=species_scores,PCoA2=species_scores2)->spp_scores
ants[match(rownames(spp_scores),ants$Species),]$Habitat_guild->guildas
guildas[is.na(guildas)]<-"Morphospecies"
data.frame(spp_scores,guildas=guildas)->spp_scores # Position of species in the PCoA plot
#
pcoaplot<-ggplot()+
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0,linetype="dashed",alpha=0.5)+
  geom_point(aes(x=eixo1,y=eixo2, colour=data_fig$HA3000*100, size=data_fig$Size,shape=data_fig$Farm))+
  scale_size_continuous(breaks = c(1, 10, 100), range = c(3, 9)) +
  labs(x="PCoA 1 - 27.25%",y="PCoA 2 - 14.31%")+
  labs(size= "Fragment size (ha)",shape="Farm",color = "Habitat amount (%)")+
  scale_colour_gradient(low="#f9c74f", high="#1f8d49",)+
  geom_point(data=spp_scores,aes(x=PCoA1,y=PCoA2),shape=21)+
  geom_text_repel(data=spp_scores[spp_scores$guildas==c("Forest","Generalist"),],aes(x=PCoA1,y=PCoA2,label = gsub("_", " ",rownames(spp_scores[spp_scores$guildas==c("Forest","Generalist"),]))),family="sans",fontface="italic",size=2.5,max.overlaps =10,force= 10,box.padding = 0.5,point.padding = 0.5 )+
  theme_classic()+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold") ,legend.justification = "top" );pcoaplot
# Save
tiff("PCoA.tiff",width = 8,height = 6,units = "in",res = 600)
pcoaplot
dev.off()
#
#### MAPAS ####
read.csv("Codigos-da-legenda-colecao-9.csv",h=T,sep="\t")->cols
#
st_read("Rodovias_Amazonas.shp")->road
st_transform(road,32621)->road
st_crop(road,ant_rast_utm)->roads
#
pdbff
frags
names(ant_rast_utm)[1]<-"brasil_coverage_1993.tif"
#
## BDFFP fragments 
mapbdffp<-ggplot() +
  geom_stars(data = ant_rast_utm,aes(fill=factor(brasil_coverage_1993.tif)), downsample=0) +
  geom_sf(data = frags, fill = "black", color = "black", alpha= 0.5) +
  #geom_sf(data = pdbff, fill = "transparent", color = "red", linetype='dashed') +
  scale_fill_manual(values = cols$Color,labels = cols$Description, breaks = cols$Class_ID, name="Land cover", na.value = "transparent") +
  labs(x = "", y = "", color = "Land Cover")+
  geom_sf(data=roads, color="black",linetype = "solid",linewidth=0.3,alpha=0.75)+
  theme_minimal()+
  theme(strip.text = element_text(family="sans", size = 12, face="bold", hjust = 0), strip.background = element_rect(fill="white", colour = "white"), axis.text.x = element_blank(),axis.text.y= element_text(family="sans", size=12, color = "black"), axis.title = element_text(family = "sans", size = 14, face="bold"), legend.position = "right", legend.text = element_text(family = "sans", size = 12, hjust = 0), legend.title = element_text(family="sans", size = 12, face = "bold", hjust = 0), legend.justification = "top")+
          annotation_scale(location = "bl", width_hint = 0.065,style="ticks", text_face="italic", text_family = "sans", text_cex = 0.75,pad_x = unit(0.85, "cm"), pad_y = unit(0.65, "cm"))+
          annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_nautical(line_width = 1,line_col = "black",fill = c("black","transparent"), text_size = 5, text_face = "italic", text_family = "sans",text_col = "black"), height = unit(1, "cm"), width = unit(1, "cm"), pad_x = unit(0.7, "cm"), pad_y = unit(1, "cm"))+
  annotate("text", x = Inf, y = Inf, label = "Dimona Farm", hjust =5.25, vjust = 11, size = 4, family = "sans",color = "black",fontface = "italic")+
  annotate("text", x = Inf, y = Inf, label = "Porto Alegre Farm", hjust = 2.5, vjust = 8, size = 4, family = "sans",color = "black",fontface = "italic")+
  annotate("text", x = Inf, y = Inf, label = "Esteio Farm", hjust = 2, vjust = 12, size = 4, family = "sans",color = "black",fontface = "italic")

tiff("test.tiff",width = 10,height = 4,units = "in",res = 600)
mapbdffp
dev.off()
#
## Landscape scales
all_buf<-lapply(1:nrow(pdbff),function(x){
  all_buf_temp<-lapply(c(1000,1500,2000,2500,3000),function(y){
    st_buffer(pdbff[x,], y)->a
    return(a)
  })
  do.call("rbind",all_buf_temp)->all_buf_temp
  return(all_buf_temp)
})
do.call("rbind",all_buf)->all_buf
#
pdballbuff<-ggplot() +
  geom_stars(data = ant_rast_utm,aes(fill=factor(brasil_coverage_1993.tif)), downsample=0) +
  geom_sf(data = frags, fill="black", color = "black", alpha=0.5) +
  geom_sf(data=roads, color="black",linetype = "solid",linewidth=0.3,alpha=0.75)+
  scale_fill_manual(values = cols$Color, 
                    labels = cols$Description, 
                    breaks = cols$Class_ID,
                    name="Land Cover",
                    na.value = "transparent") +
  theme_minimal()+
  scale_y_continuous(scales::pretty_breaks(n = 6))+
  theme(axis.title = element_blank())+
  labs(x = "", y = "", color = "Land Cover")+
  geom_sf(data=all_buf, fill = "transparent", color = "black", linetype = "dashed",linewidth=0.3)+
  geom_sf(data = pdbff, color = "red", size=2) +
  theme(plot.margin = margin(0,0,0,0, unit = "pt"),strip.text = element_text(family="sans", size = 12, face="bold", hjust = 0), strip.background = element_rect(fill="white", colour = "white"), axis.text.x = element_blank(), axis.text.y = element_text(family = "sans", size = 12, face="bold",color = "black"), legend.position = "right", legend.text = element_text(family = "sans", size = 12, hjust = 0), legend.title = element_text(family="sans", size = 12, face = "bold", hjust = 0), legend.justification = "top")+
  annotation_scale(location = "bl", width_hint = 0.065,style="ticks", text_face="italic", text_family = "sans", text_cex = 0.75,pad_x = unit(0.85, "cm"), pad_y = unit(0.65, "cm"))+
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_nautical(line_width = 1,line_col = "black",fill = c("black","transparent"), text_size = 5, text_face = "italic", text_family = "sans",text_col = "black"), height = unit(1, "cm"), width = unit(1, "cm"), pad_x = unit(0.7, "cm"), pad_y = unit(1, "cm")); pdballbuff

## Only landscape size 2500 m
all_buf<-lapply(1:nrow(pdbff),function(x){
  all_buf_temp<-lapply(c(2500),function(y){
    st_buffer(pdbff[x,], y)->a
    return(a)
  })
  do.call("rbind",all_buf_temp)->all_buf_temp
  return(all_buf_temp)
})
do.call("rbind",all_buf)->all_buf
#
st_crop(ant_rast_utm,all_buf[[length(all_buf)]]) -> ant_landscape
names(ant_landscape)[1]<-"brasil_coverage_1993.tif"

land_use_map<-ggplot() +
  geom_stars(data = ant_landscape,aes(fill=factor(brasil_coverage_1993.tif))) +
  geom_sf(data = all_buf, fill = "transparent", color = "black",linetype="dashed", linewidth= 0.5) +
  geom_sf(data = pdbff, color = "red", size=2) +
  scale_fill_manual(values = cols$Color, 
                    labels = cols$Description, 
                    breaks = cols$Class_ID,
                    name="Cobertura do solo",
                    na.value = "transparent") +
  labs(x = "", y = "", color = "Land Cover")+
  geom_sf(data=roads, color="black",linetype = "solid",linewidth=0.3,alpha=1)+
  theme_minimal()+
  theme(strip.text = element_text(family="sans", size = 12, face="bold", hjust = 0), strip.background = element_rect(fill="white", colour = "white"), axis.text = element_text(family="sans", size=12, color = "black"), axis.title = element_text(family = "sans", size = 14, face="bold"), legend.position = "right", legend.text = element_text(family = "sans", size = 12, hjust = 0), legend.title = element_text(family="sans", size = 12, face = "bold", hjust = 0), legend.justification = "top")+
  annotation_scale(location = "bl", width_hint = 0.065,style="ticks", text_face="italic", text_family = "sans", text_cex = 0.75,pad_x = unit(0.85, "cm"), pad_y = unit(0.65, "cm"))+
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_nautical(line_width = 1,line_col = "black",fill = c("black","transparent"), text_size = 5, text_face = "italic", text_family = "sans",text_col = "black"), height = unit(1, "cm"), width = unit(1, "cm"), pad_x = unit(0.7, "cm"), pad_y = unit(1, "cm"))+
  annotation_scale(location = "bl", width_hint = 0.065,style="ticks", text_face="italic", text_family = "sans", text_cex = 0.75,pad_x = unit(0.85, "cm"), pad_y = unit(0.65, "cm"))+
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_nautical(line_width = 1,line_col = "black",fill = c("black","transparent"), text_size = 5, text_face = "italic", text_family = "sans",text_col = "black"), height = unit(1, "cm"), width = unit(1, "cm"), pad_x = unit(0.7, "cm"), pad_y = unit(1, "cm")) ; land_use_map

mapa_fin<-ggarrange(mapbdffp,pdballbuff,land_use_map, 
                    common.legend = T,legend = "bottom",
                    ncol= 1, nrow=3,
                    align = "hv", labels = c("a)","b)","c)"),
                    label.x = 0.89,
                    label.y =0.94,
                    font.label = list(size = 14, color = "black", face = "bold", family = "sans",
                                      position = "top"
                    ));mapa_fin

tiff("mapafinal.tiff",width = 7,height = 11,units = "in",res = 600)
mapa_fin
dev.off()


# # # # # # # # #
#### Generic ####
# # # # # # # # #
generico<-function(tabela,gradiente,at,grad,eixoY,eixoX,excluirspp){ 
  tabela<-as.matrix(tabela)
  gradiente<-as.matrix(gradiente)
  media.pond<-colSums(tabela*gradiente[,1])/colSums(tabela) 
  sub.orden<-tabela[order(gradiente[,1],decreasing=F),]	
  sub.orde<-sub.orden[,order(media.pond,decreasing=T)]
  dados.pa<-matrix(0,nrow(tabela),ncol(tabela))
  dados.pa[tabela>0]<-1
  ordenado<-sub.orde[,which(colSums(dados.pa)>excluirspp)] 
  par(mfrow=c(ncol(ordenado)+1,1),mar=c(0,4,0.2,10),oma=c(3,1,1,6))
  layout(matrix(1:(ncol(ordenado)+1)),heights=c(3,rep(1,ncol(ordenado))))
  plot(sort(gradiente[,1]),axes=F,ylab="",mfg=c(21,1),lwd=10,las=2,lend="butt",frame.plot=F,xaxt="n",type="h",col="black",ylim=c(min(gradiente),max(gradiente)))
  axis(side=2,at=c(0,max(gradiente)),las=2)
  mtext(grad,5,outer=T,font=2,line=-10,padj=-18.5,las=2)
  for(i in 1:ncol(ordenado)){
    if(i == ncol(ordenado)){
      barplot(ordenado[,i],bty="l",axisnames=T,col.axis="white",axis.lty = 1,names.arg = 1:nrow(ordenado),cex.names=000.1,axes=FALSE,col="black")
    }
    else{barplot(ordenado[,i],bty="l",axisnames=F,axes=FALSE,col="black")}
    mtext(colnames(ordenado)[i],3,line=-1.0,adj=0,at=at,cex=.8,font=3)
  }
  mtext(eixoX,1,outer=T,font=2,line=1.2)
  mtext(eixoY,2,font=2,outer=T,line=-2)
}
#
# Plots 
generico(abund_plot,patch_size,11,"Size","Species","Sampling plots oreder by size",5)
generico(abund_plot,pred_table$HA2500,11,"Habitat Amount","Species","Sampling plots oreder by Habitat Amount",5)
#
generico(forest_abun,patch_size,11,"Size","Species","Sampling plots oreder by size",1)
generico(forest_abun,pred_table$HA3000,11,"Habitat Amount","Species","Sampling plots oreder by Habitat Amount",1)
# # # # # # # # # # # # # # # # # # #