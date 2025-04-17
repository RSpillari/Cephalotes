# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Script: Landscape habitat amount x Patch size ####
# # # # # # # # # by Ricardo Spillari # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

### Library required packages ###
library(sp)
library(sf)
library(raster)
library(vegan)
library(glmmTMB)
library(betapart)
library(bipartite)
library(ggplot2)
library(ggpubr)
library(TITAN2)
library(ggrepel)
library(reshape2)
library(MuMIn)
library(indicspecies)
# # # # # # # # # # # # # # #
#### Response Variables #####
# # # # # # # # # # # # # # #

# Ant data Vasconcelos & Delabie 2000  Updated by Ricardo Spillari #

## Richness ##
ants<-read.csv("Cephalotes_data.csv",h=T, sep = ",")
head(ants)
# Species Abundance/Frequency Table
abund_table<-table(ants$Fragment,ants$Species)
write.csv(abund_table,"Abundance_table.csv",row.names = F) # Save abundance table #
# Species Occurrence Table
occ <- decostand(abund_table, "pa")
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
rowSums(occ_op)->Richness_opens # Species richness of open habitat speciealist #
Richness_open<-c(1,3,3,2,1,1,1,0,1)

## Fragment ID ##
ID<-rownames(abund_table)

## Composition - PCoA ##
vegdist(occ,"jaccard")->dist_jac  # Using occurrence data 
betadisper((dist_jac),ID)->pcoa
pcoa$vectors[,1]->eixo1
pcoa$vectors[,2]->eixo2

## Species in the axis
wascores(eixo1, occ)->species_scores
wascores(eixo2, occ)->species_scores2
data.frame(PCoA1=species_scores,PCoA2=species_scores2)->spp_scores

# abund<-table(ants$Fragment,ants$Species,ants$Method)
# melt(abund)->abund_new
# abund_new[abund_new$value!=0,]->abund_new
# abund_new$value<-rep(1,length(abund_new$value))
# table(paste0(abund_new$Var1,"_",abund_new$Var3),abund_new$Var2)->final_mat
# substr(rownames(final_mat),6,6)->meth
# multipatt(final_mat,meth)->scor
# summary(scor)
ants[match(rownames(spp_scores),ants$Species),]$Habitat_guild->guildas
guildas[is.na(guildas)]<-"Morphospecies"

data.frame(spp_scores,guildas=guildas)->spp_scores


## Response Variable Table ##
respvars_table<-cbind(ID,Richness,Richness_forest,Richness_generalist,Richness_open,eixo1,eixo2)
resp_table<-as.data.frame(respvars_table,row.names = F)
colnames(resp_table)<-c("Fragment","Richness","Richness_forest","Richness_generalist","Richness_open","PCoA1","PCoA2")
head(resp_table) # Response variables table #
write.csv(resp_table,"resp.csv",row.names = F) # Save #
as.numeric(resp_table$Richness)->resp_table$Richness
as.numeric(resp_table$Richness_forest)->resp_table$Richness_forest
as.numeric(resp_table$Richness_generalist)->resp_table$Richness_generalist
as.numeric(resp_table$Richness_open)->resp_table$Richness_open
as.numeric(resp_table$PCoA1)->resp_table$PCoA1
as.numeric(resp_table$PCoA2)->resp_table$PCoA2

# # # # # # # # # # # # # # #
#### Predictor variables ####
# # # # # # # # # # # # # # #

## Landscape habitat amount ##
# Sample points coordinates 
coordenadas<-read.csv("points.csv",h=T, sep = ",")
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
buffer_distance_meters <- c(250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000)
# Landscape covers for each buffer size
porcentage<-{}
for (x in 1:9){
  for (y in 1:16){
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
patch<-rep(c(rep(1,16),rep(10,16),rep(100,16)),3)# Size of each fragment #

## Farm ##
farm<-c(rep("Esteio",48),rep("Dimona",48),rep("Porto Alegre",48)) # Farm where each fragment is located #
as.data.frame(farm)->farm1
as.factor(farm)->farm2

## Predictor Variables Table ##
pred_table<-cbind(habitat$frag,farm1,habitat$scale,habitat$forest,patch)
colnames(pred_table)<-c("Fragment","Farm","Scale","Forest","Size")
head(pred_table) # Predictor variables table #
write.csv(pred_table,"preds.csv",row.names = F) # Save #
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
# Part 1. Richness ####
# # # # # # # # # # # # 
# Effects of landscape habitat amount and patch size over total species richness and forest specialist richness #

#### Scale of effect ####

## Species richness ##
values<-{}
for (h in buffer_distance_meters ) {
  paisagem<-pred_table_pad[pred_table_pad$Scale==h,]
  paisagem$Richness<-resp_table$Richness
  glmmTMB(Richness~Forest+Size+(1|Farm), family = "gaussian",data=paisagem)->result
  AIC(result)->akai
  summary(result)[[6]]->coefs
  coefs$cond->coef
  data.frame(Model=coef,rsquare=MuMIn::r.squaredGLMM(result)[1],aic=akai,scale=h)->values[[length(values)+1]]
}
do.call(rbind,values)->valuesR
valuesR$Model.Pr...z..<0.05->valuesR$p ;valuesR 
write.csv(valuesR,"valuesRichness.csv",row.names = T) # Save #
paisagem_rich<-pred_table_pad[pred_table_pad$Scale=="2250",] # Scale of effect = 2000 

## Forest specialist richness ##
values2<-{}
for (h in buffer_distance_meters ) {
  paisagem<-pred_table_pad[pred_table_pad$Scale==h,]
  paisagem$Richness_Forest<-resp_table$Richness_forest
  glmmTMB(Richness_Forest~Forest+Size+(1|Farm), family = "gaussian",data=paisagem)->result
  AIC(result)->akai
  summary(result)[[6]]->coefs
  coefs$cond->coef
  data.frame(Model=coef,rsquare=MuMIn::r.squaredGLMM(result)[1],aic=akai,scale=h)->values2[[length(values2)+1]]
}
do.call(rbind,values2)->valuesRfor
valuesRfor$Model.Pr...z..<0.05->valuesRfor$p ;valuesRfor
write.csv(valuesRfor,"valuesForestRichness.csv",row.names = T) # Save #
paisagem_richfor<-pred_table_pad[pred_table_pad$Scale=="2500",] # Scale of effect = 2500

#### Final Models ####

## Species richness ##
cbind(paisagem_rich,Richness=resp_table$Richness)->table_richall
model_rich<-glmmTMB(Richness ~ Forest + Size  + (1|Farm), family = "gaussian", data=table_richall) 
summary(model_rich)
r.squaredGLMM(model_rich)
# Correlations
# Fores amount x patch size
cor(paisagem_rich$Forest,paisagem_rich$Size)
# Forest amount x Farm
aov(paisagem_rich$Forest~paisagem_rich$Farm)->anova
summary(anova) 
TukeyHSD(anova)
# Richness x Farm
aov(table_richall$Richness~paisagem_rich$Farm)->resu
summary(resu)
TukeyHSD(resu)

## Forest specialist richness ##
cbind(paisagem_richfor,Richness=resp_table$Richness_forest)->table_richfor
model_richfor<-glmmTMB(Richness ~ Forest + Size  + (1|Farm), family = "gaussian", data=table_richfor) 
summary(model_richfor)
r.squaredGLMM(model_richfor)
# Correlations
# Fores amount x patch size
cor(paisagem_richfor$Forest,paisagem_richfor$Size)
# Forest amount x Farm
aov(paisagem_richfor$Forest~paisagem_richfor$Farm)->anova
summary(anova) 
TukeyHSD(anova)
# Richness x Farm
aov(table_richfor$Richness~paisagem_richfor$Farm)->resu
summary(resu)
TukeyHSD(resu)

# # # # # # # # # # # # # #
## Part 2. Composition ####
# # # # # # # # # # # # # #

#### Scale of effect ####
values3<-{}
for (h in buffer_distance_meters ) {
  paisagem<-pred_table_pad[pred_table_pad$Scale==h,]
  paisagem$Composition<-resp_table$PCoA1
  glmmTMB(Composition~Forest+Size+(1|Farm), family = "gaussian",data=paisagem)->result
  AIC(result)->akai
  summary(result)[[6]]->coefs
  coefs$cond->coef
  data.frame(Model=coef,rsquare=MuMIn::r.squaredGLMM(result)[1],aic=akai,scale=h)->values3[[length(values3)+1]]
}
do.call(rbind,values3)->valuesCOM 
valuesCOM$Model.Pr...z..<0.05->valuesCOM$p; valuesCOM
write.csv(valuesCOM,"valuesComposition.csv",row.names = T) # Save 
paisagem_com<-pred_table_pad[pred_table_pad$Scale=="1750",] # Scale of effect = 1500  

#### Final Models ####
cbind(paisagem_com,Composition=resp_table$PCoA1)->table_comp
model_comp<-glmmTMB(Composition ~ Forest + Size  + (1|Farm), family = "gaussian", data=table_comp) 
summary(model_comp)
r.squaredGLMM(model_comp)
# Correlations
# Fores amount x patch size
cor(paisagem_com$Forest,paisagem_com$Size)
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

## Figure scale of effects ##
rsquare<-c(unique(valuesR$rsquare),unique(valuesRfor$rsquare),unique(valuesCOM$rsquare))
scale<-rep(buffer_distance_meters,3)
cbind(rsquare,scale)->datar
as.data.frame(datar)->datarsquare
datarsquare$cores<-c("black","black","black","black","black","black","black","black","red","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","red","black","black","black","black","black","black","black","black","black","black","black","black","red","black","black","black","black","black","black","black","black","black")
datarsquare$plot<-c(rep("(a) Species Richness",16),rep("(b) Forest Specialist Richness",16),rep("(c) Species Composition",16))
figscale<-ggplot(data=datarsquare,mapping=aes(x=scale,y=rsquare))+
  geom_line(linetype="dashed")+
  geom_point(size=3,aes(colour = cores))+
  scale_color_manual(values = c("black","#ff3131"))+
  guides(color="none")+
  labs(x="Landscape extent", y="Adj. R²")+
  scale_x_continuous(n.breaks = 8,)+
  facet_wrap(.~plot,scales="free_y",nrow=3,ncol=1)+
  theme_minimal()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold"),legend.position = "right",legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1,legend.justification = "top"); figscale
# Save
tiff("Scales_of_effect_all.tif",width = 5,height = 6, units = "in",res = 600)
figscale
dev.off()

## Figure Richness x Habitat amount X Patch size ##
# Plot Species Richness
pred_table[pred_table$Scale==2250,]->pred_data_rich
pred_data_rich$Richness<-resp_table$Richness
richplot<-ggplot(data=pred_data_rich,mapping = aes(x = Forest, y = Richness, color=Farm, size = Size))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  labs(x="", y="Species\nRichness", color= "Farm", size= "Patch size")+
  theme_classic()+
  scale_color_manual("Farm",values = c("Esteio"="#f8b511","Dimona"="#ff3131", "Porto Alegre" ="#148bc0"))+
  scale_y_continuous("Species\nRichness",n.breaks = 8)+
  scale_size(breaks =unique(pred_data_rich$Size))+
  geom_point()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "bottom" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "center",legend.box = "vertical"); richplot
# Plot Forest Specialist Richness
pred_table[pred_table$Scale==2500,]->pred_data_richfor
pred_data_richfor$Richness<-resp_table$Richness_forest
richfeorestplot<-ggplot(data=pred_data_richfor,mapping = aes(x = Forest, y = Richness, color=Farm, size = Size))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  labs(x="", y="Frorest Specialist \nRichness", color= "Farm", size= "Patch size")+
  theme_classic()+
  scale_color_manual("Farm",values = c("Esteio"="#f8b511","Dimona"="#ff3131", "Porto Alegre" ="#148bc0"))+
  scale_y_continuous("Forest Specialist\nRichness",n.breaks = 8)+
  scale_size(breaks =unique(pred_data_rich$Size))+
  geom_point()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.position = "bottom" ,legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "center",legend.box = "vertical"); richfeorestplot
# Plot Composition
pred_table[pred_table$Scale==1750,]->pred_data_comp
pred_data_comp$Composition<-resp_table$PCoA1
pred_data_comp$Composition2<-resp_table$PCoA2
compplot<-ggplot(data=pred_data_comp,mapping = aes(x = Forest, y = Composition*-1, color= Farm, size = Size))+
  geom_smooth(method=lm, color='black',size=1,linewidth=0.6, se=F)+
  labs(x="Landscape habitat amount", y="Composition", color= "Farm", size= "Patch size")+
  theme_classic()+
  scale_color_manual("Farm",values = c("Esteio"="#f8b511","Dimona"="#ff3131", "Porto Alegre" ="#148bc0"))+
  scale_y_continuous("Species \nComposition",n.breaks = 8)+
  scale_size(breaks =unique(pred_data_comp$Size))+
  geom_point()+
  theme(strip.text=element_text(family="sans",size=12,face="bold",hjust = 0),strip.background = element_rect(fill="white",colour = "white"),axis.text =element_text(family="sans",size=12),axis.title = element_text(family="sans",size=14,face="bold",),legend.text =element_text(family="sans",size=12,hjust = 0),legend.title =element_text(family="sans",size=12,face="bold",hjust = 0),legend.text.align = 1, legend.justification = "center",legend.position = "bottom" ,legend.box = "vertical"); compplot
# Plot Final
figure<-ggarrange(richplot,richfeorestplot,compplot, 
          common.legend = T, legend="bottom",
          labels= c("(a)","(b)","(c)"),
          ncol= 1, nrow=3,
          align = "hv", 
          font.label = list(size = 14, color = "black", face = "bold", family = "sans", position = "top")); figure
# Save #
tiff("Habitat amount.tiff",width = 4,height = 9, units = "in",res = 600)
figure
dev.off()

## Figure PCoA axis ##
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

## Figure Generico ##
# Create the generico function
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
###
# Plot generico
tiff("generico.tif",width = 8,height = 20,units = "in",res = 300)
generico(occ,pred_data_comp$Forest,11,"Forest cover","Species","Sampling plots order by habitat amount",1)
dev.off()





