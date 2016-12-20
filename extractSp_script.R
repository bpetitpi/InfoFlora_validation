library("RODBC")
#library(RPostgreSQL)
library(parallel)
library(rgdal)
library(rgeos)
library(raster)
library(FNN)
#library(sf)
#library(cleangeo)

source("extractSp_functions.R")

#my.mask<-raster("S:/Users/Petitpierre/Data/env/bio6100.tif")
#my.mask[!is.na(my.mask)]<-0
#dataType(my.mask)<-"INT1S"
#writeRaster(my.mask,"S:/Users/Petitpierre/Data/env/mask.tif", format="GTiff")
my.mask<-raster("S:/Users/Petitpierre/Data/env/mask.tif")
#
#"D:/InfoFlora/InfoFlora_EvaluationTool/BP/"

ch<-odbcConnect("PostgreSQL35W")

sp.list<-sqlQuery(ch,"SELECT public.observations.v_taxon, public.observations.v_accepted_taxon_id
FROM public.observations;") #generate the species list

sp.list<-unique(sp.list)
write.table(sp.list,file=paste(my.dir,"data/sp/sp_list.txt",sep=""),sep="/t",quote=F,row.names=F)
sp.list<-read.delim(paste(my.dir,"data/sp/sp_list.txt",sep=""),h=TRUE,sep="/t")
detectCores()

SELECT public_genre.nom_genre, public_nom.nom_espece, public_nom.auteur_espece, public_nom.rang, public_nom.nom_infra_sp, public_nom.auteur_infra_sp, public_nom.nom_complet, public_nom.no_genre, public_nom.no_nom
FROM public_genre INNER JOIN public_nom ON public_genre.no_genre = public_nom.no_genre;

