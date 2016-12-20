return.accurate<-function(x,data){
  prec<-data[x]
  return(x[which(prec==min(prec))][1])
}

wkb2geo<-function(x){
  x<-st_as_sfc(structure(list(as.character(x)),class="WKB"),EWKB=TRUE)
  return(as(x,"Spatial"))
}

sf2sp<-function(x){
  return(as(x,"Spatial"))
}



sp.extract<-function(x,exclude=500,acc=50,my.mask){

  ch<-odbcConnect("PostgreSQL35W")
  my.sp<-sqlQuery(ch,paste("
    SELECT public.observations.obs_id, public.observations.locality_id, public.observations.project_id, 
    public.observations.input_id, public.observations.releve_type, public.observations.obs_date, 
    public.observations.date_precision, public.observations.date_expert, public.observations.introduced, 
    public.observations.introduced_expert, public.observations.x_swiss, public.observations.y_swiss, 
    public.observations.xy_type, public.observations.xy_precision, public.observations.count_unit,
    public.observations.abundance_code, public.observations.v_taxon, public.observations.v_accepted_taxon_id, 
    public.observations.v_validation_status, public.observations.v_observers, public.observations.v_native_status, 
    public.observations.v_presence_status, public.observations.v_introduction_status, public.observations.v_doubt_status, 
    public.observations.v_rarity_index, public.observations.v_xy_radius, ST_AsText(ST_Simplify(ST_TRANSFORM(public.observations.geom_storage, 21781),25,TRUE)),
    ST_AsText(ST_Simplify(ST_TRANSFORM(public.observations.v_geom_buffer, 21781),25,TRUE))
    FROM public.observations 
    WHERE (((public.observations.v_accepted_taxon_id)=",x,"));
    ",sep=""))  #extract the database for a given species id_number

  odbcCloseAll()
  xcol<-11
  ycol<-12
  rcol<-26
  vcol<-19
  shpcol<-27
  bufcol<-28
  
  unvalid<-which(my.sp[,vcol]<200 | my.sp[,vcol]==250) 
  if (length(unvalid)>0){
    my.sp<-my.sp[-unvalid,]
  }
  my.sp<-my.sp[my.sp[,rcol]<exclude,]
  my.sp<-my.sp[-which(is.na(my.sp[,28])),]

  my.sp.pts<-my.sp[my.sp[,rcol]<=acc,]
  my.id<-c()
  my.data<-my.mask
  if (length(which(my.sp[,rcol]<=acc))>0){
    pts<-remove.duplicates(SpatialPointsDataFrame(coords=my.sp.pts[,c(xcol,ycol)],data=my.sp.pts,
      proj4string = CRS("+init=epsg:21781")),zero=acc)
    pts.r<-rasterize(pts,my.mask,100, background=0)
    pts.r.id<- -1*rasterize(pts,my.mask,my.sp.pts[,1])
    my.data<-my.data+pts.r
    my.id<-c(my.id,unique(pts.r.id))
  }
  
  if (nrow(my.sp.pts)<=nrow(my.sp)){
    my.sp.lyr<-my.sp[my.sp[,rcol]>acc,]
    lyr<-do.call(bind,lapply(my.sp.lyr[,bufcol],readWKT))
#    lyr<-list()
#    for (i in 1:length(my.sp[-index.pts,bufcol])){
#      lyr[[i]]<-wkb2geo(my.sp[-index.pts,bufcol][i])
#    }
#    lyr<-do.call(bind,lyr)
#    to.rm<-which(clgeo_CollectionReport(lyr)[,2]==FALSE)
#    lyr<-lyr[-to.rm]
    crs(lyr)<-CRS("+init=epsg:21781")
    lyr.clean<-over(lyr,lyr,returnList = TRUE)
    to.keep<-which(sapply(lyr.clean,length)==1)    
    to.keep<-c(to.keep,unique(sapply(lyr.clean[-to.keep],return.accurate,data=my.sp.lyr[,26]))) ### to check
    lyr<-lyr[to.keep]
    my.sp.lyr<-my.sp.lyr[to.keep,]
    lyr.r.id<-min(rasterize(as(lyr,"SpatialLines"),my.mask,field=my.sp.lyr[,1], fun='max'),
                            rasterize(lyr,my.mask,my.sp.lyr[,1], fun='max'),na.rm=T)
   
    lyr.warea<-round(10000/area(lyr)*100)
    lyr.warea[lyr.warea[]>100]<-100
    lyr.r<-min(rasterize(as(lyr,"SpatialLines"),my.mask,field=lyr.warea, fun='max'),
               rasterize(lyr,my.mask,lyr.warea, fun='max'),na.rm=T)
    lyr.r[is.na(lyr.r)]<-0
    
    my.data<-my.data+lyr.r
    my.id<-c(my.id,unique(lyr.r.id))
  }
  
  if(exists("pts")&exists("lyr.clean")){
    lyr.id.rm<-unique(lyr.r.id*(pts.r.id/pts.r.id))
    my.id<-my.id[-match(lyr.id.rm,my.id)]
    lyr.r.id<-reclassify(lyr.r.id,cbind(lyr.id.rm,rep(0,length(lyr.id.rm))))
    my.sp.lyr<-my.sp.lyr[-match(lyr.id.rm,my.sp.lyr[,1]),]    
    my.data[my.data[]>100]<-100
  }
  
  my.data[my.data[]==0]<-NA
  dist<-rasterToPoints(my.data)
  proj.dist<-knnx.dist(my.sp[match(my.id,my.sp[,1]),c(xcol,ycol)],rasterToPoints(my.mask)[,1:2],k=1)
  return(list(dist=dist,my.sp.pol=my.sp.lyr,my.sp.pts=my.sp.pts,proj.dist=proj.dist))
}  

  
writeOGR(pts,dsn="S:/Users/Petitpierre/Data/tmp",layer="pts2",driver="ESRI Shapefile")

lyr2<-SpatialPolygonsDataFrame(Sr=lyr,data=my.sp.lyr,match.ID = F)
writeOGR(lyr2,dsn="S:/Users/Petitpierre/Data/tmp",layer="poly2",driver="ESRI Shapefile")

