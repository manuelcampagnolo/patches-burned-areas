# author: Manuel campagnolo, ISA/ULisboa, mlc@isa.ulisboa.pt
# last change: june 26, 2019

########################################################################## R packages
library(raster)
library(data.table)
library(igraph)
library(RANN)

######################################################################### test data set
# data: largest single tif file from SFD
r<-raster("D:\\BA-patches\\20160801-ESACCI-L3S_FIRE-BA-MSI-AREA_h41v20-fv1.1-JD.tif") # 774508900 pixels

# for tests: 
#r<-crop(r,extent(r)/50) # 50^2 times fewer pixels

################################################################################### parameters
N<-10 # number of "horizontal" blocks (to be able to process large rasters)
DELTAD<-1 # temporal gap: if 1, then two BA pixels will belong to the same patch if they are neighbors in space and dates differ by at most 2*D

# convert lon/lat into sinusoidal x/y
R<-6371007
latlong2sin<-function(lon,lat) list(x= R* (pi * lon / 180) * cos (pi * lat /180), y=R * (pi * lat /180))

# determine approximate distance (in meters) between to adjacent diagonal grid points
res<-res(r)
ext<-as.vector(extent(r))
deltalon<-res[1]; deltalat<-res[2]
aux<-latlong2sin(lon=c(mean(ext[1:2]),mean(ext[1:2])+deltalon),lat=c(mean(ext[3:4]),mean(ext[3:4])+deltalat))  
d<-sqrt(diff(aux$x)^2+diff(aux$y)^2) # d~30 meters

# define distance to determine patches
d<-25 # if one considers only 4 neighbors
K<-5 # neighbors for nn2

############################################################################ auxiliary functions
# input: xydate=data.table with columns x, y, dtmin and dtmax: one row per pixel
# K=number of neighbors to explore
# DIST=max distance in meters to be a neighbor
create.st.patches<-function(xydate=xyd,K,DIST)
{
  # determine neighbors
  knn<-RANN::nn2(data=cbind(xydate$x,xydate$y,xydate$dtmin), k=K, searchtype = "radius", radius=DIST) 
  knn$nn.dists<-NULL
  neigh<-cbind(1:nrow(knn$nn.idx),knn$nn.idx)
  rm(knn)
  # create adjency matrix 
  edges<-data.table::data.table(v1=rep(neigh[,1],ncol(neigh)-1),v2=as.vector(neigh[,-1]))
  rm(neigh)
  edges<-edges[v2!=0]
  edges[,dtumin1:=xydate$dtmin[v1]]
  edges[,dtumax1:=xydate$dtmax[v1]]
  edges[,dtumin2:=xydate$dtmin[v2]]
  edges[,dtumax2:=xydate$dtmax[v2]]
  # remove edges which time span does not intersect
  edges<-edges[dtumax1>=dtumin2 & dtumin1 <= dtumax2]
  edges[,dtumin1:=NULL]
  edges[,dtumax1:=NULL]
  edges[,dtumin2:=NULL]
  edges[,dtumax2:=NULL]
  # clusters
  Gall<-igraph::graph_from_data_frame(d=cbind(edges$v1,edges$v2),directed=FALSE)
  G<-igraph::simplify(Gall)
  
  CL<-igraph::cluster_louvain(G) # extract patches (max modularity algorithm)
  
  # indices of patches
  idxs<-1:nrow(xydate)
  return(igraph::membership(CL)[as.character(idxs)])
}

######################################################
# main function: to be applied to the i-th  block 
# it returns a data.table xyd with column "membership" of patch indices
f<-function(i)
{
  print(paste(i, "/", N))
  if (i<N) {ROW=tr$row[i];NROWS=tr$nrows[i]}
  if (i==N) {ROW=tr$row[i];NROWS=nrow(r)-tr$row[i]+1}  
  
  #  convert raster to data.data.table & remove unburned and unclassified (dt>0)
  vr<-getValues(r,row=ROW,nrows=NROWS)
  cells <- cellFromRowCol(r, c(ROW, ROW + NROWS - 1), c(1, ncol(r)))
  coordsr <- xyFromCell(r, cell = cells[1]:cells[2])xyd<-data.table(x=coordsr[,1], y=coordsr[,2], dt=vr) 
  xyd<-xyd[dt>0] # & x==x1 & y==y1]
  colnames(xyd)<-c("lon","lat","dt")
  
  # convert coordinates to meters
  xy<-latlong2sin(xyd$lon,xyd$lat)
  xyd$x<-xy$x
  xyd$y<-xy$y
  rm(xy)
  
  # create dtmin and dtmax
  xyd$dtmin<-xyd$dt-DELTAD
  xyd$dtmax<-xyd$dt+DELTAD
  
  # eliminate non necessary columns
  xyd[,dt:=NULL]
  xyd[,lon:=NULL]
  xyd[,lat:=NULL]
  
  # determine patches
  xyd$membership<-create.st.patches(xydate=xyd,K=K,DIST=d) # 
  return(xyd)
}

#################################################################################### main

# use function "blockSize" to break up raster r in N blocks
tr <- blockSize(r, minrows=floor(dim(r)[1]/N))

# compute patches for each block
# output: list of data.tables xyd with columns  x y dtmin dtmax membership
lista<-lapply(1:N, f)

# re-number patch indices in the N blocks so that patch indices do not overlap
Npatches<-0
for (k in 1:length(lista))
{
  lista[[k]]$membership<-lista[[k]]$membership+Npatches
  Npatches<-max(lista[[k]]$membership)
}

xyd <- do.call(rbind,lista) # 185754834 rows
rm(lista)
# to save results as a txt file
# fwrite(xyd,file="D:\\BA-patches\\20160801-ESACCI-L3S_FIRE-BA-MSI-AREA_h41v20-fv1.1-JD-patches.txt")

####################################################################### explore results
# explore patch sizes
uniqueN(xyd$membership) # 6429599 patches
freqs<-xyd[,.N,by=membership]
setorder(freqs,N)
freqs$cumul<-cumsum(freqs$N)
quantile(freqs$N,prob=seq(0,1,0.1))
quantile(freqs$cumul,prob=seq(0,1,0.1))


################################################################################# plot results
# plot a small area ~10*10 km at the center of r
x1<-mean(xyd$x)
x2<-x1+10000
y1<-mean(xyd$y)
y2<-y1+10000
xydate<-xyd[x>x1 & x<x2 & y>y1 & y<y2] # change condition if needed
tbl<-table(xydate$membership)
head(tbl)
xydate$sizecomp<-tbl[as.character(xydate$membership)]
xydate<-xydate[sizecomp>=10]
Ncomp<-uniqueN(xydate$membership)

plot(xydate$x, xydate$y, pch="")
for (m in unique(xydate$membership)) 
  points(xydate[membership==m]$x,
         xydate[membership==m]$y,
         col=randomcoloR::randomColor(Ncomp)[sample(Ncomp,1)],
         pch=".",cex=2)

