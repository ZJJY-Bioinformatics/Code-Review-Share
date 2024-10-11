#----------------------
# Code form  published datasets
# Cite :Deterministic transition of enterotypes shapes the infant gut microbiome at an early age
#----------------------
setwd("MES")
#Uncomment next two lines if R packages are already installed
#install.packages("cluster")
#install.packages("clusterSim")
#install.packages("ade4")
rm(list = ls())
library(cluster)
library(clusterSim)
library(ade4)

#Download the example data and set the working directory
#setwd('<path_to_working_directory>')
data=read.table("exp1.tx", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

require(clusterSim)

data.dist=dist.JSD(data)
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")

data.cluster=pam.clustering(data.dist, k=3)
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
cat(obs.silhouette) #0.1899451

#data=noise.removal(data, percent=0.01)

## plot 1
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
dev.new()
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")


#plot 2
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
dev.new()
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")


cls_data = data.frame(sample = colnames(data),cluster = data.cluster)
write.csv(cls_data,"sample_cluster.csv")
#Uncomment next two lines if R packages are already installed
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#install.packages("igraph")
#install.packages("markovchain")
library(phyloseq)
library(igraph)
library(markovchain)

# Markov-------------
rm(list = ls())

#set function to filter samples
samdat.prune_prev <- function(samdat) {
  GAP_MIN <- 120
  GAP_MAX <- 1000
  samdf <- data.frame(samdat)
  subjects <- unique(samdf$SubjectID)
  csub <- split(samdf, samdf$SubjectID)
  for(sub in subjects) {
    cc <- csub[[sub]]
    cc <- cc[order(cc$GDColl),]
    cc$PrevID <- c(NA, cc$SampleID[-nrow(cc)])
    del <- cc$GDColl - c(-999, cc$GDColl[-nrow(cc)])
    keep <- del>=GAP_MIN & del<=GAP_MAX
    if(sum(keep) == 0) {
      csub[[sub]] <- NULL
    } else {
      cc <- cc[keep,]
      csub[[sub]] <- cc
    }
  }
  return(do.call(rbind, csub))
}

##calculate transition probability
samdf = read.table("data1.txt",header = T,sep = '\t',row.names = 1)

rownames(samdf) -> samdf$SampleID
samdf$type = as.factor(samdf$type)
CSTs <- levels(samdf$type)
nstates <- nlevels(samdf$type)

samdf_prev <- samdat.prune_prev(samdf)
rownames(samdf_prev) <- samdf_prev$SampleID
samdf_prev$PrevCST <- data.frame(samdf)[samdf_prev$PrevID,"type"]
samdf_prev$CurCST <- samdf_prev$type
ttab <- table(samdf_prev$PrevCST, samdf_prev$CurCST) # prevstate=row, curstate=col
trans <- matrix(ttab, nrow=nstates)
trans <- trans/rowSums(trans)  # Normalize row sums to 1
CSTtrans <- trans
colnames(CSTtrans) <- CSTs
rownames(CSTtrans) <- CSTs
t_persist <- -1/log(diag(CSTtrans))

##plot markov chain
mcPreg <- new("markovchain", states=CSTs,
              transitionMatrix = trans, name="PregCST")
netMC <- markovchain:::.getNet(mcPreg, round = TRUE)
wts <- E(netMC)$weight/100
edgel <- get.edgelist(netMC)
elcat <- paste(edgel[,1], edgel[,2])
elrev <- paste(edgel[,2], edgel[,1])
edge.curved <- sapply(elcat, function(x) x %in% elrev)
default.par <- par(no.readonly = TRUE)
plotMC <- function(object, ...) {
  netMC <- markovchain:::.getNet(object, round = TRUE)
  plot.igraph(x = netMC, ...)  
}
vert.sz <- 2*sapply(states(mcPreg), 
                    function(x) nrow(unique(sample_data(samdf)[sample_data(samdf)$type==x,"SubjectID"])))
vert.sz <- log(vert.sz)*5
vert.font.clrs <- c("white", "white", "white", "white")
edge.loop.angle = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.14, 0, 0, 0, 0, 0)-0.45
layout <- matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3), nrow=4, ncol=2, byrow=T)
#layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
par(mar=c(0,1,1,1)+0.1)
edge.arrow.size=1
edge.arrow.width=1
edge.width = (6*wts + 0.1)
edge.labels <- as.character(E(netMC)$weight/100)
edge.labels[edge.labels<0.2] <- NA  # labels only for self-loops

##plot
plotMC(mcPreg, edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
       edge.label = edge.labels, edge.label.font=2, edge.label.cex=1.3, edge.label.color="black",
       edge.width=edge.width, edge.curved=edge.curved, 
       vertex.size=(vert.sz),
       vertex.label.font = 2, vertex.label.cex = 2,
       vertex.label.color = vert.font.clrs, vertex.frame.color = NA, 
       vertex.color = c("#CDBE6A","#89CFBE","#86A4CF","#DFAD9F"),
       layout = layout, edge.loop.angle = edge.loop.angle)

