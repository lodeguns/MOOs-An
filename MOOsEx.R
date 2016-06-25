#########################################################
#  Multi Omic Oscillations - - - An Example
#
#  Novel algorithms to detect oscillatory patterns in 
#  multi omic metabolic networks.
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################


library(classInt)
library(KernSmooth)
library(PST)
library(cluster)


# Load 14 data-subset.

load('./ktable.eco.sel.RData')
load('./osc.list.w.in.Rdata')
load('./osc.list.w.out.Rdata')
load('./osc.eco.poly.ms.Rdata')
load('./osc.eco.operon.ms.Rdata')
load('./osc.paths.list.Rdata')
load('./osc.paths.list.cai.poly.Rdata')
load('./osc.paths.list.cai.operons.Rdata')
load('./osc.paths.list.pa.poly.Rdata')
load('./osc.paths.list.pa.operons.Rdata')
load('./osc.paths.list.caipa.poly.Rdata')
load('./osc.paths.list.caipa.operons.Rdata')
load('./osc.paths.list.suzukietall.Rdata')
load('./osc.paths.list.norfloxacin.Rdata')

#Load Functions
source_https <- function(url, ...) {
  require(RCurl)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://github.com/lodeguns/NetBasAdj/raw/master/NetBasAdj.R",
             "https://github.com/lodeguns/MOOs-An/raw/master/MOOsFun.R")

#If the are some erros in the automatic-download of the functions, 
#please dowload the files manually from the two repository url above.


osc.eco.poly.ms.and <- osc.eco.poly.ms

for(j in 1:length(osc.eco.poly.ms.and))
{
  osc.eco.poly.ms.and[[j]] <- (osc.eco.operon.ms[[j]] & osc.eco.poly.ms[[j]])*1
  
}


#################################################################
# Paper's Figures 5 (a) and (b) - Multi Omic Oscillations
#
#################################################################
path.index <- 1

osc.path.list       <-  osc.paths.list.caipa.operons[path.index]
osc.path.list[[1]]  <-osc.path.list[[1]][osc.path.list[[1]][3] != "kmeans",]
osc.path.list[[1]]  <-osc.path.list[[1]][osc.path.list[[1]][3] != "hclust",]

osc.eco.mat         <-  osc.eco.operon.ms[[path.index]]
osc.w.out           <-  osc.list.w.out[path.index]
osc.w.in            <-  osc.list.w.in[path.index]
treat.path          <-  osc.paths.list.suzukietall[path.index]
omic.val            <- "CAIPA"
h                   <- 2
n.seq               <- 2


tryCatch( seq.osc.caipa.caipa.st2 <- osc.pst.pred.diss(osc.path.list, 
                                                       osc.eco.mat, osc.w.in, osc.w.out, treat.path,
                                                       omic.val, h, n.seq, TRUE, FALSE),
          error = function(e) {
            print("No way...")
            return(NULL)})

##################################################################
# Between the generated files see:
# (a) CAIPAtoCAIPA22 normal indeg.png
# (b) CAIPAtoCAIPA22 operons compression max outdeg.png
#################################################################
osc.print.clusters(seq.osc.caipa.caipa.st2, treatlabtarget, "CAIPAtoCAIPA22")

#################################################################

#################################################################
# Paper's Figures 5 (c) Single omic oscillations
#################################################################
path.index <- 1

osc.path.list       <-  osc.paths.list.cai.operons[path.index]
osc.path.list[[1]]  <-osc.path.list[[1]][osc.path.list[[1]][3] != "kmeans",]
osc.path.list[[1]]  <-osc.path.list[[1]][osc.path.list[[1]][3] != "hclust",]

osc.eco.mat         <-  osc.eco.operon.ms[[path.index]]
osc.w.out           <-  osc.list.w.out[path.index]
osc.w.in            <-  osc.list.w.in[path.index]
treat.path          <-  osc.paths.list.suzukietall[path.index]
omic.val            <- "PA"
h                   <- 2
n.seq               <- 2


tryCatch( seq.osc.cai.pa.st2 <- osc.pst.pred.diss(osc.path.list, 
                                                  osc.eco.mat, osc.w.in, osc.w.out, treat.path,
                                                  omic.val, h, n.seq, TRUE, FALSE),
          error = function(e) {
            print("No way...")
            return(NULL)})

##################################################################
# Between the generated files see:
# (c) CAItoPA22 operons compression min outdeg.png
#################################################################
osc.print.clusters(seq.osc.cai.pa.st2, treatlabtarget, "CAItoPAst22")
#################################################################

#################################################################
# Paper's Figures 5 (d) Single omic oscillations
#################################################################


path.index <- 1

osc.path.list       <- osc.paths.list.pa.operons[path.index]
osc.path.list[[1]]  <- osc.path.list[[1]][osc.path.list[[1]][3] != "kmeans",]
osc.path.list[[1]]  <- osc.path.list[[1]][osc.path.list[[1]][3] != "hclust",]

osc.eco.mat         <-  osc.eco.operon.ms[[path.index]]
osc.w.out           <-  osc.list.w.out[path.index]
osc.w.in            <-  osc.list.w.in[path.index]
treat.path          <-  osc.paths.list.suzukietall[path.index]
omic.val            <- "PA"
h                   <- 2
n.seq               <- 7


tryCatch( seq.osc.pa.pa.st7 <- osc.pst.pred.diss(osc.path.list, 
                                                 osc.eco.mat, osc.w.in, osc.w.out, treat.path,
                                                 omic.val, h, n.seq, TRUE, FALSE),
          error = function(e) {
            print("No way...")
            return(NULL)})




##################################################################
# Between the generated files see:
# (d) PAtoPA27 operons compression min outdeg.png
#################################################################
osc.print.clusters(seq.osc.pa.pa.st7, treatlabtarget, "PAtoPAst27")
#################################################################



###########################################################################################################
##################################################################################################
#   Figure 6 and Figure Online - Cophenetic Correlogram 
#   Warning the computing requires several hours (approximately 20h)
#
################################################################################################

library(classInt)
library(KernSmooth)
library(PST)
library(cluster)

#function
#clust.hc.diss <- function(seq.osc.caipa.pa.list, j){
#  seq.osc.l.in  <- seq.osc.caipa.pa.list[[j]][[1]]
#  seq.osc.l.out <- seq.osc.caipa.pa.list[[j]][[2]]
#  type.training <- paste(kegg.code, seq.osc.likh[[j]][[3]], sep = "  ")
#  
#  df.out   <- osc.get.diss.mat(seq.osc.l.out)
#  df.in    <- osc.get.diss.mat(seq.osc.l.in )
#  
#  diss.out <- daisy(df.out, metric = c("euclidean", "manhattan", "gower"),
#                    stand = FALSE, type = list())
#  diss.in  <- daisy(df.in, metric = c("euclidean", "manhattan", "gower"),
#                    stand = FALSE, type = list())
#  
#  hc1<-hclust(diss.out, method="complete")
#  hc2<-hclust(diss.in, method="complete")
#  return(list(hc1,hc2))
#}



osc.caipa.caipa.corr <- list()

for(j in 1:length(ktable.eco.sel[,1]))   ## i.e You can try on a small set 1:3 
{ 
  
  path.index <- j
  
  osc.path.list       <-  osc.paths.list.caipa.operons[path.index]
  osc.path.list[[1]]  <-  osc.path.list[[1]][osc.path.list[[1]][3] != "kmeans",]
  osc.eco.mat         <-  osc.eco.operon.ms[[path.index]]
  osc.w.out           <-  osc.list.w.out[path.index]
  osc.w.in            <-  osc.list.w.in[path.index]
  treat.path          <-  osc.paths.list.norfloxacin[path.index]
  omic.val            <- "CAIPA"
  h                   <- 2
  n.seq               <- 3
  
  ptm <- proc.time()
  
  
  ret <- tryCatch( val <- osc.pst.pred.diss(osc.path.list,
                                            osc.eco.mat,
                                            osc.w.in, 
                                            osc.w.out,treat.path, 
                                            omic.val, h, n.seq,
                                            FALSE, TRUE) ,
                   error = function(e) {
                     print("Ops...")
                     return(NULL)})
  
  
  
  if(is.null(ret))
  {
    osc.caipa.caipa.corr[length(osc.caipa.caipa.corr)+1] <- NULL
  } else {
    
    osc.caipa.caipa.corr[length(osc.caipa.caipa.corr)+1] <- val 
    }
} #end for


names(osc.caipa.caipa.corr)    <-  ktable.eco.sel$KEGGpath #small set $KEGGpath[c(1,2,3)]

m.t <- osc.caipa.caipa.corr

clean.t <- list()
str<-c()
for(j in 1:length(m.t))
{
  if(!is.null(m.t[j][[1]][[2]][[1]]))
  { if(m.t[j][[1]][[2]][[1]][[1]]$l.s[1] != 0){
    x<- length(clean.t)+1
    clean.t[x] <- m.t[j]
    clean.t[x][[1]][[3]] <- names(m.t[j])
    str <- c(str, names(m.t[j]))
  } else
  {
    print(j)
    print(names(m.t[j]))
  }
  }
  
  
}

names(clean.t) <- str


correlogram.caipa.caipa.d <- matrix(data = 0, 
                                    nrow = length(clean.t), 
                                    ncol = length(clean.t), byrow = FALSE,
                                    dimnames = list(names(clean.t), 
                                                    names(clean.t)))


#n*n
for(j in 1:length(clean.t)){
  diss.j <- osc.get.diss.mat(clean.t[[j]][[2]])
  
  diss.j.cor <- diss.j
  diss.j.cor[is.na(diss.j.cor)] <- 0
  
  diss.j.d <- daisy(diss.j, metric = c("euclidean", "manhattan", "gower"),
                    stand = FALSE, type = list())
  
  ret <- tryCatch( hc1<-hclust(diss.j.d, method="complete") ,
                   error = function(e) {
                     print("1")
                     return(NULL)})
  
  for(j1 in 1:length(clean.t))
  { 
    print(j1)
    diss.j1 <- osc.get.diss.mat(clean.t[[j1]][[2]])
    
    diss.j1.cor <- diss.j1
    diss.j1.cor[is.na(diss.j1.cor)] <- 0
    
    diss.j1.d <- daisy(diss.j1, metric = c("euclidean", "manhattan", "gower"),
                       stand = FALSE, type = list())
    ret1<-tryCatch( hc2<-hclust(diss.j1.d, method="complete") ,
                    error = function(e) {
                      print("2")
                      return(NULL)})
    if(is.null(ret) || is.null(ret1)){
      correlogram.caipa.caipa.d[j,j1] <- NA
    } else
    {
      correlogram.caipa.caipa.d[j,j1] <- round(cor(cophenetic(as.dendrogram(hc1)), cophenetic(as.dendrogram(hc2))), digits= 7)
      
    }
    
  }
  
}



ind <- apply(correlogram.caipa.caipa.d, 1, function(x) all(is.na(x)))



correlogram.caipa.caipa.d <- correlogram.caipa.caipa.d[
  colSums(is.na(correlogram.caipa.caipa.d)) != nrow(correlogram.caipa.caipa.d),
  colSums(is.na(correlogram.caipa.caipa.d)) != nrow(correlogram.caipa.caipa.d)]



library(d3heatmap)

#Heatmap on the Paper
d3heatmap(correlogram.caipa.caipa.d , 
          scale = "none", 
          colors = "Spectral",
          dendrogram = "both", k_col = 9, k_row=9)

save(correlogram.caipa.caipa.d,    file="correlogram.caipa.caipa.d.RData")



correlogram.caipa.caipa.d.names <- correlogram.caipa.caipa.d
colnames.d <- colnames(correlogram.caipa.caipa.d.names)

for(j in 1:length(colnames.d))
{
  grep(colnames.d[j], ktable.eco.sel$KEGGpath)
  colnames.d[j]<- paste(colnames.d[j],gsub( "-.-.*$", "", ktable.eco.sel[grep(colnames.d[j], ktable.eco.sel$KEGGpath),]$descr ), sep="-")
  
}
colnames(correlogram.caipa.caipa.d.names)<-colnames.d
rownames(correlogram.caipa.caipa.d.names)<-colnames.d

#Heatmap online with complete descriptions
d3heatmap(correlogram.caipa.caipa.d.names , 
          scale = "none", 
          colors = "Spectral",
          dendrogram = "both", k_col = 9, k_row=9)


save(correlogram.caipa.caipa.d.names,    file="correlogram.caipa.caipa.d.names.RData")


########################################################################################



























######################################################################################################################################
#  Example of clustering of the PSTlikelihoods.
#  In the training step is changed the number of M.O. patterns given in input in 1,3 and 7.
#  Here we are considering only Multi Omic patterns (CAI and PA) and operons compression.
#  Warning the whole process could take several time.
######################################################################################################################################

osc.path.list       <-  osc.paths.list.caipa.operons[1] # Al the possible multi omi patterns standard and operons compressed.
osc.meth.table.rel  <-  osc.get.rel.paths(osc.path.list)
osc.eco.mat         <-  osc.eco.operon.ms[[1]]          # Adjacency Matrix of 0 and 1, if 1 occurs an operon.
osc.w.out           <-  osc.list.w.out[1]               # Reactoma weights indegree with Network based Adjacency Algorithm
osc.w.in            <-  osc.list.w.in[1]                # Reactoma weights indegree with Network based Adjacency Algorithm
treat.path          <-  osc.paths.list.suzukietall[1]   #Treatments from Suzuki et All. http://www.nature.com/ncomms/2014/141217/ncomms6792/pdf/ncomms6792.pdf
omic.val            <- "CAIPA"                          #Codon Adaptation Index and Protein Abundance

h                   <- 2                                #Oscillation units - Short Memory

n.seq               <- 1                                # Number of Multi Omic Patterns considered with the best Osc.Score.
n.seq2              <- 3
n.seq3              <- 7
w                   <- TRUE                             # Boolean value, PST likelihoods with osc.w.out or osc.w.in?
nocompr             <- FALSE                            # Boolean value, do you want to consider only patterns without compression?





tryCatch( val1 <- osc.pst.pred.diss(osc.path.list,osc.eco.mat,
                                    osc.w.in, 
                                    osc.w.out,treat.path, 
                                    omic.val, h, n.seq,
                                    w, nocompr) ,
          error = function(e) {
            print("No way...")
            return(NULL)})

tryCatch( val2 <- osc.pst.pred.diss(osc.path.list,
                                    osc.eco.mat,
                                    osc.w.in, 
                                    osc.w.out,treat.path, 
                                    omic.val, h, n.seq2,
                                    w, nocompr) ,
          error = function(e) {
            print("No way....")
            return(NULL)})


tryCatch( val3 <- osc.pst.pred.diss(osc.path.list,
                                    osc.eco.mat,
                                    osc.w.in, 
                                    osc.w.out,treat.path, 
                                    omic.val, h, n.seq3,
                                    w, nocompr) ,
          error = function(e) {
            print("No way....")
            return(NULL)})


treatlabtarget <- c("CellWall-01", "CellWall-02", "Rib3-03", "Rib3-04", "Rib3-05", "Rib5-06", "Rib5-07", "FolicAcid-08", "DNAgy-09", "DNAgy-10")
kegg.code      <- "eco00010"

osc.print.clusters(val1, treatlabtarget, "CAIPAtoCAIPA1")
osc.print.clusters(val2, treatlabtarget, "CAIPAtoCAIPA3")
osc.print.clusters(val3, treatlabtarget, "CAIPAtoCAIPA7")


