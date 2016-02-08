#########################################################
#  Multi Omic Oscillations Networks
#
#  Analysis of Omic and Multi Omic Patterns - An Example
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################


library(classInt)
library(KernSmooth)

#Load the NBA weights for indegree and outdegree
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
             "https://raw.githubusercontent.com/lodeguns/MOOs-An/master/MOOsFun.R")


osc.eco.poly.ms.and <- osc.eco.poly.ms

for(j in 1:length(osc.eco.poly.ms.and))
{
  osc.eco.poly.ms.and[[j]] <- (osc.eco.operon.ms[[j]] & osc.eco.poly.ms[[j]])*1
  
}



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

tryCatch( val1 <- osc.pst.pred.diss(osc.path.list,
                                   osc.eco.mat,
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
