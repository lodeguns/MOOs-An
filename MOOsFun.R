#########################################################
#  Multi Omic Oscillations - - - Functions 
#
#  Novel algorithms to detect oscillatory 
#  patterns in multi omic metabolic networks. 
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################
library(classInt)
library(KernSmooth)
library(PST)
library(cluster)


#Paragraph 2.5 - Identification of the Oscillatory Best Patterns
#  In particular for each i-th metabolic pathway are selected the ones
#  with the highest score bigger than $ median(score * P_u)$, where $P_u$
#  is the probability that the oscillation units are different from 0. 
osc.get.method <- function(osc.meth.table, path.name)
{
      pp.score <- as.numeric(paste(osc.meth.table$pvalueseq)) * 
                  as.numeric(paste(osc.meth.table$osc.score))
  
      df <- osc.meth.table[pp.score >= median(pp.score),]
      if(nrow(df)!= 0)
      {  df["path.name"]<- path.name
         return(df)
      }
  return(NULL)
}

# Coupled with the previous function.
osc.get.rel.paths <- function( osc.paths.list.cross)
{
  df1 <- data.frame()
  for( j in 1:length(osc.paths.list.cross))
  {
    osc.meth.table <-  osc.paths.list.cross[[j]]
    path.name      <-  names(osc.paths.list.cross[j])
    
    df <- osc.get.method(osc.meth.table, path.name)
    #print(df)
    if(!is.null(df))
    {
      
      df1 <- rbind(df1, df)
      
    }
  }
  
  return(df1)
}


#  This function is of central importance, represents a sort of interface 
#  with the specific options enumerated in the front-end example:
#  https://github.com/lodeguns/MOOs-An/blob/master/MOOsEx.R

osc.pst.pred.diss <- function(osc.path.list, osc.eco.mat, osc.w.in, osc.w.out, treat.path, omic.val, h, n.seq, w, nocomp)
{  
  osc.meth.table.rel  <-  osc.get.rel.paths(osc.path.list)
  osc.meth.table.rel  <-  data.frame(lapply(osc.meth.table.rel, as.character), stringsAsFactors=FALSE)
  osc.meth.table.rel   <- osc.meth.table.rel[which(names(treat.path)== osc.meth.table.rel$path.name),]
  
  # If it isn't selected the w_cmp option only the no-compressed patterns are considered.
  if(nocomp == TRUE)
  { 
    osc.types.list <- "normal" 
  } else
  {
    osc.types.list      <-  unique(osc.meth.table.rel$type)
  }
  
  
  # In the tables select all the columns that represent a protein abundance variation for the
  # specific pathway.
  treats.idx          <-  grep("PF",names(treat.path[[1]]))
  
  list.in.out.types <- list()
  
  for(j1 in 1:length(osc.types.list))
  {
    list.out <- list()
    list.in  <- list()
    osc.types.sel <-  osc.types.list[j1]
    
    for(j in 1:length(treats.idx))
    {
      osc.treat.lab <- names(treat.path[[1]][treats.idx[j]])
      
      # Make the training step and prediction step of the PSTs on
      # the selected pathways of a specific metabolic pathway.
      # Save a list, also considering if there is an interactome layer 
      # (in going edges and out going edges).
      ll <- osc.pst.train.pred(treat.path, 
                               osc.treat.lab, 
                               osc.eco.mat, 
                               osc.w.out, 
                               osc.w.in, 
                               osc.meth.table.rel, 
                               osc.types.sel,
                               omic.val,
                               h,
                               n.seq, w)

      if(is.atomic(ll))
      {
        return(ll)
      }
      
      if(!is.null(ll))
      {df.in  <- ll$osc.pst.pred.in}
      else
      {df.in <- NULL}
      
      if(!is.null(ll))
      {df.out  <- ll$osc.pst.pred.out}
      else
      {df.out <- NULL}
      
      if(w==FALSE)
      {df.out<-ll$osc.pst.pred}
      
      
      list.in[length(list.in)+1]   <- list(df.in,  osc.treat.lab)
      list.out[length(list.out)+1] <- list(df.out, osc.treat.lab)
      
    }
    
    list.in.out.types[[length(list.in.out.types)+1]] <- list(list.in, list.out, osc.types.sel)
    
    
  }
  
  return(list.in.out.types)
  
}



# Make a PST's prediction step and training step.
# In the same function, (absolutely not optimized), is done the discretisation
# of the patterns, the compression, the training of the PST, the prediction
# of the PST and the weighting with the two types of weights described in
# the paper. Some functions are coupled with the first repository of functions:
# https://github.com/lodeguns/NetBasAdj/blob/master/NetBasAdj.R
osc.pst.train.pred <- function(treat.path, 
                               osc.treat.lab,
                               osc.eco.mat, 
                               osc.w.out, 
                               osc.w.in, osc.meth.table.rel, osc.types.sel, omic.val, h = 1, n.seq = 1, w = FALSE)
{
  
    osc.sel.meth.tab <- osc.meth.table.rel[which(osc.meth.table.rel$type == osc.types.sel ),]
    osc.sel.meth.tab <- osc.sel.meth.tab[order(as.numeric(paste(osc.sel.meth.tab$compr)), 
                                               rev(as.numeric(paste(osc.sel.meth.tab$osc.score))),
                                               decreasing = FALSE), ]
  
    omic.pa.st  <- as.numeric(treat.path[[1]]$PA)
    omic.cai.st <- as.numeric(treat.path[[1]]$CAI)
    # Single Omic  Protein Abundance Selection
    if(omic.val == "PA")
    {
      pf.j <- scale(as.numeric(treat.path[[1]][[osc.treat.lab]]))
    } else
    {
    # Multi Omic  Protein Abundance and CAI Selection
      pf.j <- scale(as.numeric(treat.path[[1]][[osc.treat.lab]]))
      pf.cai <- scale(omic.cai.st)
    
      pf.j <- scale((pf.j+pf.cai)/2)
    
    }
  
    nc <- NC(pf.j)                               #Number of Classes (Discretisation Levels) treatment
    n.q.lev          <- osc.sel.meth.tab$n.quant # standard conditions
    names(n.q.lev)   <- osc.sel.meth.tab$m.quant
  
    tab.bind <- data.frame()
    
    # Select the number of discretisation levels in common between the patterns
    # in standard conditions and after a treatment.
    for(j in 1:length(nc))
    {
       if(length(which(nc[[j]] == unique(as.numeric(n.q.lev)))) != 0)
       {
          tab.bind <- rbind(tab.bind, osc.sel.meth.tab[which(nc[[j]] == as.numeric(n.q.lev)),])
       }
     }
  
  
  osc.sel.meth.tab <- tab.bind
  
  osc.sel.meth.tab <- osc.sel.meth.tab[order(as.numeric(paste(osc.sel.meth.tab$pvalueseq)), 
                                             rev(as.numeric(paste(osc.sel.meth.tab$n.quant))),
                                             decreasing = TRUE), ]
  osc.sel.meth.tab <- osc.sel.meth.tab[which(as.numeric(paste(osc.sel.meth.tab$pvalueseq))>=max(as.numeric(paste(osc.sel.meth.tab$pvalueseq))) - 0.1),]
  
  seq.to.pst       <- osc.sel.meth.tab$seq   #ex seq.std
  seq.w.to.pst     <- osc.sel.meth.tab$compr
  
  if(length(seq.to.pst) < n.seq)
  {
    print(paste("Val max of seq: ", length(seq.to.pst)))
    
    n.seq <- length(seq.to.pst)
  }
  
  
  
  seq.to.pst       <- osc.sel.meth.tab[c(1:n.seq),]$seq   #ex seq.std
  seq.w.to.pst     <- as.numeric(osc.sel.meth.tab[c(1:n.seq),]$compr)
  
  n.q.lev          <- osc.sel.meth.tab[c(1:n.seq),]$n.quant
  names(n.q.lev)   <- osc.sel.meth.tab[c(1:n.seq),]$m.quant
  
  m.sep.lev        <- osc.sel.meth.tab[c(1:n.seq),]$m.sep
  names(m.sep.lev) <- osc.sel.meth.tab[c(1:n.seq),]$m.quant 
  
  if(n.seq != 0){
    df <- data.frame()
    mat1 <- osc.eco.mat
    
    for(i in 1:length(mat1[1,]))
    {
      mat1[i,i]<-1
      mat1<-as.matrix(mat1)
    }
    
    for(j in 1:length(n.q.lev))
    { 
      if(!is.na(n.q.lev[[j]])){
        
        class.treat <-classIntervals(pf.j, 
                                     as.numeric(paste(n.q.lev[[j]])) , 
                                     style = paste(m.sep.lev[[j]]), 
                                     rtimes = 10,
                                     intervalClosure = c("left", "right"))
        
        q.treat <-findCols(class.treat)
        seq.treat <- paste(q.treat, collapse="-")
        q.treat.c.std <-  QT(q.treat, q.treat, 0)
        seq.treat.c.std <- paste(q.treat.c.std, collapse="-")
        seq.wlev.out <- osc.select.w.all(osc.w.out, NULL)
        seq.wlev.in  <- osc.select.w.all(osc.w.in,  NULL)
        q.treat.std.out <-  QT(q.treat, as.numeric(seq.wlev.out), 4)
        q.treat.std.in  <-  QT(q.treat, as.numeric(seq.wlev.in),  4)
        
        q.treat.cc<-q.treat
        na <-99
        
        for(k in 1:length(mat1[1,]))
        {
          
          if(length(mat1[mat1[k,]!=0,])!=length(treat.path[[1]][,1]))
          {
            m<-mat1[mat1[k,]!=0,]
            m<-m[!duplicated(m==0),]
            q.treat.cc[m!=0]<-na
            na<-na+1
          }
        }
        
        q.treat.c.max <-  QT(q.treat.cc, q.treat, 0)
        q.treat.c.min <-  QT(q.treat.cc, q.treat, 1)
        q.treat.c.med <-  QT(q.treat.cc, q.treat, 2)
        
        q.treat.cc.out <-  QT(q.treat.cc, as.numeric(seq.wlev.out), 3)
        q.treat.cc.in  <-  QT(q.treat.cc, as.numeric(seq.wlev.in),  3)
        
        seq.treat.c.max <- paste(q.treat.c.max, collapse="-")
        seq.treat.c.min <- paste(q.treat.c.min, collapse="-")
        seq.treat.c.med <- paste(q.treat.c.med, collapse="-")
        
        df <- rbind(df, data.frame(seq.std   = seq.to.pst[j],
                                   seq.treat = seq.treat, 
                                   std.comp  = seq.treat.c.std, 
                                   min.comp  = seq.treat.c.min,
                                   max.comp  = seq.treat.c.max,
                                   med.comp  = seq.treat.c.med))
      }
      
    }
    
    #  Remember that if on one hand the predictions on the PSTs are done on
    #  not compressed patterns, on the other hand the PSTs could be trained 
    #  also by compressed patterns. Obviously is a project choise and for this
    #  reason there are some parts of the source commented.

    osc.pst.pred.out <- data.frame()
    #print("df$seq.std:")
    #print(df$seq.std)
    
    for(j in 1:length(df$seq.std))
    {
      if(w==FALSE)
      {
        pred.standard <- PSTpred(seq.to.pst[c(1:n.seq)], df$seq.treat[[j]] ,  h, NULL, seq.w.to.pst[c(1:n.seq)], length(unique(q.treat)))
        
      } else
      {
        
        
        pred.standard <- PSTpred(seq.to.pst[c(1:n.seq)], df$seq.treat[[j]],   h, as.numeric(c(seq.wlev.out)), seq.w.to.pst[c(1:n.seq)], length(unique(q.treat)))}
      #                  pred.comp.std <- PSTpred(seq.to.pst[c(1:n.seq)], df$std.comp[[j]],    h, q.treat.std.out, seq.w.to.pst[c(1:n.seq)])
      #                  pred.comp.max <- PSTpred(seq.to.pst[c(1:n.seq)], df$max.comp[[j]],    h, q.treat.cc.out,  seq.w.to.pst[c(1:n.seq)])
      #                  pred.comp.min <- PSTpred(seq.to.pst[c(1:n.seq)], df$min.comp[[j]],    h, q.treat.cc.out,  seq.w.to.pst[c(1:n.seq)])
      #                  pred.comp.med <- PSTpred(seq.to.pst[c(1:n.seq)], df$med.comp[[j]],    h, q.treat.cc.out,  seq.w.to.pst[c(1:n.seq)])
      
      #print("pred standard")
      #print(pred.standard)
      
      
      pred.standard.dash <- paste(pred.standard, collapse = "z")
      #                  pred.comp.std.dash <- paste(pred.comp.std, collapse = "z")
      #                  pred.comp.max.dash <- paste(pred.comp.max, collapse = "z")
      #                  pred.comp.min.dash <- paste(pred.comp.min, collapse = "z")
      #                  pred.comp.med.dash <- paste(pred.comp.med, collapse = "z")
      
      
      osc.pst.pred.out <- rbind(osc.pst.pred.out, data.frame(pred.standard = pred.standard.dash,
                                                             # pred.comp.std = pred.comp.std.dash,
                                                             # pred.comp.max = pred.comp.max.dash,
                                                             # pred.comp.min = pred.comp.min.dash,
                                                             # pred.comp.med = pred.comp.med.dash,
                                                             l.s =length(pred.standard)
                                                             # l.c.s =length(pred.comp.std),
                                                             # l.c.mx =length(pred.comp.max),
                                                             # l.c.mi = length(pred.comp.min),
                                                             # l.c.me = length(pred.comp.med)
      )
      )
      
      
      
      
    }
    
    osc.pst.pred.in <- data.frame()
    
    if(w!=FALSE){    
      for(j in 1:length(df$seq.std))
      {
        if(w==FALSE)
        {
          pred.standard <- PSTpred(seq.to.pst[c(1:n.seq)], df$seq.treat[[j]] ,  h, NULL, seq.w.to.pst[c(1:n.seq)], length(unique(q.treat)))
          
        } else
        {
          pred.standard <- PSTpred(seq.to.pst[c(1:n.seq)], df$seq.treat[[j]] ,  h, as.numeric(c(seq.wlev.in)), seq.w.to.pst[c(1:n.seq)], length(unique(q.treat))) }
        #                pred.comp.std <- PSTpred(seq.to.pst[c(1:n.seq)], df$std.comp[[j]],    h, q.treat.std.in, seq.w.to.pst[c(1:n.seq)])
        #                pred.comp.max <- PSTpred(seq.to.pst[c(1:n.seq)], df$max.comp[[j]],    h, q.treat.cc.in,  seq.w.to.pst[c(1:n.seq)])
        #                pred.comp.min <- PSTpred(seq.to.pst[c(1:n.seq)], df$min.comp[[j]],    h, q.treat.cc.in,  seq.w.to.pst[c(1:n.seq)])
        #                pred.comp.med <- PSTpred(seq.to.pst[c(1:n.seq)], df$med.comp[[j]],    h, q.treat.cc.in,  seq.w.to.pst[c(1:n.seq)])
        
        #print("pred standard")
        #print(pred.standard)              
        
        
        pred.standard.dash <- paste(pred.standard, collapse = "z")
        #                pred.comp.std.dash <- paste(pred.comp.std, collapse = "z")
        #                pred.comp.max.dash <- paste(pred.comp.max, collapse = "z")
        #                pred.comp.min.dash <- paste(pred.comp.min, collapse = "z")
        #                pred.comp.med.dash <- paste(pred.comp.med, collapse = "z")
        
        
        osc.pst.pred.in <- rbind(osc.pst.pred.in, data.frame(pred.standard = pred.standard.dash,
                                                             # pred.comp.std = pred.comp.std.dash,
                                                             # pred.comp.max = pred.comp.max.dash,
                                                             # pred.comp.min = pred.comp.min.dash,
                                                             # pred.comp.med = pred.comp.med.dash,
                                                             l.s    = length(pred.standard)
                                                             # l.c.s  = length(pred.comp.std),
                                                             # l.c.mx = length(pred.comp.max),
                                                             # l.c.mi = length(pred.comp.min),
                                                             # l.c.me = length(pred.comp.med)
        )
        )
      }}
    
    
    if(w==FALSE)
    {
      #print("osc.pst.pred.out:")
      #print(osc.pst.pred.out)
      return(list(osc.pst.pred = osc.pst.pred.out))
    }
    return(list(osc.pst.pred.out = osc.pst.pred.out, osc.pst.pred.in = osc.pst.pred.in))}
  else
  {
    return(list(osc.pst.pred.out = data.frame(), osc.pst.pred.in = data.frame()))
  }
  
}



# Number of discretisation levels for
# a specific discretised pattern.

NC <- function(x) {
  require(KernSmooth)
  require(grDevices)
  h <- dpih(x)
  bins <- seq(min(x)-h, max(x)+h, by=h)
  
  c("Sturges26" = nclass.Sturges(x),
    "Scott92" = nclass.scott(x), 
    "FreedmanDiaconis81" = nclass.FD(x),
    "Wand05"=length(bins)
  ) 
  
}


# 3 functions for the different patterns compressions.
# is an utility function

comp.group <- function(data, output, group)
{ 
  
  if(output[group] == 0)
  {return(NULL)} 
  
  if(group != 1)
  {if(output[group-1] ==1)
  { return(NULL) }
  }
  
  for(i in group:(length(data)-1))
  {
    
    if(output[i] == 1)
    { 
      if(output[i+1] == 0)
      {
        group <- unique(c(group, i , i+1))
        #print(group)
        break
      }
      
      if(output[i+1] != 0)
      {
        group <- unique(c(group, i, i+1))
        #print(group)
        
      }
    }
  }
  return(group)
}

# is an utility function
comp.opt <- function(data)
{
  output <- integer(length(data))
  for(i in 2:length(output))
  {
    if(data[[i-1L]]==data[[i]])
      output[[i-1L]] <- 1L
  }
  return(output)
}

# Paragraph 2.3 and 2.7
# takes in input a discretised omic pattern or a pattern of interactome weights, and basically 
# returns a different type of compression relative to a specific code (0,1,2,3,4) selection:
# 0 - Operons or Operons AND Protein Polymer compression max.
# 1 - Operons or Operons AND Protein Polymer compression min.
# 2 - Operons or Operons AND Protein Polymer compression med.
# 3 -If is given a pattern of interactome weights, with respect the omic pattern, 
#    is done the compression at the average value, in the case of adjacent omic 
#    values on the pattern equal.
# 4 - Standard Compression.

QT <- function(ccat,ccat1,n)
{
  red<-c()
  
  if(length(ccat1[ccat>99])==0){
    
    red<-c(ccat[1])
    for(i in 2:length(ccat))
    {
      if(ccat[(i-1)] != ccat[i])
      {
        red<-c(red,ccat[i])    
      }
    }
    
    if( n == 4)
    { 
      
      output<-comp.opt(ccat)
      for(j in 1:(length(ccat)-1))
      {
        gg1<- comp.group(ccat, output, j)
        if(!is.null(gg1) || length(gg1)!=0){
          ccat1[gg1]<- sum(ccat1[gg1])/length(ccat1[gg1])
        }
      }
      
      red<-c(ccat1[1])
      for(i in 2:length(ccat))
      {
        if(ccat[(i-1)] != ccat[i])
        {
          red<-c(red,ccat1[i])    
        }
      }
      
      
    }
    
  }
  else
  {
    
    st1<-ccat[ccat>99]
    st1<-st1[!duplicated(st1)]
    
    
    
    for(f in 1:length(st1))
    { 
      st2<-ccat1[grep(st1[f], ccat)]
      if(n==0)
      {
        
        ccat[grep(st1[f],ccat)[1]]<-max(st2)
        #  print(ccat)
        #  print(f)
        ccat<-ccat[-grep(st1[f],ccat)]
        #  print(ccat)
        #  print(max(st2))
        
      }
      
      if(n==1)
      {
        
        ccat[grep(st1[f],ccat)[1]]<-min(st2)
        
        
        
        ccat<-ccat[-grep(st1[f],ccat)]
        
        
      }
      if(n==2)
      {
        
        
        ccat[grep(st1[f],ccat)[1]]<-round(sum(st2)/length(st2))
        ccat<-ccat[-grep(st1[f],ccat)]
        
        
        #
      }
      
      if(n==3)
      {
        ccat1[grep(st1[f],ccat)[1]]<-sum(st2)/length(st2)
        ccat1<-ccat1[-grep(st1[f],ccat)[-c(1)]]
        
      }
      
      
    }
    if(n!=3){
      red<-ccat}
    else
    {red <- ccat1}
  }
  
  
  
  
  
  return(red)
  
}

#Utility Method fot the compression
SeqId<-function(nc,ccat)
{
  if(ccat[1]>= nc/2){
    x<-paste(rep(nc:1, length(ccat)), collapse="-")
  }
  else
  {
    x<-paste(rep(1:nc, length(ccat)), collapse="-")
  }
  
  ideale<-substr(x, 1, length(ccat)*2-1)
  
  return(ideale)
}




# Probabilistic suffic trees for the prediction. Paragraph 2.6
# This function returns '1 - likelihood'
PSTpred<-function(seq.std, seq.treat, h, seq.w=NULL, seq.w.c=NULL, A=1){

    require(PST)
  
  if(!is.null(seq.w.c))
  {
    seq.w.c <- 1-seq.w.c
    seq.std <- seqdef(seq.std, weights = seq.w.c)
    
    
    out <- tryCatch( S1 <- pstree(seq.std, L = h, weighted=TRUE, ymin=0.001) ,
                     error = function(e) {
                       print("The prediction in these conditions is not possible.")
                       return(NULL)})
    
    if(is.null(out))
    {return(NULL)}
    
  } else
  {
    seq.std <- seqdef(seq.std)
    
    out<- tryCatch( S1 <- pstree(seq.std, L = h, ymin=0.001)  ,
                    error = function(e) {
                      print("The prediction in these conditions is not possible.")
                      return(NULL)})
    if(is.null(out))
    {return(NULL)}
  }
  
  A <- h-1
  
  #print(paste("ecco la A:",A))
  
  C99 <- qchisq(0.99,h-1)/2
  tryCatch( S1 <- prune(S1, gain="G2", C=C99, delete=FALSE)  ,
            error = function(e) {
              print("The prediction in these conditions is not possible.")
              return(NULL)})
  
  if(is.null(out))
  {return(NULL)}
  
  seq.treat<-seqdef(seq.treat)
  
  
  
  tryCatch( pp1<-predict(S1, seq.treat, decomp=TRUE, p1=1, L=NULL) ,
            error = function(e) {
              print("The prediction in these conditions is not possible.")
              return(NULL)})
  if(is.null(out))
  {return(NULL)}
  
  
  pp1<-c(pp1)
  
  
  if(!is.null(seq.w))
  { 
    if(length(pp1) != length(seq.w))
    {pp1<- 1-pp1 #could happen if we use kmeans or other methods non deterministic
    }
    else
    {pp1<- 1-(pp1 ^ seq.w)}
    
  }
  
  return(1-pp1)
  
}

osc.pst.split.likeh <- function(str)
{
  str    <- paste(str)
  str.om <- suppressWarnings(as.numeric(strsplit(str, "z")[[1]]))
  
  return(str.om)
  
}

#  2 functions for the study of the hierarchical clusters and printing.

osc.get.diss.mat <- function(seq.osc.l)
{
  df <- data.frame()
  osc.t <- c(0)
  
  for(j in 1:length(seq.osc.l))
  {
    if(!is.null(seq.osc.l[[j]]$pred.standard[1])){
      if( seq.osc.l[[j]]$l.s[1] !=0){
        osc.t <- osc.pst.split.likeh(seq.osc.l[[j]]$pred.standard[1])}}
    
    for(i in 2:length(seq.osc.l[[j]]$pred.standard))
    {
      if(!is.null(seq.osc.l[[j]]$pred.standard[i]) ){
        if(seq.osc.l[[j]]$l.s[1]!=0){
          osc.t <- osc.t + osc.pst.split.likeh(seq.osc.l[[j]]$pred.standard[i])}}
    }
    
    if(length(osc.t)!=1){
      osc.t <- osc.t /length(seq.osc.l[[j]]$pred.standard)
      #print(df)
      df <- rbind(df, osc.t)
      names(df) <- c(1:length(osc.t))
    }
    
    
  }
  return(df)  
  
}


osc.print.clusters<- function(seq.osc.likh, treatlabtarget, kegg.code)
{
  for(j in 1:length(seq.osc.likh)){
    if(seq.osc.likh[[j]][[1]][[1]][1,2] != 0) {
      seq.osc.l.in  <- seq.osc.likh[[j]][[1]]
      seq.osc.l.out <- seq.osc.likh[[j]][[2]]
      type.training <- paste(kegg.code, seq.osc.likh[[j]][[3]], sep = "  ")
      
      
      df.out   <- osc.get.diss.mat(seq.osc.l.out)
      df.in    <- osc.get.diss.mat(seq.osc.l.in )
      
      diss.out <- daisy(df.out, metric = c("euclidean", "manhattan", "gower"),
                        stand = FALSE, type = list())
      diss.in  <- daisy(df.in, metric = c("euclidean", "manhattan", "gower"),
                        stand = FALSE, type = list())
      
      label    <- "outdeg"
      png(file=paste(type.training,label,".png", collapse=""), width=1480, height=1240, res=300)
      
      plot(hclust(diss.out, method="complete"), main=paste(type.training, label, sep=" "), labels=treatlabtarget )
      dev.off()
      
      label    <- "indeg"
      png(file=paste(type.training, label,".png", collapse=""), width=1480, height=1240, res=300)
      
      plot(hclust(diss.in, method="complete"), main=paste(type.training, label, sep=" "), labels=treatlabtarget )
      dev.off()
      
      diss.out <- daisy(df.out, metric = c("euclidean", "manhattan", "gower"),
                        stand = TRUE, type = list())
      diss.in  <- daisy(df.in, metric = c("euclidean", "manhattan", "gower"),
                        stand = TRUE, type = list())
      
      label    <- "outdegS"
      png(file=paste(type.training,label,".png", collapse=""), width=1480, height=1240, res=300)
      
      plot(hclust(diss.out, method="complete"), main=paste(type.training, label, sep=" "), labels=treatlabtarget )
      dev.off()
      
      label    <- "indegS"
      png(file=paste(type.training, label,".png", collapse=""), width=1480, height=1240, res=300)
      
      plot(hclust(diss.in, method="complete"), main=paste(type.training, label, sep=" "), labels=treatlabtarget )
      dev.off()
      
    }
    
  }
  
}
