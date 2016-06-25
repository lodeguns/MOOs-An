# MOOs-An
##**Multi Omic Oscillations Patterns**

The procedure requires to download all the files and run the script **MOOsEX.R**, all functions are gathered in **MOOsFun.R**

The data are represented in lists saved in R data format, also all the considered pathways and metabolic pathways. 

All the data are disposed **in order** as they are listed in the supplementary material and are referred to the organism **Escherichia coli K-12 MG1655**, KEGG identifier: "eco". 

These data are calculated with the algorithm NBA and are related to the interactome's weights of the metabolic pathways.

> osc.list.w.in.Rdata

> osc.list.w.out.Rdata

These data are referred to the binary matrices. that for each
metabolic network the operons ( at a genomic level ) and protein complexes (poly)
(at the proteomic level) describe their presence.

> osc.eco.poly.ms.Rdata

> osc.eco.operon.ms.Rdata

This is a list of pathways where for each one
is associated the CAI, protein abundances, gene sequences, gene symbols,
protein names, gene lengths and  the genes' 5'-3' (+)  or 3'-5' (-) .

> osc.paths.list.Rdata

For each pathway considered are present all the possible SOPs and MOPs 
discretised.

Operons AND protein complexes compressions
> osc.paths.list.cai.poly.Rdata      #SOPs of CAI

> osc.paths.list.pa.poly.Rdata       #SOPs of Protein Abundance

> osc.paths.list.caipa.poly.Rdata    #MOPs of CAI and Protein Abundance


Operons Compressions
> osc.paths.list.pa.operons.Rdata    #SOPs of CAI

> osc.paths.list.cai.operons.Rdata   #SOPs of CAI

> osc.paths.list.caipa.operons.Rdata #SOPs of CAI



For each pathway considered are present all the PFs 
as is described in the paper.
> osc.paths.list.suzukietall.Rdata  #Suzuki et All.

> osc.paths.list.norfloxacin.Rdata  #Faith et All.



##**Oscillatory patterns identification**
In the folder **score of omic patterns - toy examples** are present some tests with all the possible combinations of 3,4,5 discretisation levels of patterns of length 9. As described in the paper is computed: P_u (labelled prob in these files), the oscillation score (labelled score) and its normalisation (nrm.sc).



# Interactome weights with Network-based Adjacency - Algorithm 1

The pseudo-code of the Algorithm 1 is described in the paper.

A small example of the Algorithm 1 is linked here: https://github.com/lodeguns/NetBasAdj and in integrated in this source code.


--------------------------------------------------------------------------------------------------------
Please note that the software and the packages linked are built on R version 3.2.5 (2016-04-14)
> Platform: x86_64-w64-mingw32/x64 (64-bit)
> Running under: Windows >= 8 x64 (build 9200)

locale:
>  LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
>  LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
>  LC_TIME=English_United Kingdom.1252    

attached base packages:
>  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
> d3heatmap_0.6.1.1  cluster_2.0.4      PST_0.88           RColorBrewer_1.1-2 TraMineR_1.8-11.1 
> KernSmooth_2.23-15 classInt_0.1-23   

loaded via a namespace (and not attached):
 > Rcpp_0.12.5         DEoptimR_1.0-4      plyr_1.8.3          prabclus_2.2-6      base64enc_0.1-3    
 > class_7.3-14        tools_3.2.5         mclust_5.2          boot_1.3-18         rpart_4.1-10       
 > dendextend_1.2.0    digest_0.6.9        jsonlite_0.9.22     gtable_0.2.0        lattice_0.20-33    
 > png_0.1-7           Matrix_1.2-6        yaml_2.1.13         mvtnorm_1.0-5       gridExtra_2.2.1    
 > e1071_1.6-7         trimcluster_0.1-2   fpc_2.1-10          htmlwidgets_0.6     diptest_0.75-7     
 > stats4_3.2.5        grid_3.2.5          nnet_7.3-12         robustbase_0.92-6   data.table_1.9.6   
 > flexmix_2.3-13      survival_2.39-4     foreign_0.8-66      latticeExtra_0.6-28 Formula_1.2-1      
 > kernlab_0.9-24      whisker_0.3-2       ggplot2_2.1.0       magrittr_1.5        modeltools_0.2-21  
 > MASS_7.3-45         Hmisc_3.17-4        scales_0.4.0        htmltools_0.3.5     splines_3.2.5      
 > colorspace_1.2-6    acepack_1.3-3.3     munsell_0.4.3       chron_2.3-47       


