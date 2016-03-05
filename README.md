# MOOs-An
##**Multi Omic Oscillations Patterns**

The procedure requires to download all the files and run the script **MOOsEX.R**, all functions are gathered in **MOOsFun.R**

The data descripted are presented in lists saved in R data format, also all the pathways and metabolic networks considered. 

All the data are disposed **in order** as are listed in the supplementary material and are referred to the organism **Escherichia coli K-12 MG1655**, KEGG identifier: "eco". 

These data are calculated with the algorithm NBA and are related to the interactome's weights of the metabolic networks.

> osc.list.w.in.Rdata

> osc.list.w.out.Rdata

These data are referred to the binary matrices that describe for each
metabolic network the operons ( at a genetic level ) and protein complexes (poly)
(at the proteic level).

> osc.eco.poly.ms.Rdata

> osc.eco.operon.ms.Rdata

This is a list of metabolic networks where for each metabolic network
is associated the CAI, protein abundances, gene sequences, gene symbols,
protein names, gene lengths and  the genes' 5'-3' (+)  or 3'-5' (-) .

> osc.paths.list.Rdata

For each pathway considered are present all the possible SOPs and MOPs.

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


