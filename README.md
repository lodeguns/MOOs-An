# MOOs-An
##**Multi Omic Oscillations Patterns Analyses **

The procedure requires to download all the files and run the script 
**MOOsEX.R** while all functions are gathered in MOOsFun.R

All the data descripted, in R format, are provided for the organism
Escherichia coli K-12 MG1655, KEGG id: "eco" and the pathways considered
in order as are listed in the files are: 

[1]  "00010" "00020" "00030" "00040" "00051" "00052" "00053" "00061" 
     "00071" "00130" "00190" "00220" "00230" "00240" "00250"
     
[16] "00260" "00261" "00270" "00280" "00290" "00300" "00310" "00330"
     "00340" "00350" "00360" "00362" "00380" "00400" "00410"
     
[31] "00440" "00450" "00480" "00500" "00520" "00521" "00540" "00550" 
     "00561" "00564" "00620" "00630" "00640" "00650" "00660"
     
[46] "00670" "00680" "00730" "00740" "00750" "00760" "00770" "00780"
     "00790" "00860" "00900" "00910" "00920" "02010" "02020"
     
[61] "02030" "02040" "02060" "03010" "03018" "03030" "03060" "03070"
     "03410" "03420" "03430" "03440" "04122"

These data are calculated with the algorithm NBA and
are related to the indegree and outdegree weights of the
metabolic networks considered.

> osc.list.w.in.Rdata

> osc.list.w.out.Rdata

These data are the binary matrices that describes for each
metabolic network the operons at a genetic level and protein polymer
at the proteic level.

> osc.eco.poly.ms.Rdata

> osc.eco.operon.ms.Rdata

This is a list of metabolic networks where for each metabolic network
is furnished the CAI, Protein Abundance, Gene Sequences, Gene Symbol,
Protein Name, Gene Length, 5'-3' (+)  or 3'-5' (-)

> osc.paths.list.Rdata

For each pathway considered are present all the possible SOPs and MOPs.

Operons AND Protein Polymers Compressions
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


