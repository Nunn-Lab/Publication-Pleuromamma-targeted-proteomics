[![OCE-1829318](https://img.shields.io/badge/NSF-1829318-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1829318)

# Publication-Pleuromamma-targeted-proteomics
Supporting documentation for Timmins-Schiffman et al. publication.<br/>
The documents included in this repository will allow you to recreate the figures in the paper.

Figure 1. Nonmetric multidimensional scaling plot of the incubation DDA dataset.<br/>
input file: ABACUS_copepods_metazoan_output.csv<br/>
R script: Figure 1 NMDS.R

Figure 2. Log fold change of differentially abundant DDA proteins.<br/>
input files: circadian DDA proteins.csv; DDA DAPs RK with annot copy.csv<br/>
R script: Figure 2 script.R

Figure 3. PCoA plots of targeted proteomics (in situ and incubation)<br/>
input file: file for PCoA.csv<br/>
R script: Figure 3 pcoa.R

Figure 4. Full cluster plots used to make this figure.<br/>
input files: file for PCoA.csv<br/>
R script: Figure 4 cluster.R

Other files included:<br/>
biostats.R and documentation (source package used for some of the plots)<br/>
comet.params file is the parameters file used in the Comet search of the DDA data
