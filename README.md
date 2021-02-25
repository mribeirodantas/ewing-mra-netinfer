# Gene Regulatory Networks Generation using RTN R/Bioconductor package.

RNA-seq data were processed according to the Tuxedo protocol. Briefly, we used fastqdump to convert SRA compressed files into FASTQ format files and fastqc to assess the quality of them. Then we run TopHat to align reads to the hg19 reference transcriptome.	For the comparison of gene expression of the master regulators among solid pediatric tumors and healthy tissues, the datasets GSE34620, GSE16476, GSE53224, GSE75271, GSE87437, GSE29684, GSE66533, and GSE3526 were normalized together. Looking for samples on the same platform type was important here since different microarray platforms have a different set of probes. With the list of differentially expressed genes (required for the Master Regulator Analysis), we proceeded with the analyses described in the codes in this reposiiotry.

Data used for this analysis are publicly available in Sequence Read Archive (SRA-NCBI) and/or Gene Expression Omnibus (GEO-NCBI).

Links to data specified on the scripts: 
* [GSE34620](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34620)
* [GSE17618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17618)
* [GSE63157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63157)
* [GSE73610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73610)
* [GSE67073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67073)
* [GSE16476](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16476)
* [GSE53224](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53224)
* [GSE75271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75271)
* [GSE87437](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87437)
* [GSE29684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29684)
* [GSE66533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66533)
* [GSE35263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35263)

For package installing, there is two sources: The Comprehensive R Archive Network - [CRAN](https://cran.r-project.org/) - and [Bioconductor](https://bioconductor.org). To install CRAN packages, on the R console, just type:

```{r}
install.packages("package")
```

For Bioconductor packages, the first installation requires the package BiocManager (from CRAN). After installing `BiocManager`, you can install Bioconductor packages using the following command on R console:

```{r}
BiocManager::install("package")
```
The parameter `package` can be substituted by an character vector (`c()`) containing all packages to install. 

R Packages used in this analysis:
1. From CRAN
* BiocManager
* data.table
* classInt 
* RColorBrewer 
* ggplot2 
* dplyr
* stringr
* FactoMineR

2. From Bioconductor
* affy
* limma
* biomaRt
* Fletcher2013b 
* RTN
* RedeR 
* enrichR
* edgeR
* DESeq2

More details will be added here.

