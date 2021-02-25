# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#
#                        Copyright (c) 2017-2021
#            Marcel Ribeiro-Dantas <marcel.ribeiro-dantas@curie.fr>
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

# Getting intersection between DEGs  ------------------------------------------------

library(readr)
library(dplyr)

hMSC_signature <- read_tsv('gene_expdiff_2hMSC_3ES.tsv')
hNCC_signature <- read_tsv('gene_expdiff_2hNCC_3ES.tsv')

hMSC_signature %>%
  filter(significant == 'yes') %>%
  select(gene) -> hMSC_signature

hNCC_signature %>%
  filter(significant == 'yes') %>%
  select(gene) -> hNCC_signature

dplyr::intersect(hMSC_signature, hNCC_signature) -> hMSC_hNCC_sig_intersection

write_tsv(hMSC_signature, path = 'hMSC_signature.tsv')
write_tsv(hNCC_signature, path = 'hNCC_signature.tsv')
write_tsv(hMSC_hNCC_sig_intersection, path = 'hMSC_hNCC_sig_intersection.tsv')

# Quality Control of MicroArray data -------------------------------------------

source("https://bioconductor.org/biocLite.R")
biocLite(c("affy", "hgu133plus2cdf", "simpleaffy", "affyPLM", "ape"))

# Reading the Affymetrix MicroArray Data
library(affy)
library(hgu133plus2cdf)
library(ape)

CELfiles <- list.celfiles()

# CELfiles

rawdata <- ReadAffy()
n_samples <- c("ES1", "ES2", "ES3", "ES4", "hMSC")

# 1. Intensity Distribution
library("RColorBrewer")
usr.col=brewer.pal(9, "Set1")
mycols=rep(usr.col, each=1)
hist(rawdata, lty=rep(1,length(CELfiles)), col=mycols)
boxplot(rawdata,las=3,cex.axis=0.5, names=n_samples, col=mycols)
legend(11,180, cex=0.7, n_samples, lty=1, col=mycols)

# 2. Probe-intensity Images
par(mfrow = c(3,2))
image(rawdata[,])

# 3. Quality control metrics (They make use of Affymetrix control probes)
# They are used to check if the arrays hybridized correctly, and whether the sample
# quality was acceptable.
# Average Background must be similar across all arrays.
# Se der paw com XML, faz: conda install -c r r-xml=3.98_1.5
# Se pedir pra atualizar, diz n
library("simpleaffy")
dat.qc<-qc(rawdata)
# average background (arrays must have similar avbgs)
if (var(avbg(dat.qc)) <= 5) {
  print("Average background OK. Variance < 5")
} else {
  print("Average background BAD. Variance >= 5")
}

# scaling needed to achieve the same mean by the Affymetrix analysis software
# One scaling should not be more than three fold to another array
arraysscaling <- sfs(dat.qc)
if (max(arraysscaling)/min(arraysscaling) < 3) {
  print("Scaling OK.")
} else {
  print("Scaling BAD")
}
#  Variance percent present genes
var(percent.present(dat.qc))

# beta-actin GAPDH
# ratios(dat.qc)[, "actin3/actin5"] > 3
ratios(dat.qc)[, "actin3/actin5"]

# RNA Degradation Curves
library("RColorBrewer")
dat.deg=AffyRNAdeg(rawdata)
summaryAffyRNAdeg(dat.deg)
par(mar=c(5,4,4,7) + .1, xpd=T)
colourCount = 117
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
plotAffyRNAdeg(dat.deg, cols=getPalette(colourCount))
legend(11,180, cex=0.7, n_samples, lty=1, col=getPalette(colourCount))
# plot qc
plot(dat.qc)

# From now on, we won't work on rawdata anymore. Buckle your seat belt!
library("affyPLM")
dat.PLM = fitPLM(rawdata, background=F, normalize=F)
par(mfrow=c(2,2))
image(dat.PLM, type="weights", which=5)
image(dat.PLM, type="resids", which=5)
image(dat.PLM, type="weights", which=5)
image(dat.PLM, type="resids", which=5)
RLE(dat.PLM, ylim=c(-0.5,7), col=getPalette(colourCount), las=3, cex.axis=0.6, names=n_samples)
NUSE(dat.PLM, ylim=c(0.95,1.35), col=getPalette(colourCount), las=3, cex.axis=0.6, names=sampleNames(rawdata))

# PCA
color = c(replicate(4, "red"), replicate(1, "blue"))
symbol = c(replicate(4, 1), replicate(1, 2))
data.PC = prcomp(t(expression)) # por paciente
# data.PC2 = prcomp(expression, scale.=TRUE) # por gene
plot(data.PC$x,col=color, pch=symbol)
legend("topleft", legend=c("ES", "hMSC"), pch = c(1, 2))
plot(data.PC)
barplot(data.PC$sdev/sum(data.PC$sdev))

# PCA-based gene filtering
library(affy)
library(pvac)
expression_set <- rma(rawdata)
ft <- pvacFilter(rawdata)
expression_set.filtered <- expression_set[ft$aset,]
plot(density(ft$pvac[ft$nullset]),xlab="PVAC score",main="", col="gray",cex.lab=0.5,xlim=c(0,1))
lines(density(ft$pvac),col=1)
abline(v=ft$cutoff,lty=2,col="gray")

# Clustering Tree for Samples
library(ape) # devido ao plot.phylo
rawdata <- ReadAffy()
expression_set <- rma(rawdata)
expression <- exprs(expression_set)
d <- cor(2^expression, method="pearson")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=1, show.node.label=TRUE, no.margin=TRUE, cex=0.3)

# GSE34620 Network inference ----------------------------------------------

source("http://bioconductor.org/biocLite.R")
biocLite(c("affy", "RTN", "hgu133plus2.db", "Fletcher2013b"))

library(affy)
library(RTN)
library(hgu133plus2.db)
library(Fletcher2013b)

# Lendo os dados de MicroArray
rawdata <- ReadAffy()
sample_names <- sampleNames(rawdata)
# Deixando mais legivel o nome das amostras
sample_names <- gsub(".CEL.gz", "", sample_names)

# Robust Multi-Array Average expression measure
# This function converts an AffyBatch object into an ExpressionSet object using the robust
# multi-array average (RMA) expression measure.
expression_set <- rma(rawdata)
expression <- exprs(expression_set)
colnames(expression) <- sample_names

# Carregando symbols dos probe_sets
symbol_map <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(symbol_map)
genesym.probeid <- as.data.frame(symbol_map[mapped_probes])
rm(symbol_map, mapped_probes)

# Preparando os TFs e pegando algumas informaçoes dos dados levantados no trecho anterior
data(miscellaneous)
rm(ESR1bdsites, consensus, fimoESR1, fimoFOXA1, fimoGATA3, FOXA1bdsites, GATA3bdsites, metaPCNA, SPDEFbdsites, chromlen, randsites, risksites)
tfs_names <- names(tfs)
tfs_df <- subset(genesym.probeid, symbol %in% tfs_names)
tfs <- tfs_df$probe_id
names(tfs) <- tfs_df$symbol
rm(tfs_names, genesym.probeid, tfs_df)

# Criando arquivo de anotaçoes com ENTREZID
# E opcional no tni.preprocess mas e recomendado, entao coloquei
probes <- c(do.call("cbind", as.list(rownames(expression))))
gene_annotations <- select(hgu133plus2.db, probes, c("SYMBOL", "ENTREZID"))
# Removing duplicated probeids
gene_annotations <- gene_annotations[!duplicated(gene_annotations$PROBEID),]

# Comecando o RTN de fato
dataRTN <- list(gexp = expression, tfs = tfs, gene_annotations=gene_annotations)
rtni <- new("TNI", gexp=dataRTN$gexp, transcriptionFactors=dataRTN$tfs)

# parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(3, "SOCK"))
rtni <- tni.preprocess(rtni, gexpIDs=dataRTN$gene_annotations)
rtni<-tni.permutation(rtni)

rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni, eps=100)
rtni <- tni.conditional(rtni, tfs=dataRTN$tfs)

saveRDS(rtni, "rtni_permutation3_34620_2.2.RData")

tni.get(rtni, what="summary")
refnet <- tni.get(rtni, what="refnet")
tnet <- tni.get(rtni, what="tnet")
stopCluster(getOption("cluster"))
library(RedeR)
rdp <- RedPort()
calld(rdp)


library(RedeR)
rdp <- RedPort()
calld(rdp)

g_rmap <- tni.graph(rtni_34620, tnet="dpi", gtype="rmap", tfs=c("ARNT2", "PAX7"))
g_treeleaf <- tni.graph(rtni_34620, tnet="dpi", gtype = "amapDend")

addGraph(rdp, g_rmap)
relax(rdp, p8=100, p4=1, p3=30, p1=40)

# GSE63157 Network inference --------

source("https://bioconductor.org/biocLite.R")
# biocLite("oligo")
library(oligo)
library(limma)
# library(affy)
library(GEOquery)

detach("package:RTN", unload=TRUE)
library(RTN, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/RTNnewest/")

setwd('/home/mribeirodantas/Documentos/Mestrado/2018/Exon Array/GSE63157')
# gds <- getGEO("GSE63157")

celFiles <- list.celfiles()
rawData <- read.celfiles(celFiles)
eset <- rma(rawData)
write.exprs(eset,file="expr_values_by_oligo.txt", sep="\t")

# Add annotation
# This assumes you already normalized the data, and the object "eset" has the data in it (from above)
# Load annotation library
library(huex10sttranscriptcluster.db)

# Strategy is to create data frame objects and merge them together - put expression info into a data frame
my_frame <- data.frame(exprs(eset))

# Put annotation information in a data frame.  To get specific fields, use packageNameSYMBOL, where the caps part names the type of data you're after
# To get a list of available annotation information, run the packagename with () at the end, i.e. mogene20sttranscriptcluster()
Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "),
                    SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "),
                    DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "))
unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}
Annot <- unfactorize(Annot)
Annot[rownames(Annot) == '2999755',]$SYMBOL <- 'AEBP1'
# Merge data frames together (like a database table join)
all <- merge(Annot, my_frame, by="row.names", all=T)

# Write out to a file:
write.table(all,file="data.ann_2.txt",sep="\t")

library(Fletcher2013b)

# Preparando os TFs e pegando algumas informaçoes dos dados levantados no trecho anterior
data(miscellaneous)
rm(ESR1bdsites, consensus, fimoESR1, fimoFOXA1, fimoGATA3, FOXA1bdsites, GATA3bdsites, metaPCNA, SPDEFbdsites, chromlen, randsites, risksites)
tfs_names <- names(tfs)
rm(tfs)
tfs_names <- tfs_names[-which(tfs_names == 'NA')]

# Quantos dos TFs que o Mauro sugeriu não constam nos arquivos de anotação?
length(which(tfs_names %in% Annot$SYMBOL == FALSE)) # 208

# Criar lista de TFs com todos os IDs existentes no Annot
tfs <- data.frame(tfs_names)
colnames(tfs) <- 'SYMBOL'
Annot$identificador <- rownames(Annot)
tfs_a <- merge(Annot,tfs,by=c("SYMBOL"))
tfs_ids <- cbind(as.character(tfs_a$SYMBOL), tfs_a$identificador)
tfs_ids <- tfs_ids[!duplicated(tfs_ids), ]
expression_data <- as.matrix(na.omit(cbind(all[1], all[5:89])))
rownames(expression_data) <- expression_data[,1]
expression_data <- expression_data[,-1]
mode(expression_data) <- 'numeric'

tfs <- tfs_ids[,2]
names(tfs) <- tfs_ids[,1]

library(RTN)
dataRTN <- list(gexp = expression_data, tfs = tfs)
rtni <- new("TNI", gexp=dataRTN$gexp, regulatoryElements=dataRTN$tfs)

setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE63157/')
save.image(file = 'antes_preprocess.RData')
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE63157/')
load('antes_preprocess.RData')
detach("package:RTN", unload=TRUE)
library(RTN, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/RTNnewest/")

# parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(7, "SOCK"))
rtni <- tni.preprocess(rtni)
rtni<-tni.permutation(rtni)
save.image(file = "rtni_permutation_63157_2.2.RData")
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni, eps=0)
save.image(file = "rtni_permutation2_63157_2.2.RData")
rtni <- tni.conditional(rtni, tfs=dataRTN$tfs)
save.image(file = "rtni_permutation3_63157_2.2.RData")

##
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE63157/')
load('rtni_permutation3_63157_2.2.RData')
##

library(affy)
library(RTN)
library(hgu133plus2.db)
library(Fletcher2013b)

# c("CREB3L1", "AEBP1", "PBX3", "PAX7", "RUNX3", "MEF2C", "GLI3", "ARNT2") %in% colnames(RegScores$dif)
common_TFs <- c("CREB3L1", "PBX3", "PAX7", "RUNX3", "MEF2C", "GLI3", "ARNT2", "AEBP1") # da análise anterior de comuns

RegScores <- tni.gsea2(rtni)
# RegScores_log <- tni.gsea2(rtni, log=TRUE)

heatmap(RegScores$dif[, common_TFs])
# heatmap(RegScores_log$dif[, common_TFs])

# Comparing the following regulons: GSE63175 and GSE34620 -----------------
# It is required to load regulons.analysis.R

setwd("/home/mribeirodantas/Documentos/Mestrado/2018/Construção Rede Regulatória/")
load('rtni_permutation3_34620_2.2.RData')
detach("package:RTN", unload=TRUE)
library(RTN, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/RTNnewest/")
g_34620 <- tni.graph(rtni, tnet="dpi", gtype = "rmap")
rm(dataRTN, expression, expression_set, gene_annotations, rawdata, rtni, symbol_map, mapped_probes, probes, sample_names, tfs)
g_63157 <- tni.graph(rtni_63157, tnet="dpi", gtype = "rmap")

# Genes regulados pelo TF probeid
# Qual o probeid?
# PAX7
select(huex10sttranscriptcluster.db, 'PAX7', keytype='SYMBOL', c("PROBEID"))$PROBEID
genes_in_63157_reg_by <- select(huex10sttranscriptcluster.db, as.vector(neighbors(g_63157, '2323347')$nodeAlias), c("SYMBOL"))$SYMBOL
genes_in_63157_reg_by <- genes_in_63157_reg_by[!is.na(genes_in_63157_reg_by)]
genes_in_63157_reg_by <- c(genes_in_63157_reg_by_PAX7[1:7], genes_in_63157_reg_by)
genes_in_34620_reg_by <- get_regulated_by_tf(g_34620, 'PAX7')
print(100*length(Reduce(intersect, list(genes_in_63157_reg_by, genes_in_34620_reg_by)))/max(length(genes_in_63157_reg_by), length(genes_in_34620_reg_by)))
# 13.43% de genes em comum
sum(genes_in_63157_reg_by %in% genes_in_34620_reg_by)

# RUNX3
select(huex10sttranscriptcluster.db, 'RUNX3', keytype='SYMBOL', c("PROBEID"))$PROBEID
genes_in_63157_reg_by <- select(huex10sttranscriptcluster.db, as.vector(neighbors(g_63157, '2401994')$nodeAlias), c("SYMBOL"))$SYMBOL
genes_in_63157_reg_by <- genes_in_63157_reg_by[!is.na(genes_in_63157_reg_by)]
# genes_in_63157_reg_by <- c(genes_in_63157_reg_by_PAX7[1:7], genes_in_63157_reg_by)
genes_in_34620_reg_by <- get_regulated_by_tf(g_34620, 'RUNX3')
print(100*length(Reduce(intersect, list(genes_in_63157_reg_by, genes_in_34620_reg_by)))/max(length(genes_in_63157_reg_by), length(genes_in_34620_reg_by)))
# 7.22% de genes em comum
sum(genes_in_63157_reg_by %in% genes_in_34620_reg_by)


# Plotting Tree and Leaf of GSE63157 network
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE63157/')
load('rtni_permutation3_63157_2.2.RData')
g_treeleaf <- tni.graph(rtni, tnet="dpi", gtype = "amapDend")

library(RedeR)
rdp <- RedPort()
calld(rdp)
addGraph(rdp, g_treeleaf$g)
relax(rdp, p8=100, p4=1, p3=30, p1=40)


# GSE17618 Network Inference ----------------------------------------------
# Not really inference, since we will use the skeleton of GSE34620

## GSE17618
## GPL570 N=51
## Excluí manualmente os últimos 7, que são linhagens celulares.
## É importante que a normalização seja feita apenas com as amostras de pacientes
## N final = 44
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17618

##
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/Construção Rede Regulatória/')
# GSE34620
load('rtni_permutation3_34620_2.2.RData')
rtni_34620 <- rtni
rm(dataRTN, expression, expression_set, gene_annotations, rawdata, symbol_map, mapped_probes, probes, sample_names, tfs, rtni)
##

setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE17618_RAW/')

library(affy)
library(RTN)
library(hgu133plus2.db)
library(Fletcher2013b)

# Lendo os dados de MicroArray
rawdata <- ReadAffy()
sample_names <- sampleNames(rawdata)
# Deixando mais legivel o nome das amostras
sample_names <- substr(sample_names, 0, 9)

# Robust Multi-Array Average expression measure
# This function converts an AffyBatch object into an ExpressionSet object using the robust
# multi-array average (RMA) expression measure.
expression_set <- rma(rawdata)
expression <- exprs(expression_set)
colnames(expression) <- sample_names

# Carregando symbols dos probe_sets
symbol_map <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(symbol_map)
genesym.probeid <- as.data.frame(symbol_map[mapped_probes])
# rm(symbol_map, mapped_probes)

# Preparando os TFs e pegando algumas informaçoes dos dados levantados no trecho anterior
data(miscellaneous)
rm(ESR1bdsites, consensus, fimoESR1, fimoFOXA1, fimoGATA3, FOXA1bdsites, GATA3bdsites, metaPCNA, SPDEFbdsites, chromlen, randsites, risksites)
tfs_names <- names(tfs)
tfs_df <- subset(genesym.probeid, symbol %in% tfs_names)
tfs <- tfs_df$probe_id
names(tfs) <- tfs_df$symbol
rm(tfs_names, genesym.probeid, tfs_df)

# Criando arquivo de anotaçoes com ENTREZID
# E opcional no tni.preprocess mas e recomendado, entao coloquei
probes <- c(do.call("cbind", as.list(rownames(expression))))
gene_annotations <- select(hgu133plus2.db, probes, c("SYMBOL", "ENTREZID"))
# Removing duplicated probeids
gene_annotations <- gene_annotations[!duplicated(gene_annotations$PROBEID),]

# Comecando o RTN de fato
dataRTN <- list(gexp = expression, tfs = tfs, gene_annotations=gene_annotations)
rtni <- new("TNI", gexp=dataRTN$gexp, regulatoryElements=dataRTN$tfs)

# parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(7, "SOCK"))
rownames(dataRTN$gene_annotations) <- dataRTN$gene_annotations$PROBEID
rtni <- tni.preprocess(rtni, rowAnnotation=dataRTN$gene_annotations)
rtni<-tni.permutation(rtni)

# save.image(file = "rtni_permutation_17618_2.2_08_2018.RData")
# parei aqui a inferência da rede
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/artigo/GSE17618_RAW/')
load('rtni_permutation_17618_2.2_08_2018.RData')

rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)

rtni <- tni.conditional(rtni, tfs=dataRTN$tfs)
save.image(file = "rtni_final_17618_2.2_08_2018.RData")
load('rtni_final_17618_2.2_08_2018.RData')
####
common_TFs <- c("CREB3L1", "PBX3", "PAX7", "RUNX3", "MEF2C", "GLI3", "ARNT2") # da análise anterior de comuns

rtni_34620@gexp <- rtni@gexp
RegScores_inserted_in_skeleton <- tni.gsea2(rtni_34620)
# RegScores_inserted_in_skeleton_log <- tni.gsea2(rtni_34620, log=TRUE)

heatmap(RegScores_inserted_in_skeleton$dif[, common_TFs])
# heatmap(RegScores_inserted_in_skeleton_log$dif[, common_TFs])

# Plot Regulatory Units
##
detach("package:RTN", unload=TRUE)
library(RTN, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/RTNnewest/")
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/Construção Rede Regulatória/')
# GSE34620
load('rtni_permutation3_34620_2.2.RData')
rtni_34620 <- rtni
rm(dataRTN, expression, expression_set, gene_annotations, rawdata, symbol_map, mapped_probes, probes, sample_names, tfs, rtni)
##

library(RedeR)
rdp <- RedPort()
calld(rdp)

# RMAP
g_rmap_8 <- tni.graph(rtni_34620, tnet="dpi", gtype="rmap", tfs=c("ARNT2", "AEBP1", "PBX3", "RUNX3", "PAX7", "GLI3", "CREB3L1", "MEF2C"))
g_rmap_sig <- tni.graph(rtni_34620, tnet="dpi", gtype="rmap", tfs=c("ARNT2", "CREB3L1", "GLI3", "PAX7", "PBX3", "RUNX3"))
g_rmap_sig_17618 <- tni.graph(rtni_34620, tnet="dpi", gtype="rmap", tfs=c("RUNX3", "PAX7", "ARNT2"))
g_rmap_sig_63157 <- tni.graph(rtni_34620, tnet="dpi", gtype="rmap", tfs=c("CREB3L1", "GLI3", "PAX7", "PBX3", "RUNX3"))
# g_treeleaf <- tni.graph(rtni_34620, tnet="dpi", gtype = "amapDend")

# AMAP
g_amap_8 <- tni.graph(rtni_34620, tnet="dpi", gtype="amap", amapFilter = 'phyper', amapCutoff = 1, tfs=c("ARNT2", "AEBP1", "PBX3", "RUNX3", "PAX7", "GLI3", "CREB3L1", "MEF2C"))
g_amap_sig <- tni.graph(rtni_34620, tnet="dpi", gtype="amap", amapFilter = 'phyper', amapCutoff = 1, tfs=c("ARNT2", "CREB3L1", "GLI3", "PAX7", "PBX3", "RUNX3"))
g_amap_sig_17618 <- tni.graph(rtni_34620, tnet="dpi", gtype="amap", amapFilter = 'phyper', amapCutoff = 1, tfs=c("RUNX3", "PAX7", "ARNT2"))
g_amap_sig_63157 <- tni.graph(rtni_34620, tnet="dpi", gtype="amap", amapFilter = 'phyper', amapCutoff = 1, tfs=c("CREB3L1", "GLI3", "PAX7", "PBX3", "RUNX3"))

addGraph(rdp, g_amap_sig_63157)
relax(rdp, p8=100, p4=1, p3=30, p1=40)

# GSE34620/GSE63157 MRA intersection --------------------------------------

# GSE34620

library(RTN, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/RTNnewest/")

setwd('/home/mribeirodantas/Dropbox/Marcel/artigo mestrado novo 2018/resultados/RDatas/')
load('rtni_permutation3_34620_2.2.RData')
# load('rtni_permutation3_GSE17679.RData')
rm(rawdata, probes, sample_names, tfs, expression, expression_set, dataRTN)

# Pegando os hits
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/RNAseq/MRA/')
deg_output_hMSC <- read.table(file='gene_exp_2hMSC_3ES.diff', header=TRUE, sep="\t", row.names = NULL)
deg_output_hNCC <- read.table(file='gene_exp_2hNCC_3ES.diff', header=TRUE, sep="\t", row.names = NULL)
deg_output_hMSC <- c(deg_output_hMSC[3], deg_output_hMSC[10], deg_output_hMSC[14])
deg_output_hMSC <- as.data.frame(deg_output_hMSC)
deg_output_hNCC <- c(deg_output_hNCC[3], deg_output_hNCC[10], deg_output_hNCC[14])
deg_output_hNCC <- as.data.frame(deg_output_hNCC)
# deg_output <- na.omit(deg_output)

# install.packages('splitstackshape')
library(splitstackshape)
# Separando linhas com mais de um gene e repetindo valor de expressão
deg_output_hMSC <- cSplit(deg_output_hMSC, "gene", ",", "long")
deg_output_hNCC <- cSplit(deg_output_hNCC, "gene", ",", "long")

deg_output_hMSC <- deg_output_hMSC$gene[deg_output_hMSC$significant == 'yes']
deg_output_hNCC <- deg_output_hNCC$gene[deg_output_hNCC$significant == 'yes']

# Convertendo gene symbol para probe id
library(hgu133plus2.db)
library(annotate)
deg_output_hMSC <- select(hgu133plus2.db, as.vector(deg_output_hMSC), keytype="ALIAS", c("SYMBOL", "PROBEID"))
deg_output_hNCC <- select(hgu133plus2.db, as.vector(deg_output_hNCC), keytype="ALIAS", c("SYMBOL", "PROBEID"))

deg_output_hMSC <- deg_output_hMSC[(which(!is.na(deg_output_hMSC$PROBEID))), ]
deg_output_hNCC <- deg_output_hNCC[(which(!is.na(deg_output_hNCC$PROBEID))), ]

hits_hMSC <- unique(deg_output_hMSC$PROBEID)
hits_hNCC <- unique(deg_output_hNCC$PROBEID)

library(snow)
options(cluster=makeCluster(7, "SOCK"))
rtna_hMSC <- tni2tna.preprocess(object=rtni, hits = hits_hMSC)
rtna_hNCC <- tni2tna.preprocess(object=rtni, hits = hits_hNCC)

rtna_hMSC <- tna.mra(rtna_hMSC)
rtna_hNCC <- tna.mra(rtna_hNCC)

MRs_hMSC <- tna.get(rtna_hMSC, what="mra")$Regulon
MRs_hNCC <- tna.get(rtna_hNCC, what="mra")$Regulon

MRs_hMSC_hNCC <- intersect(MRs_hMSC, MRs_hNCC)

# GSE63157
setwd('/home/mribeirodantas/Dropbox/Marcel/artigo mestrado novo 2018/resultados/RDatas/')
load('rtni_permutation3_63157_2.2.RData')
rm(all, dataRTN)

# Pegando os hits
setwd('/home/mribeirodantas/Documentos/Mestrado/2018/RNAseq/MRA/')
deg_output_hMSC <- read.table(file='gene_exp_2hMSC_3ES.diff', header=TRUE, sep="\t", row.names = NULL)
deg_output_hNCC <- read.table(file='gene_exp_2hNCC_3ES.diff', header=TRUE, sep="\t", row.names = NULL)
deg_output_hMSC <- c(deg_output_hMSC[3], deg_output_hMSC[10], deg_output_hMSC[14])
deg_output_hMSC <- as.data.frame(deg_output_hMSC)
deg_output_hNCC <- c(deg_output_hNCC[3], deg_output_hNCC[10], deg_output_hNCC[14])
deg_output_hNCC <- as.data.frame(deg_output_hNCC)
# deg_output <- na.omit(deg_output)

# install.packages('splitstackshape')
library(splitstackshape)
# Separando linhas com mais de um gene e repetindo valor de expressão
deg_output_hMSC <- cSplit(deg_output_hMSC, "gene", ",", "long")
deg_output_hNCC <- cSplit(deg_output_hNCC, "gene", ",", "long")

deg_output_hMSC <- deg_output_hMSC$gene[deg_output_hMSC$significant == 'yes']
deg_output_hNCC <- deg_output_hNCC$gene[deg_output_hNCC$significant == 'yes']

# Convertendo gene symbol para probe id
library(huex10sttranscriptcluster.db)
deg_output_hMSC <- select(huex10sttranscriptcluster.db, as.vector(deg_output_hMSC), keytype='ALIAS', c("SYMBOL","PROBEID"))
deg_output_hNCC <- select(huex10sttranscriptcluster.db, as.vector(deg_output_hNCC), keytype='ALIAS', c("SYMBOL","PROBEID"))

deg_output_hMSC <- deg_output_hMSC[(which(!is.na(deg_output_hMSC$PROBEID))), ]
deg_output_hNCC <- deg_output_hNCC[(which(!is.na(deg_output_hNCC$PROBEID))), ]

hits_hMSC <- unique(deg_output_hMSC$PROBEID)
hits_hNCC <- unique(deg_output_hNCC$PROBEID)

library(snow)
options(cluster=makeCluster(7, "SOCK"))
rtna_hMSC <- tni2tna.preprocess(object=rtni, hits = hits_hMSC)
rtna_hNCC <- tni2tna.preprocess(object=rtni, hits = hits_hNCC)

rtna_hMSC <- tna.mra(rtna_hMSC)
rtna_hNCC <- tna.mra(rtna_hNCC)

MRs_hMSC <- tna.get(rtna_hMSC, what="mra")$Regulon
MRs_hNCC <- tna.get(rtna_hNCC, what="mra")$Regulon

MRs_hMSC_hNCC <- intersect(MRs_hMSC, MRs_hNCC)

