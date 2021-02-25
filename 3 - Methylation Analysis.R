# https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
library(minfi)
library(IlluminaHumanMethylation27kmanifest)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
library(RColorBrewer)
library(limma)


# Loading the data --------------------------------------------------------

path <- '/home/mribeirodantas/Dropbox/Paper Mestrado/MethylationAnalysis/GSE118872/samplesheet'
targets <- read.metharray.sheet(path)
# sub("_Grn\\.idat.*", "", targets$Basename)
RGset <- read.metharray.exp(targets = targets)
RGset

# get the 27k annotation data
ann27k <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
head(ann27k)

# Quality Control ---------------------------------------------------------


# calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
RGset <- RGset[,keep]
RGset

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# Normalization -----------------------------------------------------------


# normalize the data; this results in a GenomicRatioSet object
# mSetSq <- preprocessQuantile(RGset)
# mSetSq <- preprocessIllumina(RGset, normalize = c('controls'))
mSetSq <- preprocessNoob(RGset)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(RGset)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(RGset, sampGroups=targets$Sample_Group, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))

# Data exploration --------------------------------------------------------


# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$status)])
legend("top", legend=levels(factor(targets$status)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")


# Filtering ---------------------------------------------------------------

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.05) == ncol(mSetSq)
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$status)])
legend("right", legend=levels(factor(targets$status)), text.col=pal,
       cex=0.7, bg="white")

# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))


# Probe-wise differential methylation analysis ----------------------------

# this is the factor of interest
cellType <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
individual <- factor(targets$Sample_Name)

# use the above to create a design matrix
design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- c(levels(cellType))

# fit the linear model
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(GroupA-GroupB,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2, adjust.method = 'bonferroni', p.value = 0.05))

# Store and sort by P-value
top.table <- topTable(fit2, sort.by = "P", n = Inf)

# get the table of results for the first contrast (groupA - groupB)
ann27kSub <- ann27k[match(rownames(mVals),ann27k$Name),
                      c(1:4,12:19,24:ncol(ann27k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann27k)
head(DMPs)

library(dplyr)
DMPs %>%
  filter(adj.P.Val <= 0.05) %>%
  filter(Symbol %in% c('PAX7', 'RUNX3', 'MEF2C', 'GLI3', 'CREB3L1',
                       'ARNT2', 'PBX3')) %>%
  View
# CREB3L1 and ARNT2 not there

DMPs %>%
        filter(adj.P.Val <= 0.05) %>%
        filter(Symbol %in% c('PAX7', 'RUNX3', 'MEF2C', 'GLI3', 'CREB3L1',
                             'ARNT2', 'PBX3')) %>% select(Symbol, logFC) %>%
        View

# write.table(DMPs, file="GSE118872/DMPs.csv", sep=",", row.names=FALSE)

# plot the top 4 most significantly differentially methylated CpGs
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
        plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

# Differential methylation analysis of regions ----------------------------
# This region analysis can't be done apparently for 27K arrays (in this pipeline)

# library(DMRcate)
# myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M",
#                              analysis.type = "differential", design = design,
#                              contrasts = TRUE, cont.matrix = contMatrix,
#                              coef = "GroupA - GroupB")
# str(myAnnotation)
# #endif /* NEWSTUFF */
# DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(DMRs)
# results.ranges
# # set up the grouping variables and colours
# groups <- pal[1:length(unique(targets$Sample_Group))]
# names(groups) <- levels(factor(targets$Sample_Group))
# cols <- groups[as.character(factor(targets$Sample_Group))]
# # draw the plot for the top DMR
# par(mfrow=c(1,1))
# DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols,
#          what = "Beta", arraytype = "450K", genome = "hg19")



# Gene ontology testing ---------------------------------------------------
# FROM THIS POINT ON, THE ANALYSIS IS VERY 450K-SPECIFIC SO IM NOT SURE IF WE
# CAN TRUST THE OUTPUTS
# Get the significant CpG sites at less than 5% FDR
# sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
# # First 10 significant CpGs
# sigCpGs[1:10]
# # Total number of significant CpGs at 5% FDR
# length(sigCpGs)
# # Get all the CpG sites used in the analysis to form the background
# all <- DMPs$Name
# # Total number of CpG sites tested
# length(all)
#
# library(missMethyl)
# par(mfrow=c(1,1))
# gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
#
# # Top 10 GO categories
# topGSA(gst, number=10)

# Master Regulator Analysis -----------------------------------------------
library(readr)
library(dplyr)

MRs_regulons_34620 <- read_delim('GSE118872/MRs/MRs_Regulons_34620.csv', delim = '\t')
MRs_regulons_63157 <- read_delim('GSE118872/MRs/MRs_Regulons_63157.csv', delim = '\t')
MRs_regulons_34620 <- lapply(MRs_regulons_34620, function(x) x[!is.na(x)])
MRs_regulons_63157 <- lapply(MRs_regulons_63157, function(x) x[!is.na(x)])

hMSC_signature <- read_tsv('GSE118872/MRs/gene_expdiff_2hMSC_3ES.tsv')
hNCC_signature <- read_tsv('GSE118872/MRs/gene_expdiff_2hNCC_3ES.tsv')

hMSC_signature %>%
  # filter(significant == 'yes') %>% # q-value < 0.05
  unique %>%
  select(gene) -> hMSC_signature

hNCC_signature %>%
  # filter(significant == 'yes') %>%
  unique %>%
  select(gene) -> hNCC_signature

human_genes_count_hMSC <- nrow(hMSC_signature)
human_genes_count_hNCC <- nrow(hNCC_signature)

cont_table <- data.frame(matrix(0, nrow = 2, ncol = 2),
                         row.names = c("in_regulon", "out_regulon"))
colnames(cont_table) = c("met", "not_met")

DMPs <- DMPs %>%
  # 0.05 é 7, 7, 6, 5 mas 0.01 é 7 7 7 7  :O
  filter(adj.P.Val < 0.05) %>%
  select(Symbol) %>%
  unique

tables_creation <- function(regulons, gene_count) {
  DMPs <- DMPs$Symbol
  unlist(regulons)

  cont_table["in_regulon", "met"] <- sum(DMPs %in% regulons)

  cont_table["in_regulon", "not_met"] <- sum(!is.na(regulons)) - sum(DMPs %in% regulons)

  cont_table["out_regulon", "met"] <- sum(!(DMPs %in% regulons))
  cont_table["out_regulon", "not_met"] <- abs(gene_count - sum(cont_table))

  return(cont_table)
}

#---- Significance test ----

# 34620 hMSC
tables_34620_hMSC <- lapply(MRs_regulons_34620, tables_creation,
                            gene_count = human_genes_count_hMSC)
tables <- tables_34620_hMSC

p_values <- lapply(tables, fisher.test)
p_values <- sapply(p_values, function(x) return(x$p.value))
adjusted_p <- p.adjust(p_values, method = "bonferroni")
hm_regulons_34620_hMSC <- names(adjusted_p)[adjusted_p < 0.05]

# 34620 hNCC
tables_34620_hNCC <- lapply(MRs_regulons_34620, tables_creation,
                            gene_count = human_genes_count_hNCC)
tables <- tables_34620_hNCC

p_values <- lapply(tables, fisher.test)
p_values <- sapply(p_values, function(x) return(x$p.value))
adjusted_p <- p.adjust(p_values, method = "bonferroni")
hm_regulons_34620_hNCC <- names(adjusted_p)[adjusted_p < 0.05]

# 63157 hMSC
tables_63157_hMSC <- lapply(MRs_regulons_63157, tables_creation,
                            gene_count = human_genes_count_hMSC)
tables <- tables_63157_hMSC

p_values <- lapply(tables, fisher.test)
p_values <- sapply(p_values, function(x) return(x$p.value))
adjusted_p <- p.adjust(p_values, method = "bonferroni")
hm_regulons_63157_hMSC <- names(adjusted_p)[adjusted_p < 0.05]

# 63157 hNCC
tables_63157_hNCC <- lapply(MRs_regulons_63157, tables_creation,
                            gene_count = human_genes_count_hNCC)
tables <- tables_63157_hNCC

p_values <- lapply(tables, fisher.test)
p_values <- sapply(p_values, function(x) return(x$p.value))
adjusted_p <- p.adjust(p_values, method = "bonferroni")
hm_regulons_63157_hNCC <- names(adjusted_p)[adjusted_p < 0.05]

hm_regulons_34620_hMSC
hm_regulons_34620_hNCC
hm_regulons_63157_hMSC
hm_regulons_63157_hNCC

