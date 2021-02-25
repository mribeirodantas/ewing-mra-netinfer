# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(magrittr)
library(dplyr)
library(topGO)
library(enrichplot)

orgDBName = "org.Hs.eg.db"

setwd('/home/mribeirodantas/Dropbox/Paper Mestrado/Functional Analysis/')
# Network 1 (GSE34620) ----------------------------------------------------


# ARNT2 -------------------------------------------------------------------

ARNT2_1 <- c('HLF', 'ZNF365', 'GATA2', 'ZBTB25', 'BACH2', 'IKZF4', 'LMO4', 'ITGB8', 'SCN3A',
           'TUBB2B', 'FRRS1L', 'TCEAL7', 'TCEAL2', 'TRIM36', 'PTPN22', 'MIAT', 'FAM226B',
           'TMEM200A', 'RBFOX1', 'LRCH2', 'DNAJC6', 'FAXC', 'PAPLN', 'RPS6KA6', 'CNKSR2',
           'CXXC4', 'TMEM71', 'BICD1', 'PHF21B', 'JPH4', 'NFASC', 'KIF26B', 'RAB33A',
           'KIF5C', 'VASH2', 'BTBD11', 'IGLON5', 'CEP70', 'SERP2', 'NEDD4L', 'JADE3',
           'DMRTC1', 'SYT11', 'SPINT2', 'KRT17', 'SRGAP1', 'PTGIS', 'SKIDA1', 'C1GALT1C1L',
           'ADRB3', 'IGSF9', 'NTNG1', 'NOVA2', 'LRRC10B', 'EPOR', 'IL17D', 'TRIL',
           'IL17RD', 'RALGPS1', 'GNG2', 'NFIB', 'SESTD1', 'POC1B', 'LINC01122', 'TNFRSF4',
           'B3GALT2', 'ZSWIM5', 'PCNX2', 'SYNJ2', 'GPR173', 'ZNF670', 'JUP', 'ASPHD2',
           'PPOX', 'F7', 'SLC35F2', 'OSGIN2', 'DLL3', 'FXYD6', 'TVP23A', 'ADCY7',
           'DBH-AS1', 'LOC100506497', 'SCD5', 'HEBP2', 'LRRC49', 'TRIM35', 'LRRC6',
           'HSPA12B', 'WNT2B', 'FSD1', 'HIP1R', 'PIK3R3', 'GPR161', 'PALM3', 'MDK',
           'MICAL1', 'TMEM151B', 'ZNRF1', 'ABHD17C', 'PRR36', 'SEPT3', 'VANGL1',
           'LOC100506885', 'RHBDD2', 'RNF144A', 'FERMT3', 'ZNF33A', 'AKAP5', 'EZH2',
           'REC8', 'PIH1D2', 'C16orf45', 'MLKL', 'CACNB3', 'LPGAT1', 'ADO', 'CRIP2',
           'PRMT6', 'TRIM17', 'MYH10', 'CCDC136', 'GOLGA2P7', 'CAMK2G', 'PIPOX', 'ARL9',
           'LOC100506125', 'LY6E', 'KIF5A', 'SH2B3', 'GLMP', 'PCYOX1L', 'PIANP', 'RFPL2',
           'GOLGA2P10', 'MIPEP', 'SLC39A11', 'SNN', 'GOLGA7', 'DCXR', 'TRIM46', 'FSTL3',
           'LZTS3', 'LIMK2', 'EXOG', 'AMOTL1', 'UBXN8', 'CTIF', 'TMEM99', 'FAM241B',
           'DLG3', 'NQO2', 'LYRM1', 'GTF2H4', 'A1BG-AS1', 'TMEM59L', 'NUDT10', 'TM9SF4',
           'TMEM216', 'MAST3', 'SPRYD4', 'PTTG1IP', 'SVBP', 'TMEM41B', 'BRF2', 'CERK',
           'PLPPR2', 'ZNF821', 'HSBP1', 'PIGO', 'GPI', 'CHKA', 'TEX264', 'TMEM205',
           'SMIM30', 'ALDH16A1', 'AMPD2', 'ACADS', 'LENG1', 'SMIM15', 'PEX11G', 'NPTN',
           'GFOD2', 'COASY', 'CSK', 'CAPNS1', 'PPP6R1', 'LDHA', 'VPS52')

eg = bitr(ARNT2_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "CC", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

ARNT2_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(ARNT2_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of ARNT2 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/ARNT2_GSE34620.png')

# CREB3L1 -----------------------------------------------------------------

CREB3L1_1 <- c('AEBP1', 'PBX3', 'VDR', 'ESR1', 'SIX2', 'TFEC', 'POU3F1', 'SHOX2',
             'HOXA11', 'SATB2', 'RUNX2', 'TSPY26P', 'IBSP', 'PTX3', 'MMP13',
             'MAMDC2', 'OMD', 'SCIN', 'TNFRSF19', 'SGMS2', 'ADAMTS18', 'HEPH',
             'VCAN', 'LRRC15', 'PTPRD', 'PLA2G4A', 'HS3ST3B1', 'ZFHX4-AS1',
             'WFDC1', 'PTGFR', 'ABCC9', 'SLC16A7', 'SPATS2L', 'HS3ST3A1',
             'ALPL', 'PDGFRL', 'COL24A1', 'NRP2', 'FAM69A', 'PAPSS2', 'WIPI1',
             'PLD1', 'NT5E', 'RHOBTB1', 'OLFML2B', 'KDELR3', 'GALNT5', 'SMCO4',
             'FAP', 'PPIC', 'RBP4', 'FNDC3B', 'SLC41A2', 'SDC2', 'IGF2',
             'TENM4', 'SEC24D', 'KCTD12', 'LINC02544', 'RAB27B', 'PLA2R1',
             'ALPK2', 'ACAN', 'CAMK2D', 'ADAMTS9', 'ARSB', 'FAM114A1', 'EOGT',
             'PPFIBP2', 'ADAM19', 'FAM169A', 'CHST6', 'PCOLCE', 'ARHGAP42',
             'MELTF', 'S100A4', 'WNT5B', 'MARCKS', 'TEX9', 'DTX4', 'S100A11',
             'SNHG18', 'IFITM10', 'BMP1', 'SP7', 'DYNC1I1', 'C1QTNF5', 'COPZ2',
             'HRH1', 'TMEM2', 'SLC8A3', 'TMEM106A', 'JDP2', 'LSP1', 'CTSK',
             'ATP10A', 'SORD', 'MPZL1', 'PLOD2', 'CADM3', 'SMTN', 'ARHGAP22',
             'FKBP10', 'COL6A1', 'FAM162A', 'ITGA10', 'CSF1R', 'UNC5C',
             'KIAA1958', 'COL5A1', 'RGS3', 'GBE1', 'MRAS', 'KCNK2', 'FBXL22',
             'CDR2L', 'CBLN4', 'METRN', 'FRMD6-AS1', 'GSTM5', 'GLCE', 'S100A13',
             'MRC2', 'CRABP2', 'ERMN', 'FKBP11', 'BCAR3', 'KAZALD1', 'KCNK6',
             'FAM84B', 'RRBP1', 'TOM1L1', 'CDC42EP1', 'VASN', 'C2orf88',
             'MMP14', 'NAF1', 'PAMR1', 'CPZ', 'TNFSF4', 'FBLIM1', 'MMP23B',
             'NDRG4', 'C1QTNF6', 'ZNF469', 'S100A6', 'GPR68', 'KLHL15',
             'FCHSD1', 'TMEM116', 'FBXO28', 'ACTN1', 'CARHSP1', 'TP53I11',
             'P3H1', 'FOPNL', 'CALHM5', 'SPNS2', 'LRP1', 'IRS2', 'HTRA3',
             'ZXDB', 'HHIPL1', 'IFT74', 'CMBL', 'CHPF', 'TXNDC5', 'CFHR2',
             'RIBC1', 'NKIRAS1', 'SPIN1', 'USP27X-AS1', 'ARHGAP23', 'PLBD2',
             'OARD1', 'ARF4', 'TMEM214', 'CHPF2', 'SLC17A9', 'PPIB', 'GBF1',
             'HDAC7', 'NT5C3A', 'CD63', 'HYAL2', 'WBP4', 'TSEN34', 'PRPSAP2',
             'RTL8C', 'YIPF2', 'DCP1A', 'SLC39A13', 'KIAA1191', 'SEC31A')

eg = bitr(CREB3L1_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

CREB3L1_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(CREB3L1_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of CREB3L1 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/CREB3L1_GSE34620.png')

# GLI3 --------------------------------------------------------------------

GLI3_1 <- c('AEBP1', 'TCF7L1', 'BACH2', 'ETV5', 'COL12A1', 'FST', 'RERG',
          'IGDCC4', 'LOX', 'PDGFRA', 'PDZRN3', 'FRMD6', 'PRR16', 'PCSK5',
          'RASSF8-AS1', 'MSRB3', 'FAM69A', 'DKK3', 'TSPAN12', 'MTCL1',
          'RHOBTB1', 'RASSF8', 'KIF26B', 'FAP', 'FBN1', 'PPIC', 'ST6GALNAC3',
          'MMP2', 'ADAMTS9', 'TAGLN', 'EPS8', 'SH3RF3', 'NUAK1', 'FAM110B',
          'COL6A2', 'PCOLCE', 'CPXM1', 'MFAP2', 'GATM', 'PDGFRB', 'BEX1',
          'FAM57A', 'TFB1M', 'MICAL1', 'RNF144A', 'PLPP7', 'GOLM1', 'PTGFRN',
          'TMEM44', 'NT5DC2', 'FSCN1', 'TXNDC5', 'ALDH18A1', 'CELA1', 'KLK1',
          'CD99L2', 'COPS7B', 'LINC00244', 'SELENON', 'AKR1A1')

eg = bitr(GLI3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

GLI3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(GLI3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of GLI3 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/GLI3_GSE34620.png')

# MEF2C -------------------------------------------------------------------

MEF2C_1 <- c('NFE2L1', 'AEBP1', 'PBX3', 'MEOX2', 'CITED2', 'RUNX1T1', 'HEY1',
           'CA3', 'MAMDC2', 'CAVIN4', 'TNNI1', 'PLN', 'VCAN', 'TNNC2', 'UNC45B',
           'SRPX', 'ZNF521', 'COL3A1', 'DCLK1', 'DMD', 'CTHRC1', 'COL1A1',
           'PDGFRA', 'SORBS1', 'NEXN', 'FABP4', 'PRR16', 'MSRB3', 'GPX3',
           'ADAMTS5', 'FHL1', 'TRPC1', 'TM4SF1', 'WIPI1', 'MECOM', 'MXRA5',
           'PRUNE2', 'TM4SF18', 'PTPRB', 'TEK', 'ADGRL2', 'RHOBTB1', 'TACC2',
           'SHISA2', 'SORBS2', 'GUCY1A3', 'ZNF827', 'FAP', 'PPIC', 'RBMS3',
           'H19', 'TMTC1', 'SGCE', 'NRP1', 'PALLD', 'MMRN2', 'ITGA9', 'IL32',
           'PDK4', 'RCAN2', 'NID2', 'CNN3', 'TIE1', 'LPL', 'AFAP1L1', 'EMCN',
           'COL15A1', 'CDH5', 'FHOD3', 'LIMCH1', 'CD200', 'PECAM1', 'GIMAP6',
           'AOC3', 'LAMA4', 'MAN1A1', 'CARMN', 'APLNR', 'ZEB1', 'NID1', 'PALMD',
           'GIMAP8', 'ADGRF5', 'F13A1', 'DOCK9', 'TSPAN7', 'FRY', 'FILIP1L',
           'VWF', 'ADGRL4', 'SPRY1', 'MYL9', 'ARHGAP6', 'PPP1R3A', 'SCHIP1',
           'MMP2', 'CAMK2D', 'TMCC3', 'PDLIM1', 'ADAMTS9', 'FAM198B', 'STON2',
           'HOTS', 'WWTR1', 'RCSD1', 'RPS6KA2', 'TAGLN', 'RASSF4', 'CARD10',
           'ARHGAP29', 'INPP4B', 'PLA2G16', 'HTRA1', 'SYTL2', 'KDR', 'AQP1',
           'NUAK1', 'PTPRM', 'TNS3', 'PRKG1', 'COX6A2', 'ITGA6', 'HACD1',
           'IGFN1', 'RFTN1', 'GMPR', 'KANK1', 'ECSCR', 'DOCK4', 'C1QTNF5',
           'FRMD4A', 'ROBO4', 'SYNM', 'CD93', 'NES', 'PDGFRB', 'A2M', 'ALPK3',
           'GRB10', 'C1RL', 'SPARCL1', 'SMTN', 'SEC61B', 'ACTA2', 'PODXL',
           'ARAP3', 'STAB1', 'GYG1', 'EXOC6B', 'PLVAP', 'BVES', 'HIST1H2BH',
           'CRACR2A', 'LAMC1', 'ESAM', 'GBE1', 'SPARC', 'H2BFS', 'ADD3-AS1',
           'TMEM204', 'GSN', 'S100A13', 'HSPG2', 'DYSF', 'PXDN', 'AK1', 'LAMA5',
           'HIST1H2AM', 'SHISA4', 'TMEM70', 'CD34', 'TTC32', 'LGALS1', 'GTF2F1',
           'KIAA1671', 'DUS4L', 'NPAP1', 'ZNF570', 'CCDC124', 'NDUFA6-AS1',
           'RPRD2', 'CDK2AP1', 'SPATA2L', 'PSD', 'TRHR', 'MOAP1', 'FAM83C',
           'CCDC174', 'POLR3C', 'ZNF134', 'LOC102546294', 'SELENOW', 'PRNT',
           'TP73-AS1', 'UBA2', 'TIMM17B', 'TUBGCP2', 'ZSCAN25', 'PWP1',
           'TOPORS', 'FAM187A', 'WDR83', 'HNRNPK')

eg = bitr(MEF2C_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

MEF2C_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(MEF2C_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of MEF2C in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/MEF2C_GSE34620.png')

# PAX7 --------------------------------------------------------------------

PAX7_1 <- c('HOXC4', 'NR0B1', 'HOXD10', 'PPARG', 'HOXB13', 'ARNTL', 'ZNF302',
          'SOX5', 'CHI3L1', 'FMO3', 'CPA3', 'MED12L', 'IL7', 'IL13RA1', 'EGFR',
          'ENPP3', 'SIAH3', 'NPY5R', 'PLS1', 'IQCA1', 'MTCL1', 'MAPT-AS1',
          'MS4A2', 'PHKA1', 'CORO2B', 'PAPPA', 'CDKL2', 'TPD52L1', 'CCL11',
          'REPS2', 'AKR1C3', 'PIK3CG', 'MLC1', 'PCDH8', 'IGSF21', 'FAM105A',
          'S100A2', 'JAKMIP2', 'NPY1R', 'CNMD', 'GYG2', 'CMTM8', 'HOOK1',
          'ARTN', 'ZNF385D', 'CYP2E1', 'AMER2', 'TMEM178B', 'DPF3', 'TOX2',
          'SMIM10', 'EPB41L4A-AS2', 'WWC1', 'TUBB4A', 'RIPOR3', 'CAPRIN1',
          'AKR1C2', 'PRKAG2-AS1', 'PDS5B', 'ST8SIA5', 'CYP26B1', 'ADCYAP1',
          'C9orf66', 'PLEKHG3', 'HS3ST4', 'UBXN7-AS1', 'LCN2', 'TESC',
          'TMEM246', 'LRRC75A', 'CMA1', 'TXNRD3', 'CMTM1', 'L3MBTL3', 'ARHGEF9',
          'CHP2', 'CNN2', 'THEM6', 'FAT4', 'SLC45A3', 'H2AFY2', 'ITIH3',
          'MYORG', 'LINC01607', 'CCNJL', 'KIAA1522', 'CIB2', 'COPRS', 'ANO9',
          'ZSWIM6', 'CHD3', 'UBXN2A', 'PHF10', 'APOL4', 'BAHCC1', 'FCGRT',
          'GINM1', 'NECTIN4', 'SPINT1-AS1', 'CST4', 'PHF13', 'IL25', 'RPL39')

eg = bitr(PAX7_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PAX7_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PAX7_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PAX7 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/PAX7_GSE34620.png')

# PBX3 --------------------------------------------------------------------

PBX3_1 <- c('CEBPZ', 'MEIS2', 'MEF2C', 'BTAF1', 'CREB3L1', 'ZNF423', 'ZFHX4',
          'ZNF395', 'TMEM255A', 'PDZRN4', 'SRPX', 'AIG1', 'GULP1', 'CTSC',
          'TMEM200A', 'PDGFRA', 'TMEM133', 'RHOJ', 'DPYD', 'MSRB3', 'NRP2',
          'GHR', 'NR2F2', 'FAM69A', 'PTN', 'TM4SF1', 'DNAH9', 'SLC8A1',
          'NIPAL2', 'GPR34', 'NPW', 'ADGRL2', 'RHOBTB1', 'NHSL1', 'SMCO4',
          'FBN1', 'PPIC', 'SELENOM', 'ZC2HC1A', 'NRP1', 'MSR1', 'TIE1',
          'SLC46A3', 'PDE5A', 'CORO2B', 'VAMP8', 'SNX7', 'VEGFC', 'CD200',
          'CDK6', 'LAMA4', 'KITLG', 'CPED1', 'NHS', 'MAN1A1', 'MS4A7', 'PLA2R1',
          'PIEZO2', 'ZNF503', 'FEZ1', 'PTGIS', 'C1QB', 'PLK2', 'LOC729970',
          'DOCK11', 'SCHIP1', 'DOCK10', 'CERKL', 'C1R', 'RASSF4', 'C1QA',
          'CARD10', 'NTNG1', 'NOVA2', 'GAS6', 'SIPA1L2', 'MAPT', 'ZNF608',
          'PPFIBP2', 'UBE2E2', 'NUAK1', 'GGTA1P', 'DPYSL3', 'TNS3', 'ASPHD1',
          'PCOLCE', 'TGFBR2', 'GNG2', 'XG', 'TM6SF1', 'RFTN1', 'C1QC', 'ECSCR',
          'MGLL', 'HLA-DMA', 'NR2F1-AS1', 'C1QTNF5', 'LPP', 'MYLK', 'PCNX2',
          'COPZ2', 'KATNAL2', 'MID1', 'SLC40A1', 'PDGFRB', 'A2M', 'SPTBN1',
          'TSPAN18', 'FXYD6', 'CHN1', 'HEBP2', 'LHFPL2', 'C9orf3', 'CMKLR1',
          'GPX8', 'SCPEP1', 'SLCO2B1', 'ESAM', 'RAB40B', 'BNIP3', 'ING5',
          'SEZ6L2', 'MICAL1', 'MRC2', 'FOLR2', 'LYSMD2', 'NR2F1', 'TBK1',
          'EHD3', 'LINC02361', 'LGMN', 'MAPRE2', 'RBFA', 'OPLAH', 'HPS3',
          'ISYNA1', 'SUSD1', 'METTL8', 'LINC00667', 'TCP10L', 'RCC1L', 'FNDC11',
          'PRRT3-AS1', 'NIFK-AS1', 'GAL3ST4', 'SNHG5', 'LIMK2', 'DAP', 'ZNF703',
          'GNPTAB', 'GPAT4', 'TAX1BP3', 'PEAK1', 'TECPR1', 'TXNDC5', 'A1BG-AS1',
          'CYHR1', 'LOC101927610', 'TMEM14A', 'PPM1J', 'ROCK2', 'TMEM41B',
          'EXOSC5', 'CYB5D1', 'PPT1', 'RPL10A', 'ST13', 'NCK2', 'TMEM205',
          'PLA2G6', 'CLBA1', 'EBNA1BP2', 'YTHDF3-AS1', 'CCDC92', 'LRRC56',
          'PMPCA', 'USP35', 'MRPL24', 'SUPV3L1', 'NR2C2AP', 'OTUD5', 'HADHB')

eg = bitr(PBX3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PBX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PBX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PBX3 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/PBX3_GSE34620.png')

# RUNX3 -------------------------------------------------------------------

RUNX3_1 <- c('TFAP2B', 'ERF', 'NOTCH2', 'GTF2IRD1', 'CCDC155', 'PTPN22', 'NRK',
           'CLDN1', 'MPP7', 'SIAH3', 'CD7', 'GRIK3', 'ST3GAL6', 'DCBLD1',
           'CHPT1', 'PRICKLE2', 'TJP1', 'PYGL', 'IGSF21', 'RADIL', 'SLC17A7',
           'ARTN', 'CCDC30', 'MGC16275', 'SMIM17', 'CLEC2L', 'CADPS2', 'MPZL1',
           'SMIM10L2A', 'TOX2', 'RARRES2', 'ADCY10P1', 'WWC1', 'TRANK1', 'WHRN',
           'DNAL4', 'DUSP19', 'ASIC2', 'REEP6', 'DISC1', 'ZNF571-AS1', 'GARNL3',
           'NR2F2-AS1', 'LINC02361', 'TH', 'FABP5', 'SLC16A8', 'WDR31',
           'LINC00471', 'ZNF416', 'TRAM1L1', 'SECTM1', 'SMIM8', 'ERI2',
           'RTN4RL1', 'RNASEH1-AS1', 'CACNA1H', 'CLDN16', 'SLC43A3', 'PDZD4',
           'KLHL35', 'CCDC191', 'SWT1', 'B3GNTL1', 'DEAF1', 'LINC01431', 'UMPS',
           'C11orf74', 'KIF7', 'ZNF502', 'ABCD1', 'C2orf76', 'QSOX2', 'KLF16',
           'ATP5L', 'TXNDC15', 'ASPSCR1', 'LRRC56', 'DUSP28', 'MOB4', 'ATG13',
           'IDH3G', 'RPL35')

eg = bitr(RUNX3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

RUNX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(RUNX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of RUNX3 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 1 (34620)/RUNX3_GSE34620.png')

# Enrichment
eg_PAX7_1 = bitr(PAX7_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_RUNX3_1 = bitr(RUNX3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_PBX3_1 = bitr(PBX3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_MEF2C_1 = bitr(MEF2C_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_CREB3L1_1 = bitr(CREB3L1_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_GLI3_1 = bitr(GLI3_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_ARNT2_1 = bitr(ARNT2_1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
               OrgDb = orgDBName, drop = T)

clusters1 <- list(eg_PAX7_1$ENTREZID, eg_ARNT2_1$ENTREZID, eg_CREB3L1_1$ENTREZID,
                 eg_GLI3_1$ENTREZID, eg_RUNX3_1$ENTREZID, eg_PBX3_1$ENTREZID,
                 eg_MEF2C_1$ENTREZID)
names(clusters1) <- c('PAX7', 'ARNT2', 'CREB3L1', 'GLI3', 'RUNX3', 'PBX3',
                     'MEF2C')
CompareGO_BP = compareCluster(clusters1, fun="enrichGO", pvalueCutoff=0.05,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=10, title="GO - BP 7 MRs Network 1")

# Network 2 (GSE63157) -----------------------------------------------------------

# ARNT2 -------------------------------------------------------------------

ARNT2_2 <- c('CTBP2', 'PKNOX2', 'CAMTA1', 'DNALI1', 'RIMKLA', 'HHAT', 'VASH2',
             'GSTM3', 'KCNA3', 'PCNX2', 'EFR3B', 'DNAJC27-AS1', 'HADHB', 'TTL',
             'DPP10', 'LYPD6', 'CASP8', 'ADD2', 'SEPT10', 'QTRT2', 'GAP43',
             'PLXNB1', 'TNIK', 'PCYT1A', 'SRD5A3', 'SNCA', 'ANXA5', 'RAB3C',
             'REEP2', 'PCYOX1L', 'CPLX2', 'ACOT13', 'TUBB2B', 'TUBB2A', 'CGAS',
             'ADCYAP1R1', 'AUTS2', 'LOC105375347', 'DYNC1I1', 'CCDC136',
             'MAGI2', 'EZH2', 'PTPRN2', 'MIR153-2', 'LOC105375615',
             'LOC105379583', 'OSTF1', 'ZNF483', 'RASEF', 'FRMD3', 'TACC2',
             'GRID1', 'SORCS1', 'NAV2', 'ELP4', 'MAPK8IP1', 'MAPK8IP1P2', 'MDK',
             'NCAM1', 'TMEM25', 'TTC36', 'TECTA', 'SORL1', 'STK33', 'ATL3',
             'ATL3', 'NRXN2', 'MAP6', 'B4GALNT3', 'BTBD11', 'PIANP', 'TESC',
             'LOC105370010', 'MSI1', 'REC8', 'KIAA0391', 'SLC22A17', 'C14orf37',
             'RTN1', 'MAP3K9', 'CGNL1', 'CORO2B', 'MEX3B', 'SCNN1G', 'NDRG4',
             'JPH3', 'FOPNL', 'USP31', 'MLKL', 'BCAR1', 'CHST6', 'TTC25',
             'SLFN11', 'SLFN13', 'CACNB1', 'AXIN2', 'MTCL1', 'SETBP1', 'CDH2',
             'KLHL14', 'TMEM59L', 'SPINT2', 'CPT1C', 'SERTAD3', 'MMP24',
             'GDAP1L1', 'ADAM33', 'GGT7', 'KCNQ2', 'SLC37A1', 'TIAM1', 'GNAZ',
             'KCTD17', 'SH3BP1', 'JOSD1', 'PHF21B', 'NHS', 'GPR173', 'TCEAL2',
             'UBE2A', 'PCYT1B', 'AMOT')
eg = bitr(ARNT2_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

ARNT2_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(ARNT2_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of ARNT2 in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/ARNT2_GSE63157.png')

# CREB3L1 -----------------------------------------------------------------

CREB3L1_2 <- c('HHEX', 'RARB', 'ZNF10', 'ZNF235', 'RNF207', 'PLOD1', 'AGO4',
               'HOOK1', 'VPS45', 'NPR1', 'P3H1', 'PTGER3', 'NES', 'NRBP1',
               'FBLN7', 'ZC3H6', 'GYPC', 'PTPN18', 'LOC105373618', 'CPS1',
               'ALS2', 'CHPF', 'ARL4C', 'PLXNA1', 'IGSF10', 'PLD1', 'CPEB2',
               'ABRAXAS1', 'CAMK2D', 'FAT1', 'FBXL7', 'SEMA5A', 'MCC',
               'HIST1H2BD', 'PXDC1', 'TRAM2', 'RPS6KA2', 'CHN2', 'COL1A2',
               'PLAG1', 'CUX1', 'TNS3', 'EZH2', 'PLAT', 'ADAMTSL1', 'ADAMTSL1',
               'SECISBP2', 'NEK6', 'LOC100129034', 'SIRT1', 'H2AFY2', 'KAZALD1',
               'ADRA2A', 'FRMD4A', 'MKX', 'CAPRIN1', 'LDLRAD3', 'IGF2',
               'INS-IGF2', 'ST5', 'DKK3', 'C11orf24', 'SMCO4', 'THY1', 'NUP58',
               'GAS6', 'LINC00454', 'KCTD12', 'RASA3', 'VRK1', 'LINC00618',
               'EVL', 'PYGL', 'NID2', 'BMP4', 'SLC8A3', 'POMT2', 'GLCE', 'CA12',
               'GLIS2', 'CDYL2', 'COTL1', 'FKBP10', 'MRC2', 'NXN', 'COPZ2',
               'MIR152', 'KCTD15', 'NOTCH3', 'HSD17B14', 'DLGAP4', 'SLC2A10',
               'SRXN1', 'TMEM74B', 'JAG1', 'RRBP1', 'TGM2', 'COL6A1', 'COL6A2',
               'TIAM1', 'SMTN', 'CDC42EP1', 'ARHGAP8', 'PRR5', 'PRR5-ARHGAP8',
               'ARHGAP6')

eg = bitr(CREB3L1_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

CREB3L1_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(CREB3L1_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of CREB3L1 in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/CREB3L1_GSE63157.png')

# GLI3 --------------------------------------------------------------------

GLI3_2 <- c('CREB3L2', 'NKX2-2', 'EPHB2', 'TMEM60', 'ATP1A1', 'BOLA1', 'JAK1',
            'LINC01359', 'BCAR3', 'GSTM3', 'NOTCH2', 'NOTCH2NL', 'LOC100996763',
            'LOC101929796', 'ID2', 'LOC105373412', 'BCL2L11', 'RTN4', 'CMTM7',
            'WWTR1', 'PKD2', 'CYP2U1', 'GAB1', 'RPSAP36', 'CAMK2D', 'CXXC5',
            'DAP', 'CTD-2154B17.1', 'UGT3A2', 'ARSB', 'MCC', 'MARCH3', 'SLIT3',
            'MIR218-2', 'ALDH5A1', 'HIST1H2AC', 'C6orf48', 'SNORD52', 'LRRC1',
            'DCDC2', 'GPLD1', 'AUTS2', 'LOC105375347', 'FLNC', 'GIMAP2',
            'ARHGEF10', 'PREX2', 'CLU', 'MIR6843', 'FOCAD', 'LOC105375989',
            'TLE4', 'TGFBR1', 'MSANTD3', 'MSANTD3-TMEFF1', 'AIF1L', 'TLE1',
            'DENND1A', 'MIR601', 'LOC105376266', 'LOC105376267', 'H2AFY2',
            'DNAJC12', 'LDLRAD3', 'IGF2', 'INS-IGF2', 'RNF141', 'FXYD6',
            'FXYD6-FXYD2', 'BORCS5', 'CMAS', 'RASSF8', 'TSPAN11', 'ETFBKMT',
            'GLIPR1L2', 'C12orf66', 'IPO5', 'MTMR6', 'LOC105370121', 'DOCK9',
            'AIDA', 'GMPR2', 'HAUS4', 'GLCE', 'PAQR5', 'MPI', 'MYO5C', 'PRKCB',
            'MIR1273H', 'NETO2', 'GLG1', 'CYB5D1', 'NBR1', 'NXN', 'SLFN11',
            'SLFN13', 'SLC16A6', 'RAB40B', 'SETBP1', 'PHLPP1', 'FBXO15',
            'SMIM7', 'NNAT', 'SULF2', 'XPNPEP3', 'FBLN1', 'LARGE1', 'DNAL4',
            'MID1', 'KCND1', 'AMOT')

eg = bitr(GLI3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

GLI3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(GLI3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of GLI3 in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/GLI3_GSE63157.png')

# MEF2C -------------------------------------------------------------------

MEF2C_2 <- c('ETS1', 'RUNX1', 'MACF1', 'KIAA0754', 'NEXN', 'PALMD', 'VCAM1',
             'CD53', 'CASQ1', 'PEA15', 'FCER1G', 'FCGR2A', 'LOC105371473',
             'FCGR2C', 'FCGR2B', 'LAMC1', 'PLA2G4A', 'CFH', 'CFHR1', 'PTPRC',
             'PLPP3', 'FAM69A', 'ARHGAP29', 'LOC105378859', 'CNN3', 'CASQ2',
             'NHLH2', 'NOTCH2', 'NOTCH2NL', 'LOC100996763', 'LOC101929796',
             'S100A10', 'TPM3', 'OLFML2B', 'DPT', 'CD34', 'AKT3',
             'LOC105373255', 'LTBP1', 'RASGRP3', 'GALNT5', 'ITGA6', 'ITGAV',
             'COL3A1', 'SPATS2L', 'CDC42EP3', 'CYP1B1', 'NEB', 'FAP', 'TTN',
             'FUNDC2', 'COL5A2', 'SLC40A1', 'MYL1', 'RBMS3', 'TGFBR2', 'GYG1',
             'MME', 'IQCJ', 'SCHIP1', 'IQCJ-SCHIP1', 'FNDC3B', 'TMEM110-MUSTN1',
             'MUSTN1', 'ADAMTS9', 'EOGT', 'GBE1', 'CCDC80', 'HEG1', 'SLC9A9',
             'ST13', 'WWTR1', 'GOLIM4', 'MECOM', 'PLD1', 'TNFSF10', 'SLIT2',
             'MIR218-1', 'FAM114A1', 'LIMCH1', 'FIP1L1', 'PDGFRA', 'LNX1-AS2',
             'SGMS2', 'GAB1', 'RPSAP36', 'GUCY1B3', 'PALLD', 'EMCN', 'CAMK2D',
             'SEC24D', 'SFRP2', 'FAM198B', 'ASB5', 'PDLIM3', 'FAT1', 'RAI14',
             'CMYA5', 'VCAN', 'GRAMD2B', 'CHSY3', 'TGFBI', 'DAB2', 'MCTP1',
             'ELL2', 'LOX', 'PPIC', 'SPARC', 'HLA-DRA', 'HLA-DQA1', 'HLA-DQA1',
             'DCBLD1', 'PLN', 'MYCT1', 'F13A1', 'HLA-DPA1', 'CLIC5', 'DST',
             'COL12A1', 'FILIP1', 'BVES', 'LAMA4', 'MAN1A1', 'CTGF', 'CPED1',
             'GIMAP4', 'ELMO1', 'PODXL', 'C7orf13', 'LINC01006', 'PDGFRL',
             'PREX2', 'SULF1', 'NOV', 'COL14A1', 'LOC105375730', 'MSR1',
             'NCALD', 'ENPP2', 'TEK', 'MAMDC2', 'RPL24', 'TGFBR1', 'PTPRD',
             'OGN', 'OMD', 'ASPN', 'ANGPTL2', 'ZEB1', 'PRKG1', 'ANKRD1',
             'PPP1R3C', 'HECTD2-AS1', 'TNNT3', 'OLFML1', 'PARVA', 'SERPING1',
             'MS4A6A', 'MS4A4E', 'LOC105369319', 'SLN', 'EMP1', 'MSRB3',
             'LOC105369809', 'IRAK3', 'LYZ', 'C1R', 'A2M', 'PTPRB', 'LUM',
             'GLT8D2', 'LMO7', 'POSTN', 'DOCK9', 'AIDA', 'NID2', 'FERMT2',
             'LGMN', 'TPM1', 'ACTC1', 'FBN1', 'MYO1E', 'CDH13', 'PECAM1',
             'WIPI1', 'MIR635', 'PTPRM', 'GALNT1', 'PDCD2L', 'FPR3', 'STK35',
             'CD93', 'SMIM11B', 'SMIM11A', 'NRIP1', 'KDELR3', 'MCAT', 'CYBB',
             'IL13RA1', 'LOC101928242', 'MXRA5', 'SMPX', 'DMD')

eg = bitr(MEF2C_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

MEF2C_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(MEF2C_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of MEF2C in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/MEF2C_GSE63157.png')

# PAX7 --------------------------------------------------------------------

PAX7_2 <- c('CEBPB', 'NKX2-2', 'NR0B1', 'RUNX3', 'STAT6', 'PRDM2', 'KAZN',
            'IGSF21', 'ARTN', 'TMEM56', 'TMEM56-RWDD3', 'ATP1A1', 'JAK1',
            'LINC01359', 'GNG12', 'DENND1B', 'CERS6', 'CYP26B1', 'C1QL2',
            'SESTD1', 'IRS1', 'PXYLP1', 'TRPC1', 'GBE1', 'PARM1',
            'LOC105377280', 'UGDH', 'USP46', 'CDKL2', 'CAMK2D', 'PDE5A',
            'NPY1R', 'FAM105A', 'SNX18', 'WDR41', 'ARSB', 'SAPCD1', 'MSH5',
            'MSH5-SAPCD1', 'DSE', 'ENPP3', 'STXBP5', 'ECI2', 'FAM184A', 'FUCA2',
            'TSPAN13', 'LRRN3', 'ETV1', 'SLC25A13', 'RARRES2', 'KIAA1456',
            'ZDHHC2', 'EPHX2', 'ZNF703', 'TDRP', 'ALDH1B1', 'TTC39B', 'GAS1',
            'TRIM14', 'REEP3', 'LOC101929846', 'H2AFY2', 'ADRB1', 'MPP7',
            'CCND1', 'STK33', 'FADS1', 'ORAOV1', 'CADM1', 'TMTC2', 'ACACB',
            'PARP11', 'DUSP16', 'ITPR2', 'LOC101928554', 'TMTC1', 'SLC41A2',
            'TNFRSF19', 'AMER2', 'SIAH3', 'CNMD', 'FGF14', 'FGF14-IT1',
            'KDELC1', 'RCOR1', 'DPF3', 'BCL11B', 'SPRED1', 'PAQR5', 'MAN2A2',
            'HDDC3', 'CACNA1H', 'SCNN1G', 'PRKCB', 'MIR1273H', 'HS3ST4',
            'CARHSP1', 'ARHGAP17', 'HIRIP3', 'CDH8', 'GLG1', 'CHD3', 'SCARNA21',
            'MAPT', 'MAPT-IT1', 'BAHCC1', 'KCNAB3', 'GALNT1', 'SETBP1', 'MALT1',
            'ZNF521', 'STXBP2', 'PET100', 'FCGRT', 'PRR12', 'RCN3', 'AKAP8L',
            'SHANK1', 'SIRPA', 'SIRPB1', 'SLC24A3', 'PARD6B', 'RIPOR3', 'SEPT5',
            'SEPT5-GP1BB', 'FAM19A5', 'KCND1', 'RPS6KA6', 'PCDH18')

eg = bitr(PAX7_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PAX7_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PAX7_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PAX7 in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/PAX7_GSE63157.png')

# PBX3 --------------------------------------------------------------------

PBX3_2 <- c('MEIS2', 'MEOX2', 'NKX2-2', 'SMAD9', 'ZBTB25', 'ADGRL2', 'GSTM2',
            'ERVK11-1', 'GSTM1', 'STRIP1', 'ATP1A1', 'PEA15', 'LAMC1', 'JAK1',
            'LINC01359', 'NES', 'CRABP2', 'AKT3', 'LOC105373255', 'CLIP4',
            'VAMP8', 'KYNU', 'NRP2', 'TP53I3', 'MALL', 'LOC100288570',
            'LIMS3-LOC440895', 'PLA2R1', 'CHN1', 'RBMS3', 'CMTM7', 'CD200',
            'TRPC1', 'IQCJ', 'SCHIP1', 'IQCJ-SCHIP1', 'LMLN', 'CCR1', 'IL17RD',
            'LOC105377101', 'ADAMTS9', 'ROBO1', 'SLCO2A1', 'TM4SF1', 'NCEH1',
            'ANKRD18DP', 'SHROOM3', 'SPRY1', 'ADAMTS3', 'FAM105A', 'RAI14',
            'VCAN', 'RANBP17', 'USP12', 'ENC1', 'LHFPL2', 'MBLAC2', 'PPIC',
            'FBN2', 'LOC105379167', 'HAVCR2', 'ADAM19', 'AIG1', 'CLIC5',
            'TRAM2', 'LAMA4', 'MAN1A1', 'FUCA2', 'ZSCAN25', 'ZNF398', 'ETV1',
            'CDK6', 'LOC105375397', 'SGCE', 'DPYSL2', 'PLAT', 'TOX', 'NCALD',
            'TOP1MT', 'TLE4', 'CTSL', 'CTSLP8', 'TLR4', 'AIF1L', 'PTPRD',
            'PLIN2', 'TMEM2', 'CELF2', 'ZEB1', 'RASSF4', 'LRMDA', 'PAPSS2',
            'ENTPD7', 'ACSL5', 'NRP1', 'RHOBTB1', 'PPFIBP2', 'CCND1', 'TMEM133',
            'ARHGAP42', 'CD59', 'C11orf91', 'DLG2', 'CTSC', 'NOX4', 'EMP1',
            'MSRB3', 'LOC105369809', 'IRAK3', 'PLXNC1', 'ITPR2', 'LOC101928554',
            'PRICKLE1', 'R3HDM2', 'C12orf66', 'PHLDA1', 'LOC105369847',
            'SLC41A2', 'TNFRSF19', 'DCLK1', 'EPSTI1', 'DOCK9', 'AIDA', 'KDELC1',
            'TPM1', 'GLCE', 'MYO5C', 'RAB27A', 'SHISA6', 'FKBP10', 'MRC2',
            'WIPI1', 'MIR635', 'ZNF521', 'TYROBP', 'SAMHD1', 'NRIP1', 'ADAMTS1',
            'ADAMTS5', 'KDELR3', 'CYBB', 'MAGED1', 'HEPH', 'MAP3K15', 'SH3KBP1',
            'GPC3', 'PCDH18')

eg = bitr(PBX3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PBX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PBX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PBX3 in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/Network 2 (63157)/PBX3_GSE63157.png')

# RUNX3 -------------------------------------------------------------------

RUNX3_2 <- c('ESRRB', 'PAX7', 'TEAD1', 'TFAP2B', 'ERMAP', 'ST3GAL3', 'MIR6079',
             'ARTN', 'TRIM67', 'CASZ1', 'C1orf127', 'SLC6A9', 'NOTCH2',
             'NOTCH2NL', 'LOC100996763', 'LOC101929796', 'MINDY1', 'FAM89A',
             'WNT10A', 'ACSL3', 'LIMS2', 'SP3', 'SP3P', 'AMT', 'NICN1', 'NICN1',
             'IFT57', 'PKD2', 'SGCB', 'ROPN1L', 'LOC105369214', 'IQGAP2',
             'ADRA1B', 'FAM8A1', 'ITPR3', 'GPSM3', 'REV3L', 'GPR146', 'PTPN12',
             'LOC105375363', 'GOT1L1', 'AIFM2', 'TYSND1', 'SLC25A45', 'CPT1A',
             'WNK1', 'LOC101929432', 'TBK1', 'APAF1', 'SUPT20H', 'UGGT2',
             'TTC8', 'LOC105370615', 'DPF3', 'SPRED1', 'DUOX1', 'DUOX2', 'TJP1',
             'LOC105370743', 'DUOXA1', 'AAGAB', 'CACNA1H', 'KIFC3', 'MYH10',
             'CD7', 'SECTM1', 'INSR', 'TMEM143', 'BTBD3', 'RIPOR3', 'ZNF217',
             'GTPBP1', 'PACSIN2', 'HAUS7', 'TREX2', 'ADRB3', 'GOT1L1', 'LOC100130417')

eg = bitr(RUNX3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

RUNX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(RUNX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of RUNX3 in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()

ggsave('IMGs/Network 2 (63157)/RUNX3_GSE63157.png')

# Enrichment
eg_PAX7_2 = bitr(PAX7_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
               OrgDb = orgDBName, drop = T)
eg_RUNX3_2 = bitr(RUNX3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_PBX3_2 = bitr(PBX3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
               OrgDb = orgDBName, drop = T)
eg_MEF2C_2 = bitr(MEF2C_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)
eg_CREB3L1_2 = bitr(CREB3L1_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                  OrgDb = orgDBName, drop = T)
eg_GLI3_2 = bitr(GLI3_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
               OrgDb = orgDBName, drop = T)
eg_ARNT2_2 = bitr(ARNT2_2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = orgDBName, drop = T)

clusters2 <- list(eg_PAX7_2$ENTREZID, eg_ARNT2_2$ENTREZID, eg_CREB3L1_2$ENTREZID,
                 eg_GLI3_2$ENTREZID, eg_RUNX3_2$ENTREZID, eg_PBX3_2$ENTREZID,
                 eg_MEF2C_2$ENTREZID)
names(clusters2) <- c('PAX7', 'ARNT2', 'CREB3L1', 'GLI3', 'RUNX3', 'PBX3',
                     'MEF2C')
CompareGO_BP = compareCluster(clusters2, fun="enrichGO", pvalueCutoff=0.05,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=10, title="GO - BP 7 MRs Network 2")

# Total Network 1 ---------------------------------------------------------

network1 <- unique(c(ARNT2_1, MEF2C_1, PBX3_1, CREB3L1_1, PAX7_1, RUNX3_1, GLI3_1))

# x=compareCluster(network1, OrgDb=orgDBName)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dotplot(x, showCategory=5)

eg = bitr(network1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

net1_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(net1_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of 7 MRs in GSE34620") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/7MRs_GSE34620.png')

# Enrichment
clusters <- list(c(eg_PAX7_1$ENTREZID, eg_ARNT2_1$ENTREZID, eg_CREB3L1_1$ENTREZID,
                 eg_GLI3_1$ENTREZID, eg_RUNX3_1$ENTREZID, eg_PBX3_1$ENTREZID,
                 eg_MEF2C_1$ENTREZID),
                 c(eg_PAX7_2$ENTREZID, eg_ARNT2_2$ENTREZID, eg_CREB3L1_2$ENTREZID,
                   eg_GLI3_2$ENTREZID, eg_RUNX3_2$ENTREZID, eg_PBX3_2$ENTREZID,
                   eg_MEF2C_2$ENTREZID))
names(clusters) <- c('Network 1', 'Network 2')
CompareGO_BP = compareCluster(clusters, fun="enrichGO", pvalueCutoff=0.01,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=10, title="GO - Biological Process")

# Total Network 2 ---------------------------------------------------------------

network2 <- unique(c(ARNT2_2, MEF2C_2, PBX3_2, CREB3L1_2, PAX7_2, RUNX3_2, GLI3_2))

eg = bitr(network2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

net2_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(net2_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of 7 MRs in GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/7MRs_GSE63157.png')

# Both networks -----------------------------------------------------------

# PAX7 -----------------------------------------------------------

PAX7 <- unique(c(PAX7_1, PAX7_2))

eg = bitr(PAX7, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PAX7_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PAX7_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PAX7 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/PAX7_both_networks.png')

# RUNX3 -----------------------------------------------------------

RUNX3 <- unique(c(RUNX3_1, RUNX3_2))

eg = bitr(RUNX3, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

RUNX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(RUNX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of RUNX3 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/RUNX3_both_networks.png')

# PBX3 -----------------------------------------------------------

PBX3 <- unique(c(PBX3_1, PBX3_2))

eg = bitr(PBX3, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

PBX3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(PBX3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of PBX3 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/PBX3_both_networks.png')

# CREB3L1 -----------------------------------------------------------

CREB3L1 <- unique(c(CREB3L1_1, CREB3L1_2))

eg = bitr(CREB3L1, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

CREB3L1_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(CREB3L1_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of CREB3L1 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/CREB3L1_both_networks.png')

# ARNT2 -----------------------------------------------------------

ARNT2 <- unique(c(ARNT2_1, ARNT2_2))

eg = bitr(ARNT2, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

ARNT2_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(ARNT2_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of ARNT2 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/ARNT2_both_networks.png')

# MEF2C -----------------------------------------------------------

MEF2C <- unique(c(MEF2C_1, MEF2C_2))

eg = bitr(MEF2C, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

MEF2C_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(MEF2C_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of MEF2C in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/MEF2C_both_networks.png')

# GLI3 -----------------------------------------------------------

GLI3 <- unique(c(GLI3_1, GLI3_2))

eg = bitr(GLI3, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"), OrgDb = orgDBName, drop = T)

egoBP = groupGO(gene = unique(eg$ENTREZID), OrgDb = orgDBName, ont = "BP", readable = T, level = 2)
egoBP@result$ONTOLOGY = "Biological Process"

ego = egoBP@result
ego = ego[ego$Count > 0,]

GLI3_go <- ego %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = Description))

ggplot(GLI3_go, aes(x=Description, y=Count)) +
  ggtitle("Functional enrichment of GLI3 in GSE34620 and GSE63157") +
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
ggsave('IMGs/GLI3_both_networks.png')


# Total networks enrichment ----
# Enrichment
clusters_t <- list(unique(c(eg_PAX7_2$ENTREZID, eg_PAX7_1$ENTREIZ)),
                  unique(c(eg_ARNT2_2$ENTREZID, eg_ARNT2_1$ENTREZID)),
                  unique(c(eg_CREB3L1_2$ENTREZID, eg_CREB3L1_1$ENTREZID)),
                  unique(c(eg_GLI3_2$ENTREZID, eg_GLI3_1$ENTREZID)),
                  unique(c(eg_RUNX3_2$ENTREZID, eg_RUNX3_1$ENTREZID)),
                  unique(c(eg_PBX3_2$ENTREZID, eg_PBX3_1$ENTREZID)),
                  unique(c(eg_MEF2C_2$ENTREZID, eg_MEF2C_1$ENTREZID)))
names(clusters_t) <- c('PAX7', 'ARNT2', 'CREB3L1', 'GLI3', 'RUNX3', 'PBX3',
                       'MEF2C')
CompareGO_BP = compareCluster(clusters_t, fun="enrichGO", pvalueCutoff=0.05,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=10, title="GO - BP 7 MRs 2 Networks")


# With fake RUNX3
clusters_t <- list(unique(c(eg_PAX7_2$ENTREZID, eg_PAX7_1$ENTREIZ)),
                   unique(c(eg_PAX7_2$ENTREZID, eg_PAX7_1$ENTREIZ)),
                   unique(c(eg_ARNT2_2$ENTREZID, eg_ARNT2_1$ENTREZID)),
                   unique(c(eg_CREB3L1_2$ENTREZID, eg_CREB3L1_1$ENTREZID)),
                   unique(c(eg_GLI3_2$ENTREZID, eg_GLI3_1$ENTREZID)),
                   unique(c(eg_RUNX3_2$ENTREZID, eg_RUNX3_1$ENTREZID)),
                   unique(c(eg_PBX3_2$ENTREZID, eg_PBX3_1$ENTREZID)),
                   unique(c(eg_MEF2C_2$ENTREZID, eg_MEF2C_1$ENTREZID)))
names(clusters_t) <- c('PAX7', 'RUNX3', 'ARNT2', 'CREB3L1', 'GLI3', 'RUNX3V', 'PBX3',
                       'MEF2C')
CompareGO_BP = compareCluster(clusters_t, fun="enrichGO", pvalueCutoff=0.05,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=10, title="GO - BP 7 MRs 2 Networks")
