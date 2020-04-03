# Bulk-RNA-seq all in one, and use edgeR analyze the data!
# 
BiocManager::install("DESeq2")
BiocManager::install("Rsubread")

library(DESeq2)
library(edgeR)
library(Rsubread)
library(RColorBrewer)
library(scatterplot3d)
library(pheatmap)
library("BiocParallel")
register(MulticoreParam(30))
#
# remove(list = ls())
getwd()
setwd("/Users/kerencheng/Desktop")
setwd("/home-new/ssm497/My_epi_data_part2/RNA-seq")
setwd("/Volumes/Backup/UTSA-experiment/Project-ID4GFP_SSC/Bulk.RNA-seq.comb/")
# save.image("Bulk.comb_edger.Rdata")
load("Bulk.comb_edger.Rdata")
###
# use build-in annotation
bam.0 = list.files(path = "/home-new/ssm497/My_epi_data_part2/RNA-seq", 
                   pattern= "^P6_.*\\d.bam$",
                   full.names = TRUE) # 
bam.1 = list.files(path = "/home-new/ssm497/My_epi_data_part2/BMC.Kubo", 
                   pattern= "*.BAM$",
                   full.names = TRUE) # 
bam.2 = list.files(path = "/home-new/ssm497/My_epi_data_part2/GSE102954.ATAC", 
                   pattern= "*.BAM$",
                   full.names = TRUE) # 
bam.3 = list.files(path = "/home-new/ssm497/My_epi_data_part2/GSE49624.Hammoud.CSC", 
                   pattern= "*.BAM$",
                   full.names = TRUE) # 
bam.4 = list.files(path = "/home-new/ssm497/My_epi_data_part2/GSE62355.Hammoud.GD", 
                   pattern= "*.BAM$",
                   full.names = TRUE) # 
bam.5 = list.files(path = "/home-new/ssm497/My_epi_data_part2/RNA-seq", 
                   pattern= "^P8_.*\\d.bam$",
                   full.names = TRUE) # 

# ~!___________________________________________
# use the same gtf with cellRanger?

Spg.fc <- featureCounts(files =c(bam.0, bam.1, bam.2, bam.3, bam.4, bam.5) ,
                            annot.inbuilt="mm10",
                            useMetaFeatures = TRUE,
                            isPairedEnd = TRUE,
                            nthreads = 30)

dim(Spg.fc$counts) # 27179 * 317 cells
saveRDS(Spg.fc, "Spg.fc.rds")
########################################################################################
# clean the colnames
# ,[\s\S]*$ or 
# ,.*$ to match everything after the first comma (see explanation for which one to use); or
# [^,]*$ to match everything after the last comma (which is probably what you want).
# x = colnames(E_MTAB_2600.counts$counts)
# y = gsub( "\\..*", "", x) # replace anything after first dot with nothing.把第一个dot后面的所有字符替换为空字符。

#
colnames(Spg.fc$counts)
countTable <- data.frame(Spg.fc$counts)
colnames(countTable)
# countTable[,13] = NULL  # remove Pro
colnames(countTable) = gsub( "\\.fastq.subjunc.BAM*", "", colnames(countTable)) # clean colnames
# Count matrix input
colData <- data.frame(bamnames = colnames(countTable), condition = c("P6.Id4.Bright", "P6.Id4.Bright","P6.Id4.Bright",
                                                                     "P6.Id4.Dim", "P6.Id4.Dim", "P6.Id4.Dim",
                                                                     "P6.Tspan8.high","P6.Tspan8.high","P6.Tspan8.high",
                                                                     "P6.Tspan8.low","P6.Tspan8.low","P6.Tspan8.low",
                                                                     "P0.5_Oct4", "P7.5_Kit.neg", "P7.5_Kit.pos",
                                                                     "P0.5_Oct4","P7.5_Kit.neg", "P7.5_Kit.pos",
                                                                     "P7.Thy1", "P7.Thy1", "P7.Kit", "P7.Kit",
                                                                     "Adult.Kit", "Adult.Thy1",
                                                                     "P0.Thy1", "P0.Thy1", "P7.Thy1", "P7.Thy1",
                                                                     "P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC",
                                                                     "P12.Thy1", "P12.Thy1","P12.Kit", "P12.Kit", "P14.Thy1", "P14.Kit",
                                                                     "P0.Oct4", "P7.Oct4", "P7.Id4", "P7.VasaThy1", "P7.VasaThy1", "P7.VasaThy1",
                                                                     "P8.Id4.Bright","P8.Id4.Bright","P8.Id4.Bright",
                                                                     "P8.Id4.Dim", "P8.Id4.Dim", "P8.Id4.Dim"
))


#
my.colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
              "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe", 
              "#469990", "#e6beff", "#9A6324", "#fffac8", "#800000", 
              "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", 
              "#FFF5EE", "#696969", "#000000")
my.colors <- colors[as.numeric(pca.plot$condition)]
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### heatmap ? heatmap of top 100 variable genes?
# use the log transform on the data set
# rld <- rlog(dds, blind=F)
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T), 1000)
matrix.hm <- assay(rld)[topVarianceGenes, ]
matrix.hm 
match("Pou5f1",rownames(matrix.hm))
#
quantile(matrix.hm, probs = seq(0, 1, 0.01))
q01 = quantile(matrix.hm, 0.01)
q99 = quantile(matrix.hm, 0.99)
matrix.hm[matrix.hm < q01] = q01
matrix.hm[matrix.hm > q99] = q99
#
df.anno.row <- as.data.frame(colData(dds)[,"condition"])
rownames(df.anno.row) = colnames(assay(rld))
colnames(df.anno.row) = "condition"
x = levels(df.anno.row$condition)
sort(x)
ann_colors = list( condition = c(Adult.Kit = "#e6194B",  
                                 Adult.Thy1 = "#3cb44b", 
                                 P0.5_Oct4 = "#ffe119", 
                                 P0.Oct4 ="#4363d8", 
                                 P0.Thy1 = "#f58231",
                                 P12.Kit = "#911eb4", 
                                 P12.Thy1 = "#42d4f4",
                                 P14.Kit = "#f032e6", 
                                 P14.Thy1 = "#bfef45", 
                                 P6.Id4.Bright = "#fabebe", 
                                 P6.Id4.Dim = "#469990", 
                                 P6.Tspan8.high = "#e6beff", 
                                 P6.Tspan8.low = "#9A6324", 
                                 P7.5_Kit.neg = "#fffac8",
                                 P7.5_Kit.pos = "#800000",
                                 P7.Id4 = "#aaffc3",
                                 P7.Kit = "#808000",
                                 P7.Kit_CSC = "#ffd8b1",
                                 P7.Oct4 = "#000075",
                                 P7.Thy1 = "#a9a9a9",
                                 P7.VasaThy1 = "#FFF5EE", 
                                 P8.Id4.Bright = "#696969",
                                 P8.Id4.Dim = "#000000"))


# select the 'contrast' you want
pdf("Bulk.comb.heatmap.pdf", width = 9, height = 9, useDingbats = F)
pheatmap(matrix.hm, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_colors = ann_colors,
         annotation_col = df.anno.row )
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### correlation coefficient of RNA-seq datasets ? 
library(corrplot)
library(RColorBrewer)

###########################################################################################################################################################################
# reanalyze using edgeR
library(edgeR)
library(org.Mm.eg.db)
ls()
#
colnames(Spg.fc$counts)
head(rownames(Spg.fc$counts))
group = c("P6.Id4.Bright", "P6.Id4.Bright","P6.Id4.Bright",
          "P6.Id4.Dim", "P6.Id4.Dim", "P6.Id4.Dim",
          "P6.Tspan8.high","P6.Tspan8.high","P6.Tspan8.high",
          "P6.Tspan8.low","P6.Tspan8.low","P6.Tspan8.low", "E16.5",
          "P0.5_Oct4", "P7.5_Kit.neg", "P7.5_Kit.pos",
          "P0.5_Oct4","P7.5_Kit.neg", "P7.5_Kit.pos",
          "P7.Thy1", "P7.Thy1", "P7.Kit", "P7.Kit",
          "Adult.Kit", "Adult.Thy1",
          "P0.Thy1", "P0.Thy1", "P7.Thy1", "P7.Thy1",
          "P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC","P7.Kit_CSC",
          "P12.Thy1", "P12.Thy1","P12.Kit", "P12.Kit", "P14.Thy1", "P14.Kit",
          "P0.Oct4", "P7.Oct4", "P7.Id4", "P7.VasaThy1", "P7.VasaThy1", "P7.VasaThy1",
          "P8.Id4.Bright","P8.Id4.Bright","P8.Id4.Bright",
          "P8.Id4.Dim", "P8.Id4.Dim", "P8.Id4.Dim"
)
length(group)
group = as.factor(group)
edger.count = Spg.fc$counts
dim(edger.count)
# edger.count = edger.count[,-13]  # remove E16.5_prospermatogonia 
# 
Spg.DEG.list <- DGEList(edger.count, group=group)
Spg.DEG.list$samples
# change the ENTREZID to official gene symbol using org.MM.eg.db
require(org.Mm.eg.db)
SYMBOL <- mapIds(org.Mm.eg.db, keys = rownames(Spg.DEG.list), keytype="ENTREZID", column="SYMBOL")
Spg.DEG.list$genes <- data.frame(ENTREZID=rownames(Spg.DEG.list), SYMBOL = SYMBOL, Length = Spg.fc$annotation$Length) # get the gene length from the featureCount!!

# step 3 : Filtering and normalization
keep <- filterByExpr(Spg.DEG.list)
summary(keep)
#
Spg.DEG.list <- Spg.DEG.list[keep, , keep.lib.sizes=FALSE]
Spg.DEG.list <- calcNormFactors(Spg.DEG.list) # this is TMM normalization
Spg.DEG.list$samples
# Exploring differences between libraries
pdf("Spg.DEG.plotMD_1.pdf")
plotMD(cpm(Spg.DEG.list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

# step 4 :Data exploration
points <- c(0:23)
my.colors = rep(brewer.pal(8, "Dark2"), 3)
#
pdf("Spg.DEG.plotMDS_1.pdf")
plotMDS(Spg.DEG.list,pch=points[group] , col=my.colors[group])
legend("topleft", legend=levels(group), pch=points, col=my.colors, ncol=2, cex=0.8)
dev.off()
# Spg.DEG.list
pdf("Spg.DEG.plotMDS_2.pdf")
plotMDS(Spg.DEG.list, pch=points[group], col=my.colors[group])
dev.off()

# step 5 : The design matrix : to determine your control and treament?
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

# step 6 Estimating the dispersion
Spg.DEG.list <- estimateDisp(Spg.DEG.list, design, robust=TRUE)
Spg.DEG.list$common.dispersion
# The vertical axis of the plotBCV plot shows square-root dispersion, 
# also known as biological coefficient of variation (BCV)
pdf("Spg.DEG.plotBCV_1.pdf") 
plotBCV(Spg.DEG.list)
dev.off()
#
fit <- glmQLFit(Spg.DEG.list, design, robust=TRUE)
head(fit$coefficients)
summary(fit$df.prior)
#
pdf("Spg.DEG.QLDisp_1.pdf")
plotQLDisp(fit)
dev.off()

# step 7 : Differential expression
P6.Id4.BD.con <- makeContrasts(P6.Id4.Bright - P6.Id4.Dim, levels=design)  # the left is treatment group, right one is control.
# We will use QL F-tests instead of the more usual likelihood ratio tests (LRT) 
# as they give stricter error rate control by accounting for the uncertainty in dispersion estimation:
P6.Id4.res <- glmQLFTest(fit, contrast = P6.Id4.BD.con)
P6.Id4.de <- decideTestsDGE(P6.Id4.res)
summary(P6.Id4.de)
topTags(P6.Id4.res)
#
pdf("Spg.DEG.plotMD.P6.Id4.BD.pdf")
plotMD(P6.Id4.res, status=P6.Id4.de, values=c(1,-1), col=c("red","blue"), legend="topright")
dev.off()
# set the log2 fold change to 1.5
P6.Id4.res.2 <- glmTreat(fit, contrast = P6.Id4.BD.con, lfc=log2(1.5))
str(P6.Id4.res.2$table)
topTags(P6.Id4.res.2)
summary(decideTestsDGE(P6.Id4.res.2))
P6.Id4.DEGs.results = P6.Id4.res.2$table
P6.Id4.DEGs.results$SYMBOL = rownames(P6.Id4.DEGs.results)
#
P6.Id4.up.edgr = dplyr::filter(P6.Id4.DEGs.results, logFC > 1.5 & PValue < 0.05)
# P6.Id4.PGCM.results$ENTREZID = rownames(P6.Id4.PGCM.results)
# P6.Id4.PGCM.results = P6.Id4.PGCM.results %>% left_join(Spg.DEG.list$gene , by = "ENTREZID")
# DE with a FC significantly above 1.5 at an FDR cut-off of 5%.
# make a one to one contrast. for jake. we stop at the GO analysis.
###########################################################################################################################################################################
# make a heatmap for the genes expression.
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("AnnotationDbi")
library(pheatmap)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)

# claculate the fpkm now or cpm using edgeR cpm() function.
Spg.DEG.cpm <- cpm(Spg.DEG.list, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
# rownames(Spg.DEG.cpm) <- Spg.DEG.list$genes$SYMBOL 
quantile(Spg.DEG.cpm, probs = c(0.05, 0.95)) # check the extremum.
q95 = quantile(Spg.DEG.cpm, 0.95)
q05 = quantile(Spg.DEG.cpm, 0.05)
Spg.DEG.cpm[Spg.DEG.cpm > q95] =  q95
Spg.DEG.cpm[Spg.DEG.cpm < q05] = q05

######################## 
# plot
pdf("Spg.DEG.HM.cpm_1.pdf", width = 8.5, height = 11)
pheatmap(Spg.DEG.cpm, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         scale = "row", 
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, main = NA,
)
dev.off()

#
pdf("Spg.DEG.HM.cpm_2.pdf", width = 8.5, height = 11)
hm.2 = pheatmap(Spg.DEG.cpm, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, main = NA,
)
dev.off()
##################################
# fpkm
# colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-") # some just has one replicates
# genes.x = c("Id4",  "Gfra1", "Tcl1","Zbtb16","Pou5f1","Stra8")

# gene.length.df = data.frame(ENTREZID = rownames(Spg.DEG.list)) %>%  left_join(select(gene.legnth.df.2, c("length", "ENTREZID")), by = "ENTREZID")
Spg.DEG.fpkm <- rpkm(Spg.DEG.list, normalized.lib.sizes = FALSE, log = TRUE, prior.count = 1)
Spg.DEG.fpkm[1:10, 1:10]
# remove NA
# Spg.DEG.fpkm = na.omit(Spg.DEG.fpkm)
# dim(Spg.DEG.fpkm)

pdf("Spg.DEG.HM.fpkm_2.pdf", width = 8.5, height = 15 )
pheatmap(Spg.DEG.fpkm, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         scale = "row", 
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, main = NA,
)
dev.off()
#
gene.a = c("Dazl", "Ddx4", "Zbtb16", "Pou5f1", "Nanog", "Sall4", "Foxo1", "Lin28a", "Klf4", "Sox2", "Glis1","Utf1", 
           "Id4","Gfra1", "Etv5", "Tcl1", "Foxc2", "Eomes", "Tspan8", "Pax7", "Thy1", "Pdx1", "Bmi1", 
           "Dmrt1","Kit", "Sox3", "Neurog3", "Sohlh1", "Stra8", "Sycp3", "Dmrtb1")
Spg.DEG.fpkm.subset = as.data.frame(Spg.DEG.fpkm)
Spg.DEG.fpkm.subset$symbol = mapIds(org.Mm.eg.db, keys = rownames(Spg.DEG.fpkm.subset), keytype="ENTREZID", column="SYMBOL")
Spg.DEG.fpkm.subset = dplyr::filter(Spg.DEG.fpkm.subset, symbol %in% gene.a)
Spg.DEG.fpkm.subset = column_to_rownames(Spg.DEG.fpkm.subset, var = "symbol") %>% as.matrix()
#
pdf("Spg.HM.fpkm.gene.a_1.pdf", width = 17, height = 25, useDingbats = FALSE )
pheatmap(t(Spg.DEG.fpkm.subset), 
         color = colorRampPalette(rev(brewer.pal(n = 5, name ="RdBu")))(20),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         display_numbers = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, main = NA,
)
dev.off()

#### just compare P6 Id4 Bright with Dim using edger

ls()
head(colnames(Spg.fc$counts))

P6.Id4.mex = Spg.fc$counts[,1:6]
saveRDS(P6.Id4.mex, "P6.Id4.mex.rds")



