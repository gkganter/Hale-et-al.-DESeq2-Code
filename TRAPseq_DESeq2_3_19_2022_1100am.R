rm(list = ls())
par(mfrow=c(1,1))

#following: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-quality-assessment-by-sample-clustering-and-visualization
#stopped at variations to the standard workflow

setwd("C:/Users/cmhal/Desktop/R/Counts_files_trap")

library("BiocManager")
library("DESeq2")
library("dplyr")
library("AnnotationDbi")
library("org.Dm.eg.db")
library("ggplot2")


########STEP ONE: LOAD IN DATA###################################################
C1 <- read.table("Control-s1.counts.txt",header=TRUE, quote = "\"")
C2<- read.table("Control-s2.counts.txt",header=TRUE, quote = "\"")
C3<- read.table("Control-s3.counts.txt",header=TRUE, quote = "\"")
E1<- read.table("UV-s1.counts.txt",header=TRUE, quote = "\"")
E2<- read.table("UV-s2.counts.txt",header=TRUE, quote = "\"")
E3<- read.table("UV-s3.counts.txt",header=TRUE, quote = "\"")


data2<- bind_cols(C1,C2,C3,E1,E2,E3)
data3 <-data2[!duplicated(as.list(data2))]
all_trap_count_data <-data3[, c(1,7:12)]
colnames(all_trap_count_data)<- c("GeneID", "Control-s1", "Control-s2", "Control-s3", "UV-s1", "UV-s2", "UV-s3")
countdata<- all_trap_count_data

sampleData <-read.table("testCondition.txt", header=TRUE, sep="\t", row.names= NULL, stringsAsFactors=FALSE)
geneID <- countdata$GeneID
sampleIndex <- countdata[,(2:7)] 
countdata <- as.matrix(sampleIndex)
rownames(countdata) <- geneID

head(sampleData)
rownames(sampleData) <- sampleData$SampleID
sampleData <- sampleData[,c("condition", "batch")]
sampleData$condition <- factor(sampleData$condition)
sampleData$batch <- factor(sampleData$batch)
head(sampleData)

#CMH:Put the columns of the count data in the same order as rows names of the sampleData
countdata <- countdata[,unique(rownames(sampleData))]
all(colnames(countdata) == rownames(sampleData))#doublecheck

#Create factors for future plotting
condition<-factor(sampleData$condition)
condition





############STEP 2:Quality of assessment of Raw Data################################################




#box plots of raw read counts (log10)
png(file="Raw_read_counts_per_gene.boxplot.png")
#must add one to raw counts that are "0" before log transform or get "Inf"
boxplot(log2(countdata+0.1), main="", xlab="", ylab="Raw read counts per gene (log(2))",axes=FALSE, col = c("lightblue", "lightblue", "lightblue", "red", "red", "red") )
axis(2)
axis(1,at=c(length(samples):6),labels=colnames(countdata),las=2,cex.axis=0.8)
dev.off()




############STEP 3: MAKE THE DDS OBJECT and PREFILTER (from normalized data)#############################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=sampleData, design = ~ condition, tidy= FALSE)

#normalize data before applying pre-threshold
dds <- estimateSizeFactors(dds)

dds <-estimateDispersions(dds)


dds <- dds[ rowSums(counts(dds)) >= 1, ]


#Dispersion plot and fitting alternatives
par(mfrow=c(1,2))
plotDispEsts(dds, ylim= c(1e-05, 1e+02),xlab="mean")
plotDispEsts(dds, ylim= c(1e-05, 1e+01),xlim= c(1e+00, 1e+06),xlab = "mean of normalized counts", ylab = "dispersion", main= "Dispersion plot")
plotDispEsts(dds, xlab= "mean of normalized counts", main= "Dispersion plot")

#look at histo of data before pre-thresholding counts
par(mfrow=c(1,1))

ddsdf<- as.data.frame(counts(dds))




keep <- rowSums(ddsdf)
keep<-rowMeans(ddsdf[,4:6])

hist((log2(keep+1)), breaks=100, col="blue", border= "white",
     main= expression(Log[2]~"transformed counts per gene"),
     xlab=expression(Log[2]~" counts across all samples"), ylab= "Number of genes")
abline(v=log2(120), col="red", lwd=2, lty=2)


epsilon <- 1
hist(log10(keep2+0.1), breaks=100, col="blue", border="white", 
     main="Log2-transformed counts per gene", xlab ="log2(counts+1)", ylab="Number of genes",
     las=1, cex.axis=0.7)
abline(v=log10(120), col="red", lwd=2, lty=2)

boxplot(log2(countdata + epsilon), col=sampleData$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(Counts +1)")


#scatterplots 
#https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical.pdf
#look at page 8
plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());
}


# Plot the scatter plot for a few pairs of variables selected at random
set.seed(123) # forces the random number generator to produce fixed results
pairs(log10(countdata[,sample(ncol(countdata),6)] + 1), 
      panel=plotFun, lower.panel = NULL)


pairs(log10(countdata+1))
pairs(log10(countdata[,1:3]+1))
pairs(log10(countdata[,4:6]+1))


ggplot(data=ddsdf, 
       aes(ddsdf, group= condition, fill= cut)) +
  geom_density(adjust=1.5)
  facet_wrap(~cut)


  plot.new()
  hist(as.matrix(countdata), col = "blue", border = "white", breaks = 100)  
  

pairs(assay(rld), 
        panel=plotFun, lower.panel = NULL)  
  
    
  
#look at mean-ddsdf#look at mean-variance relationship beforehand
#Mean-variance relationship
mean.counts <- rowMeans(counts(dds))
variance.counts <- apply(counts(dds), 1, var)
plot(x=mean.counts, y=variance.counts, pch=16, cex=0.3, main="Mean-variance relationship", log="xy")
abline(a=0, b=1)


#pre-threshold
keep <- rowSums(counts(dds, normalized=TRUE) >= 20 ) >=6
table(keep)

dds <- dds[keep,]

ddsbck <- counts(dds)
ddsbck <- as.data.frame(ddsbck)
ddsbck$symbol <- mapIds(org.Dm.eg.db, 
                        keys=row.names(ddsbck), 
                        column="SYMBOL", 
                        keytype="FLYBASE",
                        multiVals="first")

ddsbck$entrez <- mapIds(org.Dm.eg.db, 
                         keys=row.names(ddsbck), 
                         column="ENTREZID", 
                         keytype="FLYBASE",
                         multiVals="first")

ddsbck$name =   mapIds(org.Dm.eg.db,
                        keys=row.names(ddsbck), 
                        column="GENENAME",
                        keytype="FLYBASE",
                        multiVals="first")


write.table(ddsbck, file = "backgroundgenes_norm_20c6samp_.txt", sep= '\t', col.names= NA)

#look at histo of data after pre-thresholding counts



keep3 <- rowSums(counts(dds))

hist(log2(keep3), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))


hist(log2(rowMeans(counts(dds)+1)), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))






############STEP 4: RUN DESEQ2######################################################################
dds <- DESeq(dds)
dds <- dds[ rowSums(counts(dds)) >= 1, ]
resultsNames(dds) # lists the coefficients
res <- results(dds, name = "condition_UV_vs_Control", alpha = 0.05)
summary(res)

resLfc <- lfcShrink(dds, coef="condition_UV_vs_Control", type = "apeglm")


rld <- rlog(dds, blind=FALSE)
rlogcounts <- data.frame(assay(rld))
rownames(rlogcounts) <- rownames(rld)

vsd <- vst(dds, blind=FALSE)

mcols(res)$description

######STEP 5:QC plots after DESeq2####################################################################

##Sample Relationships
#PCA plot
plotPCA(rld,intgroup=c("condition"), ntop=100)#based on the top 100 genes


distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <-  colData(vsd)$condition
colnames(mat) <-  colData(vsd)$sampleNO





library(stats)
#Pearson correlation plot
corr_coeff <-cor((rlogcounts), method= "pearson")
corr_coeff
as.dist(1-corr_coeff, upper = TRUE) %>%
  as.matrix %>%
  pheatmap::pheatmap(., main = "Pearson correlation")

cor.test(rlogcounts$Control.s1, rlogcounts$UV.s2)


library(Hmisc)
res4<-rcorr(as.matrix(counts(dds)))
round(res4$P, 10)

install.packages("correlation")
library(correlation)

correlation::correlation(rlogcounts,
                         include_factors = TRUE, method = "auto"
)

#heatmap of count matrix
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
library("vsn")

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[("condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


#Heatmap of the sample to sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

par(mfrow=c(1,2))
##Normalization of data
###plotMA
###the function plotMA shows the log2 fold changes attributable to a 
#given variable over the mean of normalized counts for all the samples in 
#the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.05. Points which fall out of the window are plotted as open triangles pointing either up or down.
DESeq2::plotMA(res,alpha=0.05,ylim=c(-5,5), colSig="red", ylab= Log[2]~"fold change")
#It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated
#with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
DESeq2::plotMA(resLfc, alpha=0.05, ylim=c(-5,5), colSig="red",ylab= Log[2]~"fold change") ####has to have DESeq2 in front of the plotMA command bc limma interferes



#Sample normalization Boxplots
#make boxplots of log transformed normalized and raw data side by side (DEGs will be from normalized data), this is post removal of low read counts
pdf("boxplot_norm_compare.pdf")
par(mfrow=c(1,2))# Create a 2 x 2 plotting matrix
raw_box<- boxplot(log2(counts(dds, normalized=FALSE)+1), main="Raw counts", ylab=Log[2]~"Counts",las=2, col=rep(c("lightblue","lightblue", "lightblue", "red", "red", "red")))
norm_box<- boxplot(log2(counts(dds, normalized=TRUE)+1), main="Normalized counts",ylab= Log[2]~"Counts", las=2, col=rep(c("lightblue","lightblue", "lightblue", "red", "red", "red")))
dev.off()

#Mean-variance relationship
mean.counts <- rowMeans(counts(dds))
variance.counts <- apply(counts(dds), 1, var)
plot(x=mean.counts, y=variance.counts, pch=16, cex=0.3, main="Mean-variance relationship", log="xy")
abline(a=0, b=1)

#outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#Dispersion plot and fitting alternatives
par(mfrow=c(1,1))
plotDispEsts(dds, ylim= c(1e-05, 1e+02))
plotDispEsts(dds)

par(mfrow=c(1,1))
plotDispEsts(dds, ylim= c(1e-05, 1e+02),xlab="mean")
plotDispEsts(dds, ylim= c(1e-05, 1e+01),xlim= c(1e+01, 1e+06),xlab = "mean of normalized counts", ylab = "dispersion", main= "Dispersion plot")
plotDispEsts(dds, xlab= "mean of normalized counts", main= "Dispersion plot")

W <- res$statddsW <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

#Library sizes boxplot (QUALITY ASSESSMENT)
#First, we can plot how many reads we have for each sample. 
#Whilst normalisation can account for imbalance in coverage across the samples, extreme differences may be indicative of underlying problems in the samples.

par(mfrow=c(1,2))
librarySizes <- colSums(counts(dds, normalized=FALSE))
barplot(librarySizes, 
        names=names(librarySizes), 
        las=3, ylab = "Reads/Sample",
        main="Barplot of Library Sizes")
abline(h=1.58e7, lty=2)


librarySizesnorm2<-colSums(counts(dds, normalized=TRUE))
barplot(librarySizesnorm2,
        names=names(librarySizesnorm2),
        las=3, ylab= "Reads/Sample",
        main= "Barplot of Library Sizes (Normalized)")



#########STEP 6:Make contrasts for DEGs#################################################################
ddsResults <- results(dds, contrast=c("condition", "UV", "Control" ))
ddsResults_df <-as.data.frame(ddsResults)
summary(ddsResults)



#padj <0.05 number of genes
sum(ddsResults$padj < 0.05, na.rm=TRUE)

ddsResults_df <- na.omit(ddsResults)
ddsResults_df <- as.data.frame(ddsResults_df)
ddsResults_df$symbol <- mapIds(org.Dm.eg.db, 
                         keys=row.names(ddsResults_df), 
                         column="SYMBOL", 
                         keytype="FLYBASE",
                         multiVals="first")

write.table(ddsResults_df, file = "DEG_TRAP2019_counts20andabove_3_16_22.txt", sep= '\t', col.name= NA)


summary(ddsResults_df)


#CMH:make a separate file for all DEG adjpvalue < or = 0.05

SigAll <- ddsResults_df[abs(ddsResults_df$padj) <= 0.05, ]
#excluding rows with 0 counts
SigAll <- na.omit(SigAll)
write.table(SigAll, file = "All_adjpvalue20c6samp_3_2_22.txt", sep= '\t', col.names= NA)

#add gene symbols to the table
SigAll2 <- as.data.frame(SigAll)
SigAll2$symbol <- mapIds(org.Dm.eg.db, 
                      keys=row.names(SigAll2), 
                      column="SYMBOL", 
                      keytype="FLYBASE",
                      multiVals="first")

SigAll2$entrez <- mapIds(org.Dm.eg.db, 
                      keys=row.names(SigAll2), 
                      column="ENTREZID", 
                      keytype="FLYBASE",
                      multiVals="first")

SigAll2$name =   mapIds(org.Dm.eg.db,
                     keys=row.names(SigAll2), 
                     column="GENENAME",
                     keytype="FLYBASE",
                     multiVals="first")
write.table(SigAll2, file = "SigAll2__20c6samp_GO.txt",sep="\t",col.names=NA)
SigAll2





#CMH:make a separate file for all upregulated adjpvalue < or = 0.05
SigUp <- SigAll2[(SigAll2$log2FoldChange) > 0.000, ]
write.table(SigUp, file = "UP_adjpvalue20c6samp_3_10_22.txt", sep= '\t', col.names= NA)

#CMH: make a separate file for all the downregulated adjpvalue < or = 0.05
SigDown <-SigAll2[(SigAll2$log2FoldChange) < 0.000,]
write.table(SigDown, file = "Down_adjpvalue20c6samp_3_7_22.txt", sep= '\t', col.names= NA)


##Plot DEGs

#Make Volcano plot

series_matrix_data <- read.table( "DEG_TRAP2019_counts20andabove_3_16_22.txt",sep="\t",quote="\"",header=TRUE)
row.names(series_matrix_data) <- series_matrix_data[,1] 
padj=series_matrix_data$padj
FC=series_matrix_data$FC  
log2FC=series_matrix_data$log2FoldChange

library(EnhancedVolcano)#making fancy volcano plots
#CMH:Make volcano plot with Enhanced Volcano package
EVolcano_TRAP <- EnhancedVolcano(series_matrix_data, lab = series_matrix_data$symbol, x= 'log2FoldChange', y= 'padj', ylab = bquote(~Log[10]~adjusted~italic(P)),
                                 pCutoff = 0.05, FCcutoff = 0.5, labSize = 4.0, title= "Larval Nociceptors 24 hr post UV-injury",xlim= c(-7, 5), ylim = c(0, -log10(10e-12)))

EVolcano_TRAP



#Make a heatmap of top 50 DEG genes by padj
#heatmap of top rld
#https://stackoverflow.com/questions/32040195/edit-row-and-col-names-in-pheatmap
#https://towardsdatascience.com/pheatmap-draws-pretty-heatmaps-483dab9a3cc

significant_results_sorted <- SigAll[order(SigAll$padj), ]
significant_genes_50 <- rownames(significant_results_sorted[1:50, ])


significant_genes_50

rld_counts <- assay(rld)


rld_counts_sig <- rld_counts[significant_genes_50, ]

rld_counts_sig

labs.row <- significant_results_sorted$symbol

pheatmap(rld_counts_sig,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         labels_row=labs.row,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 5,
         legend=TRUE,
         legend_labels = TRUE,
         annotation_legend = TRUE,
         height = 20)




#heatmap of top Upregulated genes
library(pheatmap)
significant_results_sorted <- SigAll2[order(SigAll2$padj), ]
significant_genes_50 <- rownames(significant_results_sorted[1:50, ])


significant_genes_50

rld_counts <- assay(rld)


rld_counts_sig <- rld_counts[significant_genes_50, ]

rld_counts_sig

labs.row <- significant_results_sorted$symbol

pheatmap(rld_counts_sig,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         labels_row=labs.row,
         border_color = NA,
         fontsize = 12,
         scale = "row",
         fontsize_row = 8,
         legend=TRUE,
         legend_labels = TRUE,
         annotation_legend = TRUE,
         main = "Top 50 DEGs",
         height = 20)

labs.row <- significant_results_sorted$symbol


#Plot individuadds#Plot individual genes for group comparison
#we can use plotCounts fxn to compare the normalized counts (log2 counts)
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

gene<-series_matrix_data$symbol

plotCounts(dds, gene = "FBgn0264753", intgroup="condition", main = "Rgk1")
plotCounts(dds, gene ="FBgn0060296", intgroup="condition", main= "pain")
plotCounts(dds, gene="FBgn0005640", intgroup="condition", main= "Eip63E")#pval:0.001586451	padj:0.050464255
plotCounts(dds, gene="FBgn0030749", intgroup="condition", main="AnxB11")#pval:0.045733323	padj:0.322142635
plotCounts(dds, gene ="FBgn0036566", intgroup="condition", main = "ClC-c")

#venndiagram examples
https://sbc.shef.ac.uk/prostate-bioinformatics/rna-seq-de.nb.html#

#Library sizes boxplot (QUALITY ASSESSMENT)
#First, we can plot how many reads we have for each sample. Whilst normalisation can account for imbalance in coverage across the samples, extreme differences may be indicative of underlying problems in the samples.
librarySizes <- colSums(counts(dds))
barplot(librarySizes, 
        names=names(librarySizes), 
        las=3, ylab = "Reads/Sample",
        main="Barplot of Library Sizes")
abline(h=1.58e7, lty=2)



## The affy library has a density plotting function
library(affy)

## Create a list of 4 colors to use which are the same used throughout this chapter 
library(scales)
myColors <- hue_pal()(2)

## Plot the log2-transformed data with a 0.1 pseudocount
plotDensity(log2(countdata(dds)+0.1), col=c(" blue", " blue", "blue", "red", "red", "red"), 
            lty=c(1:ncol(countdata)),lwd=2, xlab=expression(Log[2]~"raw counts"),
            main='Expression Distribution')

## Add a legend and vertical line
legend('topright',legend=(colnames(countdata)), lty=c(1:ncol(countdata)),
       col=c(" blue", " blue", "blue", "red", "red", "red"))
abline(v=log2(20) , lwd=1, col='gray', lty=2)



