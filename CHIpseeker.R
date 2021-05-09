BiocManager::install("ChIPseeker")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("ReactomePA")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- readPeakFile("bed_file")

#Annotation

peakAnno <- annotatePeak(peak, tssRegion=c(-5000, 5000),TxDb=txdb, annoDb="org.Mm.eg.db")
annotation <- as.data.frame(peakAnno)
write.table(annotation,file="annotation.txt",sep="\t")

#plots
png("pie.png", width = 1500,  height = 1500, res = 300)
plotAnnoPie(peakAnno)
dev.off()

png(file="venn.png", width = 1500,  height = 1500, res = 300)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

#TSS regions plot
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

png(file="TSS_coverage.png")
plotAvgProf(tagMatrix, xlim=c(-5000, 5000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

#functional enrichment
library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId,organism = "mouse")
gene <- seq2gene(peak, tssRegion = c(-5000, 5000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene,organism = "mouse")
png("functional_enrichemnt.png", width = 1500,  height = 1500, res = 300)
dotplot(pathway2)
dev.off()
