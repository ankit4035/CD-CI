library(dada2)
set.seed(100)
### Download data files from https://drive.google.com/drive/folders/1ihRgMjM97CqYEJRXS8ESIXUT7SDook2w?usp=sharing  and extract #####


#import filepath
path <- "/data1/ankit/16s/cultured/run1" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#plot and check quality distribution
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

#filter and trim sequences. Also, trim primer sequences from R1 and R2
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE, maxEE = 2, truncLen = c(249,240), trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#plot and check quality of filtered data
plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

#learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#denoise data using predicted error rates
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

dadaFs[[1]]

#merge forward-reverse and make sequence table
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

#based on merged length, only keep sequences with length 400-430 (expected target fragment length) and save the final table.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:430]
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution", xlab = "Length", ylab = "Number of reads", col = "skyblue1")
saveRDS(seqtab2, "seqtab-run1.rds")

#import filepath for second run and proceed as above
path2 <- "/data1/ankit/16s/cultured/run2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path2)

fnFs2 <- sort(list.files(path2, pattern="_R1.fastq", full.names = TRUE))
fnRs2 <- sort(list.files(path2, pattern="_R2.fastq", full.names = TRUE))
sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)

filtFs2 <- file.path(path2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path2, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names2
names(filtRs2) <- sample.names2

plotQualityProfile(fnFs2[1:10])
plotQualityProfile(fnRs2[1:10])

out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, maxN=0, rm.phix=TRUE, maxEE = 2, truncLen = c(249,240), trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out2)

plotQualityProfile(filtFs2[1:10])
plotQualityProfile(filtRs2[1:10])

errF2 <- learnErrors(filtFs2, multithread = TRUE)
errR2 <- learnErrors(filtRs2, multithread = TRUE)

plotErrors(errF2, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)

dadaFs2 <- dada(filtFs2, err=errF2, multithread = TRUE)
dadaRs2 <- dada(filtRs2, err=errR2, multithread = TRUE)

dadaFs2[[1]]

mergers2 <- mergePairs(dadaFs2, filtFs2, dadaRs2, filtRs2, verbose=TRUE)
seqtab.2 <- makeSequenceTable(mergers2)
table(nchar(getSequences(seqtab.2)))

#save the table
seqtab2.2 <- seqtab.2[,nchar(colnames(seqtab.2)) %in% 400:430]
table(nchar(getSequences(seqtab2.2)))
hist(nchar(getSequences(seqtab2.2)), main = "Histogram of merged read length distribution", xlab = "Length", ylab = "Number of reads", col = "skyblue1")
saveRDS(seqtab2.2, "seqtab-run2.rds")


# import (if needed) and merge multiple runs
#seqtab2 <- readRDS("seqtab-run1.rds")
#seqtab2.2 <- readRDS("seqtab-run2.rds")
st.all <- mergeSequenceTables(seqtab2, seqtab2.2)

# Remove chimeras
seqtab.all <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
sum(seqtab.all)/(sum(st.all))
rowSums(seqtab.all)

# Save final sequence table
saveRDS(seqtab.all, "seqtab.all.rds")

# Assign taxonomy using SILVA v132
taxSilva <- assignTaxonomy(seqtab.all, "/data1/ankit/16s/DADA2-db/silva_nr_v132_train_set.fa", multithread=TRUE, minBoot = 80, tryRC=TRUE)
# Write to disk
saveRDS(taxSilva, "taxSilva.rds")
#ADD SPECIES
speciesSilva <- addSpecies(taxSilva, "/data1/ankit/16s/DADA2-db/silva_species_assignment_v132.fa.gz", tryRC = TRUE, n = 10000)
saveRDS(speciesSilva, "speciesSilva.rds")
write.table(speciesSilva, file = "speciesSilva.tsv", sep = "\t")

#import metadata file
metadata <- read.delim("metadata.tsv")
rownames(metadata) <- metadata[,1]


#loading sequence table, taxonomy assigned and metadata table in the phyloseq object, ps
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.all, taxa_are_rows=FALSE), tax_table(speciesSilva), sample_data(metadata))

#By default name of all ASVs is their sequence. Replace all the names with suffix "ASV" followed by 1,2,.... And save the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
saveRDS(object = ps, file = "ps.RDS")

sort(sample_sums(ps))
table(taxa_sums(ps))

ls <- seq(0, 8500, by = 5000)
hist(sample_sums(ps), freq = TRUE, ylab = "Frequency of samples")

#removing ASVs with toal count less than 6
ps1 <- prune_taxa(taxa_sums(ps) > 5, ps)
ps1
saveRDS(object = ps1, file = "ps1.RDS")


#load required packages
library(ggplot2)
library(plyr)
library(ggpubr)
library(data.table)
library(randomcoloR)
library(tidyr)
library(microbiome)
library(dplyr)
library(UpSetR)
width = 210
height = 297

#######Alpha diversity
#plot
alpha <- plot_richness(ps1, measures = c("Observed", "Shannon"), x = "Media", color = "Description")
#save
ggsave(filename = "alpha-plot.pdf", plot = alpha, device = "pdf", width = width, height = height, units = "mm")
#extract values for statistical comparison
div <- microbiome::alpha(x = ps1, index = "all")
metadata$ShannonDiversity = div$diversity_shannon
metadata$ObservedASV = div$observed
##Comparisons
#comparison between type of approach
compare_means(formula = ObservedASV ~ Type, data = metadata, method = "wilcox.test", p.adjust.method = "BH")
compare_means(formula = ShannonDiversity ~ Type, data = metadata, method = "wilcox.test", p.adjust.method = "BH")
#comparison between collection
compare_means(formula = ObservedASV ~ Description, data = metadata, method = "wilcox.test", p.adjust.method = "BH")
compare_means(formula = ShannonDiversity ~ Description, data = metadata, method = "wilcox.test", p.adjust.method = "BH")
###

###Taxonomy exploration
#summarizing phyla level assignments
all <- as.data.frame(table(tax_table(ps1)[,2], useNA = "ifany"), stringsAsFactors = FALSE)
#summarizing phyla from CD samples only; include samples from CD only, remove taxa with sum(taxa) =0, make phylum level count and summarize
cd <- as.data.frame(table(tax_table(prune_taxa(taxa_sums(subset_samples(physeq = ps1, Type == "CD")) > 0, subset_samples(physeq = ps1, Type == "CD")))[,2], useNA = "ifany"), stringsAsFactors = FALSE)
#summarizing and plotting
colnames(all) <- c("Phylum", "All samples")
all$Phylum[is.na(all$Phylum)] <- "Unclassified phylum"
colnames(cd) <- c("Phylum", "CD samples")
ASV.count <- gather(data = merge(x = all, y = cd, all.x = TRUE), key = "Samples", value = "Frequency", -Phylum)
count.plot <- ggplot(data=ASV.count, aes(x = Phylum, y = Frequency*10, fill = Samples)) +   geom_bar(stat="identity", position=position_dodge()) + scale_y_log10(labels=function(x)x/10) + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + theme(legend.position = c(0.5, 0.9), legend.background = element_rect(fill = NA), legend.direction = "horizontal") + geom_text(aes(label=Frequency), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) + scale_fill_manual(values=c('#999999','#E69F00')) + labs(x = "Phyla", y = "Count of unique ASVs")
ggsave(filename = "ASV-count.pdf", plot = count.plot, device = "pdf", width = height*1.6, height = width, units = "mm")


###### Plot phyla level distribution
#Agglomerate taxonomy to phylum level including not assigned ASVs on relative abundances of ps1
ps1.Phylum.rel <- tax_glom(physeq = microbiome::transform(ps1, "compositional"), taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified phylum
taxa_names(ps1.Phylum.rel) <- tax_table(ps1.Phylum.rel)[, 2]
taxa_names(ps1.Phylum.rel)[is.na(taxa_names(ps1.Phylum.rel))] <- "Unclassified Phylum"
##plot box plot for top 15 phyla
TopPhyla <- names(sort(taxa_sums(ps1.Phylum.rel), TRUE)[1:15])
psdf <- data.table(psmelt(prune_taxa(TopPhyla, ps1.Phylum.rel)))
#colcodes.phylum <- distinctColorPalette(length(unique(psdf$OTU))+5) #commented to avoid selecting new colors everytime
Phylum.box <- ggplot(psdf, aes(x=OTU, y=Abundance*100, fill = Media))  + geom_boxplot(position = position_dodge(preserve = "single"), size = 0.2) + scale_fill_manual(values = colcodes.phylum) + scale_y_log10() + theme(legend.position = "left") + annotation_logticks(sides = "l") + scale_x_discrete(labels = function(labels) {sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n\n'), labels[i]))}) + ylab("Relative abundance of Phylum") + xlab("Phyla") + theme(plot.margin=unit(c(1,0.25,0.25,0.75), "cm"), axis.text.x = element_text(size=10)) + labs(fill = "Media") + geom_vline(xintercept = seq(1.5, 14.5, by = 1))
Phylum.box
ggsave(filename = "Phylum-box-plot.pdf", plot = Phylum.box, device = "pdf", width = height, height = width/2, units = "mm")



######## Converting to binary(presence-absence)
ps1.binary <- transform_sample_counts(ps1, function(x) ifelse(x>0,1,0))

#PCoA on bray-curtis distance
binary.bray <- distance(physeq = ps1.binary, method = "bray")
binary.ordination <- ordinate(ps1.binary, "PCoA", binary.bray)
ordinationplot <- phyloseq::plot_ordination(physeq = ps1.binary, ordination = binary.ordination, type = "samples", color = "Media", shape = "Description") + geom_point(size=3) + theme(aspect.ratio=0.8)
ggsave(filename = "binary-ordination-plot.pdf", plot = ordinationplot, device = "pdf", width = height, height =width, units = "mm")
#checking PERMANOVA for significant differences
vegan::adonis(formula = binary.bray ~ Type, data = metadata, method = "bray")
vegan::adonis(formula = binary.bray ~ Description, data = metadata, method = "bray")


#prepare data frame for comparison. Add group wise values
taxa_table <- as.data.frame(tax_table(ps1.binary))
binary.table <- merge(taxa_table, (t(ps1.binary@otu_table@.Data)), all = TRUE, by = "row.names")
#adding/summarizing data group-wise
binary.table = binary.table %>% mutate(CI = select(., R31:R45) %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(CD = select(., 'S3-AIA':'S4-YMA') %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(R3CI = select(., R31:R35) %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(R3CD = select(., 'S3-AIA':'S3-YMA') %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(R4CI = select(., R41:R45) %>% rowSums(na.rm = TRUE))
binary.table = binary.table %>% mutate(R4CD = select(., 'S4-AIA':'S4-YMA') %>% rowSums(na.rm = TRUE))

#prepare list and plot venn in UpsetR
R3CI <- as.character(binary.table[binary.table$R3CI > 0, "Row.names"])
R3CD <- as.character(binary.table[binary.table$R3CD > 0, "Row.names"])
R4CI <- as.character(binary.table[binary.table$R4CI > 0, "Row.names"])
R4CD <- as.character(binary.table[binary.table$R4CD > 0, "Row.names"])
group.list = list(R3CI = R3CI, R4CI = R4CI, R3CD = R3CD, R4CD = R4CD)
#plot upset plot
upset(fromList(group.list), sets = c("R3CI", "R4CI", "R3CD", "R4CD"), keep.order = TRUE,  order.by = "degree", nsets = 4, query.legend = "none", queries = list(list(query = intersects, params = list("R3CI", "R4CI", "R3CD", "R4CD"), color = "orange", active = T, query.name = "Present in all samples"), list(query = intersects, params = list("R3CD", "R4CD"),  color = "red", active = T, query.name = "CD"),list(query = intersects, params = list("R3CD"),  color = "red", active = T, query.name = "CD3"), list(query = intersects, params = list("R4CD"),  color = "red", active = T, query.name = "CD4"), list(query = intersects, params = list("R3CI", "R4CI"), active = T, color = "blue", query.name = "CI"), list(query = intersects, params = list("R3CI"), active = T, color = "blue", query.name = "3CI"), list(query = intersects, params = list("R4CI"), active = T, color = "blue", query.name = "4CI")))
#save the plot manually

#export CD and CI exclusive ASV table to plot Krona externally and compare taxonomy
binary.table.CD <- binary.table[binary.table$CI == 0,]
write.table(x = binary.table.CD, file = "binary.table.CD.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
binary.table.CI <- binary.table[binary.table$CI > 0,]
write.table(x = binary.table.CI, file = "binary.table.CI.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


#### Heatmap
#taxa ordering to keep CD exclusive ASVs on top
names.taxa <- taxa_names(ps1.binary)
names.taxa <- names.taxa[!names.taxa %in% binary.table.CD$Row.names]
taxa.name.order <- c(binary.table.CD$Row.names, names.taxa)
#plot heatmap
ASV.pa.heatmap <- plot_heatmap(ps1.binary, method = NULL, trans = NULL, sample.order = NULL, taxa.order = rev(taxa.name.order), max.label = 6000, low = "white", na.value = "white", high = "#000033") + xlab("Samples") + ylab("ASVs")
ggsave(filename = "ASV-pa-heatmap.pdf", plot = ASV.pa.heatmap, device = "pdf", width = width, height = height, units = "mm")


########## Genus level
#agglomerate to genus level
ps1.Genus <- tax_glom(physeq = ps1, taxrank = "Genus", NArm = FALSE)
#extract taxonomy file
genus.taxonomy <- as.data.frame(tax_table(ps1.Genus)[,1:6], stringsAsFactors = FALSE)
#Many of the taxa are assigned to higher levels than genus. "_X" is added to the end of such taxa. Repitition of suffix means assignment at higher level. For example, if taxa is assigned till family level its genus name will have single suffix; if assigned till class level its genus name will have 3 suffixes.
genus.taxonomy[is.na(genus.taxonomy$Phylum),]$Phylum <- paste0(genus.taxonomy[is.na(genus.taxonomy$Phylum),]$Kingdom, "_X")
genus.taxonomy[is.na(genus.taxonomy$Class),]$Class <- paste0(genus.taxonomy[is.na(genus.taxonomy$Class),]$Phylum, "_X")
genus.taxonomy[is.na(genus.taxonomy$Order),]$Order <- paste0(genus.taxonomy[is.na(genus.taxonomy$Order),]$Class, "_X")
genus.taxonomy[is.na(genus.taxonomy$Family),]$Family <- paste0(genus.taxonomy[is.na(genus.taxonomy$Family),]$Order, "_X")
genus.taxonomy[is.na(genus.taxonomy$Genus),]$Genus <- paste0(genus.taxonomy[is.na(genus.taxonomy$Genus),]$Family, "_X")
#Convert phyloseq object to presence-absence
ps1.Genus.binary <- transform_sample_counts(ps1.Genus, function(abund) 1*(abund>0))
#adding modified taxonomy information to phyloseq object
taxa_names(ps1.Genus.binary) <- genus.taxonomy[,6]
tax_table(ps1.Genus.binary)[,6] <- genus.taxonomy[,6]

#extracting taxa information with presence-absence information
genus.binary.table <- merge(ps1.Genus.binary@tax_table@.Data, (t(ps1.Genus.binary@otu_table@.Data)), all = TRUE, by = "row.names")
#adding/summarizing group-wise counts
genus.binary.table = genus.binary.table %>% mutate(CI = select(., R31:R45) %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(CD = select(., 'S3-AIA':'S4-YMA') %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(R3CI = select(., R31:R35) %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(R3CD = select(., 'S3-AIA':'S3-YMA') %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(R4CI = select(., R41:R45) %>% rowSums(na.rm = TRUE))
genus.binary.table = genus.binary.table %>% mutate(R4CD = select(., 'S4-AIA':'S4-YMA') %>% rowSums(na.rm = TRUE))
#prepare list and plot venn in UpsetR
genus.R3CI <- as.character(genus.binary.table[genus.binary.table$R3CI > 0, "Row.names"])
genus.R3CD <- as.character(genus.binary.table[genus.binary.table$R3CD > 0, "Row.names"])
genus.R4CI <- as.character(genus.binary.table[genus.binary.table$R4CI > 0, "Row.names"])
genus.R4CD <- as.character(genus.binary.table[genus.binary.table$R4CD > 0, "Row.names"])
genus.group.list = list(R3CI = genus.R3CI, R4CI = genus.R4CI, R3CD = genus.R3CD, R4CD = genus.R4CD)
#plot upset plot
upset(fromList(genus.group.list), sets = c("R3CI", "R4CI", "R3CD", "R4CD"), keep.order = TRUE,  order.by = "degree", nsets = 4, query.legend = "none", queries = list(list(query = intersects, params = list("R3CI", "R4CI", "R3CD", "R4CD"), color = "orange", active = T, query.name = "Present in all samples"), list(query = intersects, params = list("R3CD", "R4CD"),  color = "red", active = T, query.name = "CD"),list(query = intersects, params = list("R3CD"),  color = "red", active = T, query.name = "CD3"), list(query = intersects, params = list("R4CD"),  color = "red", active = T, query.name = "CD4"), list(query = intersects, params = list("R3CI", "R4CI"), active = T, color = "blue", query.name = "CI"), list(query = intersects, params = list("R3CI"), active = T, color = "blue", query.name = "3CI"), list(query = intersects, params = list("R4CI"), active = T, color = "blue", query.name = "4CI")))
#save plot manually

#extract CD exclusive genera for observation and checking media-wise detection
genus.binary.table.CD <- genus.binary.table[genus.binary.table$CI == 0,]

################THE END#############
