# Based on this excellent tutorial (accessed 2018-10-29):
# https://benjjneb.github.io/dada2/tutorial.html
library(dada2)

path = "../1.raw"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# ===============
# Check quality scores of reads. 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280, 280),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

# Dereplicate fastq files. 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Manual inspection
dadaFs[[1]]


# Merge denoised pairs of forward and reverse reads. Those which do not sufficiently overlap or contain too many mismatches are rejected. 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Distribution of lengths is wider than I would have expected. Does this represent some error upstream?
#
# To filter on lengths, do this: seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


# What proportion retained after chimera removal?
sum(seqtab.nochim)/sum(seqtab)
# >  0.8112675

# -----------------------
# Rename samples (rows in the table, currently). 
# THIS IS DATASET-SPECIFIC. Done using the mapping file. 
new.sample.names = row.names(seqtab.nochim)
mapping_file = read.csv(file="../mapping_file.txt", sep="\t", header=TRUE)

# Probabily a more efficient way to do this using the "map" function. 
for(i in 1:nrow(mapping_file)) {
  illumina_name = mapping_file[i, "IlluminaID"]
  sample_name = as.character(mapping_file[i, "Sample"])
  index = which(new.sample.names == illumina_name)
  new.sample.names[index] = sample_name
}
row.names(seqtab.nochim) = new.sample.names


# -----------------------
# Assess DADA2 performance. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#     input filtered denoisedF denoisedR merged nonchim
# 1  140868    36613     36315     35965  33369   28014
# 10 135881    44606     44215     43894  40497   34747
# 11 119837    45453     45100     44793  40768   31416
# 12 120347    44862     44662     44497  42267   30272
# 13 151873    56803     56571     56498  54038   42891
# 14 143059    49120     48553     48253  44487   36968
#
# MR: lots of reads lost at the filtering stage. Why?


# -----------------------
taxa <- assignTaxonomy(seqtab.nochim, "referenceDBs/dada2/gg_13_8_train_set_97.fa.gz", multithread=TRUE)



# Let's inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# -----------------------
# Save to filesystem. 
# Equiv to an OTU table. 
write.table(t(seqtab.nochim), "seqtab-nochim.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(taxa, file="seqs_taxonomy.txt", sep=",", quote=FALSE)



# theme_set(theme_bw())
# library(phyloseq); packageVersion("phyloseq")
# library(ggplot2); packageVersion("ggplot2")
# theme_set(theme_bw())

# biom_file = import_biom(BIOMfilename='/Users/markread/dropbox_usyd/projects/henstridge/4.otus/otu_table_filtered_json.biom')

# treefile = read_tree_greengenes("/Users/markread/Dropbox (Sydney Uni)/projects/micro_ecology_resources/referenceDBs/gg_13_8_otus/trees/97_otus.tree")
# refseqfile = ("/Users/markread/Dropbox (Sydney Uni)/projects/micro_ecology_resources/referenceDBs/gg_13_8_otus/rep_set/97_otus.fasta")
# samplefile = import_qiime_sample_data("/Users/markread/dropbox_usyd/projects/henstridge/mapping_file.txt")
