#This code is meant to be run in three code blocks with user input given at the beginning of blocks 2 and 3.
#See `https://benjjneb.github.io/dada2/tutorial.html` for more details on dada2 functions.

#Begins Block 1

library(dada2) #Necessary for analyzing sequencing runs.
library(pdftools) #Helps to concatenate pdfs of quality plots.
library(RPushbullet) #Used to notify the user when runs are completed. It is optional.

options(error = recover)

path <- "./"

list.files(path)

#change the pattern to whatever represents the forward and reverse reads.
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plot_quality_profiles <- function(plot = TRUE) {
  #This begins a block to concatenate the profile plots for every sample.
  maxF <- dim(array(fnFs))
  maxR <- dim(array(fnRs))

  for(i in 1:maxF) {
    profile_plot_F <- dada2::plotQualityProfile(fnFs[i])
    pdf(paste("profile_plot_F", i, ".pdf", sep = ""))
    plot(profile_plot_F)
    dev.off()
  }

  for(i in 1:maxR) {
    profile_plot_R <- dada2::plotQualityProfile(fnRs[i])
    pdf(paste("profile_plot_R", i, ".pdf", sep = ""))
    plot(profile_plot_R)
    dev.off()
  }

  #be very careful to run this in a pdf free folder, otherwise it will add extraneous files to the profile_plots pdf.
  pdfs <- sort(list.files(path, pattern = ".pdf", full.names = TRUE))
  #pdf_combine() will fail if it has to concatenate more than 500 pdfs at a time. Subset the data accordingly with:
  #pdf_combine(c(pdfs[1:500]), output = "profile_plots_1-500.pdf")
  pdf_combine(c(pdfs), output = "profile_plots.pdf")
  if(file.exists(c(pdfs)))
    unlink (c(pdfs))
}

#Evaluate quality profiles of forward reads and adjust truncLen parameters accordingly.
#Forward reads are usually higher quality than reverse reads.
#Generally, samples with forward reads which have quality greater than 30 through most of the read should be kept.
plot_quality_profiles(plot = TRUE)

#End of Block 1

#Eliminate samples with quality profiles which are below a specified threshold.
#A strict threshold would be to include samples with mean quality score above 30 for roughly 200 bases.
#One way to eliminate samples from the pipeline is to move them into a separate directory.

#Begining of Block 2
#change the pattern to whatever represents the forward and reverse reads.
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#adjust the truncLen according to what seems to be good from looking at the profile plots and take note of the truncation lengths used.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
              maxN=0, maxEE=c(2,7), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)

filtered_reads <- out

write.csv(filtered_reads , file = "filtered_reads.csv")

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#This is a pushbullet command. Completely optional.
pbPost("note", "Srt Err Lrn", "Now")

errF <- learnErrors(derepFs, multithread=TRUE)

errR <- learnErrors(derepRs, multithread=TRUE)

Errors_plot <- plotErrors(errF, nominalQ=TRUE)

pdf('Errors_plot.pdf')
plot(Errors_plot)
dev.off()

#Before the DADA2 main function it is a good idea to save an image just in case something goes wrong during analysis.
save.image(file = "Example_Image_before_dada2.Rda")
savehistory(file = "Example_History_before_dada2.Rhistory")

#This is the main dada2 function. Method can be "independent", "pseudopooled", or "pooled."
dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect the distribution of sequence lengths and decide if and where to have merged sequence length cutoffs.
table(nchar(getSequences(seqtab)))

#After the DADA2 main function it is also a good idea to save an image just in case something goes wrong during analysis.
save.image(file = "Example_Image_after_dada2.Rda")
savehistory(file = "Example_History_after_dada2.Rhistory")

#End of Block 2

#Beginning of Block 3

#trim sequences around target length and adjust if necessary.
#If there is a wide range of merged sequence lengths it would be good to trim them.
#This is currently set for V1-2 anticipated length.
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 252:287]
table(nchar(getSequences(seqtab2)))

#Remove chimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method=consensus, multithread=FALSE, verbose=TRUE)

#transpose seqtab to make easier to copy over
seqtab.nochim.t <- t(seqtab.nochim)
write.csv(seqtab.nochim.t, file = "seqtab.nochim.csv")

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/Path_To_Files/Dada2_formatted_silva/silva_nr_v138_train_set.fa.gz", minBoot=80, multithread=TRUE)

write.csv(taxa , file = "taxa.csv")

taxa.species <- addSpecies(taxa, "~/Path_To_Files/Dada2_formatted_silva/silva_species_assignment_v138.fa.gz")

write.csv(taxa.species , file = "taxa_with_species.csv")

merged <- cbind(taxa.species, seqtab.nochim.t)
write.csv(merged, file = "merged.csv")

save.image(file= "Example_Image_Final.rda")
savehistory(file= "Example_History_Final.Rhistory"); pbPost("note", "Dada2", "Completed")
#end of third chunk

#Proceed to install package dada2tools and follow along for post-processing analysis.
