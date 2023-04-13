

## CROPS_DROUGHT EXPERIMENT - Yudi M. Lozano_ 2023 ##
## Bioinfo ##
## Dada2 bioinformatic pipeline for ITS sequences based in tutorial version 1.8 ## 
## https://benjjneb.github.io/dada2/ITS_workflow.html ##

# https://www.youtube.com/watch?v=t08I1uaim8k  # nice tutorial step by step 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

install.packages("ShortRead");install.packages("Rcpp");install.packages("rlang")
install.packages("vctrs");install.packages("cli")

library(Rcpp);library(dada2);library(ShortRead)
library(tidyverse);library(dplyr)

######## Adapter, primer removal and filtering steps ########
path <- "C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/sequences_1_2"

list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#p5fits7 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGARTCATCGAATCTTTG"
#p7its4 <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTCCTCCGCTTATTGATATGC"

#Remove ambiguous bases -> "N" 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE) # on windows, set multithread = FALSE

#sanity check
test <-letterFrequency(sread(readFastq(fnFs[2])), letters = "N")
sum(test)
#[2] 3902
test <-letterFrequency(sread(readFastq(fnFs.filtN[2])), letters = "N")
sum(test)
#[2] 0  

# Adapter sequences and orientations
p5 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
p7 <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

allOrients_ad <- function(adapter) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(adapter)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
p5.orients.ad <- allOrients_ad(p5) #all possible orientations of p5
p7.orients.ad <- allOrients_ad(p7) #all possible orientations of p7

# find out whether adapters are present among sequences
adapterHits <- function(adapter, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(adapter, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(p5.ForwardReads = sapply(p5.orients.ad, adapterHits, fn = fnFs.filtN[[2]]), # change number within brackets to examine other samples
      p5.ReverseReads = sapply(p5.orients.ad, adapterHits, fn = fnRs.filtN[[2]]), 
      p7.ForwardReads = sapply(p7.orients.ad, adapterHits, fn = fnFs.filtN[[2]]), 
      p7.ReverseReads = sapply(p7.orients.ad, adapterHits, fn = fnRs.filtN[[2]]))

# Reverese compliments of adapters are present

# Check if Cutadapt can be run form R
cutadapt <- "C:/Users/Admin/Desktop/BIBS/bioinfotrial/cutadapt.exe"
system2(cutadapt, args = "--version")

# Path to files with adapter free reads
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnf.cut <- file.path(path.cut, basename(fnFs.filtN))
fnr.cut <- file.path(path.cut, basename(fnRs.filtN))

p7.RC <- dada2:::rc(p7)
p5.RC <- dada2:::rc(p5)

# Arguments to trim the reverse-complement of p7 off of R1 (Fowards)
R1.flags.ad <- paste("-b", p7.RC, "--minimum-length 10", "-O 10") #sequences with less than 10 bp were removed
# Arguments to trim the reverse-complement of p5 off of R2 (Reverse)
R2.flags.ad <- paste("-B", p5.RC, "--minimum-length 10", "-O 10") #sequences with less than 10 bp were removed

# Trim adapters with Cutadapt 
# If running on windows, set multiple cores are not supported. So set multiple core option as FALSE

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags.ad, R2.flags.ad, "-n", 1, # -n 1 required to remove adapter only once
                             "-o", fnf.cut[i], "-p", fnr.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

## check for adapter presence in files
rbind(p5.ForwardReads = sapply(p5.orients.ad, adapterHits, fn = fnf.cut[[1]]), #Juan tenia [2], pero el video [1]
      p5.ReverseReads = sapply(p5.orients.ad, adapterHits, fn = fnr.cut[[1]]), 
      p7.ForwardReads = sapply(p7.orients.ad, adapterHits, fn = fnf.cut[[1]]), 
      p7.ReverseReads = sapply(p7.orients.ad, adapterHits, fn = fnr.cut[[1]]))

# all clean yupiiii!!

# check presence of Ns
Ntest <-letterFrequency(sread(readFastq(fnf.cut[1])), letters = "N")
sum(test)
#[1] 0

# Now remove primers
fITS7 <- "GTGARTCATCGAATCTTTG"
ITS4 <- "TCCTCCGCTTATTGATATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(fITS7) #all possible orientations of forward prim
REV.orients <- allOrients(ITS4) #all possible orientations of reverse prim
FWD.orients

#find out whether primer sequences are present in sequences
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnf.cut[[1]]), #JUan tenia en todo esto [2]
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnr.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnf.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnr.cut[[1]]))

# Mostly original orientation of primers and reverse compliments. 
#There are odd ocurrences of forward primers in 5' of reverse reads

# Path to files with primer free reads
path.cut1 <- file.path(path, "cutadapt1")
if(!dir.exists(path.cut1)) dir.create(path.cut1)
fnfs.cut <- file.path(path.cut1, basename(fnf.cut))
fnrs.cut <- file.path(path.cut1, basename(fnr.cut))

# get orientations of interest
FWD.RC <- dada2:::rc(fITS7)
REV.RC <- dada2:::rc(ITS4)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", fITS7, "-a", fITS7, "-a", REV.RC, "--minimum-length 10") #sequences with less than 10 bp were removed
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", ITS4, "-A", fITS7, "-A", FWD.RC, "--minimum-length 10") #sequences with less than 10 bp were removed

# Run Cutadapt 
for(i in seq_along(fnf.cut)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnfs.cut[i], "-p", fnrs.cut[i], # output files
                             fnf.cut[i], fnr.cut[i])) # input files
}

#check for primer presence in files
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnfs.cut[[1]]), #Juan tenia [2] pero video [1]
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnrs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnfs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnrs.cut[[1]]))

# All clear # All with zeros which is perfect

# check presence of Ns
Ntest <-letterFrequency(sread(readFastq(fnfs.cut[1])), letters = "N") #Juan tenia 2...
sum(test)
#[1] 0

#concatenate files containing primer and adapter free reads into R objects

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut1, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Check if forward and reverse files match (FROM THE VIDEO)
if(length(cutFs) == length(cutRs)) print ("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print ("Forward and reverse files do match. Go back anc check")

# Extract sample names for latter (for the phyloseq):
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
saveRDS(sample.names,"C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/names_file")

#Generate quality profile plots for our reads (FROM THE VIDEO)
# In case we have more than 20 fastqfiles, the following command will randomly choose 20 files to be plotted

if (length(cutFs) <=20) {
  fwd_qual_plots <-plotQualityProfile(cutFs) +
    scale_x_continuous(breaks = seq(0,300,20)) +
    scale_y_continuous (breaks = seq(0,40,5)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <-plotQualityProfile(cutRs) +
    scale_x_continuous(breaks = seq(0,300,20)) +
    scale_y_continuous (breaks = seq(0,40,5)) +
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size=20, 1:length(cutFs))#grab 20 random samples to plot
  fwd_qual_plots <-plotQualityProfile(cutFs[rand_samples]) +
    scale_x_continuous(breaks = seq(0,300,20)) +
    scale_y_continuous (breaks = seq(0,40,5)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <-plotQualityProfile(cutRs[rand_samples]) +
    scale_x_continuous(breaks = seq(0,300,20)) +
    scale_y_continuous (breaks = seq(0,40,5)) +
    geom_hline(yintercept = 30)
}

fwd_qual_plots
rev_qual_plots

#explanation quality plots
#All of our samples have at least 280 bp (reads length)
# we defined quality at 30, and later we need to trimm some ends here. Bt maybe just put at 20?
# 30 maybe would remove a lot
#In the rev plots que quality drops a bit earlier which is is expected in Miseq
#Green line is the mean quality 
#gray colors is a heatmap of the quality at that position 

#just to see but no need of it
plotQualityProfile(fnFs) # before cutadapt. the primers still are 

#In case you wanto to plot reads of a particular sample

fwd_qual_plots <- plotQualityProfile(cutFs[1]) + #first sample. Also can use cutFs[1:5],cutFs[c(2,4,7)] etc
  scale_x_continuous(breaks = seq(0,300,20)) +
  scale_y_continuous(breaks = seq(0,40,5)) +
  geom_hline(yintercept = 30)
rev_qual_plots<-plotQualityProfile(cutRs[1])+
  scale_x_continuous(breaks = seq(0,300,20)) +
  scale_y_continuous(breaks = seq(0,40,5)) +
  geom_hline(yintercept = 30)


####### FILTER AND TRIM

#create folder "filtered" within cutadapt folder
filtFs <- file.path(path.cut1, "filtered", basename(cutFs))
filtRs <- file.path(path.cut1, "filtered", basename(cutRs))

# get sample names object
sample.names <- readRDS("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/names_file")

# Filter and trim reads
# I chose not to truncate reads at the same lenghts because ITS2 varies in lenght across fungal taxa.
#In the video they agree. But if we want to use it it will be as follows
#truncLen=c(200, 140) # following the plotquality where the quality dropped 
# However, I specified a truncate threshold when mean quality score of base calls drops below 10. 
# Reads less than 50 bp in lenght are discarded
# MaxEE. Maximum expected error set to 2 bp . reads higher than this will be discarded. You should
#set a higher MaxEE for reverse reads du to their lower quality (2,4).
#There is an equation (video) but somehow for 300 reads we would expect 3% error 
#maxN=0 . No allowing reads with ambiguous bases
#rm.phix=TRUE. Removes any reads that match the phix bacteriphage genome which is typically added
#to Illumina sequencing
#If we leave multihread TRUE is not problem. In windows it doesnot work.
#notes from Juan and Video min 29 


#filter command  
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN=0, maxEE = c(2,4), 
                     truncQ = 10, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE

# table tracking the loss of samples during initial filtering steps
head(out)
tail(out) # looks good, did not lose too much reads.
saveRDS(out,"C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/filter_and_trim_out_FUN_run1.rds") 

# Extract sample names for latter:
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
saveRDS(sample.names,"C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/names_file")

#Generate sample names with out the file extensions and check if the file names match

sample.names <- sapply(strsplit(basename(filtFs), "_L001"),'[',1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_L001"),'[',1) 
if (identical(sample.names, sample.namesR)) {print("Files are still matching...congratulations")
} else {stop("Forward and reverse files do not match..sorry")}

# some stringent filtering routines eliminated sample S059, so filtFs objects need to be updated or commands below will throw an error
#filtFs <- filtFs[-17] 
#filtRs <- filtRs[-17] 

#Error model generation  (use machine learning algorithm to "learn" error nd then see how error rate relates to consesus quality score)

#dada 2 learns the specific error-signature of our data set
#This is why files from separate sequencing runs have to be processed separately till later on
#this later with help us to know which are real sequences, biological seequences, error sequences.
set.seed (123)
errF <- learnErrors(filtFs, multithread = F)
errR <- learnErrors(filtRs, multithread = F)
#?learnErrors

#In the case the filter and trim step has removed all reads from a file (e.g in blanck samples).Erro
#stimation does not work, Then run the following 
#Define filtFs and filtRs and sample.names again
#path_new <-"-/miseq_run1_fastq_files/cutadapt/filtered
# there are more lines see video min 37 

#save error calculation as RDs files
setwd("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou") # check the directory where you are working 
saveRDS (errF, "errF_FUN_run1.rds")
saveRDS (errR, "errR_FUN_run1.rds")

#Visualize the errors
plotErrors(errF, nominalQ = TRUE) # error frequency in base calls decrease with higer quality score - looks good
plotErrors(errR, nominalQ = TRUE) #same

#it shows the error rates for each possible transition (A=>C, A=>G, etc)  
#Points are the observed error rates for each consensus quality scores.
#Error rates should drop with increasing quality score

###DEREPLICATION

#### The dada2 tutorial still implements a dereplication at this point 
#This does not seem to be necesary any more with the newer dada2 version

#Apply the dada2 core sequence-variant inference algoritm 

#Set pool: TRUE to allow information to be shared across samples
#This makes it easier to resolve rare variants whcih occur just once or twice
#but increase computation time (problem with larga data sets)and can give you false positives 
#if you dont pool you will not take into accout some rare sequences with e.g.,
#low abundance may be identified as sequencing error..you can use ". "pseudo is something in between

#Find 'real' variants with DADA2 algorithm = otherwise known as 'sample inference'
dadaFs <-dada(filtFs, err=errF, multithread = F, pool = "pseudo") 
dadaRs <-dada(filtRs, err=errR, multithread = F, pool = "pseudo")

#Save sequence-variant inference as RDS files which my be uploaded 
getwd()
saveRDS (dadaFs, "dataFs_FUN_run1.rds")
saveRDS (dadaRs, "dataRs_FUN_run1.rds")

# Inspecting the returned dada-class object of the first sample
dadaFs[[1]] # 116 sequences variants from the total 9887 input unique sequences
# This wouls be the true biological sequences 

### There are mor of dereplication in Mps_drought_fungal 

######## MERGE READS ######

#Note. If you have data from 2 different pools sequencing ..check video aroun min 48

#lets continue with our data 
#Merge forward and reverse together to obtain the full denoised
#Adjust the minium overlap (default =12) and the maximum mismatch allowed
#e.g., mismatch of one nucleotide

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 10, maxMismatch = 1, verbose=TRUE)
saveRDS (mergers, "mergers_FUN_run1.rds")
head(mergers[[1]])

#Construct an amplicon variance sequence table (ASV table)
#if maxmitmach >0 has been allowed in the mergePairs step,
#"duplicate sequences detected and merged" may appear as output

ASV_table <- makeSequenceTable(mergers)

#How many seequences variants were inferred?
dim(ASV_table) #2 samples, 2 sequences in total
#Merging discards a lot of reads! Means there is not much overlap between R1 and R2

ASV_tab_Forw <- makeSequenceTable(dadaFs) # Try again using only forward reads
dim(ASV_tab_Forw) 
#[1]   2 159 
#Forward reads are good quality and long enough. 
# I will use only the forward reads table and sequencess from now on according to the recomendations in
# Pauvert et al. 2019. Fungal Ecology

#REMOVE CHIMERAS
#PCR artefacts taht may inflate diversity estimates
#Chimeric sequences

ASV_tab.nochim_Forw <- removeBimeraDenovo(ASV_tab_Forw, method="consensus", multithread=F, verbose=TRUE)
# Identified 7 bimeras out of 159 input sequences.

# inspect distribution of ASV lengths
table(nchar(getSequences(ASV_tab.nochim_Forw))) 
hist(nchar(getSequences(ASV_tab.nochim_Forw)), xlab="Lenght of reads (bp)", 
     main = "ASV length distribution")
#looks good -- most reads have between 150-200 bp

# Save ASV table
getwd()
saveRDS(ASV_tab.nochim_Forw, "ASV_tab_nochim_Forw.rds")

#calculate percentage of the reads that were non-chimeric
ASV_tab_Forw #filtered
ASV_tab.nochim_Forw #non chimeras 

sum(ASV_tab.nochim_Forw)/sum(ASV_tab_Forw) # only 0.3% non chimeras..that is very god

# inspect the number of reads removed in every step 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(ASV_tab.nochim_Forw))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "denoised","denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/track_samples1.csv", fileEncoding = "utf8")

#write a fasta file of the final, non chimeric sequences
asv_seqs <-colnames (ASV_tab.nochim_Forw)
asv_headers <-vector(dim(ASV_tab.nochim_Forw)[2],mode="character")
for (i in 1:dim(ASV_tab.nochim_Forw)[2]){
  asv_headers[i] <-paste (">ASV", i, sep="")
}
asv_fasta<-c(rbind(asv_headers,asv_seqs))
write(asv_fasta, "non_chimeras_ASV.fa")

#write an ASV count table of the final non chimeric sequences

colnames (ASV_tab.nochim_Forw) <-paste0("ASV",seq(ncol(ASV_tab.nochim_Forw)))
ASV_counts<-t(ASV_tab.nochim_Forw)
write.table(ASV_counts, file="ASV_counts_nochim.txt", sep="\t",quote = F, col.names = NA)

#Remove singletones from the non-chimeric ASVs ( si aun hay)

#First read back  your table "ASV_tab.nochim_Forw" into R again.
myASV_nonchim <- readRDS("ASV_tab_nochim_Forw.rds")
#transform they to numeric as they will likely be integer
mode (myASV_nonchim)="numeric"
#subset columns with counts of >1
myASV_nonchim_nosingle <-myASV_nonchim[,colSums(myASV_nonchim)>1]

#write a fasta file of the final, non-chimeric, non-singleton sequences
asv_seqs <-colnames(myASV_nonchim_nosingle)
asv_headers <-vector(dim(myASV_nonchim_nosingle)[2], mode="character")
for (i in 1:dim(myASV_nonchim_nosingle)[2]){
  asv_headers[i] <-paste ("ASV", i, sep="")
}

#write an ASV count table of the final, non-chimeric sequences 
colnames (myASV_nonchim_nosingle) <-paste0("ASV", seq(ncol(myASV_nonchim_nosingle)))
ASV_counts <-t(myASV_nonchim_nosingle) # transposing table
ASV_counts <-as.data.frame(ASV_counts)
ASV_counts$ASV_ID <-rownames(ASV_counts) #Add new column with ASV_ID 
write.table (ASV_counts, file="ASV_counts_no_chim_nosingle.txt", sep="\t", quote = F, col.names = NA)

#track reads
getN <- function(x) sum(getUniques(x))
track_nonchim <-cbind(sapply(mergers, 
                             getN), rowSums(ASV_tab.nochim_Forw), rowSums(myASV_nonchim_nosingle))
colnames(track_nonchim) <- c("merged","nonchim","nonsingle")
track_nonchim
write.table(track_nonchim, "track_nonchim_nonsingle.txt",sep="\t", col.names = NA)

#Joining tracka
firsttrack <- read.csv("track_samples1.csv")
trackALL<-cbind(firsttrack, track_nonchim)
write.table (trackALL, file="track_ALL_samples.txt", sep="\t", quote = F, col.names = NA)
# remember that I do NOT use the "merged samples" but the FORWARD samples

# TAXONOMIC ASSIGNATION 

#Download general release UNITE database
#https://doi.plutof.ut.ee/doi/10.15156/BIO/1280049

#the command is an implementation of RDP bayesian classifier Wang et al. 2007 doi:10.1128/AEM.00062-07
# assignation specified with minBoot confidence to 80 and repeated with default (minboot = 50)

# ran it on FUs HPC cluster with the folling parameters:
# Note that in this cluster one thread is equivalent to one core.
# mem-per-cpu refers to RAM memory 

#!/bin/bash
#SBATCH -N1 -n4 --mem-per-cpu=8192M --qos=hiprio -t00:30:00 # this setup (32 gb RAM, 4 cores) took about 11 minutes
#SBATCH -N1 -n8 --mem-per-cpu=1024M --qos=hiprio -t00:10:00 # this setup (8 gb RAM, 8 cores) took about 7 minutes
#SBATCH -N1 -n16 --mem-per-cpu=512M --qos=hiprio -t00:10:00 # this setup (8 gb RAM, 16 cores) took about 7 minutes
#module add R-bundle-Bioconductor/3.9-foss-2019a-R-3.6.0
#Rscript RDPclassifierR.r

# 8 GB RAM is more than enough for a job of about 1 million sequences
# the optimal set up was with 8 cores - with more cores the process is less efficient.

#contents of RDPclassifierR.r

setwd("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou")
library(dada2)
set.seed(123)
ASV_tab.nochim <- readRDS("ASV_tab_nochim_Forw.rds") # this is the only .RDS ASv_tab
taxa_80 <- assignTaxonomy(ASV_tab.nochim, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/unite_2021/sh_general_release_dynamic_10.05.2021.fasta", multithread = TRUE,
                          tryRC = F, minBoot = 80)
saveRDS(taxa_80, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxaunite_file_stringent.rds")


taxa_50 <- assignTaxonomy(ASV_tab.nochim, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/unite_2021/sh_general_release_dynamic_10.05.2021.fasta", multithread = TRUE,
                          tryRC = F, minBoot = 50)
saveRDS(taxa_50, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxaunite_file_default.rds")


######## CHECK DIFFERENCES IN TAXONOMIC ASSIGNATION####

library(dplyr)
setwd("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou")
taxa.un.de <- as.data.frame(readRDS(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxaunite_file_default.rds", sep="")), stringsAsFactors = F)

unid_s_de = taxa.un.de %>% summarise(sum(is.na(Species))) 
unid_g_de = taxa.un.de %>% summarise(sum(is.na(Genus)))
unid_f_de = taxa.un.de %>% summarise(sum(is.na(Family)))
unid_o_de = taxa.un.de %>% summarise(sum(is.na(Order)))
unid_p_de = taxa.un.de %>% summarise(sum(is.na(Phylum))) 
unid_k_de = taxa.un.de %>% summarise(sum(is.na(Kingdom))) # all fungal sequences

taxa.un.str <- as.data.frame(readRDS(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxaunite_file_stringent.rds", sep="")), stringsAsFactors = F)

unid_s_str = taxa.un.str %>% summarise(sum(is.na(Species)))
unid_g_str = taxa.un.str %>% summarise(sum(is.na(Genus)))
unid_f_str = taxa.un.str %>% summarise(sum(is.na(Family)))
unid_o_str = taxa.un.str %>% summarise(sum(is.na(Order)))
unid_p_str = taxa.un.str %>% summarise(sum(is.na(Phylum)))
unid_k_str = taxa.un.str %>% summarise(sum(is.na(Kingdom))) # all fungal sequences

#create a table
unid.sum <- rbind(cbind(unid_k_de, unid_p_de, unid_o_de,unid_f_de,unid_g_de,unid_s_de),
                  cbind(unid_k_str, unid_p_str, unid_o_str,unid_f_str,unid_g_str,unid_s_str))
rownames(unid.sum) <- c("RDP_Default","RDP_Stringent")
colnames(unid.sum) <- c("NAs_Kingdom","NAs_Phylum","NAs_Order","NAs_Families","NAs_Genus","NAs_Species")
unid.sum
#unid.sum["Total_ASVs"] <- c("1137","1137") # no se que es esto

write.csv(unid.sum, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxid.csv", fileEncoding = "utf8")

# The cost of being stringent in terms of identification thesholds translates in 
unid.sum$NAs_Genus[2]-unid.sum$NAs_Genus[1] # 12 ASVs or
round(((unid.sum$NAs_Genus[2]-unid.sum$NAs_Genus[1])/unid.sum$NAs_Genus[2])*100,2) # 17.65% 
# decrease in ASVs identified to the level of genus but a gain in 17.65% percent more accuracy. 

# I will use the identities assigned by the stringent identification schema.

######## Trait assignation ####

# clean up if required #
rm(list=ls())
library(dplyr)
setwd("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou")
# Load and edit taxonomic assignation table (here the 80 stringeny)
taxa.un.str <- as.data.frame(readRDS(paste("taxaunite_file_stringent.rds", sep="")), stringsAsFactors = F)

taxa.print <- taxa.un.str # create new object
rownames(taxa.print) <- NULL # Removing sequence rownames for display only
head(taxa.print)


# Edit header format; subset only ASV names along with assigned genera 
library(dplyr)
taxa.tr   = taxa.print %>%
  mutate(Kingdom=stringr::str_replace(Kingdom, c("k__"),""))%>%
  mutate(Phylum=stringr::str_replace(Phylum, c("p__"),""))%>%
  mutate(Class=stringr::str_replace(Class, c("c__"),""))%>%
  mutate(Order=stringr::str_replace(Order, c("o__"),""))%>%
  mutate(Family=stringr::str_replace(Family, c("f__"),""))%>%
  mutate(Genus=stringr::str_replace(Genus, c("g__"),""))%>%
  mutate(Species=stringr::str_replace(Species, c("s__"),""))%>%
  as_tibble() %>%
  tibble::rownames_to_column() %>% dplyr::rename(ASV_ID=rowname) %>%
  mutate(ASV_ID=paste0("ASV", 1:nrow(taxa.print))) %>%
  as.data.frame()#%>%
#select(ASV_ID,Genus)
#the above did not function if R is so loaded... so restart R
View(taxa.tr) # same as taxa print but without the ASV column
###############################

# load "FungalTraits" database Põlme et al. 2021. Fungal Diversity https://doi.org/10.1007/s13225-020-00466-2
#see this
#https://github.com/traitecoevo/fungaltraits
#https://link.springer.com/article/10.1007/s13225-020-00466-2
# I have downloaded the table GENUS from Polme

FT <- read.csv("C:/Users/Admin/Desktop/BIBS/bioinfotrial/FUNGALTRAITSDB_Genus.csv") %>%
  rename(Genus=GENUS,Comment.on.genus=COMMENT.on.genus,secondary_lifestyle=Secondary_lifestyle)%>%
  select(-jrk_template,-Phylum,-Class,-Order,-Family) #remove these columns

# left join trait and taxa tables to get traits

taxa.tr = taxa.tr %>% left_join(FT, by="Genus")
taxa.tr <- as.matrix(taxa.tr) #transform to a character matrix 
rownames(taxa.tr) <- taxa.tr[,1]#use col 1 as rownames to be able to assemple a phyloseq object

# save table 
write.csv(taxa.tr, "C:/Users/Admin/Desktop/BIBS/bioinfotrial/taxtr.csv", fileEncoding = "utf8")


## JOIN ASV and TAXA in a single table 
getwd()
taxa <-read.csv("taxtr.csv") #taxonommy table
ASV_counts_final <- read.table("ASV_counts_no_chim_nosingle.txt")
FIN_ASV_taxa <- taxa %>% left_join(ASV_counts_final, by="ASV_ID")

write.csv(FIN_ASV_taxa,"C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/FIN_ASV_taxa.csv", fileEncoding = "utf8")

str(FIN_ASV_taxa)
str(taxa)
#
#


#####ESTO NO ESTA FUNCIONANDO ..CREO Q PORQUE NO TENGO TODAS LAS READS! 

######## Merge all available data into a phyloseq object #####


# clean up if required #
rm(list=ls())

setwd("C:/Users/Admin/Desktop/BIBS/bioinfotrial")
# necesary library
library(dplyr)
library(phyloseq)

# get sample names object
sample.names <- readRDS(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/names_file"))
View(sample.names)
# read experimental metadata
meta_crop <- read.csv(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/metadata_crops_drou.csv", sep=""))
rownames(meta_crop) <- sample.names #give sample names used in original files to experimental data
meta_crop$Treatment <- paste0(meta_crop$species,":",meta_crop$treatment)
''
#load ASV tab
ASV_tab.nochim_Forw <- readRDS(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/ASV_tab_nochim_Forw.rds", sep=""))
ASV_tab_final <- ASV_tab.nochim_Forw # create new object
colnames(ASV_tab_final) <- paste0("ASV", 1:ncol(ASV_tab_final)) #edit new object colnames - remove sequences
class(ASV_tab_final) <- "numeric"

# load taxonomy as character matrix
taxa.tr <- as.matrix(read.csv(paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/taxtr.csv", sep=""), fileEncoding = "utf8", row.names=1))

# create phyloseq object
crops_dr <- phyloseq(tax_table(taxa.tr), otu_table(ASV_tab_final, taxa_are_rows=F), sample_data(meta_mp))
crops_dr

saveRDS(crops_dr, paste("C:/Users/Admin/Desktop/BIBS/bioinfotrial/crops_drou/mic_dr_trial.rds", sep=""))

