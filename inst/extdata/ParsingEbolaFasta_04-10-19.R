#### Load libraries ####
library(ape)
library(pathogenGenomicsWorkshop)

#### Read in the ebola FASTA file ####

# NOTE: FASTA file retrieved from: https://github.com/nextstrain/ebola/blob/master/example_data/ebola.fasta

# Read in the nucleotide alignment
fastaFile <- paste0("~/Desktop/pathogenGenomicsWorkshop/inst/extdata/ebola.fasta")
nucleotideAlignment <- read.dna(fastaFile, format = "fasta", as.character=TRUE)

# Convert the nucleotide characters to upper case
nucleotideAlignment <- toupper(nucleotideAlignment)

# Replace any gaps with Ns
nucleotideAlignment[nucleotideAlignment == "-"] <- "N"

#### Remake the sequence names ####

newNames <- c()
for(i in seq_len(nrow(nucleotideAlignment))){
  
  parts <- strsplit(rownames(nucleotideAlignment)[i], split="\\|")[[1]]
  parts[6] <- gsub(parts[6], pattern="_", replacement="-")
  newNames[i] <- paste0("Seq-", i, "_", parts[6], "_", parts[4])
}
rownames(nucleotideAlignment) <- newNames

#### Check for uninformative sites ####

# Count the nucleotides at each site in the alignment
nucleotideCountsAtEachSite <- countNucleotidesAtEachSite(nucleotideAlignment)

# Identify the uninformative sites
uninformativeSites <- which(nucleotideCountsAtEachSite < 2)

# Create a new nucleotide alignment without the uninformative sites
nucleotideAlignmentInformative <- nucleotideAlignment[, -uninformativeSites]

# Create a dataframe to record the positions of the variant positions that were retained
positions <- seq_len(ncol(nucleotideAlignment))
informativePositions <- data.frame("Position"=positions[-uninformativeSites])

#### Print the fasta to file ####

# Write the informative FASTA to file
outputFile <- paste0("~/Desktop/pathogenGenomicsWorkshop/inst/extdata/ebola_parsed.fasta")
write.dna(nucleotideAlignmentInformative, outputFile, format="fasta", colsep="")

# Write the gneome positions to file
outputFile <- paste0("~/Desktop/pathogenGenomicsWorkshop/inst/extdata/ebola_FASTApositions.txt")
write.table(informativePositions, file=outputFile, quote=FALSE, sep="\t", row.names=FALSE)
