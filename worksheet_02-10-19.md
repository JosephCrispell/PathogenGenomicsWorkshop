---
title: "Pathogen Genomics Workshop"
author: "Joseph Crispell"
date: "04 Oct 2019"
output: 
  html_document:
    keep_md: true
---



---

## Introduction

Pathogens threaten the health of people and animals. Understanding pathogen transmission can help us understand how to control it.

Today we are going to be working with genomic data for the pathogen *Mycobacterium bovis*. *M. bovis* causes bovine tuberculosis in cattle and many other species. It costs Ireland millions to control. How can genomic data help?

We'll use some *M. bovis* genomic data sourced from infected cattle and wildlife to try and understand the role of wildlife. Wildlife species have been shown to harbour and transmit infection to cattle, we want to know if that is the case here in Ireland.

### Learning objectives

After our workshop we hope that you have learnt about:
- Programming in R
- Loading and using a FASTA nucleotide file
- Constructing and plotting a phylogenetic tree
- Why github is really useful
- Using and building R packages

---

## Step 1: Set your working directory

We'll be creating a few different plots over the course of todays workshop. As we progress, we would also recommend that you working along with us with your R script.

We'd recommend storing today's work (plots and your script) in a single folder. We'll start that process by creating a folder and setting it up as our working directory:


```r
# Get the current date
today <- format(Sys.Date(), "%d-%m-%y")

# Create a new directory on our desktop to store todays outputs
directory <- file.path("~", "Desktop", paste0("CRT_PathogenGenomicsWorkshop_", today), "")
dir.create(directory)

# Set your working directory
setwd(directory)
```



> QUESTION:<br>
> 1. What is a working directory?

---

## Step 2: Getting started

Next, we are going to install some R packages that we'll use throughout the workshop. The packages are `ape`, `phangorn` and `PathogenGenomicsWorkshopPackage`. The first two packages are commonly used for phylogenetic analyses in R. 

The `PathogenGenomicsWorkshopPackage` is an R package that we have specifically developed for this course. It has a few functions that we'll use later on.

To install the R packages, you can use the following code:


```r
# Install the 'ape' package
install.packages("ape", repos="https://cloud.r-project.org")

# Install the 'phangorn' package
install.packages("phangorn", repos="https://cloud.r-project.org")

# Install the 'pathogenGenomicsPackage' package
install.packages("devtools", repos="https://cloud.r-project.org") # We need devtools to install packages from github
devtools::install_github("JosephCrispell/pathogenGenomicsWorkshop")
```

With the packages successfully installed, we can now load them:



```r
# Load the required libraries
library(ape) # For reading in sequence
library(phangorn) # For testing substitution models and building and plotting the phylogeny
library(pathogenGenomicsWorkshop) # Our custom package for the current course
```

> QUESTION:<br>
> 1. Why create/use R packages?

---

## Step 3: Reading in the FASTA file

A FASTA file stores one or multiple nucleotide sequences. Our FASTA file stores the nucleotides present at a subset of genomic positions in `48` different *M. bovis* genomes. We can read in our *M. bovis* FASTA file with the following code:

```r
# Read in the FASTA file
fastaFile <- system.file("extdata", "Mbovis.fasta", package = "pathogenGenomicsWorkshop")
nucleotideAlignment <- read.dna(fastaFile, format = "fasta", as.character=TRUE)

# Convert the nucleotides to upper case
nucleotideAlignment <- toupper(nucleotideAlignment)
```

Notice by default nucleotides are stored in lower case, we've chosen to convert them to uppercase.

> QUESTIONS:<br>
> 1. Can anyone tell me what class of variable we have stored the sequences in?<br>
> 2. Do we have `48` sequences?<br>
> 3. How many positions are in the FASTA file?<br>

We also have a file that tells us which position on the *M. bovis* genome that each site in the FASTA file relates to. Let's read that in:

```r
# Read in the genome positions
positionsFile <- system.file("extdata", "Mbovis_FASTApositions.txt", package = "pathogenGenomicsWorkshop")
genomePositions <- read.table(positionsFile, header=TRUE)
```

> EXERCISE:<br>
> 1. Create the plot below<br>

<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

> QUESTION:<br>
> 1. Why might there be areas of the genome with more variants?<br>

Let's also take a quick look at the FASTA file. The code below will create a plot that is saved as a PDF in your working directory:

```r
plotFASTA(nucleotideAlignment, pdfFileName=paste0("FullNucleotideAlignment_", today, ".pdf"))
```

<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

---

## Step 4: Cleaning up the FASTA file

There are a lot of sites in the *M. bovis* alignment that aren't informative. We can clean up the alignment using the code below.


```r
# Count the nucleotides at each site in the alignment
nucleotideCountsAtEachSite <- countNucleotidesAtEachSite(nucleotideAlignment)

# Identify the uninformative sites
uninformativeSites <- which(nucleotideCountsAtEachSite < 2)

# Create a new nucleotide alignment without the uninformative sites
nucleotideAlignmentInformative <- nucleotideAlignment[, -uninformativeSites]
informativeGenomePositions <- genomePositions[-uninformativeSites, ]
```

> QUESTIONS:<br>
> 1. What is an uninformative site?<br>
> 2. What does line 5 in the above code block do?<br>

Now, let's take another look at the alignment, how has it changed?


```r
plotFASTA(nucleotideAlignmentInformative, 
          pdfFileName=paste0("InformativeSitesAlignment_", today, ".pdf"))
```
<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

> QUESTIONS:<br>
> 1. Can anyone guess what the nucleotide sequence at the top of the plot is?<br>
> 2. Could removing uninformative sites have any effect?<br>

---

## Step 5: Extract the sequence metadata from the IDs

As you will have seen, the sequence labels in our alignment contain some information about our sequences. Let's extract these data and store them in a `data.frame`:


```r
# Extract metadata from sequences
sequenceInfo <- getSequenceInfoFromNames(rownames(nucleotideAlignment))

# Take a quick look at the metadata
head(sequenceInfo)
```

```
##                        Name  Species SamplingDate
## 1 Seq-1_Wildlife_2014-11-13 Wildlife   2014-11-13
## 2 Seq-2_Wildlife_2014-12-09 Wildlife   2014-12-09
## 3 Seq-3_Wildlife_2014-11-18 Wildlife   2014-11-18
## 4 Seq-4_Wildlife_2014-12-09 Wildlife   2014-12-09
## 5      Seq-5_Cow_2015-05-08      Cow   2015-05-08
## 6 Seq-6_Wildlife_2015-03-19 Wildlife   2015-03-19
```

> EXERCISE:<br>
> 1. Calculate the number of samples sourced from wildlife and the number sourced from cattle<br>

---

## Step 6: Examine the quality of the nucleotide sequences

We don't have extensive data on the quality of our nucleotide sequences available but we can learn something about their quality from the nucleotide alignment. There are some `N`s in the alignment. 

> QUESTION:<br>
> 1. What do `N`s in a nucleotide alignment mean?<br>

Let's calculate the proportion of nucleotides in each sequence that are `N`s:

```r
# Calculate the proportion of Ns for each sequence
proportionNsInInformativeSites <- calculateProportionNsOfEachSequence(nucleotideAlignmentInformative)
```

> EXERCISE:<br>
> 1. Create the plot below:<br>

<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

There are a couple of nucleotide sequences that don't have data for ~8% of the genome sites.

> QUESTION:<br>
> 1. How might these differences in sequence quality impact our analyses?<br>

---

## Step 7: Build a phylogenetic tree

To build a phylogenetic tree we need to calculate the number of differences between each of our nucleotide sequences. We need to construct a genetic distance matrix:


```r
# Build a genetic distance matrix
distances <- dist.dna(as.DNAbin(nucleotideAlignmentInformative), model="raw")
```

Note that in the code above, we had to change the class (format) that we were storing our nucleotide alignment in.

Next, we'll build an initial neighbour-joining phylogenetic tree:


```r
# Build a quick initial phylogenetic tree
initialNJTree <- nj(distances)
```

The neighbour joining algorithm is a fast method to construct a phylogenetic tree but it isn't very robust. We are now going to construct a tree using the Maximum Likelihood algorithm. In addition, we are going to use bootstrapping to investigate the robusting of the phylogenetic tree structure.

> QUESTIONS:<br>
> 1. Why the Maximum Likelihood algorithm is a more robust tree building algorithm?<br>
> 2. How does bootstrapping work?<br>

To prepare to build our phylogeny, we'll construct a likelihood object using an initial tree and our nucleotide alignment:


```r
# Convert the nucleotide sequences into the PHYDAT format
sequencesInPhyDatFormat <- phyDat(nucleotideAlignmentInformative, type="DNA")

# Compute likelihood of the initial Neighbour Joining tree given sequences
likelihoodObject <- pml(initialNJTree, sequencesInPhyDatFormat)
```

With that object, we'll first run our maximum likelihood algorithm without bootstrapping:


```r
# Run maximum likelihood
fittingOutput <- optim.pml(likelihoodObject,
                           optNni = TRUE, # Optimise topology
                           optInv = TRUE, # Optimise proportion of variable sites
                           optBf = TRUE, # Optimise the base frequencies
                           model = "HKY", # Substitution model
                           rearrangement = "NNI", # Nearest Neighbour Interchanges
                           control = pml.control(maxit=100000)) # Set the maximum number of iterations
```

Lastly, now we'll take the output of our maximum likelihood analysis and feed it into a bootstrapping analysis:


```r
# Build a bootstrapped maximum likelihood phylogeny
bootstrapResults <- bootstrap.pml(fittingOutput,
                                  bs = 100,
                                  jumble = TRUE,
                                  control = pml.control(maxit=100000)) # Set maximum iteration number
```

---

## Step 8: Plotting the phylogenetic tree

Now that we have constructed and bootstrapped our maximum likelihood phylogenetic tree, let's take a look at it. First we'll need to extract the phylogeny from our bootstrapping output:


```r
# Get phylogenetic tree with bootstrap values
# Returns phylogenetic tree with bootstrap values as node labels
mlTreeBS <- plotBS(fittingOutput$tree, bootstrapResults, type="fan")
```



With the phylogeny stored as an object, we are going to create a simple plot:


```r
# Convert the branch lengths to approximate SNPs
mlTreeBS$edge.length <- mlTreeBS$edge.length * ncol(nucleotideAlignmentInformative)

# Set the plotting margins
par(mar=c(0.1, 0.1, 0.1, 0.1))

# Plot the phylogeny 
plot.phylo(mlTreeBS, show.tip.label = TRUE, edge.width = 2, type = "phylogram", edge.color = "grey")
```

<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

> QUESTIONS:<br>
> 1. What does line 2 in the above code block do?<br>
> 2. Why does the reference stick out so far?<br>

Let's remove the reference sequence and take a closer look at the phylogenetic relationships:


```r
# Remove the reference
mlTreeBSWithoutRef <- drop.tip(mlTreeBS, tip = "Reference_Cow_1997-10-15")

# Set the plotting margins
par(mar=c(0.1, 0.1, 0.1, 0.1))

# Plot the phylogeny - add in bootstrap values and species shapes
plot.phylo(mlTreeBSWithoutRef, show.tip.label = FALSE, edge.width = 3, type = "phylogram",
           edge.color = "grey")

# Add bootstrap values
bootstrapValues <- mlTreeBSWithoutRef$node.label/100
nodelabels(pch=19, frame="none", col=rgb(0,0,0, bootstrapValues), cex=0.75)

# Add node labels
tiplabels(pch=ifelse(grepl(mlTreeBSWithoutRef$tip.label, pattern = "Wildlife"), 21, 24), 
          bg=ifelse(grepl(mlTreeBSWithoutRef$tip.label, pattern = "Wildlife"), "red", "blue"),
          col = "dimgrey", cex=1.5)

# Add a scale
addSNPScale(position="topright", lineWidth=2)

# Add a species legend
legend("right", legend=c("Wildlife", "Cow"), pch=c(19, 17), col=c("red", "blue"), bty="n", xpd=TRUE)
```

<img src="worksheet_02-10-19_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

> QUESTIONS:<br>
> 1. How similar are the *M. bovis* bacteria infected cattle and wildlife?<br>
> 2. Does that have any implications for control?<br>
> 3. What further analyses could we use on these data?<br>

---

## Step 9: Wrapping up

Today we have analyses nucleotide sequence data derived from whole genome sequence *M. bovis* data. *M. bovis* is an important bacterial pathogen. 

Our samples were sourced from infected cattle and wildlife here in Ireland. Using these data, we are hoping to learn about what role wildlife are playing in Ireland's bovine tuberculosis problem.

Alongside the analyses of the *M. bovis* data, we've introduced aspects of programming in R, using github and creating and using an R package.

### Some useful resources

To finish up, we would like to point out some helpful resources:

- [Phylogenetic tree building](https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/)
- [Programming in R](https://www.tutorialspoint.com/r/index.htm)
- [Ask and answer questions](https://stackoverflow.com/questions/tagged/r)
- [Building an R package](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/)
- [Getting started with github](https://guides.github.com/activities/hello-world/)
- [Using R markdown](https://rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf)

---

## Additional dataset

If you've raced through the analyses of the *M. bovis* data. We have added an *Zaire ebolavirus* (Ebola virus) nucleotide sequence that you could work on. Ebola virus presents a massive risk to human health. The genomic data you'll examine is a subset of that sourced from infected people in West Africa [more information here](https://www.nature.com/articles/nature22040).

There are some really nice visualisations of these data can be found [here](https://nextstrain.org/ebola). You can get started with the Ebola virus data with the following code:

```r
# Read in the FASTA file
fastaFile <- system.file("extdata", "ebola_parsed.fasta", package = "pathogenGenomicsWorkshop")
nucleotideAlignment <- read.dna(fastaFile, format = "fasta", as.character=TRUE)

# Read in the genome positions
positionsFile <- system.file("extdata", "ebola_FASTApositions.txt", package = "pathogenGenomicsWorkshop")
genomePositions <- read.table(positionsFile, header=TRUE)
```

---

## Challenge

If you're having no trouble with the workshop content, take a look at the functions in our `pathogenGenomicsWorkshop` package. 

**Can you improve them?**

You'll find the code for each of the functions [here](https://github.com/JosephCrispell/pathogenGenomicsWorkshop/blob/master/R/pathogenGenomicsWorkshop.R). The `system.time()` function might be useful to calculate how long a function takes.


## TO DO

# Remove the boostrapping?

# Remove uninformative sites

# Remove the reference?

# Remove genome positions?

# Test on MAC and WINDOWS machines

# Make talking notes
- Everyone working together
- Main hopes from the workshop
- Challenges and extra data at the end
- Guage experience levels in the group
- Move between github, R and html worksheet
- Hope for some engagement with exercises and questions
