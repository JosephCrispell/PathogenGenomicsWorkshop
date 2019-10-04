# A collection of functions for pathogen genomics workshop for the SFI Centre for Research Training in Genomics Data Science
# Date created: 25-09-19
# Author: Joseph Crispell

### Building package commands (more info here: https://josephcrispell.github.io/BlogPosts/CreatingAnRPackage_11-09-18/CreatingAnRPackage_11-09-18.html)

## Packages to load
#library("devtools")
#library("roxygen2")

## Creating package
#packageDirectory <- "/home/josephcrispell/Desktop/pathogenGenomicsWorkshop/"
#create(packageDirectory)

## Documenting changes
#setwd(packageDirectory)
#document()

## Install package
#devtools::install_github("JosephCrispell/pathogenGenomicsWorkshop")

## Remove package
#remove.packages("pathogenGenomicsWorkshop")

#' An object of class "matrix" storing an alignment fo upper-case nucleotide sequences
#'
#' @name nucleotideAlignment
#' @docType data
#' @author Joseph Crispell \email{crispelljoseph@@gmail.com}
#' @keywords data
NULL

#' An object of class "phylo" storing a boostrapped Maximum Likelihood phylogenetic tree
#'
#' @name mlTreeBS
#' @docType data
#' @author Joseph Crispell \email{crispelljoseph@@gmail.com}
#' @keywords data
NULL

#' Calculate the proportion of \code{N}s in each nucleotide sequence
#'
#' Function calculates the proportion of \code{N}s in each nucleotide sequence in an alignment
#' @param nucleotideAlignment A matrix or a list containing the DNA sequences, or an object of class "alignment"
#' @keywords nucleotide alignment DNAbin
#' @export
#' @return Returns a numeric vector reporting the porportion of \code{N}s in each sequence in the alignment
#' @examples
#' # Load the nucleotide alignment
#' data("nucleotideAlignment")
#' 
#' # Calculate the proportion of Ns in each nucleotide sequence
#' proportionNs <- calculateProportionNsOfEachSequence(nucleotideAlignment)
calculateProportionNsOfEachSequence <- function(nucleotideAlignment){
  
  # Initialise a vector to store the proportion of Ns found in each sequence
  proportionNs <- c()
  
  # Note the number of sequences in the nucleotide alignment
  nSequences <- nrow(nucleotideAlignment)
  
  # Examine each sequence
  for(sequenceIndex in seq_len(nSequences)){
    
    # Count the number of each nucleotide in the current sequence
    nucleotideCounts <- table(nucleotideAlignment[sequenceIndex, ])
    
    # Check if any Ns were found
    if("N" %in% names(nucleotideCounts)){
      
      # Calculate the proportion of Ns in the current sequence
      proportionNs[sequenceIndex] <- nucleotideCounts[["N"]] / sum(nucleotideCounts)
      
    # If no Ns found, proportion = 0
    }else{
      proportionNs[sequenceIndex] <- 0
    }
  }
  
  return(proportionNs)
}

#' Extract species and sampling date from sequence names
#'
#' A specific function to retrieve metadata from formatted character strings
#' @param names A character vector containing specifically formatted character strings containing metadata: ID_species_\%Y-\%m-\%d
#' @keywords character metadata
#' @export
#' @return Returns a data.frame containing the metadata extracted from each character string
#' @examples
#' # Load the nucleotide alignment
#' data("nucleotideAlignment")
#' 
#' # Extract metadata from sequences
#' sequenceInfo <- getSequenceInfoFromNames(rownames(nucleotideAlignment))
getSequenceInfoFromNames <- function(names){
  
  # Create a dataframe to store the sequence information
  sequenceInfo <- data.frame("Name"=names, "Species"=NA, "SamplingDate"=NA)
  
  # Examine each of the sequence names
  for(index in seq_along(names)){
    
    # Split the current sequence's name into its parts
    parts <- strsplit(names[index], split="_")[[1]]
    
    # Store the species of the current sequence
    sequenceInfo[index, "Species"] <- parts[2]
    
    # Store the sampling date of the current sequence
    sequenceInfo[index, "SamplingDate"] <- parts[3]
  }
  
  # Convert the sampling dates into a date format
  sequenceInfo$SamplingDate <- as.Date(sequenceInfo$SamplingDate, format="%Y-%m-%d")
  
  return(sequenceInfo)
}

#' Count the number of different nucleotides in a sequence alignment 
#'
#' Function counts the number of each nucleotide (A, C, G, and T) at each position in a nucleotide alignment
#' @param nucleotideAlignment A matrix or a list containing the DNA sequences, or an object of class "alignment"
#' @keywords nucleotide alignment DNAbin
#' @export
#' @return Returns a numeric vector reporting the number of different nucleotides at each position in the input alignment
#' @examples
#' # Load the nucleotide alignment
#' data("nucleotideAlignment")
#' 
#' # Count the number of each different nucleotide at each position
#' nucleotideCounts <- countNucleotidesAtEachSite(nucleotideAlignment)
countNucleotidesAtEachSite <- function(nucleotideAlignment){
  
  # Initialise a vector to count the number of different nucleotides at each position (not including 'N's)
  nNucleotidesAtEachSite <- c()
  
  # Note the number of sites in the alignment
  nSites <- ncol(nucleotideAlignment)
  
  # Examine each site in the alignment
  for(siteIndex in seq_len(nSites)){
    
    # Get the unique nucleotides at the current position in the alignment
    nucleotides <- unique(nucleotideAlignment[, siteIndex])
    
    # Remove the N, if present
    nucleotides <- nucleotides[nucleotides != 'N']
    
    # Count the number of different nucleotides observed at the current site
    nNucleotidesAtEachSite[siteIndex] <- length(nucleotides)
  }
  
  return(nNucleotidesAtEachSite)
}

# Function to plot a nucleotide sequence alignment
# NOTE: designed to work with character matrix returned from read.dna()
#' Plot the different nucleotides present in an alignment of nucleotide sequences
#'
#' Function uses coloured polygons to illustrate which nucleotides are present in each sequence at each position
#' @param nucleotideAlignment A matrix or a list containing the DNA sequences, or an object of class "alignment"
#' @param pdfFileName A character string providing the name of a PDF to send plot into. Default value is NULL (plot not sent to file)
#' @param pdfWidth The width of the graphics region in inches. Default value is 14.
#' @param pdfHeight The height of the graphics region in inches. Default value is 7.
#' @param labelSpace Numeric multiplier to change the size of the space alotted to the sequence labels
#' @param lineForSequenceNames The number of lines into the margin at which the sequence names will be plotted. Default value is 0.
#' @param labelCex Numeric multiplier to change the size of the sequence labels. Default value is 0.5.
#' @keywords nucleotide alignment DNAbin plot
#' @export
#' @examples
#' # Load the nucleotide alignment
#' data("nucleotideAlignment")
#' 
#' # Plot the nucleotide alignment
#' plotFASTA(nucleotideAlignment)
plotFASTA <- function(nucleotideAlignment, pdfFileName=NULL, pdfWidth=14, pdfHeight=7, labelSpace=1, lineForSequenceNames=0, labelCex=0.5){

  # Open a pdf if requested
  if(is.null(pdfFileName) == FALSE){
    pdf(pdfFileName, width=pdfWidth, height=pdfHeight)
  }
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(4.1, 5*labelSpace, 3.1, 0.5))
  
  # Note the colour of each nucleotide
  nucleotideColours <- list('A'="green", 'C'="blue", 'G'="black", 'T'="red", 'N'="white")
  
  # Note the number of sequences and sites
  nSequences <- nrow(nucleotideAlignment)
  nSites <- ncol(nucleotideAlignment)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=c(1, nSites), ylim=c(1, nSequences), bty="n", yaxt="n", ylab="", xlab="Position")
  
  # Examine each sequence
  for(sequenceIndex in seq_len(nSequences)){
    
    # Examine each position in the current sequence
    for(siteIndex in seq_len(nSites)){
      
      # Get the nucleotide at the current position
      nucleotide <- nucleotideAlignment[sequenceIndex, siteIndex]
      
      # Draw a polygon for the current poisition coloured by nucleotide
      polygon(x=c(siteIndex-0.5, siteIndex-0.5, siteIndex+0.5, siteIndex+0.5), 
              y=c(sequenceIndex-0.5, sequenceIndex+0.5, sequenceIndex+0.5, sequenceIndex-0.5),
              border=rgb(0,0,0,0), col=nucleotideColours[[nucleotide]])
    }
  }
  
  # Add the sequence names
  axis(side=2, at=seq_len(nSequences), labels=rownames(nucleotideAlignment), tick=FALSE, las=1, line=-1.75+lineForSequenceNames, 
       cex.axis=labelCex)
  
  # Add nucleotide legend
  legend(x=nSites/2, y=nSequences + (0.15*nSequences), horiz=TRUE, xpd=TRUE, pch=22, col="grey", bty="n", xjust=0.5,
         legend=names(nucleotideColours), pt.bg=unlist(nucleotideColours), pt.cex=1.5)
  
  # Reset the plotting margins
  par(mar=currentMar)
  
  # Close the PDF if requested
  if(is.null(pdfFileName) == FALSE){
    dev.off() 
  }
}

#' Add a SNP scale to phylogeny
#'
#' Function that adds a SNP scale onto a plotted phylogeny
#' @param x An optional numeric X coordinate for the scale
#' @param y An optional numeric Y coordinate for the scale
#' @param postion An optional character vector detailing the location of the plot (top, middle, left, right, etc.)
#' @keywords phylogeny scale nucleotide snp
#' @export
#' @examples
#' # Load the nucleotide alignment
#' data("mlTreeBS")
#' 
#' # Remove the reference
#' mlTreeBSWithoutRef <- drop.tip(mlTreeBS, tip = "Reference_Cow_1997-10-15")
#'
#' # Plot the phylogeny - add in bootstrap values and species shapes
#' plot.phylo(mlTreeWithoutRef, show.tip.label = FALSE, edge.width = 3, type = "phylogram", edge.color = "grey")
#'
#' # Add a scale
#' addSNPScale(position="topright", lineWidth=2)
addSNPScale <- function(x=NULL, y=NULL, position=NULL, size=1, lineWidth=1){
  
  # Get the axis limits of the current plot
  axisLimits <- par("usr")
  
  # Get the label location
  coords <- c(x, y)
  
  # If no coordinates then calculate location based on a string (if no strin, defaults to centre)
  if(is.null(x)){
    coords <- calculateLocation(position, axisLimits)
    x <- coords[1]
    y <- coords[2]
  }
  
  # Add a scale
  lines(x=c(x-(0.5*size), x+(0.5*size)), y=c(y, y), lwd=lineWidth, xpd=TRUE)
  text(x=x, y=y, labels=paste(size, ifelse(size > 1, "SNPs", "SNP")), pos=1, xpd=TRUE)
}

#' Get the X and Y coordinates of a string coded location on a plot
#'
#' Function used by \code{addSNPScale()}
#' @param location A character string indicating a location on a plot
#' @keywords internal
#' @return Returns a numeric vector storing the coordinates determined by string location
calculateLocation <- function(location=NULL, axisLimits){
  
  # Initialise a vector to store the coordinates of the location
  coords <- c(NA, NA)
  
  # Calculate the length of each axis
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  
  # Calculate the centre if no location provided as a default
  coords[1] <- axisLimits[1] + (xLength*0.5)
  coords[2] <- axisLimits[3] + (yLength*0.5)
  
  # Check if location text given
  if(is.null(location) == FALSE){
    
    # Check if location is one of the options
    if(location %in% c("top", "bottom", "middle", "left", "right",
                       "topleft", "topright",
                       "bottomleft", "bottomright") == FALSE){
      location <- NULL
      stop("Location string not recognised. Please try one of: c(\"top\", \"bottom\", \"middle\",", 
           "\"left\", \"right\", \"topleft\", \"topright\", \"bottomleft\", \"bottomright\")")
    }
    
    # Define the location using the string code
    if(grepl(location, pattern="top")){
      coords[2] <- axisLimits[4] - 0.1*yLength
    }
    if(grepl(location, pattern="left")){
      coords[1] <- axisLimits[1] + 0.1*xLength
    }
    if(grepl(location, pattern="right")){
      coords[1] <- axisLimits[2] - 0.1*xLength
    }
    if(grepl(location, pattern="bottom")){
      coords[2] <- axisLimits[3] + 0.1*yLength
    }
  }
  
  return(coords)
}
