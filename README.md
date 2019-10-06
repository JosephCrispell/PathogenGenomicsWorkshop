<img align="left" src="Logos/ucd.gif" width=12%>
<img align="right" src="Logos/nuig.png" width=40%>

<br/><br/><br/><br/><br/>

# Pathogen Genomics Workshop
## Author: Joseph Crispell
## Licence: GPL-3
## Depends: ape, phangorn

<br/><br/>

## Description

The SFI Centre for Research Training in Genomics Data Science ([visit the website](https://genomicsdatascience.ie/)) has organised a two day workshop on the 7<sup>th</sup> & 8<sup>th</sup> of October (2019). As part of this workshop, we are running a Pathogen Genomics session.

Our two hour session aims to introduce the 2019 cohort of PhD students to building a phylogenetic tree in R. After our workshop we hope that you have learnt about:

- Programming in R
- Loading and using a FASTA nucleotide file
- Constructing and plotting a phylogenetic tree
- Why github is really useful
- Using and building R packages

In preparation for the current workshop, we have created the current respository. This respository stores an R package `pathogenGenomicsWorkshop` that we created to act as a teaching tool in our session.

The `pathogenGenomicsWorkshop` R package contains a set of functions that were specifically designed with the workshop and datasets in mind. During the workshop, feel free to take a look around the package.

## Installation

If you have time before the workshop, you can install the `pathogenGenomicsWorkshop` R package using the following code:
```
# Install the devtools package - required to install an R package stored on github
install.packages("devtools", repos="https://cloud.r-project.org")

# Install the pathogenGenomicsWorkshop package
devtools::install_github("JosephCrispell/pathogenGenomicsWorkshop")
```

## Preparation

Alongside install our R package. It would also be helpful if you were able to install the following additional packages:
```
# Installing the 'ape' package
install.packages("ape", repos="https://cloud.r-project.org")

# Installing the 'phangorn' package
install.packages("phangorn", repos="https://cloud.r-project.org")
```
