# Script Name: mds_phylogeny.R
#
# Author: Helena Cooper
# Last edited: 10/12/2020
#
# Description: The functions in this script perform the MDS reduction of multiple phylogeny trees and prioritisation of 
#              species for ribosome structures.
#

#-----------------------------------------------------------------------------------------------------------------------

######## Specify requried packages 
packages <- c("ape", "ggplot2","ggforce", "dplyr","ggrepel","ggpubr","stringr")

######## Now install and/or load all packages.
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#-----------------------------------------------------------------------------------------------------------------------

####### Function for converting a Newick phylogeny tree to a patrisitic distance matrix.

patristic_distance <- function(outtree) {
  
  ### outtree is the generated phylogeny tree for each ribosomal sequence type.
  
  # Reads in Newick format as a "phylo" object
  tree <- read.tree(outtree)
  # Computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
  matrix <- cophenetic.phylo(tree)

  ### Remove species that have dodgy sequence, and sort matrix alphabetically.
  
  toMatch <- c("Structural", "Other", "Pathogenic")

  # Removes anything with incomplete metadata - should end with either "Structural", "Other" or "Pathogenic".
  matrix <- matrix[ grep(paste(toMatch,collapse="|"), rownames(matrix), value=TRUE), grep(paste(toMatch,collapse="|"), colnames(matrix), value=TRUE)]

  # Remove genome ID so that the data is alphabetically ordered by species name.
  i = 1
  while ( i <= nrow(matrix) ){
    rownames(matrix)[i] <- gsub("^.*?_" ,"" , rownames(matrix)[i])
    colnames(matrix)[i] <- gsub("^.*?_" ,"" , colnames(matrix)[i])
    i = i + 1
  }
  
  # Sort matrix alphabetically and return generated distance matrix to environment.
  matrix <- matrix[order(row.names(matrix)),order(colnames(matrix))] ; return(matrix)

}

#-----------------------------------------------------------------------------------------------------------------------

###### Function for adding two distance matrices together, while removing species that aren't present in both matrices.

combine_two_matrices <- function(one,two) {
  
  ### one and two are distance matrices produced by patristic_distance()
  
  # Convert rownames to a column for filtering purposes
  one.rows <- tibble::rownames_to_column(as.data.frame(one),"Name")
  two.rows <- tibble::rownames_to_column(as.data.frame(two),"Name")
  
  ### Process first matrix to contain species in second matrix.

  # Subset species/rows from matrix one that are present in matrix two.
  one.subset <- one.rows[(one.rows$Name %in% two.rows$Name),]
  # Get a list of column numbers of species/rows that are not present in matrix two.
  list <- as.integer(rownames(one.rows[!(one.rows$Name %in% two.rows$Name),]))
  # Make the species names as rownames, and remove the column containing rownames. 
  rownames(one.subset) <- one.subset[,1] ; one.subset <- one.subset[,c(2:ncol(one.subset))]
  # Delete associated columns for species aren't present in matrix two (ie: create a square matrix)
  one.subset[,list] <- NULL 
  # Convert back into a matrix.
  one.matrix <- as.matrix(one.subset)

  ### Process second matrix to contain species in first matrix (same logic as matrix one).

  two.subset <- two.rows[(two.rows$Name %in% one.rows$Name),]
  list <- as.integer(rownames(two.rows[!(two.rows$Name %in% one.rows$Name),]))
  rownames(two.subset) <- two.subset[,1] ; two.subset <- two.subset[,c(2:ncol(two.subset))]
  two.subset[,list] <- NULL 
  two.matrix <- as.matrix(two.subset)

  ### Add matrices together.
    
  sum.matrix <- one.matrix + two.matrix

  ### Return generated distance matrix to environment.
  
  return(sum.matrix)

}

#-----------------------------------------------------------------------------------------------------------------------

####### Function for the MDS reduction of the combined phylogeny distance matrix

mds_reduction <- function(m) {
  
  ### m is the final distance matrix that is being used for plotting. 
  
  # MDS reduction using classical (metric) multidimensional scaling method. 
  mds <- as.data.frame(cmdscale(as.matrix(m)),quote=FALSE)
  
  ### Format metdata columns for plotting purposes. 
  
  # Make the species metadata (rownames) as a column.
  mds$Description <- rownames(mds)
  # Extract and format species name, and target type into separate columns.
  n = 1 ; bacteria <- c() ; target <- c() ; phyla <- c() ; class <- c()
  
  while ( n <= nrow(m) ) {
  
    # Extract metadata for each row with spaces, rather than underscores.
    name <- gsub("_"," ",rownames(m)[n])
    species_name <- str_split(name, " ")[[1]][1]
    species_name <- sub("^(\\w).*", "\\1.", species_name)
    species_name <- paste(species_name,str_split(name, " ")[[1]][2])
    
    # Extracts species name, without an extra space at the end (eg: "Escherichia coli" not "Escherichia coli ") and append to list.
    bacteria <- append(bacteria, sub("^(\\S*\\s+\\S+).*", "\\1", species_name)) 
    
    # Extracts target type (ie: if it's a species of interest) and appends to list.
    phyla <- append(phyla, str_split(name, " ")[[1]][3]) 
    class <- append(class, str_split(name, " ")[[1]][4]) 
    target <- append(target, str_split(name, " ")[[1]][5])          
    
    # Update count.
    n = n + 1
  }

  # Append formatted metadata to MDS dataframe and return to environment.
  mds$Bacteria <- bacteria ; mds$Target <- target ; mds$Phyla <- phyla ; mds$Class <- class ; return(mds)

}

#-----------------------------------------------------------------------------------------------------------------------

###### Calculate minimum distance between bacteria and ribosome structure species.

min_solved_dist <- function(matrix, mds){

  ### matrix is the combined distance matrix.
  ### mds is the reduced distance matrix.

  # Identify distance matrix positions for each structural species.
  n = 1
  structures <- c()
  while ( n <= nrow(matrix) ) {
    if ( grepl("Structural",rownames(matrix)[n],fixed=TRUE) == TRUE ) {
      structures <- append(structures, n)
    }
    n = n + 1
  }

  # Calculate all distances between species and structures, then determine the minimum.
  o = 1
  min_dist <- c()   # Distance from species to structural species.
  min_name <- c()   # Name of the structural species which had the smallest phylogenetic distance.
  while ( o <= nrow(mds) ) {
    dist = 10000    # Set large maximum so first hit can be set as the "minimum", which is the threshold for all remaining loops.
    name = ''
    for ( struc in structures ) {
      d = matrix[o,struc]
      if ( d < dist ) {
        dist = d
	      name = rownames(matrix)[struc]
      }
    }
    
    min_dist <- append(min_dist,dist)
    min_name <- append(min_name,name)
    
    o = o + 1
  }

  # Append data and return datafrane to the environment.
  mds$Min_Distance <- min_dist ; mds$Min_Name <- min_name ; return(mds)
    
}

#-----------------------------------------------------------------------------------------------------------------------

###### Optional function for parsing through proposed sturctures to find an appropriate/non-pathogenic option.

select_non_pathogen <- function(options) {

  ### options is a vector of all potential structures that meet the filtering thresholds.
  
  response = n
  i = 1
  while ( response != "y" & i <= length(options) ) {
    prompt <- paste("Is this species non-pathogenic to humans (y/n):",options[i]," ",sep=" ")
    response <- (readline(prompt))
    if ( response == "y" ){
      species <- options[i]
    } else {
      species <- "Unavailable"
    }
    i = i + 1
  }
  return(species)
  
}

#-----------------------------------------------------------------------------------------------------------------------

####### Function for sampling the most representative species per phyla.

sampling_by_distance <- function(matrix, mds, automate=FALSE){   

  ### matrix is the combined distance matrix.
  ### mds is the reduced distance matrix.
  ### automate calls select_non_pathogen should specific species need to be selected.

  # Get list of unique phyla, excluding the candidate phyla.
  class <- unique(mds[(mds$Phyla != "CandidatusDependentiae" & mds$Phyla != "Proteobacteria" & mds$Phyla != "candidatedivisionNC10"),6])
  class <- append(class, unique(mds[(mds$Phyla == "Proteobacteria"),7]) )
  
  solved_struc_threshold <- summary(mds$Min_Distance)[2]  # Lowest 25% for all distance to solved structure.

  # Select a representative structure for each phyla.
  i = 1
  while ( i <= length(class) ){
  
    # Identify all species above the lower quartile threshold.
    if ( grepl("proteobacteria",class[i]) ){
      species_list <- na.omit(mds[( mds$Class == class[i] & mds$Min_Distance > solved_struc_threshold & ! grepl("sp.", mds$Bacteria) ),3])
    } else {
      species_list <- na.omit(mds[( mds$Phyla == class[i] & mds$Min_Distance > solved_struc_threshold & ! grepl("sp.", mds$Bacteria)),3])
    }

    N = length(species_list)
    j = 1
    index_list=c()  # Index and Species lists are in the same order 
    
    # If the phyla contains at least two species, then select a representative structure.
    if ( N > 1 ){       
    
      # Get list of the positions in the combined distance matrix for all potential species in the phyla. 
      while ( j <= N ){
        index_list <- append(index_list, which(grepl(species_list[j], rownames(matrix))))
        j = j + 1
      }
    
      # Subset the combined distance matrix to only include those from the same phyla.
      matrix.subset <- matrix[index_list,index_list]

      species_average <- c()   # Average phylogenetic distance to all species in phyla/class.
      j = 1
    
      # Get list of the average distance observed between individual and all other members of their phyla.
      while ( j <= N ){      
        species_average <- append(species_average,mean(matrix.subset[j,]))
        j = j + 1
      }

      # Select the most representative structure.
      if ( "Pathogenic" %in% unique(mds[index_list,5])[2] ){    # If pathogen exists, prioritise the most representative one.
        min_test=1000 ; j = 1
        # Select the most representative (has the lowest average).
        while ( j <= N ){
          if ( grepl("Pathogenic",species_list[j]) & (species_average[j] < min_test) ){
	          min_index <- which(grepl(species_list[j], mds$Description))
	          min_test <- species_average[j]
	        } 
	        j = j + 1
        }
      } else {   # If no pathogen exists, select the species with the minimum average distance recorded.
        index <- which.min(species_average)
        min_index <- which(grepl(species_list[index], mds$Description))
      }
      
      # Change target classifications for filtering and re-calculation of the minimum distances.
      
      mds$Target[min_index] <- "Solved"   
      rownames(matrix)[min_index] <- gsub("Other|Pathogenic","Structural",rownames(matrix)[min_index])   
      colnames(matrix)[min_index] <- gsub("Other|Pathogenic","Structural",colnames(matrix)[min_index])
      
      # Update count.
      i = i + 1
      
    } else {
      # If there is one or no species in the phyla, only update the count and don't select a representative.
      i = i + 1
    }
  }

  # Recalculate the distances to a sovled structure taking the proposed structures into account.
  mds <- min_solved_dist(matrix,mds)
  
  # Save the simulated data as a separate matrix and dataframe, and return these to the environment.
  simulatedData <- list("mds" = mds, "matrix" = matrix) ; return(simulatedData)

}

