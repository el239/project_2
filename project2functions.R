getX <- function(theWAM){
  returnVector <- vector()
  for(i in 1:5){
    element <- paste(Goutput$regulatoryGene[i],Goutput$targetGene[i], sep = "->")
    returnVector <- c(returnVector,element)
  } # end for
  return(returnVector)
} # end function

getY <- function(theWAM){
  returnVector <- vector()
  for(i in 1:5){
    element <- Goutput$weight[i]
    returnVector <- c(returnVector,element)
  } # end for
  return(returnVector)
} # end function

# removes columns when sum of specified index elements exceed 10%
filterMatrix <- function(theMatrix, indeces){
  toDelete <- vector()
  for (i in 1:ncol(theMatrix)){
    total <- sum(theMatrix[,i])
    spike <- sum(theMatrix[indeces,i])
    if (spike/total > .1){
      toDelete <- c(toDelete,i)
    }
  }
  if (length(toDelete) > 0){
    return(theMatrix[,-(toDelete)])
  }
  else{
    return(theMatrix) 
  }
}
# removes rows where zero count values exceed input
filterZeros <- function(theMatrix, fractionZeros){
  toDelete <- vector()
  for (i in 1:nrow(theMatrix)){
    zeros <- which(0 == theMatrix[i,])
    if (length(zeros)/ncol(theMatrix) > fractionZeros){
      toDelete <- c(toDelete,i)
    }
  }
  if (length(toDelete) > 0){
    return(theMatrix[-(toDelete),])
  }
  else{
    return(theMatrix) 
  }
}
getCPM <- function(theMatrix, theFactors){
  for (i in 1:ncol(theMatrix)){
    theMatrix[,i] <- (theMatrix[,i]*1000000)/theFactors[i]
  } # end for
  return (theMatrix)
  # end function
} # end function

# finds size factor for each sample by getting median of mean-divided values, saves to vector
getSizeFactors <- function(theMatrix){
  # pseudocount
  theMatrix <- theMatrix + 1
  returnVector <- vector()
  means <- rowMeans(theMatrix)
  # gene count avg is divident for each column
  for(i in 1:ncol(theMatrix)){
    # gets ratios for the cell, vectorwise division using single element
    halfNormalized <- theMatrix[i]/means[i]
    # gets median (size factor for column/cell)
    factor <- median(halfNormalized)
    returnVector <- c(returnVector, factor)
  }
  return (returnVector)
}
