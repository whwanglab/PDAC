#' @title run_Q3
#'
#' Run Q3 normalization
#'
#' @param dataset the full dataset
#'
#' @param the_prob a value between 0 and 1 for the quantile, default is 0.75
#'
#' @return dfs the new data object post q3 normalization, updated norm counts and properties with new q3 norm factor
#'
#' @examples
#'  run_Q3(dataset = dfs)
#'
#' @export run_Q3

run_Q3 <- function(dataset = dfs, the_prob=0.75){
  # Check that the_prob is between 0 and 1
  if(class(the_prob)!="numeric"){
    stop("the_prob needs to be numeric")
  } else if(the_prob < 0 | the_prob > 1){
    stop("the_prob is numeric but not between 0 and 1 (inclusive)")
  }
  # If the number of genes left is less than 1000, stick with original Q3 norm
  if(dim(dataset[[2]])[1] < 1000){

    print("There are less than 1000 genes, will not re-calculate Q3 normalization factors, keeping current")

  }else {

    # Drop the genes from the raw counts file that match the normalized counts file
    indx <- which(!dataset[[1]]$TargetName %in% dataset[[2]]$Gene)
    if (length(indx) > 0){
      dataset[[1]] <- dataset[[1]][-indx,]
    }

    # Do a check to make sure the dimensions match
    if(dim(dataset[[1]])[1] == dim(dataset[[2]])[1]){
      print("Successfully removed genes, commencing Q3 renormalization")# Calculate the new Q3 normalization factor
      new_q3_normFactor <- nanoNormFactors(apply(dataset[[1]][2:ncol(dataset[[1]])], 2,
                                                 quantile, probs = the_prob,  na.rm=T))

      # Add it to the properties sheet with "Post" to show that it is after dropping genes
      dataset[[3]]$NormFactorQ3_postDrop <- new_q3_normFactor[match(dataset[[3]]$Sample_ID, names(new_q3_normFactor))]

      # Renormalize the data using the new Q3 factor
      rownames(dataset[[1]]) <- dataset[[1]]$TargetName
      q3 <- data.frame(t(t(dataset[[1]][,-1])/new_q3_normFactor))

      # Add the gene names back into the normalized dataframe
      q3 <- cbind("Gene" = rownames(q3), q3)

      # Ensure name formatting preserved
      colnames(q3)[-1] <- colnames(dataset[[1]])[-1]

      # Add Q3 normalized data back to the dataset object
      dataset[[2]] <- q3
    }else{
      print("Error in dropping genes from raw count files")
      stop()
    }
  }

  return(dataset)
}
