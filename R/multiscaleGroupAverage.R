#' Get the bin resolution from the multiscale_methylation Snakemake workflow
#'
#' Imports bin/config.yaml and returns the steps
#' 
#' @param config path to the bin/config.yaml of the multiscale Snakemake
#'
#' @return vector of the steps
#'
#' @import yaml
#'
#' @export
#'
#' @examples
#'
#'
#'

getMultiscaleSteps <- function(config){
  if(!grep('\\.yaml$',config)) stop("config must be a yaml file")
  multiscale.config <- yaml::yaml.load_file(config)
  
  # Get the number of steps
  steps <- seq(multiscale.config$bounds$min,multiscale.config$bounds$max,by=multiscale.config$bounds$step)
  # nSteps <- length(steps)
  return(steps)
}


#' Helper function to get the average for a defined
#' 
#' @param means.list List of length nSamples. Each element of means.list is itself
#' a list of Granges of length nSteps
#' @param i index for which step to get the average from
#' 
#' @return a genomic ranges with the average values
#' 
#' @import GenomicRanges
#' 
#' @export
#' 
#' @examples
#'
#' 

.getMultiscaleStepGroupAverage = function(means.list,i){
  if (!is(means.list, "list")) stop("means.list needs to be a list")
  if (!is(means.list[[1]], "list")) stop("means.list[[1]] needs to be a list of lists.")
  if (!is(means.list[[1]][[1]], "GRanges")) stop("means.list[[1]][[1]] needs GRanges object.")
  step.gr <- lapply(means.list, `[[`, i) # get the ith step for all samples
  # now get the average
  sum <- rep(0,length(step.gr[[1]]))
  for(j in 1:length(step.gr)){
    sum <- sum + mcols(step.gr[[j]])$score
  }
  avg <- sum / length(step.gr) # get the average
  
  return.gr <- means.list[[1]][[i]] # take the first one for step i
  mcols(return.gr)$score <- avg # replace the score with the average
  return(return.gr)
}


#' Function to return a list of GRangesList for each steps
#' The output of this function can be passed directory to multiscaleMethylationPlot
#' @param multiscale_means_directory List of length nSamples. Each element of means.list is itself
#' a list of Granges of length nSteps
#' @param samples character vector of samples to include. These must match the sample column from bin/samples.tsv
#' in the multiscale workflow
#' @param config config.yaml from the multiscale_methylation Snakemake workflow
#' @param which GRanges for rtracklayer::import
#' 
#' @return a genomic ranges with the average values
#' 
#' @import GenomicRanges
#' 
#' @export
#' 
#' @examples
#'
#' 

multiscaleGroupAverage <- function(multiscale_means_directory,samples,config,which){
  steps <- getMultiscaleSteps(config)
  nSteps <- length(steps)
  # Get all of the files for in-group samples from the Snakemake means directory
  files <- unlist(lapply(samples,function(x){paste0(multiscale_means_directory,'/',x,'_',sprintf("%.1f", steps),'.bed.gz')}))
  
  # Now import each
  means <- lapply(
    files, function(x) {
      return(rtracklayer::import(x, format = "bedGraph", which = which))
    }
  )
  
  # Now split this to a list equal to the length of samples
  # Each element of means.list is a list of nStep granges (similar to if one sample or the data from stats dir is imported)
  means.list <- split(means, ceiling(seq_along(means)/nSteps))
  
  if (!length(samples)==length(means.list)) stop("Length of the means list should be the same as the number of samples.")
  # cat('Getting average from',length(means.list),'samples\n')
  
  # apply multiscaleGroupAverage for each step size to get a list of granges of length nStep
  group.average.mean <- lapply(seq_along(1:nSteps), .getMultiscaleStepGroupAverage, means.list=means.list)
  
  stopifnot(length(group.average.mean)==nSteps)
  
  group.average.mean.grl <- as(group.average.mean, "GRangesList")
  
  names(group.average.mean.grl) <- steps
  
  return(group.average.mean.grl)
}

