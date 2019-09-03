
#'@title parse_phylowgs
#'
#'  Parse phylowgs outputs to visualize using EvoFreq.
#'  
#'@param json_file json file output from PhyloWGS containing tree structure and cellular prevalance information.
#'@param return_all if false, return the best result based on the linearity index provided by PhyloWGS.
#'@return A list of dataframes of resolved phylogenies from PhyloWGS. 
#'
#'@examples
#'\donttest{
#'phylowgs_output="run_name.summ.json"
#'
#'tree_data <- parse_phylowgs(json_file=phylowgs_output)
#'
#'#EvoFreqPlots
#'pdf('./evofreqs.pdf', width=8, height=4, onefile = T)
#'for (i in 1:length(f$all)){
#'  clone_dynamics_df <- get_evofreq(tree_data[[i]][,c(5,length(colnames(tree_data)))], clones=tree_data[[i]]$clone, parents=tree_data[[i]]$parent, clone_cmap = "jet")
#'  p <- plot_evofreq(evofreq_df)
#'  print(p)
#'}
#'dev.off()
#'}
#'
#'@export
parse_phylowgs <- function(json_file, return_all=TRUE){
  if(substr(json_file, nchar(json_file)-2+1, nchar(json_file))=="gz"){
    theFile <- gzfile(json_file, open="rb")
    result <- rjson::fromJSON(file = theFile)
    close(theFile, type="rb")
  } else {
    result <- rjson::fromJSON(file = json_file)
  }
  
  i=1
  trees=list()
  for(item in result$trees){
    #clones
    clones <- unlist(item$structure)
    
    #clone parents
    parents <- names(item$structure)
    parentNames <- list()
    count = 1
    for(children in item$structure){
      parentNames[[count]] = rep(parents[count], length(children))
      count=count+1
    }
    parents <- as.numeric(unlist(parentNames))
    
    # Edge dataframe
    edges <- data.frame(parent=parents, clone=clones)
    root <- data.frame(parent=0, clone=0)
    edges <- rbind(root,edges)
    
    # Pull information
    node <- names(item$structure)
    count=1
    cellPrev <- list()
    snvs <- list()
    cnvs <- list()
    for(details in item$populations){
      cellPrev[[count]] = details$cellular_prevalence
      snvs[[count]] = details$num_ssms
      cnvs[[count]] = details$num_cnvs
      count=count+1
    }
    cellPrev <- do.call(rbind,cellPrev) %>% as.data.frame
    names(cellPrev) <- seq(1,length(names(cellPrev)))
    snvs <- do.call(rbind,snvs) %>% as.data.frame
    cnvs <- do.call(rbind,cnvs) %>% as.data.frame
    
    
    df <- cbind(data.frame(snv=snvs, cnv=cnvs), edges, cellPrev)
    colnames(df) <- c("snv","cnv", "parent","clone",  seq(1,length(names(cellPrev))))
    
    trees[[i]] <- list(evofreq=df, llh=item$llh, linearity_index=item$linearity_index, clustering_index=item$clustering_index)
    i=i+1
  }
  
  minVal <- Inf
  minIdx <- 0
  treeOut <- list()
  if(return_all==FALSE){
    for(i in 1:length(trees)){
      if(trees[[i]]$linearity_index<minVal){
        minVal=trees[[i]]$linearity_index
        minIdx=i
      }
    }
    treeOut[[1]] <- trees[[minIdx]]
  } else {
    treeOut <- trees
  }
  return(treeOut)
}


#'@title parse_calder
#'
#'  Parse Calder outputs to visualize using EvoFreq.
#'  
#'@param soln Clonal composition matrix file (soln csv file).
#'@param tree dot file produced by Calder.
#'@return A list of dataframes of resolved phylogenies from PhyloWGS. 
#'
#'@examples
#'\donttest{
#'theFile <- "SA501_tree1.dot"
#'theSoln <- "SA501_soln1.csv"
#'calder.data <- parse_calder(theSoln, theFile)
#'
#' ### Use the long_to_wide_size_df function to get the right data structure.
#' wide_df <- long_to_wide_size_df(long_pop_sizes_df = calder.data$sizeDf,
#'                                edges_df = calder.data$edges,
#'                                time_col_name = "time",
#'                                clone_col_name = "clone",
#'                               parent_col_name = "parent",
#'                                size_col_name = "size",
#'                                fill_gaps_in_size = T
#')
#'clones <- as.character(wide_df$clones)
#'parents <- as.character(wide_df$parents)
#'size_df <- wide_df$wide_size_df

#'clone_dynamics_df <- get_evofreq(size_df, clones, parents, clone_cmap = "jet", data_type = "size", threshold=0, test_links = T, add_origin = T, interp_method = "bezier")
#'plot_evofreq(clone_dynamics_df)
#'@export
parse_calder <- function(soln, tree){
  # Parse tree file
  lines <- readLines(tree)
  # Getting label information
  meta <- lines[grep(" [",lines,fixed=TRUE)]
  metaData <- unlist(strsplit(meta, " "))
  clones <- metaData[seq(1,length(metaData),2)]
  clones <- as.character(clones)
  cloneLabels <- metaData[seq(2,length(metaData),2)]
  cloneLabels <- unlist(strsplit(cloneLabels, "[label=\"]"))[grep(":",unlist(strsplit(cloneLabels, "[label=\"]")))]
  
  body <- lines[grep("->", lines, fixed = TRUE)]
  nodePairs <- sub("^[[:space:]]+\"", "\"", sub("\"[;[:space:]]+$","\"", unlist(strsplit(body, "->"))))
  nodeLists <- split(nodePairs, 1:length(nodePairs)%%2)
  nodes <- unique(nodePairs)
  edges <- as.matrix(data.frame(parent = sub(" ","",nodeLists[[2]]), 
                                clone = sub(" ","",sub(";","",nodeLists[[1]])), stringsAsFactors = F))
  
  lines <- readLines(soln)
  fileInfo <- lines[seq(grep("U",lines, fixed=T)+1,length(lines))]
  headerInfo = unlist(strsplit(fileInfo[1],","))
  proportion <- list()
  timepoints <- list()
  for(i in 2:length(fileInfo)){
    d <-unlist(strsplit(fileInfo[i], ','))
    timepoints[[i]] <- d[1]
    proportion[[i]] <- as.numeric(d[2:length(d)])
  }
  prop.mat <- as.data.frame((do.call(rbind, proportion[2:length(proportion)])), stringsAsFactors = F)
  colnames(prop.mat) <- clones
  dataMat <- cbind(unlist(timepoints), prop.mat)
  names(dataMat) <- c("time", names(prop.mat))
  df <- melt(dataMat, id.vars = c("time"))
  colnames(df) <- c("time","clone","size")
  df <- data.frame(time = as.character(df$time), clone = as.character(df$clone), size = df$size, stringsAsFactors = F)
  
  out <- list(clones = clones, edges = as.data.frame(edges, stringsAsFactors = F), cloneLabels=cloneLabels, sizeDf = df)
  print(out)
  print(str(out))
  wide_df <- long_to_wide_size_df(long_pop_sizes_df = out$sizeDf,
                                  edges_df = out$edges,
                                  time_col_name = "time",
                                  clone_col_name = "clone",
                                  parent_col_name = "parent",
                                  size_col_name = "size",
                                  fill_gaps_in_size = T
  )
  
  return(wide_df)
}

#' @title read.HAL
#' @param path_to_file String defining the location of the data output from HAL
#' @param fill_name Optional string defining which attribute to color by. See \code{\link{update_colors}} for more details. 
#' @param get_evofreq_df Boolean defining if a freq_frame and evofreq should be returned. 
#' @param get_evofreq_arg_list List containing additional arguments passed to \code{\link{get_evofreq}}
#' @param get_dendogram_df Boolean defining if the dataframes to plot evograms and the evogram itself should be returned.
#' @param get_evogram_arg_list List containing additional arguments passed to \code{\link{get_evogram}}
#' @return List containing information get and plot Muller plots and dendrograms
#' 
#' @examples
#' \donttest{
#' ## Default is to return the plot and all info needed to create new ones
#' hal_info <- read.HAL(path_to_hal_results)
#' print(hal_info$evofreq_plot)
#' 
#' ### Can define column to use for coloring
#' hal_info <- read.HAL(path_to_hal_results, fill_name = "Passengers")
#' print(hal_info$evofreq_plot)
#' 
#' ### The information needed to create new plots is also returned. This can be useful if you want to change colors
#' updated_plot <- update_colors(hal_info$freq_frame, hal_info$clones, hal_info$attributes$Drivers)
#' plot_evofreq(updated_plot)
#' 
#' ### Can use a list to pass additional arguments to get_evofreq
#' evofreq_args <- list("threshold"=0.0, "clone_cmap"="plasma")
#' hal_info <- read.HAL(path_to_hal_results, fill_name = "Passengers", get_evofreq_arg_list = evofreq_args)
#' print(hal_info$evofreq_plot)
#' 
#' ### Use the same approach to get dendrograms by setting  get_dendogram_df = T.  Set return_dendrogram_plot = T to get the dendrogram plot too
#' evogram_args <- list("threshold"=0.0, "clone_cmap"="plasma")
#' hal_info <- read.HAL(path_to_hal_results, fill_name = "Drivers", get_evogram_arg_list = evofreq_args, get_dendogram_df=T)
#' print(hal_info$evogram_plot)
#' }
#' @export
read.HAL <- function(path_to_file, fill_name=NULL, get_evofreq_df=TRUE, get_evofreq_arg_list=list(), get_dendogram_df=FALSE,  get_evogram_arg_list=list()){
  ## FOR TESTING ###
  # path_to_file <- f
  # fill_name <- "Drivers"
  # get_evofreq_arg_list <- list("threshold"=0.1, "clone_cmap"="jet")
  # return_evofreq_plot <- TRUE
  # get_dendogram_df <- TRUE
  # get_evofreq_arg_list = evofreq_args
  # get_evogram_arg_list <- evofreq_args
  # return_dendrogram_plot = TRUE
  ###
  
  clone_df <- read.csv(path_to_file, check.names = FALSE, stringsAsFactors = FALSE)  
  df_cols <- colnames(clone_df)
  time_col_idx <- suppressWarnings(which(!is.na(as.numeric(df_cols))))
  time_cols <- df_cols[time_col_idx]
  size_df <- clone_df[, time_cols]
  
  return_list <- list("clones"=clone_df$CloneID, "parents"=clone_df$ParentID, "size_df"=size_df)
  
  hal_args <- list("size_df" = size_df, "clones"=clone_df$CloneID, "parents"=clone_df$ParentID)
  
  attribute_col_idx <- which(!df_cols %in% c(time_cols, "CloneID", "ParentID"))
  if(length(attribute_col_idx) > 0){
    attribute_df <- cbind(data.frame("CloneID"=clone_df$CloneID),clone_df[, attribute_col_idx])
    return_list[["attributes"]] <- attribute_df
    if(!is.null(fill_name)){
      hal_args[["fill_value"]] <- clone_df[fill_name]
    }
    
  }
  
  ### Get info to plot evofreq
  if(get_evofreq_df){
    freq_df <- do.call(get_evofreq, c(hal_args, get_evofreq_arg_list))  
    return_list[["freq_frame"]] <- freq_df
    return_list[["evofreq_plot"]] <- plot_evofreq(freq_df)
  }
  
  
  ### Get info to plot evogram
  if(get_dendogram_df){
    dendro_df <- do.call(get_evogram, c(hal_args, get_evogram_arg_list))
    return_list[["dendro_links"]] <- dendro_df$links
    return_list[["dendro_pos"]] <- dendro_df$dendro_pos
    return_list[["evogram_plot"]] <- plot_evogram(dendro_df$dendro_pos, dendro_df$links)
  }
  
  return(return_list) 
}


