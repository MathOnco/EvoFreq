
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





