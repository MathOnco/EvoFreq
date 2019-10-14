get_heirarchy_level <- function(clone_id, clone_list, parent_list){
  ###FOR TESTING ###
  # clone_list <- clones_for_d
  # parent_list <- parents_for_d
  # clone_id <- 48
  ####
  ### Assigns y-position, based on number of parents. Represents number of steps back to common ancestor
  
  ### clone 48 is at level 11. It's parent, 32 should be at level 10 but it is at level 5. OR, clone 48 should be at level 6
  
  root_id <- get_root_id(clone_list, parent_list)
  if(root_id==clone_id){
    return(NULL)
  }
  clone_idx <- which(clone_list == clone_id)
  clone_parent <- parent_list[clone_idx]
  h_level <- 0
  while(clone_parent != root_id){
    h_level <- h_level + 1
    parent_idx <- which(clone_list == clone_parent)
    clone_parent <- parent_list[parent_idx]
  }
  return(h_level)
}

get_dendrogram_pos <- function(attr_df, clones_for_d, parents_for_d){
  ### FOR TESTING
  ### NOTE attr_df needs to at least have clone ids, parent_ids, and node_ids
  # clones_for_d <- to_plot_df$clones
  # parents_for_d <- to_plot_df$parents
  # attr_df <- to_plot_df$attributes
  ###
  ### All have been sorted by origin
  h_level_df <- attr_df
  h_level_df$clone_id <- clones_for_d
  h_level_df$parents <- parents_for_d
  
  h_level_df$y <- sapply(clones_for_d, function(x){get_heirarchy_level(x, clones_for_d, parents_for_d)})
  h_level_df$tree_y <- h_level_df$y
  max_level <- max(h_level_df$y)
  leaves_idx <- which(!clones_for_d %in% parents_for_d)
  h_level_df$y[leaves_idx] <- max_level ###can put all leaves at the bottom
  h_level_df$dendro_y <- h_level_df$y ### Their position if drawn so all leaves on bottom
  
  h_level_df$x <- NA
  
  ### Get position of all leaves
  n_leaves <- length(leaves_idx)
  if(n_leaves == 1){
    leaf_x <- 0.5
  }else{
    leaf_x <- seq(0, 1, length.out = n_leaves)         
  }
  
  ### Start at top, move down. 
  ### For each node, get leaves of subtree. Center node around leaves of this subtree
  ### TODO Order leaves so that they go
  ### Keep track of where in leaf positions
  root_id <- get_root_id(clones_for_d, parents_for_d)
  parent_pos_range_list <- list()
  parent_pos_range_list[[as.character(root_id)]] <- leaf_x
  traversal_idx <- get_node_idx(root_id, c_list = h_level_df$clone_id, p_list = h_level_df$parents)
  traversal_order <- order(traversal_idx$vertex)
  for(i in traversal_order){
    ### Can go through tree using traversal order. Gurantees parent's positions will have been determined already
    if(i == 1){
      ### This is the root 
      next()
    }
    
    cid <- traversal_idx$vertex[i] #h_level_df$clone_id[traversal_idx$idx[i]]
    if(as.character(cid) %in% names(parent_pos_range_list)){next()}
    clone_parent <- h_level_df$parents[h_level_df$clone_id == cid] #h_level_df$parents[traversal_idx$idx[i]]
    parents_pos_range <- parent_pos_range_list[[as.character(clone_parent)]]
    all_nodes <- h_level_df$clone_id[which(h_level_df$parents == clone_parent)] ### Get all nodes with this parent

    pos_idx <- 0
    for(cidx in seq(1, length(all_nodes))){
      subtree_idx <- get_all_idx(all_nodes[cidx], c_list = h_level_df$clone_id, p_list = h_level_df$parents) ### includes clone, so should always have at least length 1
      n_leaves <- length(which(!h_level_df$clone_id[subtree_idx] %in% h_level_df$parents[subtree_idx]))
      next_pos_idx <- n_leaves + pos_idx
      leaf_range <- parents_pos_range[seq(pos_idx + 1, next_pos_idx)]
      
      h_level_df$x[h_level_df$clone_id == all_nodes[cidx]] <- mean(leaf_range)
      
      # print(all_nodes[cidx])
      # print(c(n_leaves, length(subtree_idx)))
      
      # print(leaf_range)
      # print(mean(leaf_range))
      # print("")
      
      pos_idx <- next_pos_idx
      parent_pos_range_list[[as.character(all_nodes[cidx])]] <- leaf_range      
      
    }
    
  }
  
  return(h_level_df)
  
}

filter_edges_at_time <- function(size_df, clones, parents, time_pt=NULL, attribute_df=NULL, threshold=0.01, data_type="size",  fill_gaps_in_size = FALSE, test_links=TRUE, rescale_after_thresholding=FALSE){
  ### FOR TESTING ###
  # clones <- clone_df$CloneID
  # parents <- clone_df$ParentID
  # size_df
  
  ####
  if(any(duplicated(clones))){
    warning("Some clones have the same ID. Each clone should have a unique ID")
  }
  if(is.null(time_pt)){
    time_pt <- tail(colnames(size_df), n=1)
  }else{
    time_pt <- as.character(time_pt)
  }
  
  if(fill_gaps_in_size){
    size_df <- as.data.frame(t(apply(size_df, 1, fill_in_gaps)))
  }
  
  ### Make sure clone does not have the same id as it's parent. If true, it can cause infinite recursion
  if(test_links){
    updated_edges <- check_and_update_edges(clones, parents)
    clones <- updated_edges$updated_clones
    parents <- updated_edges$updated_parents    
  }
  
  if(data_type=="size"){
    ### data are sizes of each clone. Here, the frequency of each mutation will be determined
    freq_df <- get_mutation_df(size_df, clones = clones, parents = parents)
    freq_mat <- as.matrix(freq_df)
    freq_mat <- sweep(freq_mat, 2, colSums(size_df), "/")
    
  }else{
    check_freq_mat(size_df, clones, parents)
    freq_mat <- size_df
    updated_edges <- check_and_update_edges(clones, parents)
    clones <- updated_edges$updated_clones
    parents <- updated_edges$updated_parents
  }
  
  
  ### Subset time points to be plotted
  if(is.null(nrow(freq_mat))){
    ###Only 1 clone
    freq_mat <- t(as.matrix(freq_mat))
  }
  
  colnames(freq_mat) <- colnames(size_df)
  row.names(freq_mat) <- clones
  freq_array <- freq_mat[,time_pt]
  origin_times <- as.numeric(apply(freq_mat, 1, function(x){which(x>0)[1]}))
  
  
  ### Possible that a clone was recorded, but didn't exist at this time point
  existed_idx <- which(freq_array > 0)
  freq_array <- freq_array[existed_idx]
  clones <- clones[existed_idx]
  parents <- parents[existed_idx]
  origin_times <- origin_times[existed_idx]
  if(! is.null(attribute_df)){
    attribute_df <- attribute_df[existed_idx, ]
  }
  
  ### Filter values based on threshold
  if(threshold > 0){
    filtered_idx <- which(freq_array >= threshold)
    freq_array <- freq_array[filtered_idx]
    parents <- parents[filtered_idx]
    clones <- clones[filtered_idx]
    origin_times <- origin_times[filtered_idx]
    
    if(rescale_after_thresholding){
      freq_mat <- rescale_frequencies(size_df, clones, parents, filtered_idx, data_type)
    }
    
    if(!is.null(attribute_df)){
      attribute_df <- attribute_df[filtered_idx, ]
    }
  }
  if(is.null(attribute_df)){
    attribute_df <- data.frame("clone_id"=clones, "parents"=parents)
  }
  
  attribute_df$freq <- freq_array
  attribute_df$origin <- origin_times
  
  ordered_idx <- order(origin_times)
  clones <- clones[ordered_idx]
  parents <- parents[ordered_idx]
  attribute_df <- attribute_df[ordered_idx, ]
  attribute_df$node_id <- seq(1, length(clones))
  
  # clone_col_idx <- as.numeric(which(sapply(colnames(attribute_df), FUN = function(x){all(clones %in% as.character(unique(attribute_df[,x])))})==TRUE))
  # 
  # if(length(clone_col_idx) > 1){
  #   #### If all clones in pos df are also parents, and there are both parent and clone ids in attribute_df, then more than 1 column will considered the clone_id column in attribute_df
  #   ### Not possible for all clones in attribute df to also all be parents 
  #   ### So clone_id column is the one that has more unique elements
  #   n_unique_clones <- apply(attribute_df[clone_col_idx], 2, function(x){length(unique(x))})
  #   clone_col_name <- names(n_unique_clones)[which(n_unique_clones==max(n_unique_clones))]
  # }
  # 
  # colnames(attribute_df)[clone_col_idx] <- "clone_id"
  
  # to_plot_df$attributes
  # clone_col_idx <- which(colnames(to_plot_df$attributes)==clone_id_col_in_att_df)
  # if(length(clone_col_idx) == 0 | clone_id_col_in_att_df != "clone_id" | ! "clone_id" %in% colnames(to_plot_df$attributes){
  #   print("do not know column in attribute_df that contains clone IDs. Please check the clone_id_col_in_att_df argument or set the clone ID column in attribute_df to be 'clone_id'")
  # }
  
  # clone_col_idx <- which(colnames(attribute_df)==clone_id_col_in_att_df)
  # if(length(clone_col_idx) == 0 & clone_id_col_in_att_df != "clone_id" & ! "clone_id" %in% colnames(attribute_df)){
  #   print("do not know column in attribute_df that contains clone IDs. Please check the clone_id_col_in_att_df argument or set the clone ID column in attribute_df to be 'clone_id'")
  # }
  # colnames(attribute_df)[clone_col_idx] <- "clone_id"
  
  
  return(list("clones"=clones, "parents"=parents, "attributes"=attribute_df))
}

get_straight_links <- function(dendro_df){
  ### FOR TESTING ###
  # dendro_df <- dendro_pos
  ####
  link_df <- dendro_df
  link_df$link_type <- "straight"
  link_df$end_x <- NA
  link_df$end_y <- NA
  link_df$end_tree_y <- NA
  link_df$end_origin <- NA
  link_df$end_dendro_y <- NA
  n_clones <- length(link_df$clone_id)
  for(cidx in seq(1, n_clones)){
    
    pid <- link_df$parents[cidx]
    p_idx <- which(link_df$clone_id == pid)
    if(length(p_idx)==0){
      ## root ##
      next()
    }
    link_df$end_x[cidx] <- link_df$x[p_idx]
    link_df$end_y[cidx] <- link_df$y[p_idx]
    link_df$end_origin[cidx] <- link_df$origin[p_idx]
    link_df$end_dendro_y[cidx] <- link_df$dendro_y[p_idx]
    link_df$end_tree_y[cidx] <- link_df$tree_y[p_idx]
  }
  return(link_df)
}

get_elbow_links <- function(dendro_df){
  ### FOR TESTING ###
  # dendro_df <- dendro_pos
  ####
  info_cols <- colnames(dendro_df)[! colnames(dendro_df) %in% c("x", "dendro_y", "tree_y", "y")]
  link_df_list <- list()
  parents <- unique(dendro_df$parents)
  for(p in parents){
    # p <- parents[2]
    pidx <- which(dendro_df$clone_id == p)
    if(length(pidx)== 0){
      ### root
      next()
    }
    
    children_idx <- which(dendro_df$parents == p)
    
    mid_y <- mean(c(dendro_df$y[pidx],  min(dendro_df$y[children_idx])))
    mid_tree_y <- mean(c(dendro_df$tree_y[pidx],  min(dendro_df$tree_y[children_idx])))
    mid_origin <- mean(c(dendro_df$origin[pidx],  min(dendro_df$origin[children_idx])))
    mid_dendro_y <- mean(c(dendro_df$dendro_y[pidx],  min(dendro_df$dendro_y[children_idx])))
    
    for(cidx in children_idx){
      ### each parent has a path: 
      ### 2) its xy
      ### 2) its same x, down to y between it and child
      ### 3) child's x, and the y between it and child
      ### 4) child's xy
      link_info <- dendro_df[cidx, info_cols]
      # mid_y <- mean(c(dendro_df$y[pidx],  dendro_df$y[cidx]))
      # mid_origin <- mean(c(dendro_df$origin[pidx],  dendro_df$origin[cidx]))
      # mid_dendro_y <- mean(c(dendro_df$dendro_y[pidx],  dendro_df$dendro_y[cidx]))
      # mid_tree_y <- mean(c(dendro_df$tree_y[pidx],  dendro_df$tree_y[cidx]))
      single_link_df <- link_info[rep(1, 4), ]
      single_link_df$x <- c(dendro_df$x[pidx], dendro_df$x[pidx], dendro_df$x[cidx], dendro_df$x[cidx])
      single_link_df$y <- c(dendro_df$y[pidx], mid_y, mid_y, dendro_df$y[cidx])
      single_link_df$origin_y <- c(dendro_df$origin[pidx], mid_origin, mid_origin, dendro_df$origin[cidx])
      single_link_df$dendro_y <- c(dendro_df$dendro_y[pidx], mid_dendro_y, mid_dendro_y, dendro_df$dendro_y[cidx])
      single_link_df$tree_y <- c(dendro_df$tree_y[pidx], mid_tree_y, mid_tree_y, dendro_df$tree_y[cidx])
      link_df_list[[as.character(cidx)]] <- single_link_df
    }
  }
  
  link_df <- do.call(rbind, link_df_list)
  link_df$link_type <- "elbow"
  # link_df$tree_y <- link_df$y
  return(link_df)
  
  
}

#'@title get_evogram
#'
#'Collect information to plot dendrogram
#'@inheritParams get_evofreq
#'@param time_pt Timepoint with which the dendrgram should be drawn
#'@param link_type Defines the shape of the edges berween nodes: "straight" draws straight lines between parents and childrend, while "elbow" draws a step from parent to child
#'@return List containing two dataframes: "dendro_pos" contains the positions of the nodes, "links" contains the positions of the edges between nodes
#' @examples
#' data("example.easy.wide")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide)))))
#' size_df <- example.easy.wide[, time_col_idx]
#' parents <- example.easy.wide$parent
#' clones <- example.easy.wide$clone
#'
#' tree_info <- get_evogram(size_df, parents = parents, clones = clones)
#' tree_pos <- tree_info$dendro_pos
#' elbow_links <- tree_info$links
#' tree_p <- plot_evogram(tree_pos, elbow_links)
#' 
#' ### Can also plot with straight links
#' tree_info <- get_evogram(size_df, parents = parents, clones = clones, link_type = "straight")
#' tree_pos <- tree_info$dendro_pos
#' straight_links <- tree_info$links
#' tree_straight_p <- plot_evogram(tree_pos, straight_links)
#' 
#' ### Can also set nodes to be colored by attribute
#' data("example.easy.wide.with.attributes")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attr_size_df <- example.easy.wide.with.attributes[, time_col_idx]
#' attr_parents <- example.easy.wide.with.attributes$parent
#' attr_clones <- example.easy.wide.with.attributes$clone
#' fitness <- example.easy.wide.with.attributes$fitness
#'
#' attribute_dendro_df <- get_evogram(attr_size_df, attr_clones, attr_parents, fill_value = fitness)
#' attribute_tree_pos <- attribute_dendro_df$dendro_pos
#' attribute_elbow_links <- attribute_dendro_df$links
#' attribute_tree_elbow_p <- plot_evogram(attribute_tree_pos, attribute_elbow_links, scale_by_node_size = TRUE)
#' @export
get_evogram <- function(size_df, clones, parents, fill_value=NULL, fill_range = NULL, time_pt=NULL, clone_cmap=NULL, threshold=0.01, data_type="size", fill_gaps_in_size = FALSE, test_links=TRUE, link_type="elbow", rescale_after_thresholding=FALSE, shuffle_colors=FALSE){
  # ## FOR TESTING ###
  # Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values#
  # data("example.easy.wide")
  # ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
  # time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
  # attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide)))))
  # size_df <-size_df #example.easy.wide[, time_col_idx]
  # parents <- clone_df$ParentID  #example.easy.wide$parent
  # clones <- clone_df$CloneID #example.easy.wide$clone
  # fill_value <- clone_df$Passengers
  # threshold <- 0.02
  # clone_cmap <- "rainbow_soft"
  
  # size_df <- hal_info$size_df
  # clones <- hal_info$clones
  # parents <- hal_info$parents
  # time_pt <- NULL
  # data_type <- "size"
  # scale_by_sizes_at_time <- FALSE
  # fill_gaps_in_size <- FALSE
  # test_links <- TRUE
  # threshold <- 0.01
  # link_type <- "elbow"
  # rescale_after_thresholding <- F
  # fill_value <- NULL
  # # ###
  

  if(!is.null(fill_value)){
    fill_name <- colnames(fill_value) ### Value was passed in  using a string to get the column, e.g. df[fill_name]
    if(is.null(fill_name)){
      paresed_fill_name <- deparse(substitute(fill_value))
      fill_name <- get_argname(paresed_fill_name)
    }
    attribute_df <- data.frame("clone_id"=clones)
    attribute_df[fill_name] <- fill_value
  }else{
    attribute_df <- NULL
    fill_name <- NULL
  }
  if(!is.null(fill_name)){
    fill_value <- attribute_df[, fill_name]
    if(is.null(fill_range)){
      fill_range <- range(fill_value, na.rm = TRUE)
    }
  }

  ### Filter info ###
  to_plot_df <- filter_edges_at_time(size_df = size_df, clones = clones, parents = parents, time_pt = time_pt, attribute_df = attribute_df,  threshold = threshold, data_type = data_type,  fill_gaps_in_size = fill_gaps_in_size, test_links=test_links, rescale_after_thresholding=rescale_after_thresholding)
  ### Get dendrogram positions and attributes ###

  dendro_pos <- get_dendrogram_pos(to_plot_df$attributes, clones_for_d = to_plot_df$clones, parents_for_d = to_plot_df$parents)
  if(!is.null(attribute_df)){
    attribute_df <- to_plot_df$attributes
  }
  
  ### At this point, origin times are based on the column idx. Needs to also include actual timepoints
  actual_time_pts <- colnames(size_df)
  dendro_pos$origin_time <- actual_time_pts[dendro_pos$origin]
  any_non_numeric_time_pts <- suppressWarnings(any(is.na(as.numeric(dendro_pos$origin_time))))
  if(! any_non_numeric_time_pts){
    dendro_pos$origin <- as.numeric(dendro_pos$origin_time)
  }
  
  dendro_pos <- update_colors(evo_freq_df = dendro_pos, clones = clones, fill_value = fill_value, clone_cmap = clone_cmap, fill_range = fill_range, fill_name=fill_name, shuffle_colors=shuffle_colors)

  ### Add line segments: Elbow or diagonal ###
  if(link_type == "straight"){
    link_df <- get_straight_links(dendro_pos)
  }else{
    link_df <- get_elbow_links(dendro_pos)
  }
  
  return(list("dendro_pos"=dendro_pos, "links"=link_df))
}


#'@title plot_evogram
#'Plot dendrogram
#'@inheritParams plot_evofreq
#'@param dendro_pos_df Position of nodes returned by \code{\link{get_evogram}} 
#'@param link_df Position of edges returned by \code{\link{get_evogram}} 
#'@param node_size Determines the size of nodes if \code{scale_by_node_size} is FALSE
#'@param scale_by_node_size Boolean defining if the size of each node should be proportional to its frequency. If FALSE, then all nodes have the same size, determined by \code{node_size}
#'@param orientation Defines orientation of the tree: "td" draws it top-down, "lr" draws it left-right
#'@param depth Determines the y-position of the nodes: "level" places the nodes based on how deeply they are nested in the tree, "origin" makes the y-position the same as the clone's origin time, "bottom" positions nodes so that all leaves have the same position, furthest from the root
#'@return a ggplot object to plot the tree
#'@examples 
#' data("example.easy.wide")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide)))))
#' 
#' size_df <- example.easy.wide[, time_col_idx]
#' parents <- example.easy.wide$parent
#' clones <- example.easy.wide$clone
#' 
#' tree_info <- get_evogram(size_df, parents = parents, clones = clones)
#' tree_pos <- tree_info$dendro_pos
#' elbow_links <- tree_info$links
#' tree_p <- plot_evogram(tree_pos, elbow_links)
#' 
#' ### Can also plot with elbow links
#' tree_info <- get_evogram(size_df, parents = parents, clones = clones, link_type = "straight")
#' tree_pos <- tree_info$dendro_pos
#' straight_links <- tree_info$links
#' tree_straight_p <- plot_evogram(tree_pos, straight_links)
#' 
#' ### Default is for y to show time of clonal origin , all leaves can be at the highest level
#' elbow_bottom_p <- plot_evogram(tree_pos, elbow_links, depth = "bottom")
#' ### Can view left to right by changing orientation argument
#' lr_tree_elbow_p <- plot_evogram(tree_pos, elbow_links, orientation = "lr")
#' #' ### Can color using attributes
#' data("example.easy.wide.with.attributes")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attr_size_df <- example.easy.wide.with.attributes[, time_col_idx]
#' attr_parents <- example.easy.wide.with.attributes$parent
#' attr_clones <- example.easy.wide.with.attributes$clone
#' fitness <- example.easy.wide.with.attributes$fitness
#' #' ### Can set color using attributes. Default colormap is viridis, but can be changed to any colormap available in the colormaps packageattribute_dendro_df <- get_dendrogram(attr_size_df, attr_clones, attr_parents, attribute_df = attribute_df, fill_name = "fitness", clone_id_col_in_att_df = clone_id_col,  clone_cmap = "magma", link_type = "elbow")
#' attribute_dendro_df <- get_evogram(attr_size_df, attr_clones, attr_parents, fill_value = fitness)
#' attribute_tree_pos <- attribute_dendro_df$dendro_pos
#' attribute_elbow_links <- attribute_dendro_df$links
#' attribute_tree_elbow_p <- plot_evogram(attribute_tree_pos, attribute_elbow_links)
#' @export
plot_evogram <- function(dendro_pos_df, link_df, fill_range=NULL, node_size=5, scale_by_node_size=TRUE, orientation="td", depth="origin"){
  ### FOR TESTING ### 
  # dendro_pos_df <- hal_dendro_df$dendro_pos
  # link_df <- hal_dendro_df$links
  # scale_by_node_size = TRUE
  # orientation="td" #= top_down, "lr" left right
  # depth <- "origin"
  # fill_range <- NULL
  ###
  
  link_df <- link_df[complete.cases(link_df$x),]
  link_type <- unique(link_df$link_type)
  
  if(depth == "origin"){
    y_name <- "origin_y"
    end_y_name <- "end_origin"
    dendro_pos_df$origin_y <- dendro_pos_df$origin
    time_levels <- as.numeric(as.character(unique(dendro_pos_df$origin)))
    time_labels <- as.character(unique(dendro_pos_df$origin_time))
    if(link_type=="straight"){
      link_df$origin_y <- link_df$origin      
    }
  }else if(depth=="level"){
    y_name <- "tree_y"
    end_y_name <- "end_tree_y"
  }else if(depth=="bottom"){
    y_name <- "dendro_y"
    end_y_name <- "end_dendro_y"
  }

  color_attribute_name <- unique(dendro_pos_df$efp_color_attribute)
  if(is.na(color_attribute_name)){
    color_attribute_name <- "plot_color"
  }else{
    if(is.null(fill_range)){
      fill_range <- range(dendro_pos_df[, color_attribute_name], na.rm = TRUE)  
    }
  }
  # unique(link_df$clone_id)
  # unique(dendro_pos_df$clone_id)
  if(link_type == "straight"){
    p <- ggplot2::ggplot(dendro_pos_df, ggplot2::aes_string(x="x", y=y_name, color=color_attribute_name, size="freq")) +
      ggplot2::geom_segment(data=link_df, ggplot2::aes_string(x="x", y=y_name, xend="end_x", yend=end_y_name), inherit.aes = FALSE) #+
      # ggplot2::scale_color_identity()
    
  }else{
    p <-  ggplot2::ggplot(dendro_pos_df, ggplot2::aes_string(x="x", y=y_name, color=color_attribute_name, size="freq")) +
      ggplot2::geom_path(data=link_df, ggplot2::aes_string(x="x", y=y_name, group="clone_id"), inherit.aes = FALSE) #+
      # ggplot2::scale_color_identity()
  }
  
  if(color_attribute_name=="plot_color"){
    p <- p + ggplot2::scale_color_identity()
  }else{
    colormap_name <- unique(dendro_pos_df$cmap)
    p <- p + colormap::scale_color_colormap(color_attribute_name, colormap=colormap_name, limits=fill_range)
    # p <- p + ggplot2::scale_color_gradientn(colours = colorbar_colors)  #colormap::scale_fill_colormap(color_attribute_name, colormap=colormap_name)
  }
  
  
  if(depth=="origin"){
    y_title <- "Time"  
  }else{
    y_title <- "Level"  
  }
  
  
  if(scale_by_node_size){
    p <- p + ggplot2::geom_point() 
  }else{
    p <- p + ggplot2::geom_point(size=node_size) 
  }
  
  p <- p + 
    ggplot2::theme_classic()  + 
    ggplot2::ylab(y_title)
  
  if(orientation=="lr"){
    if(depth=="origin"){
      p <- p + ggplot2::scale_y_continuous(breaks=time_levels, labels=time_labels) 
    }
    p <- p + ggplot2::coord_flip() + 
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    
  }else{
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank())
    if(depth=="origin"){
       p <- p + ggplot2::scale_y_reverse(breaks=time_levels, labels=time_labels) 
    }else{
      p <- p + ggplot2::scale_y_reverse() 
    }
  }
  
  return(p)
}