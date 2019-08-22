#'@title animate_evogram
#'
#'View changes in the dendrogram over time. Has most of the same options as \code{\link{plot_evogram}}, but does not have option to set link type
#'@inheritParams get_evofreq
#'@inheritParams plot_evogram
#'@return List infromation need to craete the animation: "dendro_pos_df" contains all of the positions of the nodes and edges at each timestep, while "animation_plot" that can be used with gganimate to create the animation. For details on how to customize and save the animation, see \code{\link[gganimate]{gganimate}}
#'@examples
#' data("example.easy.wide.with.attributes")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' size_df <- example.easy.wide.with.attributes[, time_col_idx]
#' parents <- example.easy.wide.with.attributes$parent
#' clones <- example.easy.wide.with.attributes$clone
#' fitness <- example.easy.wide.with.attributes$fitness
#' dendro_movie_info <- animate_evogram(size_df, clones = clones, parents = parents, fill_value = fitness)
#' dendro_movie_plot <- dendro_movie_info$animation_plot
#'
#'### Add gganimate object to dendro_movie_plot to finish animation
#' \donttest{
#' library(gganimate)
#' movie_p <- dendro_movie_plot +
#'   transition_time(time) +
#'   ease_aes() +
#'   exit_shrink() +
#'   ggtitle('Time {frame_time}')
#' 
#' print(movie_p)
#'}
#'@export
animate_evogram <- function(size_df, clones, parents, fill_value=NULL, time_pts=NULL, attribute_df=NULL, clone_cmap=NULL, threshold=0.01, data_type="size", fill_gaps_in_size = F, test_links=T, fill_range = NULL, node_size=5, scale_by_node_size=T, orientation="td", depth="origin"){
  # ## FOR TESTING ###
  # data("example.easy.wide.with.attributes")
  # # Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
  # time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
  # attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
  # attribute_df <- example.easy.wide.with.attributes[, attribute_col_idx]
  # size_df <- example.easy.wide.with.attributes[, time_col_idx]
  # parents <- example.easy.wide.with.attributes$parent
  # clones <- example.easy.wide.with.attributes$clone
  # fill_value <- example.easy.wide.with.attributes$fitness
  #
  # clone_cmap <- NULL
  # size_df <- size_df
  # threshold <- 0.01
  # clones <- clones
  # parents <- parents
  # time_pts <- NULL
  # attribute_val_name <- "fitness" #"new_antigenicity"
  # attribute_df <- attribute_df
  # data_type <- "size"
  # scale_by_sizes_at_time <- F
  # interpolation_steps <- 10
  # attribute_cmap <- "viridis"
  # fill_gaps_in_size <- F
  # test_links <- T
  # threshold <- 0.01
  # attribute_val_range <- NULL
  # 
  # node_size <- 5
  # scale_by_node_size <- T
  # orientation <- "td"
  # clone_id_col_in_att_df <- "clone"
  #####
  
  if(!is.null(fill_value)){
    attribute_val_name <- deparse(substitute(fill_value))
    if(grepl("\\$", attribute_val_name)){
      ### Value was passed in as a column in a dataframe
      attribute_val_name <- strsplit(attribute_val_name, split = "$", fixed = T)[[1]][2]
    }else if(grepl("\\[.*\\]", attribute_val_name)){
      ### Value was passed in  using a string to get the column, e.g. df[attribute_val_name]
      attribute_val_name <- strsplit(attribute_val_name, '\\"')[[1]][2]
    }
    
    attribute_df <- data.frame("clone_id"=clones, "parents"=parents)
    attribute_df[attribute_val_name] <- fill_value
  }else{
    attribute_df <- NULL
    attribute_val_name <- NULL
  }
  
  if(!is.null(attribute_val_name)){
    attribute_vals <- attribute_df[,attribute_val_name]
    if(is.null(fill_range)){
      fill_range <- range(attribute_vals, na.rm = T)
    }
  }
  
  
  if(is.null(time_pts)){
    time_pts <- colnames(size_df)
  }
  
  all_time_df_list <- list()
  
  cat("\n")
  print("Making Movie Frames")
  pb <- utils::txtProgressBar(min = 1, max = length(time_pts), style = 3)
  for(time_idx in seq(length(time_pts))){
    utils::setTxtProgressBar(pb, time_idx)
    tp <- time_pts[time_idx]
    sink("/dev/null") ### Suppress output
    # to_plot_df <- filter_edges_at_time(size_df = size_df, clones = clones, parents = parents, time_pt = tp, attribute_df = attribute_df, clone_id_col_in_att_df=clone_id_col_in_att_df, threshold = threshold, data_type = data_type,  fill_gaps_in_size = fill_gaps_in_size, test_links=test_links)
    to_plot_df <- filter_edges_at_time(size_df = size_df, clones = clones, parents = parents, time_pt = tp, attribute_df = attribute_df,  threshold = threshold, data_type = data_type,  fill_gaps_in_size = fill_gaps_in_size, test_links=test_links)
    sink()
    time_dendro_pos <- get_dendrogram_pos(to_plot_df$attributes, clones_for_d = to_plot_df$clones, parents_for_d = to_plot_df$parents)
    time_dendro_pos <- get_straight_links(time_dendro_pos)
    time_dendro_pos$time <- time_idx
    
    all_time_df_list[[as.character(tp)]] <- time_dendro_pos

    if(time_idx==1){next()}
    
    prev_pos_df <- all_time_df_list[[as.character(time_pts[time_idx-1])]]
    
    new_clone_idx <- which(!time_dendro_pos$clone_id %in% prev_pos_df$clone_id)
    
    if(length(new_clone_idx)>0){
      for(cidx in new_clone_idx){
        new_clones_parent <- time_dendro_pos$parent[cidx]
        parent_idx <- which(prev_pos_df$clone_id == new_clones_parent)
        new_clone_df <- time_dendro_pos[cidx, ]
        new_clone_df$time <- prev_pos_df$time[parent_idx]
        new_clone_df$freq <- 0
        new_clone_df[c("y", "tree_y", "dendro_y", "x", "time", "origin")] <- prev_pos_df[parent_idx, c("y", "tree_y", "dendro_y", "x", "time", "origin")]
        new_clone_df$end_x <- new_clone_df$x
        new_clone_df$end_y <- new_clone_df$y
        new_clone_df$end_tree_y <- new_clone_df$tree_y
        new_clone_df$end_origin <- new_clone_df$origin
        new_clone_df$end_dendro_y <- new_clone_df$dendro_y
        prev_pos_df[parent_idx,]
        # new_clone_df[c("end_x", "end_y", "end_tree_y", "end_origin", "end_dendro_y")] <- NA
        prev_pos_df <- rbind(prev_pos_df, new_clone_df)
      }
      all_time_df_list[[as.character(time_pts[time_idx-1])]] <- prev_pos_df
    }
  }
  
  
  d_pos_df <- do.call(rbind, all_time_df_list)
  d_pos_df <- d_pos_df[order(d_pos_df$time), ]
  
  # if(is.null(attribute_val_name)){
  #   attribute_df <- NULL 
  # }
  # 
  # d_pos_df <- update_colors(evo_freq_df = d_pos_df, attribute_df = attribute_df, attribute_val_name = attribute_val_name, clone_id_col_in_att_df=clone_id_col_in_att_df, clone_cmap = clone_cmap, attribute_range = attribute_val_range)
  length(fill_value)
  d_pos_df <- update_colors(evo_freq_df = d_pos_df, clones = clones, fill_value = fill_value, clone_cmap = clone_cmap, fill_range = fill_range, attribute_val_name=attribute_val_name)
  
  
  
  if(depth == "origin"){
    y_name <- "origin_y"
    end_y_name <- "end_origin"
    d_pos_df$origin_y <- d_pos_df$origin
  }else if(depth=="level"){
    y_name <- "tree_y"
    end_y_name <- "end_tree_y"
  }else if(depth=="bottom"){
    y_name <- "dendro_y"
    end_y_name <- "end_dendro_y"
  }
  
  color_attribute_name <- unique(d_pos_df$efp_color_attribute)
  if(is.na(color_attribute_name)){
    color_attribute_name <- "plot_color"
  }else{
    if(is.null(fill_range)){
      fill_range <- range(d_pos_df[, color_attribute_name], na.rm = T)  
    }
  }
  
  
  p <- ggplot2::ggplot(d_pos_df, ggplot2::aes_string(x="x", y=y_name, color=color_attribute_name, group="clone_id", xend="end_x", yend=end_y_name, size="freq")) +
    ggplot2::geom_segment(color="black", size=1, na.rm=T)
  
  if(scale_by_node_size){
    p <- p + ggplot2::geom_point() 
  }else{
    p <- p + ggplot2::geom_point(size=node_size) 
  }
  if(orientation=="lr"){
    p <- p + ggplot2::coord_flip()
  }else{
    p <- p + ggplot2::scale_y_reverse()
  }
  
  if(color_attribute_name=="plot_color"){
    p <- p + ggplot2::scale_color_identity()
  }else{
    # p <- p + ggplot2::scale_color_gradientn(colours = colorbar_colors)  #colormap::scale_fill_colormap(color_attribute_name, colormap=colormap_name)
    colormap_name <- unique(d_pos_df$cmap)
    p <- p + colormap::scale_color_colormap(color_attribute_name, colormap=colormap_name, limits=fill_range)
  }
  
  p <- p + ggplot2::theme_void()
  
  return(list("dendro_pos_df"=d_pos_df, "animation_plot"=p))
}

add_links <- function(d_pos_list){
  ###FOR TESTING ###
  # d_pos_list <- d_pos_df
  ####
  
  for(tidx in names(d_pos_list)){
    t_df <- d_pos_list[[tidx]]
    t_df$end_x <- NA
    t_df$end_y <- NA
    t_df$end_origin <- NA
    t_df$end_dendro_y <- NA
    t_df$end_tree_y <- NA
    n_clones <- length(t_df$clone_id)
    for(cidx in seq(1, n_clones)){
      
      pid <- t_df$parents[cidx]
      p_idx <- which(t_df$clone_id == pid)
      if(length(p_idx)==0){
        ## root ##
        next()
      }
      t_df$end_x[cidx] <- t_df$x[p_idx]
      t_df$end_y[cidx] <- t_df$tree_y[p_idx]
      t_df$end_origin[cidx] <- t_df$origin[p_idx]
      t_df$end_dendro_y[cidx] <- t_df$dendro_y[p_idx]
      t_df$end_tree_y[cidx] <- t_df$tree_y[p_idx]
    }
    d_pos_list[[tidx]] <- t_df
  }
  
  return(d_pos_list)
}

set_clone_origin_d_pos <- function(d_pos_list){
  ### FOR TESTING ###
  # d_pos_list <- recentered_d_pos_list
  ###
  ### Give children parent's position in frame before origin
  time_pts <- names(d_pos_list)
  n_time_pts <- length(time_pts)
  
  for(t_idx in seq(2, n_time_pts)){
    # t_idx <- 5
    current_time <- time_pts[t_idx]
    current_d <- d_pos_list[[current_time]]
    current_d <- current_d[order(current_d$x, decreasing = T), ]
    current_d <- current_d[order(current_d$node_id), ]
    ### prev time df needs to be based on updated prev, so that positions are consistent
    prev_time <- time_pts[t_idx - 1]
    prev_pos_df <- d_pos_list[[prev_time]]
    prev_pos_df <- prev_pos_df[order(prev_pos_df$node_id), ]
    # og_prev_d <- d_pos_list[[prev_time]]
    prev_clones <- prev_pos_df$clone_id[prev_pos_df$present_at_time == T]
    
    # available_pos_df <- current_d
    n_pos <- nrow(current_d)
    n_levels <- max(current_d$tree_y)
    for(l in seq(1, n_levels)){
      ### Clone is new: Add to updated previous position list, using parent's original position. Give it a left over postion in current time point    
      new_clone_idx <- which(!current_d$clone_id %in% prev_clones & current_d$tree_y == l)
      for(ncidx in new_clone_idx){
        ### For previous timestep, give clone its parent's position 
        pid <- current_d$parents[ncidx]
        prev_pidx <- which(prev_pos_df$clone_id==pid)
        new_row <- current_d[ncidx, ]
        new_row$x <- prev_pos_df$x[prev_pidx]
        new_row$tree_y <- prev_pos_df$tree_y[prev_pidx]
        new_row$y <- prev_pos_df$y[prev_pidx]
        new_row$origin <- prev_pos_df$origin[prev_pidx]
        new_row$dendro_y <- prev_pos_df$dendro_y[prev_pidx]
        new_row$time <- prev_pos_df$time[prev_pidx]
        new_row$present_at_time <- F
        prev_pos_df <- rbind(prev_pos_df, new_row)
        
      }
    }
    d_pos_list[[prev_time]] <- prev_pos_df
  }
  return(d_pos_list)
}

recenter_parents <- function(d_pos_list){
  ##FOR TESTING ###
  # d_pos_list <- matched_d_pos_list
  ####
  
  for(tidx in names(d_pos_list)){
    t_df <- d_pos_list[[tidx]] 
    n_levels <- max(t_df$tree_y)
    rev_levels <- seq(n_levels-1, 1)
    for(l in rev_levels){
      p_in_level_idx <- which(t_df$y== l)
      for(pidx in p_in_level_idx){
        p <- t_df$clone_id[pidx]
        children_idx <- which(t_df$parents==p)
        if(length(children_idx)>0){
          children_pos <- t_df$x[children_idx]
          t_df$x[pidx] <- mean(children_pos)        
        }
      }
    }
    d_pos_list[[tidx]] <- t_df
    # tp2 <- ggplot(t_df, aes(x=x, y=-tree_y, color=as.factor(clone), group=clone, size=size)) +
    #   # geom_segment(color="black", size=1) +
    #   geom_point() +
    #   theme_void()
    # 
    # f <- paste0(t_idx, "_final.png")
    # ggsave(f, tp2)
    
  }
  return(d_pos_list)
}

match_d_pos <- function(d_pos_list){
  ### FOR TESTING ###
  # d_pos_list <- intital_d_pos_list
  ###
  
  ### Keep each node near its previous position
  
  ### 3 Scenarios:
  ### 1) Clone existed in previous time step. Update it's position based on updated previous positions
  ### 2) Clone is new: Add to updated previous position list, using parent's original position. Give it a left over postion in current time point
  updated_df_h_list <- list()
  time_pts <- names(d_pos_list)
  n_time_pts <- length(time_pts)
  updated_df_h_list[[time_pts[1]]] <- d_pos_list[[time_pts[1]]]
  # for(t_idx in seq(2, 4)){
  for(t_idx in seq(2, n_time_pts)){
    # t_idx <- 5
    current_time <- time_pts[t_idx]
    current_d <- d_pos_list[[current_time]]
    current_d <- current_d[order(current_d$x, decreasing = T), ]
    current_d <- current_d[order(current_d$node_id), ]
    ### prev time df needs to be based on updated prev, so that positions are consistent
    prev_time <- time_pts[t_idx - 1]
    updated_prev_pos_df <- updated_df_h_list[[prev_time]]
    updated_prev_pos_df <- updated_prev_pos_df[order(updated_prev_pos_df$node_id), ]
    prev_clones <- updated_prev_pos_df$clone_id[updated_prev_pos_df$present_at_time == T]
    
    available_pos_df <- current_d
    n_pos <- nrow(available_pos_df)
    n_levels <- max(current_d$tree_y)
    for(l in seq(1, n_levels)){
      # l <- 1
      #### update from top down. Will allow to adjust children based on parent's new position
      extant_clone_idx <- which(current_d$clone_id %in% prev_clones & current_d$tree_y == l)
      idx_in_prev_df <- which(updated_prev_pos_df$clone_id %in% current_d$clone_id[extant_clone_idx])
      
      ### Avoid big jumps by taking care of worst case scenarios first (largest jumps possible)
      pos_dmat <- abs(outer(updated_prev_pos_df$x[idx_in_prev_df], available_pos_df$x[available_pos_df$tree_y == l], "-"))
      max_d <- apply(pos_dmat, 1, max)
      extant_clone_idx <- extant_clone_idx[order(max_d, decreasing = T)]
      ### Clone existed in previous time step. Change its position based on its updated previous positions
      for(ecidx in extant_clone_idx){
        ecid <- current_d$clone_id[ecidx]
        
        
        idx_in_updated_prev_pos <- which(updated_prev_pos_df$clone_id==ecid)
        prev_x <- updated_prev_pos_df$x[idx_in_updated_prev_pos]
        prev_y <- updated_prev_pos_df$tree_y[idx_in_updated_prev_pos]
        dx <- rep(NA, n_pos)
        for(i in seq(1, n_pos)){
          if(available_pos_df$pos_available[i] & available_pos_df$tree_y[i]==l){
            ###Only consider available spaces in the same heirarchy level
            dx[i] <- abs(available_pos_df$x[i] - prev_x)
          }
        }
        closest_idx <- which.min(dx)[1]
        current_d$x[ecidx] <- available_pos_df$x[closest_idx]
        current_d$pos_available[closest_idx] <- F
        available_pos_df$pos_available[closest_idx] <- F
      }  
      
      ### Clone is new: Add to updated previous position list, using parent's original position. Give it a left over postion in current time point    
      # new_clone_idx <- which(!current_d$clone %in% og_prev_d$clone)
      new_clone_idx <- which(!current_d$clone_id %in% prev_clones & current_d$tree_y == l)
      for(ncidx in new_clone_idx){
        ### Give clone one of the left over positions ###
        available_pos <- which(available_pos_df$pos_available & available_pos_df$tree_y == l)[1]
        current_d$x[ncidx] <- available_pos_df$x[available_pos]
        current_d$y[ncidx] <- available_pos_df$y[available_pos]
        current_d$tree_y[ncidx] <- available_pos_df$tree_y[available_pos]
        current_d$pos_available[available_pos] <- F
        available_pos_df$pos_available[available_pos] <- F
      }
    }
    
    # any(duplicated(current_d[c("x","y")]))
    updated_df_h_list[[prev_time]] <- updated_prev_pos_df
    updated_df_h_list[[current_time]] <- current_d
    # 
    # tp2 <- ggplot(current_d, aes(x=x, y=-tree_y, color=as.factor(clone), group=clone, size=size)) +
    #   # geom_segment(color="black", size=1) +
    #   geom_point() +
    #   theme_void()
    # 
    # f <- paste0(t_idx, "_updated.png")
    # ggsave(f, tp2)
  }
  return(updated_df_h_list)
}

get_initial_d_pos <- function(mut_freq, clones, parents, attribute_df, attribute_val_name=NULL){
  ### Get positions of nodes for each dendrogram
  ### FOR TESTING ###
  # clones <- to_plot_df$clones
  # parents <- to_plot_df$parents
  # mut_freq <- to_plot_df$freq_mat
  # attribute_df <- to_plot_df$attributes
  # attribute_val_name <- "fitness"
  ###
  
  attribute_df$clone_id <- clones
  attribute_df$parents <- parents
  time_pts <- colnames(mut_freq)
  n_time_pts <- length(time_pts)
  
  first_time_dendro <- attribute_df[1, ]  ### Ordered by origin time, so root is a top
  first_time_dendro$x <- 0.5
  first_time_dendro$y <- 0
  first_time_dendro$tree_y <- 0
  first_time_dendro$dendro_y <- 0
  first_time_dendro$time <- as.numeric(time_pts[1])
  first_time_dendro$present_at_time <- T
  first_time_dendro$pos_available <- T
  first_time_dendro$freq <- mut_freq[1, 1]
  
  # first_time_dendro <- update_colors(first_time_dendro, attribute_df = first_time_dendro, attribute_val_name = attribute_val_name)
  
  attribute_df <- attribute_df[order(attribute_df$origin), ]
  
  pos_df_list <- list()
  pos_df_list[[as.character(first_time_dendro$time)]] <- first_time_dendro
  for(time_idx in seq(2, n_time_pts)){
    # time_idx <- 2
    freq_at_time <- as.numeric(mut_freq[, time_idx])
    present_idx <- which(freq_at_time > 0)
    at_time_df <- attribute_df[present_idx, ]
    at_time_df$freq <- freq_at_time[present_idx]
    at_time_dendro <- get_dendrogram_pos(at_time_df, at_time_df$clone_id, at_time_df$parents)
    at_time_dendro[1,]$x <- 0.5 ### Keep root in middle
    at_time_dendro$time <- as.numeric(time_pts[time_idx])
    at_time_dendro$present_at_time <- T
    at_time_dendro$pos_available <- T
    # at_time_dendro <- update_colors(at_time_dendro, attribute_df = at_time_dendro, attribute_val_name = attribute_val_name)
    pos_df_list[[as.character(unique(at_time_dendro$time))]] <- at_time_dendro
    
    # tp <- ggplot(at_time_dendro, aes(x=x, y=-tree_y, color=as.factor(clone), group=clone, size=size)) +
    #   # geom_segment(color="black", size=1) +
    #   geom_point() +
    #   theme_void()
    # 
    # f <- paste0(time_idx, ".png")
    # ggsave(f, tp)
  }
  
  return(pos_df_list)
  
}
