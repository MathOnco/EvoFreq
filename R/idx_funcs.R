get_all_idx <- function(c, c_list, p_list){
  all_children_idx <- c()
  add_child_idx <- function(i, c_list){
    child_idx <- which(c_list==i)
    all_children_idx <<- append(all_children_idx, child_idx)
  }
  traverse(c, c_list, p_list, add_child_idx)
  
  return(all_children_idx)
}

get_root_id <- function(clones, parents){
  ### FOR TESTING ###
  # clones <- clone_list
  # parents <- parent_list
  ###################
  ### Root is either a clone that is also its own parent, or is not in the clone list (i.e. root's root)
  parents_not_in_clone <- unique(parents[which(! parents %in%clones)])
  if(length(parents_not_in_clone)==1){
    root_id <- parents_not_in_clone
    return(root_id)
  }else if(length(parents_not_in_clone) > 1){
    warning("More that 1 parent is not in clone list. Cannot determine root id")
  }
  clone_as_parent <- unique(parents[which(parents==clones)])
  if(length(clone_as_parent) == 1){
    root_id <- clone_as_parent
    return(root_id)
  }else if(clone_as_parent > 1){
    warning("More that 1 clone is also its own parent. Cannot determine root id")
  }
  warning("Expect ID of root to be either a parent not in the clone list, or a clone with the same ID as its parent. Neither seems to be true, and so cannot determine root id")
}

get_ancestor_idx_evo_freq <- function(clone_id, clone_list, parent_list){
  ###FOR TESTING ###
  # clone_id <- cid
  # clone_list <- clones
  # parent_list <- parents
  ####
  root_id <- get_root_id(clone_list, parent_list)
  if(root_id==clone_id){
    return(NULL)
  }
  clone_idx <- which(clone_list == clone_id)
  clone_parent <- parent_list[clone_idx]
  ancestral_idx <- c()
  while(clone_parent != root_id){
    parent_idx <- which(clone_list == clone_parent)
    ancestral_idx <- append(ancestral_idx, parent_idx)
    clone_parent <- parent_list[parent_idx]
  }
  
  return(ancestral_idx)
}

get_pos <- function(clones, parents, mut_mat, og_time_pts=NULL){
  ### FOR TESTING ###
  # parents <- parents
  # clones <- clones
  # mut_mat <- freq_mat
  # og_time_pts <- NULL #as.numeric(time_pts)
  # get_pos(clones, parents, freq_mat)
  ###############
  # roots_root_idx <- which(!parents %in% clones)
  # roots_root <- parents[roots_root_idx]
  # root_id <- clones[roots_root_idx]
  
  roots_root <- get_root_id(clones, parents)
  # root_idx <- which(parents==roots_root)
  # root_id <- clones[root_idx]
  
  # g <- graph_from_data_frame(data.frame("from"=parents, "to"=clones))
  # plot(g)
  
  if(is.null(og_time_pts)){
    og_time_pts <- seq(1, ncol(mut_mat))
  }
  
  og_mut_mat <- mut_mat
  mut_mat <- interp_mut_mat(og_mut_mat, og_time_pts)
  time_pts <- as.numeric(colnames(mut_mat))
  
  origin_times <- apply(mut_mat, 1, function(x){which(x>0)[1]})
  unique_ordered_origin_times <- sort(unique(origin_times)) ### ensure earlier clones are drawn first
  y_btm_list <- list()
  y_top_list <- list()
  
  clone_id <- roots_root
  clone_str_id <- as.character(clone_id)
  n_time_pts <- ncol(mut_mat)
  
  max_size <- max(mut_mat)
  y_btm <- rep(0, n_time_pts)
  y_top <- rep(1, n_time_pts)
  y_btm_list[[clone_str_id]] <- y_btm
  y_top_list[[clone_str_id]] <- y_top
  
  x_vals <- time_pts
  
  clone_x <- c(x_vals, rev(x_vals))
  
  clone_str_id <- as.character(clone_id)
  clone_pos_list <- list()
  draw_order <- 0
  
  # print(unique_ordered_origin_times)
  for(ot in unique_ordered_origin_times){
    # ot <- unique_ordered_origin_times[1]
    parents_at_time <- unique(parents[origin_times==ot])
    n_parents <- length(parents_at_time)
    if(n_parents > 1){
      parents_at_time <- order_clones_at_time(parents_at_time, parents = parents, clones=clones)
    }
    
    for(parent in parents_at_time){
      
      parent_str_id <- as.character(parent)
      children_idx <- which(parents==parent)
      if(length(children_idx)==1){
        total_child_area <- mut_mat[children_idx,]
      }else{
        total_child_area <- as.numeric(colSums(mut_mat[children_idx,]))
      }
      
      
      if(parent==roots_root){
        parent_area <- rep(1, n_time_pts)
        n_children <- 1
      }else{
        parent_idx <- which(clones==parent)
        parent_area <- as.numeric(mut_mat[parent_idx,])
        n_children <- length(children_idx)### TODO PUT HERE???
        # children_sizes <- mut_mat[children_idx,]
        # if(length(children_idx) > 1){
        #   # print(children_sizes)
        #   n_children <- apply(children_sizes, MARGIN = 2, function(x){length(x[x>0])})
        # }else{
        #   n_children <- rep(0, length(children_sizes))
        #   n_children[children_sizes > 0] <- 1
        #   print(n_children)
        # }
        # # print(n_children)
        
      }
      # n_children <- length(children_idx)### NOTE WAS HERE
      # n_children <- length(children_idx)
      spacing <- (parent_area - total_child_area)/(n_children+1)
      # spacing[spacing < 0] <- min(spacing[spacing >= 0]) ###TODO DOES THIS WORK?
      
      ### Get bottom of parent's polygon
      current_btm <- y_btm_list[[parent_str_id]]
      loc <- spacing + current_btm
      # print(paste("parent", parent_str_id))
      # print(spacing)
      # print(current_btm)
      # print(loc)
      ### Draw points around center
      if(n_children > 3){
        children_origin_times <- as.numeric(origin_times[children_idx])
        children_idx <- distribute_x(children_idx[order(children_origin_times)])
        
        ### First sort children_idx by time or origin?
      }
      # print(parent_str_id)
      # print(parent_str_id %in% names(y_btm_list))
      # print(clones[children_idx])
      # print(current_btm)
      
      for(cidx in children_idx){
        clone_id <-clones[cidx]
        
        clone_id_str <- as.character(clone_id)
        
        if(clone_id_str %in% names(clone_pos_list)){
          next
        }
        
        child_size <- mut_mat[cidx,]
        non_zero_idx <- which(child_size > 0)
        # print(clone_id_str)
        # print(non_zero_idx)
        # print(child_size)

        bottom <- loc
        loc <- loc + child_size
        top <- loc
        loc <- loc + spacing
        
        # if(any(is.na(bottom))){
        #   print(paste("bottom positions are NA for clone", clone_id))
        #   print(paste("bottom positions are NA for clone", clone_id, ". Previous clone", clones[children_idx[z-1]], "had top y values of:"))
        #   print(temp_loc)
        # }
        # if(any(is.na(top))){
        #   print(paste("top positions are NA for clone", clone_id))
        #   print(paste("top positions are NA for clone", clone_id, ". Previous clone", clones[children_idx[z-1]], "had top y values of:"))
        #   print(temp_loc)
        # }
        
        y_btm_list[[clone_id_str]] <- bottom
        y_top_list[[clone_id_str]] <- top
        
        clone_x <- c(x_vals[non_zero_idx], rev(x_vals[non_zero_idx]))
        clone_y <- c(bottom[non_zero_idx], rev(top[non_zero_idx]))
        clone_sizes <- c(child_size[non_zero_idx], rev(child_size[non_zero_idx]))
        clone_shape_df <- data.frame(x=clone_x, y=clone_y, clone_id = clone_id, parent= parent, origin_time=ot, draw_order = draw_order, size=clone_sizes)
        draw_order <- draw_order + 1
        clone_pos_list[[clone_id_str]] <- clone_shape_df
        
      }
    }
  }
  
  clone_pos_df <- do.call(rbind, clone_pos_list)
  
  # print(any(is.na(clone_pos_df$y)))
  clone_pos_df$draw_order <- factor(clone_pos_df$draw_order, ordered=TRUE, levels = sort(unique(clone_pos_df$draw_order), decreasing = FALSE))
  # clone_pos_df$color <- as.character(clone_pos_df$color)
  clone_pos_df <- clone_pos_df[order(clone_pos_df$draw_order),]
  clone_pos_df$x <- as.numeric(clone_pos_df$x)
  clone_pos_df$y <- as.numeric(clone_pos_df$y)
  return(clone_pos_df)
  
}

get_children <- function(c, clone_list, parent_list){
  children_idx <- which(parent_list==c)
  children <- clone_list[children_idx]
  return(children)
}

get_node_idx <- function(c, c_list, p_list){
  node_idx <- c()
  c_ids <- c()
  i <- 0
  add_idx <- function(c, c_list){
    node_idx <<- append(node_idx, i)
    c_ids <<- append(c_ids, c)
    i <<- i + 1
  }
  traverse(c, c_list, p_list, add_idx)
  
  return(list("idx"=node_idx, "vertex"=c_ids))
}

get_last_idx_above_0 <- function(x){
  
  idx <- utils::tail(which(x > 0), 1)
  if(length(idx)==0){
    ###clone was never above 0
    idx <- 0
  }
  return(idx)
}

get_last_idx_at_0 <- function(x){
  idx <- utils::tail(which(x == 0), 1)
  if(length(idx)==0){
    ### clone was never 0
    idx <- 0
  }
  return(idx)
}