#'@import colormap ggplot2 dplyr bezier

#'@title get_evofreq
#'
#'  Collect information to plot frequency dynamics
#'  
#'@param size_df Dataframe in a wide format, where each row corrsponds to a single clone, and the columns are the sizes of that clone at each timepoint
#'@param clones Array containing the clone ids. The index of each clone must correspond to the same index of the row in \code{size_df} that contains the sizes of that clone over time
#'@param parents Array containing the ids of the parent of each clone in the \code{clones} array.
#'@param fill_value Array containing information that can be used to color each clone. If NULL (the default), each clone is assigned a color. If values are a clone attribute, e.g. fitness, then the colors are generated using  \code{\link[colormap]{colormaps}}. The user can also provide custom colors in 3 ways: 1) hexcode; 2) rgb values as a string, with each value being a the intensity of the color channel, each separated by commas, e.g. "255, 10, 128"; 3) Any of the named in colors in R, which can be found with \code{\link[grDevices]{colors}}
#'@param fill_range Array containing the minimum and maximum values to set the range of colors. If NULL (the default), the range is determined directly from \code{fill_value}.
#'@param time_pts Array containing the name of the timepoints. If NULL, then the name of timepoints will be a sequence from 1 to the number of columns in \code{size_df}.
#'@param clone_cmap Colormap to use for the clones. For a list of available colormaps, see \code{\link[colormap]{colormaps}}.
#'@param threshold The minimum frequency of clones to be plotted. Clones with with a frequency below this value will not be plotted
#'@param scale_by_sizes_at_time Boolean defining whether or not the plot should represent the size or frequency of each clone at each timepoint. If TRUE, the sizes are scaled by the maximum size at each timepoint, and the plot thus represents the clonal frequencies at each timepoint. If FALSE, the sizes are scaled using the maximum size in \code{size_df}, thus reflecting relative population sizes
#'@param data_type String defining what kind of information is in size_df. If "size", then the values in \code{size_df} are the population sizes. If "mutation", the values are the frequencies, between 0 and 1, of each mutation in the population over time
#'@param interpolation_steps Integer defining the number of knots to use in the spline interpolation used to fill in the gaps between observed population sizes. For sparse data, this smooths out the curves in the plot. Not recommended if the data is dense, as this is slow and may not have noticable effects
#'@param interp_method String identifying the interpolation method to use. Either "bezier", or a method used by \code{\link[stats]{splinefun}}
#'@param fill_gaps_in_size Boolean defining whether or not missing sizes should be filled in
#'@param test_links Make sure clone does not have the same id as it's parent. If true, it can cause infinite recursion. 
#'@param add_origin Boolean defining whether or not to add origin positions to founder clones, even if not present in the data. Best for sparse observed data
#'@param tm_frac Value between 0 and 1 that determines where the maximum growth rate is in the inferred origin sizes. Lower values result in earlier maximum growth
#'@param rescale_after_thresholding Boolean determining if frequencies should be rescaled after thresholding, so that frequencies are based on what was above the threshold.
#'@param shuffle_colors Boolean determining if colors should be shuffled before being assigned to each clone. Only applies when fill_value = NULL
#'@return Formatted dataframe called a "freq_frame" containing the information needed to plot the frequency dynamics over time.
#'
#'@examples
#'\donttest{
#' data("example.easy.wide")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
#' size_df <- example.easy.wide[, time_col_idx]
#' parents <- example.easy.wide$parents
#' clones <- example.easy.wide$clones
#' 
#' ### Default is to plot size
#' freq_frame <- get_evofreq(size_df, clones, parents)
#' evo_p_by_size <- plot_evofreq(freq_frame)
#' 
#' ### Can also plot frequency by setting scale_by_sizes_at_time = TRUE.
#' freq_frame <- get_evofreq(size_df, clones, parents, scale_by_sizes_at_time = TRUE)
#' evo_p_by_freq <- plot_evofreq(freq_frame)
#' 
#' ### Default is to mildly smooth corners, but this can be turned by setting interpolation_steps = 0
#' freq_frame <- get_evofreq(size_df, clones, parents, interpolation_steps = 0)
#' raw_evo_p <- plot_evofreq(freq_frame)
#' 
#' ### Several other methods to smooth corners, including using Bezier curves. However, Bezier curves dont represent the data as accurately as the methods that use splinefun, i.e. c("fmm", "periodic", "natural", "monoH.FC", "hyman")
#' freq_frame <- get_evofreq(size_df, clones, parents, interp_method = "bezier")
#' bez_evo_p <- plot_evofreq(freq_frame)
#' 
#' ### Data can also be provided as mutaiton frequencies by setting data_type = "mutation"
#' mutation_count_df <- get_mutation_df(size_df, clones, parents)
#' freq_frame <- get_evofreq(mutation_count_df, clones, parents, data_type = "mutation")
#' evo_p_from_mutation <- plot_evofreq(freq_frame)
#' 
#' ### Input needs to be in wide format, but can be converted to long format data to wide format using \code{\link{long_to_wide_freqframe}}
#' wide_df_info <- long_to_wide_freqframe(long_pop_sizes_df = example.easy.long.sizes, time_col_name = "Time", clone_col_name = "clone", parent_col_name = "parent", size_col_name = "Size", edges_df = example.easy.long.edges)
#' clones_from_long <- wide_df_info$clones
#' parents_from_long <- wide_df_info$parents
#' size_df_from_long <- wide_df_info$wide_size_df
#' freq_frame <- get_evofreq(size_df_from_long, clones_from_long, parents_from_long)
#' evo_p_from_long <- plot_evofreq(freq_frame)
#'
#' ### Setting of colors can be done when getting the freq_frame, or by updating the color later using \code{\link{update_colors}}. Available colormaps are those found in \code{\link[colormap]{colormaps}}
#' ### Default colormap is rainbow_soft, but this can be changed using the \code{clone_cmap} argument. 
#' jet_freq_frame <- get_evofreq(size_df, clones, parents, clone_cmap = "jet")
#' jet_evo_p <- plot_evofreq(jet_freq_frame)
#' 
#' ### Can color each clone by an attribute by providing a \code{fill_value}. Default colormap is viridis, but this can be changed using the \code{clone_cmap} argument
#' fitness <- runif(length(clones))
#' fitness_freq_frame <- get_evofreq(size_df, clones, parents, fill_value = fitness)
#' fitness_evo_p <- plot_evofreq(fitness_freq_frame)
#' 
#' ### The user can also provide custom colors for each clone, which will need to be passed into the \code{fill_value} argument
#' ### Custom colors can be defined using RGB values. Each color should be a string specifying the color channel values, separated by commas.
#' rgb_clone_colors <- sapply(seq(1, length(clones)), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
#' rgb_freq_frame <- get_evofreq(size_df, clones, parents, rgb_clone_colors)
#' rgb_evo_p <- plot_evofreq(rgb_freq_frame)
#' 
#' ### Custom colors can also be any of the named colors in R. A list of the colors can be found with \code{colors()}
#' named_clone_colors <- sample(colors(), length(clones), replace = FALSE)
#' named_freq_frame <- update_colors(rgb_freq_frame, clones = clones, fill_value = named_clone_colors)
#' named_evo_p <- plot_evofreq(named_freq_frame)
#' 
#' ### Custom colors can also be specified using hexcode
#' hex_clone_colors <- sample(colormap::colormap(colormap=colormaps$temperature, nshades=length(clones)))
#' hex_freq_frame <- update_colors(rgb_freq_frame, clones = clones, fill_value = hex_clone_colors)
#' hex_evo_p <- plot_evofreq(hex_freq_frame)
#'
#' ### Can revert back to original colors
#'freq_frame_default_color <- update_colors(fitness_freq_frame, clones=clones)
#'default_cmap_evo_p <- plot_evofreq(freq_frame_default_color)
#'}
#'@export
get_evofreq <- function(size_df, clones, parents, fill_value=NULL, fill_range = NULL, time_pts=NULL, clone_cmap=NULL, threshold=0.01, scale_by_sizes_at_time = FALSE, data_type="size", interpolation_steps = 20, interp_method = "monoH.FC", fill_gaps_in_size = FALSE, test_links=TRUE, add_origin=FALSE, tm_frac=0.6, rescale_after_thresholding=FALSE, shuffle_colors=FALSE){
  # # ## FOR TESTING ###
  # data("example.easy.wide.with.attributes")
  # ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
  # time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
  # size_df <- example.easy.wide[, time_col_idx]
  # parents <- example.easy.wide$parents
  # clones <- example.easy.wide$clones
  # fill_value <- clone_attribute_colors
  
  # parents <- clone_df$Parent
  # clones <- row.names(clone_df)
  # time_pts <- as.numeric(colnames(clone_df))
  # time_pts <- which(!is.na(time_pts))
  # size_df <- clone_df[time_pts]
  # fill_value <- clone_df$Drivers
  
  ### HAL
  # size_df <- hal_info$size_df
  # clones <- hal_info$clones
  # parents <- hal_info$parents
  # fill_value <- NULL
  
  # threshold <- 0.01
  # clone_cmap <- NULL
  # time_pts <- NULL
  # fill_range <-  NULL
  # scale_by_sizes_at_time <- FALSE
  # interpolation_steps <- 10
  # fill_gaps_in_size <- FALSE
  # test_links <- TRUE
  # data_type <- "size"
  # interp_method <- "monoH.FC"# "bezier"
  # add_origin <- FALSE
  # tm_frac <- 0.6
  # rescale_after_thresholding <- FALSE 
  # shuffle_colors <- FALSE
  # # # ###
  if(!is.null(fill_value)){
    fill_name <- colnames(fill_value) ### Value was passed in as a single column dataframe
    if(is.null(fill_name)){
    paresed_fill_name <- deparse(substitute(fill_value))
    print(paresed_fill_name)
    fill_name <- get_argname(paresed_fill_name)
    print(fill_name)
    }
    
    attribute_df <- data.frame("clone_id"=clones)
    attribute_df[fill_name] <- fill_value
  }else{
    attribute_df <- NULL
    fill_name <- NULL
  }

  og_time_pts <- colnames(size_df)
  to_plot_df <- filter_data(size_df = size_df, clones = clones, parents = parents, time_pts = time_pts, attribute_df = attribute_df, threshold = threshold, scale_by_sizes_at_time = scale_by_sizes_at_time, data_type = data_type,  fill_gaps_in_size = fill_gaps_in_size, test_links=test_links, add_origin=add_origin, tm_frac=tm_frac, rescale_after_thresholding=rescale_after_thresholding)
  clones <- to_plot_df$clones
  parents <- to_plot_df$parents
  freq_mat <- to_plot_df$freq_mat
  max_mutation_size <- to_plot_df$max_size
  time_pt_names <- to_plot_df$og_colnames
  time_pt_df <- data.frame("x"=as.numeric(colnames(freq_mat)), "Time_label"=time_pt_names) ### Have to make x numeric. Otherwise, float column names are converted to strings

  if(!is.null(attribute_df)){
    attribute_df <- to_plot_df$attributes
  }

  if(!is.null(fill_name)){
    fill_value <- attribute_df[,fill_name]
    if(is.null(fill_range)){
      fill_range <- range(fill_value, na.rm = TRUE)
    }
  }
  
  cat("\n")
  print("Getting Plot Positions")
  time_pts <- colnames(freq_mat)
  plot_pos_df <- get_pos(clones, parents, freq_mat, as.numeric(time_pts))
  

  if(interpolation_steps > 0){
    cat("\n")
    print("Smoothing Polygons")
    plot_pos_df <- smooth_pos(plot_pos_df, n_intermediate_steps = interpolation_steps, interp_method=interp_method)
  }
  plot_pos_df$extinction_time <- max(plot_pos_df$x)
  
  true_time_pt_idx <- which(colnames(freq_mat) %in% og_time_pts)
  for(cidx in seq(1, length(clones))){
    clone_freq <- freq_mat[cidx, true_time_pt_idx]
    zero_idx <- which(clone_freq==0)
    origin_time <- which(clone_freq!=0)[1]
    time_dif <- origin_time - zero_idx
    if(any(time_dif < 0)){
      extinction_time_idx <- zero_idx[which(time_dif < 0)[1]]
      # extinction_time <- as.numeric(names(time_dif)[extinction_time_idx])
      extinction_time <- og_time_pts[extinction_time_idx]
      cidx_in_pos_df <- which(plot_pos_df$clone_id == clones[cidx])
      plot_pos_df$extinction_time[cidx_in_pos_df] <- extinction_time
    }
  }
  
  ### Supply attribute name since using deparse inside get_evofreq will return fill_value for the name of the attribute
  plot_pos_df <- update_colors(evo_freq_df = plot_pos_df, clones = clones, fill_value = fill_value, clone_cmap = clone_cmap, fill_range = fill_range, fill_name=fill_name, shuffle_colors = shuffle_colors)

  if(!scale_by_sizes_at_time){
    plot_pos_df$y <- plot_pos_df$y*max_mutation_size
    plot_pos_df$y_label <- "Population Size"
    
  }else{
    plot_pos_df$y_label <- "Frequency"
  }
  
  
  plot_pos_df$row_id <- seq(1, nrow(plot_pos_df))
  updated_time_pts <- unique(plot_pos_df$x)
  ###Replace closest x value with original timepoint value before merge
  for(tpt in time_pt_df$x){
    closest_x_idx <- which.min(abs(as.numeric(tpt) - updated_time_pts))
    closest_x <- updated_time_pts[closest_x_idx]
    plot_pos_df$x[which(plot_pos_df$x==closest_x)] <- as.numeric(tpt)
  }
  
  plot_pos_df <- merge(plot_pos_df, time_pt_df, by="x", all=TRUE)
  plot_pos_df <- plot_pos_df[order(plot_pos_df$row_id), ]
  
  return(plot_pos_df)
}

scale_values <- function(x, out_range=c(0, 1)){
  a <- min(out_range)
  b <- max(out_range)
  in_min <- min(x, na.rm = TRUE)
  in_max <- max(x, na.rm = TRUE)
  
  scaled_x <- (b-a)*(x-in_min)/(in_max - in_min) + a
  
  return(scaled_x)
}

get_genomes <- function(c_list, p_list, out="binary"){
  ###out=binary or bases
  n_clones <- length(c_list)
  if(out=="binary"){
    gene_mat <- matrix(0, nrow=n_clones, ncol=n_clones)
  }else{
    bases <- c("A", "C", "T", "G")
    gene_mat <- matrix(sample(bases, n_clones**2, replace = TRUE), nrow=n_clones, ncol=n_clones)
  }
  
  row.names(gene_mat) <- c_list
  for(i in seq(c_list)){
    cid <- c_list[i]
    children_idx <- get_all_idx(cid, c_list, p_list)
    if(out=="binary"){
      gene_mat[children_idx, i] <-  1
    }else{
      gene_mat[children_idx, i] <- sample(bases, 1)
    }
  }
  
  return(as.matrix(gene_mat))
}

traverse <- function(c, clones, parents, fnx){
  f <- match.fun(fnx)
  f(c, clones)
  for(n in get_children(c, clones, parents)){
    traverse(n, clones, parents, fnx)
  }
}

fill_in_gaps <- function(size_array){
  #'
  #'
  ### FOR TESTING ###
  # og_size_array <- c(0,0,0,1,2,3,4,0,0,0, 1, 2, 3, 0, 0, 1, 2, 0, 0)
  # size_array <- og_size_array
  #######
  
  
  above_zero_size_idx <- which(size_array > 0) ### indices for size array
  time_between_non_zero_idx <- diff(above_zero_size_idx) ## If there are no gaps in sizes, all these should be 1. Dif is difference between idx and idx + 1 
  left_side_of_above_0_gap_idx <- above_zero_size_idx[which(time_between_non_zero_idx > 1)] ### which indices in the original array are to the left side of a gap
  
  for(lidx in left_side_of_above_0_gap_idx){
    left_size <- size_array[lidx]
    
    temp_above_zero_idx <- which(size_array>0)
    right_non_zero_idx <- temp_above_zero_idx[which(temp_above_zero_idx>lidx)][1]
    right_size <- size_array[right_non_zero_idx]
    
    m <- (left_size-right_size)/(lidx-right_non_zero_idx)
    b <- left_size - m*lidx 
    time_btwn <- seq(lidx+1, right_non_zero_idx - 1)
    sizes_btwn <- m*time_btwn + b
    size_array[time_btwn] <- sizes_btwn
  }
  return(size_array)
}

check_for_missing_links_fnx <- function(clones, parents){
  for(cid in clones){
    get_ancestor_idx_evo_freq(cid, clone_list = clones, parent_list = parents)
  }
  return("No missing links")
}

check_and_update_edges <- function(clones, parents, check_for_missing_links=TRUE){
  ##### FOR TESTING ###
  # clones <- clone_list
  # parents <- parent_list
  ######
  
  ### Make sure each clone only appears once in edge list
  if(all(as.numeric(table(clones)) != 1)){
    warning("Clone occurs more than once in edgelist")
  }
  
  ### Update so that parent of root has an id not in the clone_list. 
  clone_as_parent_idx <- which(clones==parents)
  if(length(clone_as_parent_idx)>0){
    if(length(clone_as_parent_idx)>1){
      warning("More than 1 clone is its own parent. Cannot determine root")
    }
    rand_parent_id <- stats::runif(1, min = -1, max=0)
    parents[clone_as_parent_idx] <- rand_parent_id
  }
  
  ### Determine if root can be found ###
  root_id <- get_root_id(clones, parents)
  
  return(list("updated_clones"=clones, "updated_parents"=parents))
}

check_freq_mat <- function(freq_mat, clones, parents){
  ### FOR TESTING ###
  # freq_mat <- mutation_df
  # clones <- clone_list
  # parents <- parent_list
  #########
  
  ### Make sure total number children carrying the mutation  ###
  for(i in seq(1, length(clones))){
    cln <- clones[i]
    cln_freq <- freq_mat[i, ]
    all_children_idx <- get_all_idx(cln, clones, parents)
    children <- clones[all_children_idx]
    children_idx <- all_children_idx[clones[all_children_idx] != cln]
    if(length(children == 0)){
      next
    }
    children_df <- freq_mat[children_idx, ]
    parent_greater_than_children <- apply(children_df, 1, function(x){all(cln_freq>=x)})
    if(all(parent_greater_than_children)==FALSE){
      too_big_idx <- which(parent_greater_than_children==FALSE)
      clones_too_big <- clones[too_big_idx]
      stop(paste(cln, "has descendent mutations greater than that are greater than it's size. Descendents are:", clones_too_big))
    }
  }
}

#'\code{get_mutation_df} Converts sizes to the frequency of the mutations in the population
#'@inheritParams get_evofreq
#'@param clone_size_df Dataframe in a wide format, where each row corrsponds to a single clone, and the columns are the sizes of that clone at each timepoint
#'@examples
#'data("example.easy.wide")
#'### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#'time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
#'attribute_col_idx <-suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide)))))
#'attribute_df <- example.easy.wide[, attribute_col_idx]
#'size_df <- example.easy.wide[, time_col_idx]
#'parents <- example.easy.wide$parents
#'clones <- example.easy.wide$clones
#'### Size data can be converted to mutation counts for additional analyses, like calculating mutation frequencies
#'mutation_count_df <- get_mutation_df(size_df, clones, parents)
#'total_chromosomes <- 2*colSums(size_df) ### Calculation assumes that mutations are on a single chromosome
#'allele_frequency_df <- sweep(mutation_count_df,2,total_chromosomes,"/")
#'@export
get_mutation_df <- function(clone_size_df, clones, parents){
  ### FOR TESTING ####
  # clone_size_df <- as.data.frame(size_df)
  # clones <- clone_list
  # parents <- parent_list
  ####
  cat("\n")
  print("Getting Mutation Counts")
  greater_than_one_time_pt <- ncol(clone_size_df)>1
  mutation_df <- matrix(NA, nrow = nrow(clone_size_df), ncol = ncol(clone_size_df))
  
  pb <- utils::txtProgressBar(min = 1, max = nrow(clone_size_df), style = 3)
  for(i in seq(1, length(clones))){
    utils::setTxtProgressBar(pb, i)
    cln <- clones[i]
    all_children_idx <- get_all_idx(cln, clones, parents)
    
    child_sizes_df <- clone_size_df[all_children_idx, ]
    if(greater_than_one_time_pt){
      total_cells_carrying_mutation <- colSums(child_sizes_df) ### recursion already includes adding the current clone  
    }else{
      total_cells_carrying_mutation <- sum(child_sizes_df)
    }
    
    mutation_df[i, ] <- total_cells_carrying_mutation
    
  }
  
  # check_freq_mat(mutation_df, clones, parents)
  colnames(mutation_df) <- colnames(clone_size_df)
  rownames(mutation_df) <- rownames(clone_size_df)
  return(mutation_df)
}

predict_sizes <- function(x1, y1, x2, y2, new_x){
  ### FOR Testing ###
  # y1 <- 1.00
  # y2 <- 0.02
  # x1 <- 0
  # x2 <- 30
  # new_x <- 15
  ####
  
  
  fit <- stats::lm(c(y1, y2) ~ c(x1, x2))
  new_y <- fit$coefficients[[1]] + fit$coefficients[[2]]*new_x
  return(new_y)
}

make_pointy <- function(dull_mat){
  ### FOR TESTING ###
  # dull_mat <- new_mat
  ####
  ### Make that polygons do not have abrupt starts and stops
  ### Set size to 1 for intermediate timestep before origin
  ### If clone went extinct, set intemediate time step after extinction to 1
  
  
  time_pts <- colnames(dull_mat)
  max_time <- max(time_pts)
  min_time <- min(time_pts)
  for(cidx in seq(nrow(dull_mat))){
    # cidx <- 2
    clone_origin_time_idx <- which(dull_mat[cidx, ] > 0)[1]
    clone_origin_time <- time_pts[clone_origin_time_idx]
    if(clone_origin_time > min_time){
      dull_mat[cidx, clone_origin_time_idx-1] <- 0.0001
    }
    time_zero_idx <- which(dull_mat[cidx, ] == 0)
    if(length(time_zero_idx) > 0){
      
      if(any(time_zero_idx > clone_origin_time_idx)){
        extinction_time <- names(which(time_zero_idx > clone_origin_time_idx)[1])
        extinction_time_idx <- which(time_pts == extinction_time)
        dull_mat[cidx, extinction_time_idx] <- 0.0001       
        
      }      
    }
  }
  return(dull_mat)
}

interp_mut_mat <- function(mut_mat, time_pts= NULL){
  #For each clone, estimate intermediate sizes
  ### FOR TESTING ###
  # mut_mat <- og_mut_mat
  # time_pts <- og_time_pts
  ####
  
  if(is.null(time_pts)){
    time_pts <- seq(1, ncol(mut_mat))
  }
  
  intermediate_vals <- mapply(i=seq(2, length(time_pts)), FUN=function(i){(time_pts[i]-time_pts[i-1])/2 + time_pts[i-1]})
  new_times <- rep(NA, length(time_pts)+length(intermediate_vals))
  og_time_idx <- seq(1, length(new_times),2)
  new_time_idx <- seq(2, length(new_times)-1,2)
  
  new_times[og_time_idx] <- time_pts
  new_times[new_time_idx] <- intermediate_vals
  n_new_times <- length(new_times)
  ## Increase width of matrix
  
  new_mat <- matrix(NA, nrow=nrow(mut_mat), ncol=n_new_times)
  colnames(new_mat) <- new_times
  row.names(new_mat) <- row.names(mut_mat)
  
  # max_time <- max(time_pts)
  # min_time <- min(time_pts)
  
  for(cidx in seq(nrow(new_mat))){
    # cidx <- 2
    clone_sizes <-as.numeric( mut_mat[cidx, ])
    intermediate_sizes <- mapply(i=seq(2, length(time_pts)), function(i){
      # i <- 2
      if(clone_sizes[i]== 0 | clone_sizes[i-1]== 0){
        return(0)
      }

      #
      fit <- stats::lm(c(clone_sizes[i], clone_sizes[i-1]) ~ c(time_pts[i], time_pts[i-1]))
      new_y <- fit$coefficients[[1]] + fit$coefficients[[2]]* intermediate_vals[i-1]
      # new_y <- fit$coefficients[[1]] + fit$coefficients[[2]]* intermediate_vals[i]
      return(new_y)
    })
    # intermediate_sizes <- clone_sizes[1:length(clone_sizes)-1] + diff(clone_sizes)/2
    new_mat[cidx, og_time_idx] <- clone_sizes
    new_mat[cidx, new_time_idx] <- intermediate_sizes
    
  }
  
  new_mat <- make_pointy(new_mat)
  return(new_mat)
  
}

order_clones_at_time <- function(clones_at_time, parents, clones){
  ### possible that parent and child arose within the same period. Need to ensure that parents are drawn first
  
  #### For testing ###
  # clones_at_time <- parents_at_time
  # parents <- parents
  # clones <- clones
  ####
  root_id <- get_root_id(clones, parents)
  n_clones <- length(clones_at_time)
  n_ancestors_in_t <- rep(NA, n_clones)
  
  for(ct_i in seq(n_clones)){
    cid <- clones_at_time[ct_i]
    if(cid==root_id){
      ### makes sure root is always drawn first
      n_ancestors_in_t[ct_i] <- -1
      next()
    }
    ancestor_idx <- get_ancestor_idx_evo_freq(cid, clones, parents)
    if(is.null(ancestor_idx)){
      n_ancestors_in_t[ct_i] <- 0
    }else{
      n_ancestors_in_t[ct_i] <- length(ancestor_idx)
    }
  }
  new_order <- order(n_ancestors_in_t, decreasing = FALSE)
  clones_at_time_ordered <- clones_at_time[new_order]
  
  return(clones_at_time_ordered)
}

distribute_x <- function(x){
  #### For testing ####
  # x <- seq(1, 6)
  # x <- c(1)
  # x <- c(1800, 2008, 2166, 2202, 3564, 3694)
  #####################
  
  ### Put even numbers on one side of first point, odd numbers on the other
  x_idx <- seq(2, length(x))
  even_idx <- x_idx[x_idx %% 2 == 0]
  odd_idx <- x_idx[x_idx %% 2 != 0]
  new_order_idx <- c(rev(even_idx), 1, odd_idx)
  new_order <- x[new_order_idx]
  return(new_order)
}

filter_freq_mat <- function(clones, parents, freq_mat, threshold, attr_vals=NULL){
  ### FOR TESTING ###
  # attr_vals <- attribute_vals
  ##
  max_freqs <- apply(freq_mat, 1, max)
  above_thresh_idx <- which(max_freqs >= threshold)
  
  filtered_freq_mat <- freq_mat[above_thresh_idx, ]
  filtered_clones <- clones[above_thresh_idx]
  filtered_parents <- parents[above_thresh_idx]
  if(!is.null(attr_vals)){
    filtered_attr_vals <- attr_vals[above_thresh_idx]
  }else{
    filtered_attr_vals <- NULL
  }
  
  return_vals <- list('freq_mat'=filtered_freq_mat, 'clones'= filtered_clones, 'parents'=filtered_parents, 'attr_vals'=filtered_attr_vals, 'thresh_idx'=above_thresh_idx)
}

linear_interp_point <- function(x1, x2, y1, y2){
  #### FOR TESTING ####
  # x1 <- new_x[fp_idx-1]
  # x2 <- new_x[fp_idx+1]
  # y1 <- new_forward[fp_idx-1]
  # y2 <- new_forward[fp_idx+1]
  #####
  m <- (y1-y2)/(x1-x2)
  b <- y1 - m*x1
  middle_x <- 0.5*(x1 + x2)
  new_y <- m*middle_x + b
  return(new_y)
}

get_linear_m_and_b <- function(x1, x2, y1, y2){
  #### FOR TESTING ####
  # x1 <- new_x[fp_idx-1]
  # x2 <- new_x[fp_idx+1]
  # y1 <- new_forward[fp_idx-1]
  # y2 <- new_forward[fp_idx+1]
  #####
  m <- (y1-y2)/(x1-x2)
  b <- y1 - m*x1
  return(list("m"=m, "b"=b))
}

smooth_pos <- function(sparse_pos_df, n_intermediate_steps=20, interp_method = "monoH.FC"){
  ###FOR TESTING ##
  # sparse_pos_df <- evo_pos
  # n_intermediate_steps <- 3
  # interp_method <- "bezier"
  ###
  
  interp_df_list <- list()
  unique_clones <- unique(sparse_pos_df$clone_id)
  n_clones <- length(unique_clones)
  
  all_x <- unique(sparse_pos_df$x)
  x_range <- range(all_x, na.rm = TRUE)
  all_new_x <- seq(x_range[1], x_range[2], length.out = length(all_x)*n_intermediate_steps)

  pb <- utils::txtProgressBar(min = 0, max = n_clones, style = 3)
  for(i in seq(n_clones)){
    # print(i)
    # i <- 4
    cid <- unique_clones[i]
    utils::setTxtProgressBar(pb, i)
    clone_pos_df <- subset(sparse_pos_df, clone_id==cid)
    nx <- length(unique(clone_pos_df$x))
    
    if(nx >= 4){
      
      forward_idx <- which(duplicated(clone_pos_df$x)==FALSE)
      forward_df <- clone_pos_df[forward_idx,]
      ###GO BACK TO THIS??
      # new_x <- seq(min(forward_df$x), max(forward_df$x), length.out = length(forward_df$x)*n_intermediate_steps)
      ###
      
      clone_x_range <- range(forward_df$x)
      new_x <- all_new_x[all_new_x >= clone_x_range[1] & all_new_x <= clone_x_range[2]]
      
      ### need at least 4 points for spline interpolation
      # print(interp_method)
      # print(interp_method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman"))
      if(interp_method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman")){
        # print("FW")
        # interp_method <- "monoH.FC"
        func = stats::splinefun(x=forward_df$x, y=forward_df$y, method=interp_method,  ties = mean)
        new_forward <- func(new_x)
        new_forward_x <- new_x
      
      }else if(interp_method=="bezier"){
        bezier_forward_x <- seq(0, 1, length.out = length(new_x))
        forward_bezier_points <- bezier::bezier(t=bezier_forward_x, p=forward_df[,c("x","y")])
        new_forward <- forward_bezier_points[, 2]
        new_forward_x <- forward_bezier_points[, 1]
        
      }else if(interp_method=="bezier_curve_fit"){
        bezier_forward_x <- seq(0, 1, length.out = length(new_x))
        bz_forward_fit <- suppressWarnings(bezier::bezierCurveFit(as.matrix(forward_df[,c("x","y")]), na.fill = TRUE, maxiter=500, max.rse.percent.change = 0.2, fix.start.end = TRUE))
        forward_bezier_points <- bezier::bezier(t=bezier_forward_x, p=bz_forward_fit$p)
        new_forward <- forward_bezier_points[, 2]
        new_forward_x <- forward_bezier_points[, 1]
      }else if(interp_method=="loess"){
        fit <- loess(forward_df$y ~ forward_df$x)
        new_forward <- predict(fit, new_x)
        new_forward_x <- new_x
      }
      
      
      reverse_idx <- which(duplicated(clone_pos_df$x)==TRUE)
      reverse_df <- clone_pos_df[reverse_idx,]
      
      if(interp_method %in% c("fmm", "periodic", "natural", "monoH.FC", "hyman")){
        # print("REVERSE")
        rev_func = stats::splinefun(x=rev(reverse_df$x), y=rev(reverse_df$y), method=interp_method,  ties = mean) ### x has to increase
        new_reverse <- rev_func(new_x)
        new_reverse_x <- rev(new_x)
        
      }else if(interp_method=="bezier"){
        bezier_reverse_x <- seq(0, 1, length.out = length(new_x))
        reverse_bezier_points <- bezier::bezier(t=bezier_reverse_x, p= reverse_df[seq(nrow(reverse_df), 1),c("x","y")])
        new_reverse <- reverse_bezier_points[, 2]
        new_reverse_x <- rev(reverse_bezier_points[, 1])
      }else if(interp_method=="bezier_curve_fit"){
        bezier_reverse_x <- seq(0, 1, length.out = length(new_x))
        bz_reverse_fit <- suppressWarnings(bezier::bezierCurveFit(as.matrix(reverse_df[seq(nrow(reverse_df), 1),c("x","y")]), na.fill = TRUE, maxiter=500, max.rse.percent.change = 0.2, fix.start.end = TRUE))
        reverse_bezier_points <- bezier::bezier(t=bezier_reverse_x, p=bz_reverse_fit$p)
        new_reverse <- reverse_bezier_points[, 2]
        new_reverse_x <- rev(reverse_bezier_points[, 1]) 
      }else if(interp_method=="loess"){
        rev_fit <- loess(rev(reverse_df$y) ~ rev(reverse_df$x))
        new_reverse <- predict(rev_fit, new_x)
        new_forward_x <- rev(new_x)
      }
      

      #### Make sure that there are no times where the bottom and top cross. Bottom (forward) should always be less than top (reverse)
      rev_forward_dist <- new_reverse - new_forward
      flip_idx <- which(rev_forward_dist < 0)
      n_new_x <- length(new_x)
      if(length(flip_idx) > 0){
        # fp_idx <- 2
        for(fp_idx in flip_idx){
          updated_dist <- new_reverse[fp_idx] - new_forward[fp_idx]
          if(updated_dist > 0){
            #### position was already updated
            next()
          }
          if(fp_idx > n_new_x | n_new_x <= 1 | fp_idx == 1){
            ### first and last point are the new pointy ends. This just flips them to make a triangle
            replace_forward <- new_reverse[fp_idx]
            replace_reverse <- new_forward[fp_idx]
            new_forward[fp_idx] <- replace_forward
            new_reverse[fp_idx] <- replace_reverse
            
          }else{
            
  
            ### For forward: get next positive difference (should have been negative because forward should be < reverse). Connect to that point. Will also need to update any other points along that path
            next_seq <- seq(fp_idx+1, n_new_x)
            # next_positive_idx <- fp_idx + which(rev_forward_dist[fp_idx+1:n_new_x] > 0)[1]
            next_positive_idx <- fp_idx + which(rev_forward_dist[next_seq] > 0)[1]
            if(is.na(next_positive_idx)){
              ###ALL positions after this one need to be flipped
              replace_forward <- new_reverse[next_seq]
              replace_reverse <- new_forward[next_seq]
              new_forward[next_seq] <- replace_forward
              new_reverse[next_seq] <- replace_reverse
              break()
            }
            replace_seq <- seq(fp_idx, next_positive_idx - 1)
            
            forward_lm <- get_linear_m_and_b(new_x[fp_idx-1], new_x[next_positive_idx], new_forward[fp_idx-1], new_forward[next_positive_idx])
            new_forward[replace_seq] <- forward_lm$m*new_x[replace_seq] + forward_lm$b
            # new_forward[replace_seq] <- func(forward_lm$m*new_x[replace_seq] + forward_lm$b)
            
            ### For reverse: get next negative difference. Connect to that point. Will also need to update any other points along that path
            reverse_lm <- get_linear_m_and_b(new_x[fp_idx-1], new_x[next_positive_idx], new_reverse[fp_idx-1], new_reverse[next_positive_idx])
            new_reverse[replace_seq] <- reverse_lm$m*new_x[replace_seq] + reverse_lm$b
            
            # replace_forward <- linear_interp_point(new_x[fp_idx-1], new_x[fp_idx+1], new_forward[fp_idx-1], new_forward[fp_idx+1])
            # replace_reverse <- linear_interp_point(new_x[fp_idx-1], new_x[fp_idx+1], new_reverse[fp_idx-1], new_reverse[fp_idx+1])
          }
          
          # if(replace_reverse - replace_forward ){
          # new_forward[fp_idx] <- replace_reverse
          # new_reverse[fp_idx] <- replace_forward
          # }else{
          
        }
      }
      
      # print(c(length(new_x), length(new_forward), length(new_reverse)))
      new_forward_df <- data.frame('x'=new_x, 'y'=new_forward, 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])
      new_reverse_df <- data.frame('x'=rev(new_x), 'y'=rev(new_reverse), 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])
      
      # new_forward_df <- data.frame('x'=new_forward_x, 'y'=new_forward, 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])
      # new_reverse_df <- data.frame('x'=new_reverse_x, 'y'=rev(new_reverse), 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])
      
      new_clone_pos_df <- rbind(new_forward_df, new_reverse_df)
      
      # gp <- ggplot2::ggplot(new_clone_pos_df, ggplot2::aes(x=x, y=y)) +
      # ggplot2::geom_polygon()
      # print(gp)
      
      
    }else{
      new_clone_pos_df <- clone_pos_df[colnames(clone_pos_df) != "size"]
    }
    interp_df_list[[as.character(cid)]] <- new_clone_pos_df
  }
  interp_df <- do.call(rbind, interp_df_list)
  
  
  ### interpolation may cause points to go out of bounds
  interp_df$y[interp_df$y < 0] <- 0
  interp_df$y[interp_df$y > 1] <- 1
  
  
  return(interp_df)
  
}

#'\code{plot_evofreq} Plots the frequency dynamics using ggplot. Can also use info to create animations 
#'@param freq_frame Proper dataframe returned by \code{\link{get_evofreq}}
#'@param n_time_pts Integer defining how many time points to plot, evenly spaced between \code{start_time} and \code{end_time}. If NULL, all timepoints are plotted
#'@param start_time Integer defining the timepoint at which to start plotting frequency dynamics. If NULL, plotting will begin at the first timepoint
#'@param end_time Integer defining the timepoint at which to stop plotting frequency dynamics. If NULL, plotting will end at the last timepoint
#'@param bw width of lines surrounding polygons
#'@param bc color of lines surrounding polygons
#'@param show_axes Whether to display axes
#'@param fill_range Array containing the minimum and maximum values to set the range of colors. If NULL (the default), the range is determined directly from \code{fill_value}.
#'@return ggplot of the frequency dynamics
#'@examples
#' 
#' data("example.easy.wide.with.attributes")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#' size_df <- example.easy.wide.with.attributes[, time_col_idx]
#' parents <- example.easy.wide.with.attributes$parent
#' clones <- example.easy.wide.with.attributes$clone
#' fitness <- example.easy.wide.with.attributes$fitness
#' 
#' #' 
#' ### Setting of colors can be done when getting the freq_frame, or by updating the color later using \code{\link{update_colors}}. Available colormaps are those found in \code{\link[colormap]{colormaps}}
#' ### Default colormap is rainbow_soft, but this can be changed using the \code{clone_cmap} argument. 
#' freq_frame <- get_evofreq(size_df, clones, parents)
#' evo_p <- plot_evofreq(freq_frame)
#' 
#' ### Can color each clone by an attribute by providing a \code{fill_value}. Default colormap is viridis, but this can be changed using the \code{clone_cmap} argument
#' fitness <- runif(length(clones))
#' fitness_freq_frame <- get_evofreq(size_df, clones, parents, fill_value = fitness)
#' fitness_evo_p <- plot_evofreq(fitness_freq_frame)
#' #' 
#' ### The user can also provide custom colors for each clone, which will need to be passed into the \code{fill_value} argument
#' ### Custom colors can be defined using RGB values. Each color should be a string specifying the color channel values, separated by commas.
#' rgb_clone_colors <- sapply(seq(1, length(clones)), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
#' rgb_freq_frame <- get_evofreq(size_df, clones, parents, rgb_clone_colors)
#' rgb_evo_p <- plot_evofreq(rgb_freq_frame)
#' 
#' ### Custom colors can also be any of the named colors in R. A list of the colors can be found with \code{colors()}
#' named_clone_colors <- sample(colors(), length(clones), replace = FALSE)
#' named_freq_frame <- update_colors(rgb_freq_frame, clones = clones, fill_value = named_clone_colors)
#' named_evo_p <- plot_evofreq(named_freq_frame)
#' 
#' ### Custom colors can also be specified using hexcode
#' hex_clone_colors <- sample(colormap::colormap(colormap=colormaps$temperature, nshades=length(clones)))
#' hex_freq_frame <- update_colors(rgb_freq_frame, clones = clones, fill_value = hex_clone_colors)
#' hex_evo_p <- plot_evofreq(hex_freq_frame)
#'
#' ### Can revert back to original colors
#'freq_frame_default_color <- update_colors(fitness_freq_frame, clones=clones)
#'default_cmap_evo_p <- plot_evofreq(freq_frame_default_color)
#'
#'### Can add gganimate objects to evo_p_by_size to create the animation
#'
#'\donttest{
#' library(gganimate)
#' movie_p <- fitness_evo_p +
#'   transition_reveal(x, range=range(fitness_freq_frame$x)) +
#'    view_follow()
#' # print(movie_p)
#' anim_save("evofreq_movie.gif", movie_p)
#' }
#'@export
plot_evofreq <- function(freq_frame, n_time_pts=NULL, start_time=NULL, end_time=NULL, bw=0.05, bc="grey75", show_axes=TRUE, fill_range=NULL){
  
  ### FOR TESTING ###
  # n_time_pts <- NULL
  # start_time <- NULL
  # end_time <-  NULL
  # bw <- 0.05
  # bc <- "grey75"
  # freq_frame <- hal_plot_df
  #####
  
  unique_time_pts <- unique(freq_frame$x)
  if(is.null(n_time_pts)){
    n_time_pts <- length(unique_time_pts)
  }
  if(is.null(start_time)){
    start_time <- 0
  }
  if(is.null(end_time)){
    end_time <- max(freq_frame$x)
  }

  
  # end_time <- 50
  view_df <- subset(freq_frame, x <= end_time)
  
  color_attribute_name <- unique(view_df$efp_color_attribute)
  if(is.na(color_attribute_name)){
    color_attribute_name <- "plot_color"
  }else{
    if(is.null(fill_range)){
      fill_range <- range(view_df[, color_attribute_name], na.rm = TRUE)  
    }
    
  #   color_df <-  view_df[duplicated(view_df$plot_color)==FALSE, ]
  #   color_df <- color_df[order(color_df[color_attribute_name]), ]
  #   colorbar_colors <- color_df$plot_color
  #   
  }
  
  
  # all_x_maxs_at_time <- aggregate(x ~ clone_id, data=subset(view_df, extinction_time >= end_time), FUN=max)
  # max_x <- min(all_x_maxs_at_time$x)  
  # view_df$x <- as.numeric(view_df$x)
  time_pts_df <- unique(view_df[complete.cases(view_df["Time_label"]), ][c("x", "Time_label")])
  time_pts_df$x <-  as.numeric(as.character(time_pts_df$x))
  time_pts_df <- time_pts_df[order(time_pts_df$x), ]
  time_pts_df <- subset(time_pts_df, x > 0)
  time_pts_df$Time_label <-  as.character(time_pts_df$Time_label)
  # time_pts_df$Time_label <-  as.character(as.numeric(time_pts_df$Time_label))

  # time_pts_df$Time_label[1] <- "0" ### TODO Why is this being set to 0?
  
  y_label <- unique(freq_frame$y_label)
  
  ggevodyn <- ggplot2::ggplot(view_df, ggplot2::aes_string(x="x", y="y", group="draw_order", fill=color_attribute_name)) +
    ggplot2::geom_polygon(size=bw, color=bc) +
    ggplot2::ylab(y_label) +
    ggplot2::xlab("Time") +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(breaks=as.numeric(as.character(time_pts_df$x)), labels=as.character(time_pts_df$Time_label)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if(color_attribute_name=="plot_color"){
    ggevodyn <- ggevodyn + ggplot2::scale_fill_identity()
  }else{
    colormap_name <- unique(view_df$cmap)
    # ggevodyn <- ggevodyn + ggplot2::scale_fill_gradientn(colours = colorbar_colors)  #colormap::scale_fill_colormap(color_attribute_name, colormap=colormap_name)
    ggevodyn <- ggevodyn + colormap::scale_fill_colormap(color_attribute_name, colormap=colormap_name, limits=fill_range)
  }
  
  if(! show_axes){
    ggevodyn <- ggevodyn + ggplot2::theme_void()
  }
  return(ggevodyn)
}

#'Convert data in long format to wide format
#'\code{long_to_wide_freq_ready} Converts long data frame to wide format used by evofreq.
#'@param edges_df Dataframe with 2 columns, where each row defines relation between each clone and their parent. One column contains the clone ids, and the other contains the id of that clone's parent.
#'@param long_pop_sizes_df Dataframe with 3 columns, one containing the clone id, one containing the time points, and one containing the size of the clone at that time point. There must be a size for each clone at every time point, so it did not exist, it's size is 0.
#'@param time_col_name String that specifies which column in long_pop_sizes_df contains the timepoint information
#'@param clone_col_name String that specifies which column in in edges_df contains the names of the clones. I.e. the "to" node in a graph
#'@param parent_col_name String that specifies which column in in edges_df contains the names of the parents of the descendent clone in the edges_df with colname clone_col_name. #I.e. the "from" node in a graph
#'@param size_col_name String that specifies which column in long_pop_sizes_df contains the size of each clone for the corresponding timepoint and clone in long_pop_sizes_df
#'@param fill_gaps_in_size Boolean specificing if gaps in sizes over time should be filled in, assuming linear changes. 
#'@return List containing: wide_size_df, a dataframe of the sizes in wide format, where each row contains the sizes over time for a single clone; parents, a vector containing the parent id for each clone in clones; clones, a vector of clones ids correspond to each row in wide_size_df
#'@examples
# ### Input needs to be in wide format, but can convert long format data to wide format using \code{long_to_wide_freq_ready}
# wide_df_info <- long_to_wide_freq_ready(long_pop_sizes_df = example.easy.long.sizes, time_col_name = "Time", clone_col_name = "clone", parent_col_name = "parent", size_col_name = "Size", edges_df = example.easy.long.edges)
# clones_from_long <- wide_df_info$clones
# parents_from_long <- wide_df_info$parents
# size_df_from_long <- wide_df_info$wide_size_df
# pos_from_long_df <- get_evofreq(size_df_from_long, clones_from_long, parents_from_long)
# evo_p_from_long <- plot_evofreq(pos_from_long_df)
#'@export
long_to_wide_freq_ready <- function(edges_df, long_pop_sizes_df, time_col_name, clone_col_name, parent_col_name, size_col_name, fill_gaps_in_size=FALSE){
  long_pop_sizes_df <- as.data.frame(long_pop_sizes_df)
  unique_times <- unique(long_pop_sizes_df[,time_col_name])
  unique_times <- unique_times[order(unique_times)]
  n_time_pts <- length(unique_times)
  
  clones_in_edge_df <- edges_df[, clone_col_name]
  parents_in_edge_df <- edges_df[, parent_col_name]
  clones_in_size_df <- unique(long_pop_sizes_df[,clone_col_name])
  clones_not_in_clone_list <- clones_in_size_df[!clones_in_size_df %in% clones_in_edge_df]
  
  if(length(clones_not_in_clone_list)>0){
    if(length(clones_not_in_clone_list)>1){
      stop("More thant 1 clone in size data frame but not in clone list")
    }
    warning("Clone in size data frame but not in clone list. Assuming this is the root")
    rand_parent_id <- stats::runif(1)
    clones_in_edge_df <- c(clones_in_edge_df, clones_not_in_clone_list)
    parents_in_edge_df <- c(parents_in_edge_df, rand_parent_id)
  }else{
    # root_id <- get_root_id(clones_in_edge_df, parents_in_edge_df)
    clone_as_parent_idx <- which(parents_in_edge_df==clones_in_edge_df)
    if(length(clone_as_parent_idx)>0){
      rand_parent_id <- stats::runif(1)
      parents_in_edge_df[clone_as_parent_idx] <- rand_parent_id
      # clones_in_edge_df <- c(clones_in_edge_df, clones_not_in_clone_list)
      # parents_in_edge_df <- c(parents_in_edge_df, rand_parent_id)
      if(length(clone_as_parent_idx)>1){
        warning("More than 1 clone is its own parent. Cannot determine root")
      }  
    }
  }
  
  edges_df <- data.frame(parents_in_edge_df, clones_in_edge_df)
  names(edges_df) <- c(parent_col_name, clone_col_name)
  
  n_clones <- length(clones_in_edge_df)
  size_mat <- matrix(0, nrow=n_clones, ncol=length(unique_times))
  colnames(size_mat) <- unique_times
  
  size_mat_rownames <- rep(NA, n_clones)
  cat("\n")
  print("Converting From Long to Wide Format")
  pb <- utils::txtProgressBar(min = 1, max = n_clones, style = 3)
  for(i in seq(n_clones)){
    utils::setTxtProgressBar(pb, i)
    cid <- clones_in_edge_df[i]
    cidx <- which(long_pop_sizes_df[clone_col_name]==cid)
    if(length(cidx)== 0){
      size_mat[i, ] <- 0
    }else{
      cid_sizes <- long_pop_sizes_df[cidx, size_col_name]  
      cid_times <- long_pop_sizes_df[cidx, time_col_name]
      c_time_idx <- sapply(cid_times, function(x){which(unique_times==x)})
      
      size_mat[i, c_time_idx] <- cid_sizes
    }
    size_mat_rownames[i] <- cid
  }
  rownames(size_mat) <- size_mat_rownames
  
  if(fill_gaps_in_size){
    wide_size_df <- as.data.frame(t(apply(size_mat, 1, fill_in_gaps)))
  }else{
    wide_size_df <- as.data.frame(size_mat)
  }
  
  return(list("wide_size_df"=wide_size_df, "parents"=edges_df[, parent_col_name], "clones"=edges_df[, clone_col_name]))
  
  
  
  ###Reorder edges so that they are the same as 
  # long_pop_sizes_df <- long_pop_sizes_df[order(long_pop_sizes_df[,time_col_name]),]
  # wide_size_df <- reshape(as.data.frame(long_pop_sizes_df), timevar = time_col_name, idvar = clone_col_name, direction = "wide")
  
  # parents <- edges_df[, parent_col_name]
  # clones <- edges_df[, clone_col_name]
  # clones_in_size_df <- unique(long_pop_sizes_df[,clone_col_name])
  # clones_not_in_clone_list <- clones_in_size_df[!clones_in_size_df %in% clones]
  # if(length(clones_not_in_clone_list)>0){
  #   if(length(clones_not_in_clone_list)>1){
  #     stop("More thant 1 clone in size data frame but not in clone list")
  #   }
  #   warning("Clone in size data frame but not in clone list. Assuming this is the root")
  #   rand_parent_id <- runif(1)
  #   clones <- c(clones, clones_not_in_clone_list)
  #   parents <- c(parents, rand_parent_id)
  #   edges_df <- data.frame(clones, parents)
  #   names(edges_df) <- c(clone_col_name, parent_col_name)
  # }else{
  #   
  # }
  
  # wide_size_df <- merge(wide_size_df, edges_df)
  
  # ordered_clones <- wide_size_df[,clone_col_name]
  # ordered_parents <- wide_size_df[,parent_col_name]
  # ordered_size_df <- wide_size_df[,!colnames(wide_size_df) %in% c(clone_col_name, parent_col_name)]
  # return(list("wide_size_df"=ordered_size_df, "parents"=ordered_parents, "clones"=ordered_clones))
}

#' @title get_evofreq_labels
#' @param freq_frame Properly formatted dataframe from \code{\link{get_evofreq}}.
#' @param evofreq_plot Plot returned from \code{\link{plot_evofreq}}.
#' @param apply_labels Whether to plot the labels. Must include evofreq_plot if TRUE. Default FALSE.
#' @param custom_label_text Custom labels to use (e.g. c("Mutant X","Mutant Y"))
#' @param extant_only Get labels for only the extant clones.
#' @param clone_list Only label clones with this id. Passed as a vector.
#' @param line_color Single color value. Default "black".
#' @param adj.factor Value to scale the label and lines by so that they are far enough from the plot. Default value is 10 (a 10th of the distance between max/min at x from max/min across x).
#' @return Dataframe containing the necessary coordinates to add labels in the proper locations.
#' 
#' @examples
#' data("example.easy.wide") # Load Data
#' ### Get frequency dataframe
#' freq_frame <- get_evofreq(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma")
#' evofreq_plot <- plot_evofreq(freq_frame)
#' clone_labels <- get_evofreq_labels(freq_frame, extant_only=FALSE)
#' 
#' ### Set evofreq_plot and apply_labels = TRUE to get a labeled plot back
#' get_evofreq_labels(freq_frame, extant_only=FALSE, evofreq_plot = evofreq_plot, apply_labels = TRUE)
#' @export
get_evofreq_labels <- function(freq_frame, apply_labels=FALSE, custom_label_text=NULL, evofreq_plot=NULL, clone_list=NULL, extant_only=FALSE, line_color="darkgrey", adj.factor=10){
  ### FOR TESTING ###
  # freq_frame <- freq_frame
  # apply_labels=TRUE
  # evofreq_plot=evo_freq_p
  # extant_only=FALSE
  # line_color="darkgrey"
  # adj.factor=10
  # clone_list = NULL
  ####
  
  # Get midpoint of clone
  max_ever = max(freq_frame$y)[1]
  min_ever = min(freq_frame$y)[1]
  mid <- freq_frame[freq_frame$x==min(freq_frame$x),]$y[[1]]
  spreader <- sd(freq_frame$x)/4.
  # Split by clone
  freq_frame_split <- split( freq_frame , f = freq_frame$clone_id )

  # Main values
  # position_list <- lapply(freq_frame_split, FUN=function(z){return( z[z$x==min(z$x),][1,c(1,2,3,4,7)] ) })
  info_cols <- c("clone_id",  "x", "y", "parent", "extinction_time", "plot_color")
  position_list <- lapply(freq_frame_split, FUN=function(z){return( z[z$x==min(z$x),][1,info_cols]) })
  position_df <- do.call(rbind, position_list)
  # Get the max and min value at each x regardless of clone
  position_df$max <- unlist(lapply(position_list, FUN=function(z){ return( max(freq_frame[ (freq_frame$x<z$x+spreader) & (freq_frame$x>z$x-spreader),]$y) ) }))
  position_df$min <- unlist(lapply(position_list, FUN=function(z){ return( min(freq_frame[ (freq_frame$x<z$x+spreader) & (freq_frame$x>z$x-spreader),]$y) ) }))
  
  # Assign information based on midline
  position_df$y_plotValue <- unlist(lapply(split( position_df , f = position_df$clone_id ), FUN=function(z){ if(z$y>mid){return(z$max+(abs(max_ever-z$max))/adj.factor)}else{return(z$min-(abs(z$min-min_ever))/adj.factor)} }))
  
  # Adjust if multiple have the same x:
  
  dup_x <- position_df$x[which(duplicated(position_df$x))]
  keeper_df <- subset(position_df, ! x %in% dup_x)
  for(dx in dup_x){
    x_df <- subset(position_df, x==dx)
    new_mid <- mean(x_df$y)
    x_df$y_plotValue <- unlist(lapply(split(x_df, f=x_df$clone_id), FUN=function(z){ if(z$y>new_mid){return(z$max+(abs(max_ever-z$max))/adj.factor)}else{return(z$min-(abs(z$min-min_ever))/adj.factor)} }))
    keeper_df <- rbind(keeper_df, x_df)
  }
  position_df <- keeper_df
  

  ### Get origin time, final time, and extinct/extant status
  max_time <- max(freq_frame$x)
  final_time_df <- aggregate(list("final_time"=freq_frame$x), by = list("clone_id" = freq_frame$clone_id), FUN=max)
  origin_time_df <- aggregate(list("origin_time"=freq_frame$x), by = list("clone_id" = freq_frame$clone_id), FUN=min)
  position_df <- Reduce(function(x, y) merge(x, y, all=TRUE, by="clone_id"), list(position_df, origin_time_df, final_time_df))
  position_df$extant <- TRUE
  position_df$extant[position_df$final_time < max_time] <- FALSE
  if(extant_only){
    position_df <- subset(position_df, position_df$extant==TRUE)
  }
  
  if(!is.null(clone_list)){
    position_df <- subset(position_df, position_df$clone_id %in% clone_list)
  } else {
    position_df <- position_df
  }
  
  if(apply_labels & !is.null(evofreq_plot) & is.null(custom_label_text)){
    evo_freq_p_labeled <- evofreq_plot + geom_segment(data = position_df, aes(x=position_df$x,xend=position_df$x,y=position_df$y,yend=position_df$y_plotValue), color=line_color, size=1, inherit.aes = FALSE) +
      geom_label(data=position_df, aes(x=position_df$x, y=position_df$y_plotValue, fill=position_df$plot_color, label=paste("Clone",position_df$clone_id)), color="white", inherit.aes = FALSE)
    return(evo_freq_p_labeled)
  } else if(apply_labels & !is.null(evofreq_plot) & !is.null(custom_label_text)) {
    evo_freq_p_labeled <- evofreq_plot + geom_segment(data = position_df, aes(x=position_df$x,xend=position_df$x,y=position_df$y,yend=position_df$y_plotValue), color=line_color, size=1, inherit.aes = FALSE) +
      geom_label(data=position_df, aes(x=position_df$x, y=position_df$y_plotValue, fill=position_df$plot_color, label=custom_label_text), color="white", inherit.aes = FALSE)
    return(evo_freq_p_labeled)
  } else if(apply_labels){
    # No evofreq_plot
    stop("Please provide an evofreq_plot created by plot_evofreq() or ensure that custom_label_text is provided")
  } else {
    return(position_df) # just return dataframe
  }
}

filter_data <- function(size_df, clones, parents, time_pts=NULL, attribute_df=NULL, threshold=0.01, scale_by_sizes_at_time = FALSE, data_type="size", fill_gaps_in_size = FALSE, test_links=TRUE, add_origin=FALSE, tm_frac=0.6, rescale_after_thresholding=FALSE){
  # data("example.easy.wide.with.attributes")
  # ## Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
  # time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
  # attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
  # attribute_df <- example.easy.wide.with.attributes[, attribute_col_idx]
  # size_df <- example.easy.wide.with.attributes[, time_col_idx]
  # parents <- example.easy.wide.with.attributes$parent
  # clones <- example.easy.wide.with.attributes$clone
  # 
  # clone_cmap <- "rainbow_soft"
  # # size_df <- mut_freq
  # threshold <- 0.0
  # # clones <- mut_df$clone
  # # parents <- mut_df$parent
  # time_pts <- NULL
  # fill_name <- "fitness" #"new_antigenicity"
  # attribute_df <- attribute_df
  # data_type <- "size"
  # data_type <- "size"
  # scale_by_sizes_at_time <- FALSE
  # interpolation_steps <- 10
  # fill_gaps_in_size <- FALSE
  # test_links <- TRUE
  # threshold <- 0.01
  # fill_range <- NULL
  # add_origin <- TRUE
  # tm_frac <- 0.6
  #####
  
  
  ### Possible that column names are numbers, but timepoints are unequally spaced. Only change if cannot be converted to numeric
  og_colnames <- colnames(size_df)
  if(suppressWarnings(all(is.na(as.numeric(colnames(size_df)))))){
    colnames(size_df) <- seq(1, ncol(size_df)) 
  }
  
  if(any(duplicated(clones))){
    warning("Some clones have the same ID. Each clone should have a unique ID")
  }
  if(is.null(time_pts)){
    time_pts <- colnames(size_df)
  }else{
    time_pts <- as.character(sort(as.numeric(unique(time_pts))))
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
    if(add_origin){
      og_time_pts <- colnames(freq_mat)
      freq_mat <- add_origin_mat(freq_mat, tm_frac = tm_frac)
      # neg_time_pts <- colnames(freq_mat)[! colnames(freq_mat) %in% og_time_pts]
      added_time_pts <- colnames(freq_mat)[! colnames(freq_mat) %in% og_time_pts]
      time_pts <- c(added_time_pts, time_pts)
      og_colnames <- c(rep("", length(added_time_pts)), og_colnames)
    }
    
  }else{
    check_freq_mat(size_df, clones, parents)
    freq_mat <- size_df
    freq_mat <- as.data.frame(freq_mat)
    if(add_origin){
      og_time_pts <- colnames(freq_mat)
      freq_mat <- add_origin_mat(freq_mat, tm_frac = tm_frac)
      added_time_pts <- colnames(freq_mat)[! colnames(freq_mat) %in% og_time_pts]
      time_pts <- c(added_time_pts, time_pts)
      og_colnames <- c(rep("", length(added_time_pts)), og_colnames)
      
    }
    
    # parents[which(! parents %in% clones)]
    # updated_edges <- check_and_update_edges(clones, parents)
    # clones <- updated_edges$updated_clones
    # parents <- updated_edges$updated_parents
  }
  
  ### Possible that parent has frequency of 0. Shouldn't, but possible output of bioinformatics tools
  # size_zero_idx <- which(apply(freq_mat, 1, function(x){all(x==0)}))
  # clones_with_all_size_zero <- clones[clones %in% clones[size_zero_idx]]
  # p_with_all_size_zero <- clones_with_all_size_zero[which(clones_with_all_size_zero %in% parents)]
  # for(p in p_with_all_size_zero){
  #   children_idx <- get_all_idx(p, clones, parents)
  #   children_ids <- clones[children_idx]
  #   actual_children_idx <- which(children_ids != p)
  #   children_origins <- apply(freq_mat[children_idx[actual_children_idx], ], 1, function(x){which(x>0)[1]})
  #   first_child <- names(which.min(children_origins))
  #   first_child_idx <- which(clones==first_child)
  #   
  #   pidx <- which(clones==p)
  #   freq_mat[pidx, ] <- freq_mat[first_child_idx, ]
  # }
  
  
  
  
  max_mutation_size <- max(freq_mat)
  if(scale_by_sizes_at_time){
    max_sizes_at_each_time <- apply(freq_mat, 2, max)
    freq_mat <- sweep(freq_mat, MARGIN = 2, max_sizes_at_each_time, FUN = "/")
    
  }else{
    freq_mat <- freq_mat/max_mutation_size
  }
  
  
  ### Subset time points to be plotted
  
  # freq_mat <- freq_mat[,which(colnames(size_df) %in% time_pts)]
  freq_mat <- freq_mat[, time_pts]
  if(is.null(nrow(freq_mat))){
    ###Only 1 clone
    freq_mat <- t(as.matrix(freq_mat))
  }
  colnames(freq_mat) <- time_pts
  row.names(freq_mat) <- clones
  
  ### Possible that a clone was recorded, but never had any sizes greater than 0
  existed_idx <- which(rowSums(freq_mat) > 0)
  freq_mat <- freq_mat[existed_idx, ]
  clones <- clones[existed_idx]
  parents <- parents[existed_idx]
  
  ### Get idx of clones that are extant after subsetting timepoints ###
  origin_times <- apply(freq_mat, 1, function(x){which(x>0)[1]})
  idx_after_time_thresh <- which(is.na(origin_times)==FALSE)
  freq_mat <- freq_mat[idx_after_time_thresh, ]
  
  if(is.null(nrow(freq_mat))){ 
    ###Only 1 clone
    freq_mat <- t(as.matrix(freq_mat))
  }
  
  clones <- clones[idx_after_time_thresh]
  parents <- parents[idx_after_time_thresh]
  origin_times <-origin_times[idx_after_time_thresh] ### Order by origin times at the end
  time_pts <- as.numeric(time_pts)
  
  if(!is.null(attribute_df)){
    attribute_df <- attribute_df[idx_after_time_thresh, ]
  }
  
  ### Filter values based on threshold
  if(threshold > 0){
    filtered_info <- filter_freq_mat(clones, parents, freq_mat, threshold)
    freq_mat <- filtered_info$freq_mat
    parents <- filtered_info$parents
    
    clones <- filtered_info$clones
    filtered_idx <- filtered_info$thresh_idx
    origin_times <- origin_times[filtered_idx]
    
    if(rescale_after_thresholding){
      freq_mat <- rescale_frequencies(size_df, clones, parents, filtered_idx, data_type, scale_by_sizes_at_time)
      # if(data_type=="size"){
      #   sink("/dev/null") ### Suppress output
      #   freq_df <- get_mutation_df(size_df[filtered_idx, ], clones = clones, parents = parents)
      #   sink()
      #   freq_mat <- as.matrix(freq_df)
      # }else{
      #   freq_mat <- as.matrix(freq_mat[filtered_idx,])
      # }
      # 
      # 
      # max_mutation_size <- max(freq_mat)
      # if(scale_by_sizes_at_time){
      #   max_sizes_at_each_time <- apply(freq_mat, 2, max)
      #   freq_mat <- sweep(freq_mat, MARGIN = 2, max_sizes_at_each_time, FUN = "/")
      #   
      # }else{
      #   freq_mat <- freq_mat/max_mutation_size
      # }
    }
    
    if(is.null(nrow(freq_mat))){ 
      ### Only 1 clone after thresholding
      freq_mat <- t(as.matrix(freq_mat))
      row.names(freq_mat) <- clones
      filtered_idx <- as.numeric(filtered_idx)
    }
    
    if(!is.null(attribute_df)){
      attribute_df <- attribute_df[filtered_idx, ]
    }
  }
  if(is.null(attribute_df)){
    attribute_df <- data.frame("clone_id"=clones, "parents"=parents)
  }
  
  attribute_df$origin <- origin_times
  
  ordered_idx <- order(origin_times)
  clones <- clones[ordered_idx]
  parents <- parents[ordered_idx]
  attribute_df <- attribute_df[ordered_idx, ]
  freq_mat <- freq_mat[ordered_idx, ]
  # 
  # clone_col_idx <- as.numeric(which(sapply(colnames(attribute_df), FUN = function(x){all(clones %in% as.character(unique(attribute_df[,x])))})==TRUE))
  # 
  # if(length(clone_col_idx) > 1){
  #   #### If all clones in pos df are also parents, and there are both parent and clone ids in attribute_df, then more than 1 column will considered the clone_id column in attribute_df
  #   ### Not possible for all clones in attribute df to also all be parents 
  #   ### So clone_id column is the one that has more unique elements
  #   n_unique_clones <- apply(attribute_df[clone_col_idx], 2, function(x){length(unique(x))})
  #   clone_col_name <- names(n_unique_clones)[which(n_unique_clones==max(n_unique_clones))]
  #   clone_col_idx <- which(colnames(attribute_df)==clone_col_name)
  # }
  # 
  # colnames(attribute_df)[clone_col_idx] <- "clone_id"
  
  return(list("clones"=clones, "parents"=parents, "attributes"=attribute_df, "freq_mat"=freq_mat, "max_size"=max_mutation_size, "og_colnames"=og_colnames))
  
  
}

beta_growth_fxn <- function(x, te, tm, y_max){
  ###https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4244967/
  y <- y_max*(1+ (te-x)/(te-tm))*((x/te)**(te/(te-tm)))
  return(y)
}

add_origin_mat <- function(mut_freq, tm_frac=0.6){
  ### FOR TESTING 
  # mut_freq = freq_mat
  #####
  # colnames(mut_freq) <- seq(1, ncol(mut_freq))
  time_step_size <- mean(diff(as.numeric(colnames(mut_freq))))
  ntime_steps <- ncol(mut_freq)
  te <- ntime_steps*time_step_size
  nx <- 10
  new_x <- seq(0.000001, te, length.out = nx)
  
  
  
  origin_times <- apply(mut_freq, 1, FUN = function(x){which(x>0)[1]})
  origin_at_t1_idx <- which(origin_times==1)
  size_order <- order(mut_freq[origin_at_t1_idx, 1], decreasing = TRUE)
  origin_at_t1_idx <- origin_at_t1_idx[size_order]
  
  nudge_val <- floor(nx/(length(origin_at_t1_idx)+1))
  nudge <-  0 
  origin_size_mat <- matrix(0, nrow=nrow(mut_freq), ncol = nx)
  colnames(origin_size_mat) <- as.character(new_x)
  for(idx in origin_at_t1_idx){
    initial_size <- mut_freq[idx, 1]
    
    ### temp_te gets small and smaller: 0 to nx - nudge
    temp_te <- nx - nudge
    temp_tm <- tm_frac*temp_te
    ### Also adjust tm, but do same as with above
    
    ### shift position in origin matrix by same amount temp_te is decreasing by: nude to ncol
    
    clone_new_x <- seq(0, temp_te, length.out = nx-nudge)
    origin_size_mat[idx, seq(nudge+1, nx)] <- beta_growth_fxn(clone_new_x, temp_te, temp_tm, initial_size)
    nudge <- nudge + nudge_val
  }
  
  
  colnames(origin_size_mat) <- rev(-as.numeric(colnames(origin_size_mat)))
  new_origin_mat <- origin_size_mat[, seq(1, ncol(origin_size_mat)-1)] 
  colnames(new_origin_mat) <- colnames(origin_size_mat)[seq(2, ncol(origin_size_mat))]
  new_mat <- cbind(new_origin_mat, mut_freq) #### Last column of origin_size_mat is same as first column of mut_freq
  
  return(new_mat)
  
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

get_argname <- function(paresed_name, default_name = "fill_value"){
  ### argname isn't from a dataframe

  dollar_in <- grepl("\\$", paresed_name)
  bracket_in <- grepl('\\[', paresed_name)
  if(dollar_in == FALSE & bracket_in == FALSE){
    if( grepl(',', paresed_name)){
      ## Value passed in as an array. I.e. c(1,2,3...)
      return(default_name)
    }else{
      return(paresed_name)
    }
  }
  
  ## value passed in from dataframe
  no_dollar <- tail(strsplit(paresed_name, split = "$", fixed = TRUE)[[1]], n = 1)
  xname <- strsplit(no_dollar, split = '\"', fixed = TRUE)[[1]]
  keep_idx <- which(grepl("\\[|\\]", xname) == FALSE)
  xname <- xname[keep_idx]

  return(xname)
}

rescale_frequencies <- function(freq_mat, clones, parents, filtered_idx=NULL,  data_type="size", scale_by_sizes_at_time=FALSE){
  ### Can rescale the frequencies to only show what was detected
  if(is.null(filtered_idx)){
    filtered_idx <- seq(1, nrow(freq_mat))
  }
  
  if(data_type=="size"){
    sink("/dev/null") ### Suppress output
    freq_df <- get_mutation_df(freq_mat[filtered_idx, ], clones = clones, parents = parents)
    sink()
    freq_mat <- as.matrix(freq_df)
  }else{
    freq_mat <- as.matrix(freq_mat[filtered_idx,])
  }
  
  
  max_mutation_size <- max(freq_mat)
  if(scale_by_sizes_at_time){
    max_sizes_at_each_time <- apply(freq_mat, 2, max)
    freq_mat <- sweep(freq_mat, MARGIN = 2, max_sizes_at_each_time, FUN = "/")
    
  }else{
    freq_mat <- freq_mat/max_mutation_size
  }
  
  return(freq_mat)
}
