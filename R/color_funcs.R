
data(cmap_hex_list)

get_colors <- function(n_clones, cmap_name){
  ### FOR TESTING ###
  # n_clones <- 10
  # cmap_name = "rainbow_soft"
  #####
  cmap_hex_vals <- cmap_hex_list[[cmap_name]]
  color_ramp_fxn <- colorRampPalette(cmap_hex_vals)
  cmap_vals <- color_ramp_fxn(n_clones)
  return(cmap_vals)
}

#'@title update_colors
#'
#'Update the colors of clones in the frequency dynamics dataframe or dendrogram plot dataframe
#'
#'@inheritParams get_evofreq
#'@param evo_freq_df Dataframe returned by \code{\link{get_evofreq}} or \code{\link{get_evogram}} , which contains the information to plot the frequency dynamics
#'@param clone_cmap String defining which colormap should be used to color the clones (nodes) if no attributes to color by. For a list of available colormaps, see \code{\link[colormap]{colormaps}}. If color not in  \code{\link[colormap]{colormaps}}, it is assumed all colors should be the same
#'@param fill_name String defining the name of the attribute used for coloring. If NULL, the default, the name will be inferred from the \code{fill_value} argument. 
#'@examples
#' data("example.easy.wide")
#' ### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#' time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide)))))
#' size_df <- example.easy.wide[, time_col_idx]
#' parents <- example.easy.wide$parents
#' clones <- example.easy.wide$clones
#' ### Setting of colors can be done when getting the freq_frame, or by updating the color later using \code{\link{update_colors}}. Available colormaps are those found in \code{\link[colormap]{colormaps}}
#' ### Default colormap is rainbow_soft, but this can be changed using the \code{clone_cmap} argument.
#' freq_frame <- get_evofreq(size_df, clones, parents)
#' evo_p <- plot_evofreq(freq_frame)
#'   
#' ### Can color each clone by an attribute by providing a \code{fill_value}. Default colormap is viridis, but this can be changed using the \code{clone_cmap} argument. There is also the option to set the color range using the \code{fill_range} argument
#' fitness <- runif(length(clones),0, 100)
#' fitness_freq_frame <- update_colors(freq_frame, clones = clones, fill_value = fitness, fill_range= c(0, 100))
#' fitness_evo_p <- plot_evofreq(fitness_freq_frame,  fill_range = c(0, 100))
#' 
#' ### The user can also provide custom colors for each clone, which will need to be passed into the \code{fill_value} argument
#' ### Custom colors can be defined using RGB values. Each color should be a string specifying the color channel values, separated by commas.
#' rgb_clone_colors <- sapply(seq(1, length(clones)), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
#' rgb_freq_frame <- update_colors(freq_frame, clones = clones, fill_value = rgb_clone_colors)
#' rgb_evo_p <- plot_evofreq(rgb_freq_frame)
#' 
#' ### Custom colors can also be any of the named colors in R. A list of the colors can be found with \code{colors()}
#' named_clone_colors <- sample(colors(), length(clones), replace = FALSE)
#' named_freq_frame <- update_colors(freq_frame, clones = clones, fill_value = named_clone_colors)
#' named_evo_p <- plot_evofreq(named_freq_frame)
#' 
#' ### Custom colors can also be specified using hexcode
#' hex_clone_colors <- c("#614099ff", "#1d347eff", "#94558aff", "#c96872ff", "#f1884dff", "#e8fa5bff", "#042333ff","#f9bb41ff")
#' hex_freq_frame <- update_colors(freq_frame, clones = clones, fill_value = hex_clone_colors)
#' hex_evo_p <- plot_evofreq(hex_freq_frame)
#'
#' ### Method can also be used to update node colors of the dendrograms
#' dendro_df <- get_evogram(size_df, clones, parents)
#' tree_pos <- dendro_df$dendro_pos
#' tree_links <- dendro_df$links
#' tree_p <- plot_evogram(tree_pos, tree_links)
#'
#' ### Color node by fitness
#' tree_pos_fitness_color <- update_colors(evo_freq_df = tree_pos, clones = clones, fill_value = fitness)
#' tree_p_fitness <- plot_evogram(tree_pos_fitness_color, tree_links)
#' ### Color node using custom colors in hexcode format
#' tree_info_custom_color <- update_colors(evo_freq_df = tree_pos, clones = clones, fill_value = hex_clone_colors)
#' tree_p_custom_color <- plot_evogram(tree_info_custom_color, tree_links)
#'@export
update_colors <- function(evo_freq_df, clones, fill_value=NULL, clone_cmap=NULL, fill_range=NULL, fill_name=NULL, shuffle_colors=FALSE){
  ### FOR TESTING ###
  # evo_freq_df <- d_pos_df
  # attribute_range <- NULL
  # clone_cmap <- NULL #"viridis"
  # fill_value <- fill_value
  # attribute_range <- fill_range
  # fill_name <- fill_name
  ####
  if(!is.null(fill_value)){
    ### Value is an attribute to be colored by
    if(is.null(fill_name)){
      fill_name <- colnames(fill_value) ### Value was passed in  using a string to get the column, e.g. df[fill_name]
      if(is.null(fill_name)){
        paresed_fill_name <- deparse(substitute(fill_value))
        fill_name <- get_argname(paresed_fill_name)
      }      
    }
    attribute_df <- data.frame("clone_id"=clones)
    attribute_df[fill_name] <- fill_value

    r_colors <- colors()
    str_fill_value <- as.character(fill_value)
    all_hex <- all(sapply(str_fill_value, function(x){startsWith(x, "#") & nchar(x)> 6 & nchar(x) <=9}))
    all_rgb <- all(sapply(str_fill_value, function(x){length(strsplit(x, ",")[[1]])==3}))
    all_named <- all(sapply(str_fill_value, function(x){x %in% r_colors}))
    
    user_defined_colors <- any(all_hex, all_rgb, all_named)
    if(user_defined_colors){
      fill_name <- NULL
      clone_cmap <- "custom"
      if(all_hex){
        clone_color <- fill_value
      }else if(all_rgb){
        clone_color <- rgb2hex(fill_value)
      }else if(all_named){
        clone_color <- named2hex(fill_value)
      }
      
    }else{
      ### If values are not hex or cannot be converted to hex, then assume user wants colors generated using a colormap
      
      if(is.null(clone_cmap)){
        clone_cmap <- 'viridis'      
      }
      if(is.null(fill_range)){
        fill_range <- range(fill_value, na.rm = TRUE)
      }
      clone_color <- get_attribute_colors(fill_value, min_x = fill_range[1], max_x = fill_range[2], cmap=clone_cmap)
    }
    # attribute_df <- data.frame("clone_id"=clones, "parents"=parents, "plot_color"=clone_color)
    attribute_df <- data.frame("clone_id"=clones, "plot_color"=clone_color)
    if(!user_defined_colors){
      attribute_df[fill_name] <- fill_value
    }
    
  }else{
    ### Assign unique color for each clone
    fill_name <- NULL
    if(is.null(clone_cmap)){
      clone_cmap <- 'rainbow_soft'      
    }
    clone_color <- get_clone_color(length(clones), clone_cmap)
    if(shuffle_colors){
      clone_color <- sample(clone_color)
    }
    # attribute_df <- data.frame("clone_id"=clones, "parents"=parents, "plot_color"=clone_color)
    attribute_df <- data.frame("clone_id"=clones, "plot_color"=clone_color)
  }

  attribute_df$cmap <- clone_cmap
  
  evo_freq_df <- evo_freq_df[! colnames(evo_freq_df) %in% c("efp_color_attribute", "plot_color", "cmap")] ### remove previous color and attribute names
  updated_attribute_df <- merge(evo_freq_df, attribute_df, all.x = TRUE, all.y = FALSE, by="clone_id")
  
  ### If any duplicated columns, keep those that were originally in X
  dup_x_col_idx <- grep("\\.x", colnames(updated_attribute_df))
  if(length(dup_x_col_idx)>0){
    dup_y_col_idx <- grep("\\.y", colnames(updated_attribute_df))
    new_x_colnames <- colnames(updated_attribute_df)[dup_x_col_idx]
    new_x_colnames <- gsub( "\\.x", "", new_x_colnames)
    colnames(updated_attribute_df)[dup_x_col_idx] <- new_x_colnames
    updated_attribute_df <- updated_attribute_df[-dup_y_col_idx]    
  }
  
  if(is.null(fill_name)){
    fill_name <- NA
  }
  updated_attribute_df$efp_color_attribute <- fill_name
  updated_attribute_df$plot_color <- as.character(updated_attribute_df$plot_color)
  updated_attribute_df$cmap <- clone_cmap
  return(updated_attribute_df)
}

get_clone_color <- function(n_clones, cmap="rainbow_soft"){
  ### FOR TESTING ###
  # n_clones <- length(clones)
  ####
  # pop_colors <-  colormap::colormap(colormap=colormap::colormaps[cmap][[1]], nshades=n_clones+1)[1:n_clones]
  # n_clones, cmap_name
  # pop_colors <-  get_colors(n_clones+1, cmap)[1:n_clones]
  pop_colors <-  get_colors(n_clones+1, cmap)[1:n_clones]
  return(pop_colors)
}

get_attribute_colors <- function(x, min_x =NULL, max_x = NULL, n_color_bins = 100, cmap="viridis"){
  ### FOR TESTING ###
  # x <- clone_antigenicity
  # x <- fill_value
  # min_x <- fill_range[1]
  # max_x <- fill_range[2] 
  # cmap <- clone_cmap
  # n_color_bins <- 100
  ####
  if(is.null(min_x)){
    min_x <- min(x)
  }
  if(is.null(max_x)){
    max_x <- max(x)
  }
  
  bin_breaks <- seq(min_x, max_x , length.out = n_color_bins)
  bin_number <- findInterval(x, bin_breaks, rightmost.closed = FALSE, all.inside=TRUE)
  
  # cmap_colors <- colormap::colormap(colormaps[cmap][[1]], nshades=n_color_bins)
  cmap_colors <- get_colors(n_color_bins, cmap)
  colors <- cmap_colors[bin_number]
  return(colors)
}

scale_value <- function(x, in_min, in_max, out_min, out_max){
  new_x <- ((out_max - out_min)*(x-in_min))/(in_max - in_min) + out_min
  return(new_x)
}

rgb2hex <- function(cvals){
  ### FOR TESTING ###
  # cvals <- sapply(seq(1, 10), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
  # cvals <- sapply(seq(1, 10), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=" ")})
  ### FOR TESTING ###
  delimiters <- c(" ", ",")
  delimiter <- delimiters[which(sapply(delimiters, function(x){grepl(x, cvals[1])}))]
  fill_mat <- t(sapply(cvals, function(x){as.numeric(strsplit(x, delimiter)[[1]])}))
  max_fill_val <- max(fill_mat)
  if(max_fill_val <= 1){
    max_fill_val <- 1
  }else{
    max_fill_val <- 255
  }
  
  hex_colors <- sapply(seq(nrow(fill_mat)), function(x){rgb(fill_mat[x, 1], fill_mat[x, 2], fill_mat[x, 3], maxColorValue=max_fill_val, alpha=max_fill_val)})
  return(hex_colors)
}

named2hex <- function(cvals){
  ### FOR TESTING ###
  # cvals <- sample(colors(), 10)
  ####
  
  fill_mat <- t(sapply(cvals, function(x){as.numeric(col2rgb(x))})) ### Each column is a color
  max_fill_val <- max(fill_mat)
  if(max_fill_val <= 1){
    max_fill_val <- 1
  }else{
    max_fill_val <- 255
  }
  
  hex_colors <- sapply(seq(nrow(fill_mat)), function(x){rgb(fill_mat[x, 1], fill_mat[x, 3], fill_mat[x, 3], maxColorValue=max_fill_val, alpha=max_fill_val)})
  
  return(hex_colors)
  
}