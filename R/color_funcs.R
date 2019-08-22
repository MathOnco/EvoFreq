#'@title update_colors
#'
#'Update the colors of clones in the frequency dynamics dataframe or dendrogram plot dataframe
#'
#'@inheritParams get_evofreq
#'@param evo_freq_df Dataframe returned by \code{\link{get_evofreq}} or \code{\link{get_evogram}} , which contains the information to plot the frequency dynamics
#'@param clone_cmap String defining which colormap should be used to color the clones (nodes) if no attributes to color by. For a list of available colormaps, see \code{\link[colormap]{colormaps}}. If color not in  \code{\link[colormap]{colormaps}}, it is assumed all colors should be the same
#'@param attribute_range Range of values for the attribute to color by. If NULL, then range is determined from the attribute data
#'@param clone_id_col_in_att_df Name of the column in attribute_df that contains the clone IDs
#'@examples
#'### It is also possible to color genotypes by attributes, custom colors, or colormaps avaiable in \code{\link[colormap]{colormaps}
#'data("example.easy.wide.with.attributes")
#'### Split dataframe into clone info and size info using fact timepoint column names can be converted to numeric values
#'time_col_idx <- suppressWarnings(which(! is.na(as.numeric(colnames(example.easy.wide.with.attributes))))) 
#'attribute_col_idx <- suppressWarnings(which(is.na(as.numeric(colnames(example.easy.wide.with.attributes)))))
#'attribute_df <- example.easy.wide.with.attributes[, attribute_col_idx]
#'attr_size_df <- example.easy.wide.with.attributes[, time_col_idx]
#'attr_parents <- example.easy.wide.with.attributes$parent
#'attr_clones <- example.easy.wide.with.attributes$clone
#'clone_id_col <- "clone" ### Define which column in attribute_df contains the clone_ids
#'### Can set color using attributes. Default colormap is viridis, but can be changed to any colormap available in the colormaps package
#'attribute_pos_df <- get_evofreq(attr_size_df, attr_clones, attr_parents, attribute_df = attribute_df, attribute_val_name = "fitness", clone_id_col_in_att_df = clone_id_col)
#'fitness_evo_p <- plot_evofreq(attribute_pos_df)
#'
#'###  Update color with custom color values in hexcode format
#'pos_df_custom_color <- update_colors(attribute_pos_df, attribute_df, attribute_val_name = "color", clone_id_col_in_att_df = clone_id_col)
#'custom_cmap_evo_p <- plot_evofreq(pos_df_custom_color)
#'
#'### Change clone color based on clone ids, specifying colormap. Default colormap is soft_rainbow, but can be changed to any colormap available in the colormaps package
#'pos_df_clone_color <- update_colors(attribute_pos_df, clone_cmap = "jet")
#'clone_cmap_evo_p <- plot_evofreq(pos_df_clone_color)
#'
#'### Revert to default colors based on clone ids
#'pos_df_default_color <- update_colors(attribute_pos_df)
#'default_cmap_evo_p <- plot_evofreq(pos_df_default_color)
#'
#'### Method can also be used to update node colors of the dendrograms
#'attribute_dendro_df <- get_dendrogram(attr_size_df, attr_clones, attr_parents)
#'attribute_tree_pos <- attribute_dendro_df$dendro_pos
#'attribute_links <- attribute_dendro_df$links
#'tree_p <- plot_dendro(attribute_tree_pos, attribute_links)
#'
#'### Color node by fitness
#'tree_pos_fitness_color <- update_colors(evo_freq_df = attribute_tree_pos, attribute_df = attribute_df, attribute_val_name = "fitness", clone_id_col_in_att_df=clone_id_col)
#'tree_p_fitness <- plot_dendro(tree_pos_fitness_color, attribute_links)
#'### Color node using custom colors in hexcode format
#'tree_info_custom_color <- update_colors(evo_freq_df = attribute_tree_pos, attribute_df = attribute_df, attribute_val_name = "color", clone_id_col_in_att_df=clone_id_col)
#'tree_p_custom_color <- plot_dendro(tree_info_custom_color, attribute_links)
#'@export
update_colors <- function(evo_freq_df, clones, fill_value=NULL, clone_cmap=NULL, attribute_range=NULL, attribute_val_name=NULL){
  ### FOR TESTING ###
  # evo_freq_df <- plot_pos_df
  # attribute_range <- NULL
  # clone_cmap <- NULL #"viridis"
  # fill_value <- clone_attribute_colors
  # attribute_range <- fill_range
  ####

  if(!is.null(fill_value)){
    ### Value is an attribute to be colored by
    if(is.null(attribute_val_name)){
      ### Use name of variable passed to fill_value to 
      attribute_val_name <- deparse(substitute(fill_value))
      if(grepl("\\$", attribute_val_name)){
        ### Value was passed in  using $ to get the column, e.g. df$attribute_val_name
        attribute_val_name <- strsplit(attribute_val_name, split = "$", fixed = T)[[1]][2]
      }else if(grepl("\\[.*\\]", attribute_val_name)){
        ### Value was passed in  using a string to get the column, e.g. df[attribute_val_name]
        attribute_val_name <- strsplit(attribute_val_name, '\\"')[[1]][2]
      }
    }


    attribute_df <- data.frame("clone_id"=clones, "parents"=parents)
    attribute_df[attribute_val_name] <- fill_value

    r_colors <- colors()
    str_fill_value <- as.character(fill_value)
    all_hex <- all(sapply(str_fill_value, function(x){startsWith(x, "#") & nchar(x)> 6 & nchar(x) <=9}))
    all_rgb <- all(sapply(str_fill_value, function(x){length(strsplit(x, ",")[[1]])==3}))
    all_named <- all(sapply(str_fill_value, function(x){x %in% r_colors}))
    
    user_defined_colors <- any(all_hex, all_rgb, all_named)
    if(user_defined_colors){
      attribute_val_name <- NULL
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
      if(is.null(attribute_range)){
        attribute_range <- range(fill_value, na.rm = T)
      }
      clone_color <- get_attribute_colors(fill_value, min_x = attribute_range[1], max_x = attribute_range[2], cmap=clone_cmap)  
      
    }
    attribute_df <- data.frame("clone_id"=clones, "parents"=parents, "plot_color"=clone_color)
    if(!user_defined_colors){
      attribute_df[attribute_val_name] <- fill_value
    }
    
  }else{
    ### Assign unique color for each clone
    attribute_val_name <- NULL
    if(is.null(clone_cmap)){
      clone_cmap <- 'rainbow_soft'      
    }
    clone_color <- get_clone_color(length(clones), clone_cmap)
    attribute_df <- data.frame("clone_id"=clones, "parents"=parents, "plot_color"=clone_color)
  }

  attribute_df$cmap <- clone_cmap
  
  evo_freq_df <- evo_freq_df[! colnames(evo_freq_df) %in% c("efp_color_attribute", "plot_color", "cmap")] ### remove previous color and attribute names
  updated_attribute_df <- merge(evo_freq_df, attribute_df, all.x = T, all.y = F, by="clone_id")
  
  ### If any duplicated columns, keep those that were originally in X
  dup_x_col_idx <- grep("\\.x", colnames(updated_attribute_df))
  if(length(dup_x_col_idx)>0){
    dup_y_col_idx <- grep("\\.y", colnames(updated_attribute_df))
    new_x_colnames <- colnames(updated_attribute_df)[dup_x_col_idx]
    new_x_colnames <- gsub( "\\.x", "", new_x_colnames)
    colnames(updated_attribute_df)[dup_x_col_idx] <- new_x_colnames
    updated_attribute_df <- updated_attribute_df[-dup_y_col_idx]    
  }
  
  if(is.null(attribute_val_name)){
    attribute_val_name <- NA
  }
  updated_attribute_df$efp_color_attribute <- attribute_val_name
  updated_attribute_df$plot_color <- as.character(updated_attribute_df$plot_color)
  updated_attribute_df$cmap <- clone_cmap
  return(updated_attribute_df)
}

get_clone_color <- function(n_clones, cmap="rainbow_soft"){
  ### FOR TESTING ###
  # n_clones <- length(clones)
  ####
  pop_colors <-  colormap::colormap(colormap=colormap::colormaps[cmap][[1]], nshades=n_clones+1)[1:n_clones]
  return(pop_colors)
}

get_attribute_colors <- function(x, min_x =NULL, max_x = NULL, n_color_bins = 100, cmap="viridis"){
  ### FOR TESTING ###
  # x <- clone_antigenicity
  # root_idx <- which(clone_antigenicity==0)
  ####
  if(is.null(min_x)){
    min_x <- min(x)
  }
  if(is.null(max_x)){
    max_x <- max(x)
  }
  
  bin_breaks <- seq(min_x, max_x , length.out = n_color_bins)
  bin_number <- findInterval(x, bin_breaks, rightmost.closed = F)
  cmap_colors <- colormap::colormap(colormaps[cmap][[1]], nshades=n_color_bins)
  colors <- cmap_colors[bin_number]
  return(colors)
}

scale_value <- function(x, in_min, in_max, out_min, out_max){
  new_x <- ((out_max - out_min)*(x-in_min))/(in_max - in_min) + out_min
  return(new_x)
}

rgb2hex <- function(cvals){
  ### FOR TESTING ###
  # cvals <- sapply(seq(1, length(clones)), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
  ### FOR TESTING ###
  
  fill_mat <- t(sapply(cvals, function(x){as.numeric(strsplit(x, ",")[[1]])}))
  max_fill_val <- max(fill_mat)
  if(max_fill_val <= 1){
    max_fill_val <- 1
  }else{
    max_fill_val <- 255
  }
  
  hex_colors <- sapply(seq(nrow(fill_mat)), function(x){rgb(fill_mat[x, 1], fill_mat[x, 3], fill_mat[x, 3], maxColorValue=max_fill_val, alpha=max_fill_val)})
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