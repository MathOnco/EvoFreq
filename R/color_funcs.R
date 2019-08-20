
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
update_colors <- function(evo_freq_df, attribute_df=NULL, attribute_val_name=NULL, clone_id_col_in_att_df="clone_id", clone_cmap="rainbow_soft", attribute_range=NULL){
  ### FOR TESTING ###
  # evo_freq_df <- dendro_pos
  # attribute_df <- example.easy.long.edges
  # attribute_val_name <- "fitness"
  # attribute_range <- NULL
  # clone_id_col_in_att_df="clone"
  # clone_cmap <- "viridis"
  
  ###
  if(is.null(attribute_df)){
    # attribute_val_name <- NA
    unique_clone_idx <- which(duplicated(evo_freq_df$clone_id)==F)
    unc <-evo_freq_df$clone_id[unique_clone_idx]
    
    if(clone_cmap %in% names(colormaps)){
      clone_color <- get_clone_color(length(unc), clone_cmap)
    }else{
      clone_color <- clone_cmap
    }
    attribute_df <- data.frame("clone_id"=unc, "plot_color"=clone_color)
  }else{
    
    if(!"clone_id" %in% colnames(attribute_df)){
      clone_col_idx <- which(colnames(attribute_df)==clone_id_col_in_att_df)
      if(length(clone_col_idx) == 0){
        print("do not know column in attribute_df that contains clone IDs. Please check the clone_id_col_in_att_df argument or set the clone ID column in attribute_df to be 'clone_id'")
      }
      colnames(attribute_df)[clone_col_idx] <- "clone_id"
    }

    # attribute_df$efp_color_attribute <- attribute_val_name
    attribute_df <- attribute_df[attribute_df$clone_id %in% evo_freq_df$clone_id, ]
    if(!is.null(attribute_val_name)){
      attr_vals <- attribute_df[, attribute_val_name]
      if(is.factor(attr_vals)){
        attr_vals <- as.character(attr_vals)
      }
      all_hex <- all(sapply(attr_vals, function(x){substr(x, 1, 1)=="#" & nchar(x)> 6 & nchar(x) <=9}))
      if(all_hex){
        attribute_df$plot_color <- attr_vals
        attribute_val_name <- "plot_color"
        # attribute_df$efp_color_attribute <- NA
      }else{
        if(is.null(attribute_range)){
          attribute_range <- range(attr_vals, na.rm = T)
        }
        attribute_df$plot_color <- get_attribute_colors(attr_vals, min_x = attribute_range[1], max_x = attribute_range[2], cmap=clone_cmap)  
      }
    }
  }
  
  evo_freq_df <- evo_freq_df[! colnames(evo_freq_df) %in% c("efp_color_attribute", "plot_color")] ### remove previous color
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

