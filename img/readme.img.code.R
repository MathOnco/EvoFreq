library(gridExtra)

# Note: this can be copy and pasted after installed

data("example.easy.wide") # Load a simple Data Frame example
str(example.easy.wide) # Inspect the data structure

# 'data.frame':	8 obs. of  10 variables:
#  $ parents: num  0 1 1 3 1 5 5 5
#  $ clones : num  1 2 3 4 5 6 7 8
#  $ 1      : num  1 0 0 0 0 0 0 0
#  $ 2      : num  100 5 0 0 0 0 0 0
#  $ 3      : num  200 100 5 0 0 0 0 0
#  $ 4      : num  400 5 100 1 1 0 0 0
#  $ 5      : num  0 0 200 100 100 1 0 1
#  $ 6      : num  0 0 200 125 200 10 1 15
#  $ 7      : num  0 0 300 200 300 20 10 25
#  $ 8      : num  0 0 300 300 300 25 25 100

# You have A column of parents and a column of clones then you have a column for each of the timepoints with sizes for that clone.

# Then get the frequency data. (Use ?get_freq_dynamics for options)
clone_dynamics_df <- get_freq_dynamics(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma")

# Create the plot (shown on the left below)
evo_freq_p <- plot_freq_dynamics(clone_dynamics_df)
print(evo_freq_p)

# We can also choose to update the colors or do this during the first creation. (shown on the right below)
clone_dynamics_df_jet <- update_colors(clone_dynamics_df, clone_cmap = "jet")
evo_freq_p_jet <- plot_freq_dynamics(clone_dynamics_df_jet)
print(evo_freq_p_jet)

p <- arrangeGrob(evo_freq_p, evo_freq_p_jet, ncol=2, nrow=1)
ggsave(plot=p, filename="img/easy.wide.image.png", units=c("mm"),width=317.5, height=79.375, dpi=300, device="png")

# Note: this can be copy and pasted after installed

# In this example you have two files. One is the edge list of clones and their parents
# The other file is the sizes over time for those clones
# examine the structures to see how to format
# example.easy.long.edges now has an attribute column as well
data("example.easy.long.edges")
data("example.easy.long.sizes")

# Use the long_to_wide_size_df function to get the right data structure.
wide_df <- long_to_wide_size_df(long_pop_sizes_df = example.easy.long.sizes, time_col_name = "Time", clone_col_name = "clone", parent_col_name = "parent", size_col_name = "Size", edges_df = example.easy.long.edges)
clones <- wide_df$clones
parents <- wide_df$parents
size_df <- wide_df$wide_size_df

clone_dynamics_df <- get_freq_dynamics(size_df, clones, parents, clone_cmap = "inferno")
evo_freq_p <- plot_freq_dynamics(clone_dynamics_df)

# Add custom ggplot features
evo_freq_labeled_p <- evofreq_clone_labels(clone_dynamics_df, clone_list=c(3,5,6), extant_only=F, evo_freq_p=evo_freq_p, apply_labels=T)

p <- arrangeGrob(evo_freq_p, evo_freq_labeled_p, ncol=2, nrow=1)
ggsave(plot=p, filename="img/easy.long.image.png", units=c("mm"),width=317.5, height=79.375, dpi=300, device="png")

# Network plots
data("example.wide")
time_col_idx <- which(substr(colnames(example.wide), start=1, stop = 1)=="X")
attribute_col_idx <- which(substr(colnames(example.wide), start=1, stop = 1) != "X")
colnames(example.wide)[time_col_idx] <- substring(colnames(example.wide), first=2)[time_col_idx]
attribute_df <- example.wide[, attribute_col_idx]
size_df <- example.wide[, time_col_idx]
parents <- example.wide$parent
clones <- example.wide$clone

clone_network <- get_clone_network(size_df, clones, parents)
evofreq_plot <- plot_clone_network(clone_network, layout="circlepack")


