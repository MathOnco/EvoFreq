% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvoFreq_funcs.R
\name{long_to_wide_freq_ready}
\alias{long_to_wide_freq_ready}
\title{Convert data in long format to wide format
\code{long_to_wide_freq_ready} Converts long data frame to wide format used by evofreq.}
\usage{
long_to_wide_freq_ready(edges_df, long_pop_sizes_df, time_col_name,
  clone_col_name, parent_col_name, size_col_name,
  fill_gaps_in_size = FALSE)
}
\arguments{
\item{edges_df}{Dataframe with 2 columns, where each row defines relation between each clone and their parent. One column contains the clone ids, and the other contains the id of that clone's parent.}

\item{long_pop_sizes_df}{Dataframe with 3 columns, one containing the clone id, one containing the time points, and one containing the size of the clone at that time point. There must be a size for each clone at every time point, so it did not exist, it's size is 0.}

\item{time_col_name}{String that specifies which column in long_pop_sizes_df contains the timepoint information}

\item{clone_col_name}{String that specifies which column in in edges_df contains the names of the clones. I.e. the "to" node in a graph}

\item{parent_col_name}{String that specifies which column in in edges_df contains the names of the parents of the descendent clone in the edges_df with colname clone_col_name. #I.e. the "from" node in a graph}

\item{size_col_name}{String that specifies which column in long_pop_sizes_df contains the size of each clone for the corresponding timepoint and clone in long_pop_sizes_df}

\item{fill_gaps_in_size}{Boolean specificing if gaps in sizes over time should be filled in, assuming linear changes.}
}
\value{
List containing: wide_size_df, a dataframe of the sizes in wide format, where each row contains the sizes over time for a single clone; parents, a vector containing the parent id for each clone in clones; clones, a vector of clones ids correspond to each row in wide_size_df
}
\description{
Convert data in long format to wide format
\code{long_to_wide_freq_ready} Converts long data frame to wide format used by evofreq.
}
