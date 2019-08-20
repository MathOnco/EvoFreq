context("evofreq_clone_labels testing")

test_that("evofreq_clone_labels evaluation on easy examples.",{
  data("example.easy.wide") # Load Data
  clone_dyn_df <- get_freq_dynamics(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma")
  evo_freq_p <- plot_freq_dynamics(clone_dyn_df)
  clone_labels <- evofreq_clone_labels(clone_dyn_df, extant_only=T, evo_freq_p=evo_freq_p, apply_labels=T)
  expect_equal(class(clone_labels), c("gg","ggplot"))
  clone_labels <- evofreq_clone_labels(clone_dyn_df, extant_only=T, evo_freq_p=evo_freq_p, apply_labels=F)
  expect_equal(clone_labels$clone_id, c(1,3,4,5,6,7,8))
  expect_equal(round(clone_labels$y_plotValue[4],4), 802.9992)
  clone_labels <- evofreq_clone_labels(clone_dyn_df, extant_only=F, evo_freq_p=evo_freq_p, apply_labels=F)
  expect_equal(clone_labels$clone_id, c(1,2,3,4,5,6,7,8))
  expect_equal(round(clone_labels$y_plotValue[2],4), 433.7922)
  clone_labels <- evofreq_clone_labels(clone_dyn_df, extant_only=F, evo_freq_p=evo_freq_p, apply_labels=F, clone_list = c(2))
  expect_equal(length(clone_labels$clone_id), 1)
})

