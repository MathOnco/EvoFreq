context("get_freq_dynamics testing")

test_that("get_freq_dynamics evaluation on easy examples.",{
  data("example.easy.wide")
  clone_dyn_df.wide <- get_freq_dynamics(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma",scale_by_sizes_at_time = F)
  clone_dyn_df.wide.scale <- get_freq_dynamics(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma", scale_by_sizes_at_time = T)
  val <- clone_dyn_df.wide.scale[66,3]
  val2 <- clone_dyn_df.wide[26,3]
  expect_equal(colnames(example.easy.wide), c("parents","clones","1","2","3","4","5","6","7","8"))
  expect_equal(unique(unique(as.numeric(clone_dyn_df.wide$draw_order))==c(1,2,3,5,4,6,7,8)), TRUE)
  expect_equal(unique(as.numeric(clone_dyn_df.wide$origin_time)), c(1,2,6,8))
  expect_equal(7.727264e-18*1e18, expected = val*1e18, tolerance = 0.001, scale = val)
  expect_equal(335.661975, expected = val2, tolerance = 0.001, scale = val2)
  clone_dyn_df.wide <- get_freq_dynamics(example.easy.wide[,seq(3,10)], example.easy.wide$clones, example.easy.wide$parents, clone_cmap = "magma",scale_by_sizes_at_time = T, threshold=0.1)
  expect_equal(0.0, expected = clone_dyn_df.wide[75,3])
})