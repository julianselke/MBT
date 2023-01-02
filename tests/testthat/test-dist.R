load(file.path('expected_outputs_dist.RData'))
set.seed(42)

FT <- rlnorm(1000, 2, 4) %>% round() %>% matrix(., 10, 100) %>%
  `row.names<-`(., LETTERS[1:10]) %>% otu_table(., taxa_are_rows = FALSE)
SD <- cbind(1:10,LETTERS[1:10]) %>% as.data.frame() %>%
  `row.names<-`(LETTERS[1:10]) %>% sample_data(.)
TR <- ape::rtree(100, tip.label = taxa_names(FT), rooted = TRUE) %>% phy_tree()
x <- phyloseq::merge_phyloseq(FT,SD,TR)

test_that("dist matches expectes dist", {
  expect_equal(MBT:::dist(x, dis = "euclidian"), expected_dist_euclidian)
  expect_equal(MBT:::dist(x, dis = "bray"), expected_dist_bray)
  expect_equal(MBT:::dist(x, dis = "robust.aitchison"),
               expected_dist_robust.aitchison)
  expect_equal(suppressMessages(MBT:::dist(x, dis = "philr")),
               expected_dist_philr)
  expect_equal(MBT:::dist(x, dis = "w.unifrac"), expected_dist_w.unifrac)
})






