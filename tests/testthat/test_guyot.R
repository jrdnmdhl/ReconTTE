test_that("Test that guyot implementation matches original", {
  ex1 <- readRDS(system.file("inst", "example_1", "digitized.rds", package = "ReconTTE"))
  original_guyot <- guyot_old(ex1$survival, ex1$n_at_risk)
  new_guyot <- guyot(ex1$survival, ex1$n_at_risk)
  expect_equal(new_guyot, original_guyot)
})
