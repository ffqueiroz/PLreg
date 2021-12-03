


test_that("GJS dist with mu equal 0.5 is symetric", {
  expect_equal(pPL(0.5, mu = 0.5, sigma = 2, lambda = 1, family = "NO"), 0.5)
  expect_equal(pPL(0.5, mu = 0.5, sigma = 2, lambda = 1, family = "SLASH"), 0.5)
  expect_equal(pPL(0.5, mu = 0.5, sigma = 2, lambda = 1, family = "LO"), 0.5)
})
