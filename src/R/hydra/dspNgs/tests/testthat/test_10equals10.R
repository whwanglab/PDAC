context("Hello World")

# Useful references for making your own.
# https://r-pkgs.org/tests.html
# https://testthat.r-lib.org/reference/index.html

test_that("10 equal 10 (testing workflow is working)", {
  expect_equal(-15+25, 10)
  expect_true(10==10)
  expect_false(10!=10)
})
