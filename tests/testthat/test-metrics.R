test_that("compute_rafe stub errors until implemented", {
  expect_error(compute_rafe(rnorm(3), rnorm(3), Sigma = diag(3)),
               "TODO: implement compute_rafe")
})

test_that("compute_crafe stub errors until implemented", {
  expect_error(compute_crafe(diag(3), diag(3)),
               "TODO: implement compute_crafe")
})

test_that("compute_trafe stub errors until implemented", {
  expect_error(
    compute_trafe(rnorm(3), rnorm(3), diag(3), diag(3)),
    "TODO: implement compute_trafe"
  )
})
