test_that("parameter dimensions are correct", {
  Tt = 10
  Beta = generate_parameters(Tt, range_ymeans = qlogis(c(0.01, 0.05)))

  check_dimensions <- function(Beta, Tt) {
    expect_equal(length(Beta$U), 1)
    expect_equal(dim(Beta$A), c(Tt + 1, 5))
    expect_equal(dim(Beta$W1),c(Tt + 1, 2))
    expect_equal(dim(Beta$W2),c(Tt + 1, 2))
    expect_equal(dim(Beta$Y), c(Tt + 1, 6))
  }

  check_dimensions(Beta, Tt)
})
