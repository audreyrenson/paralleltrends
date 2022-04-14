test_that("data dimensions are correct (wide)", {

  N = 100
  Tt = 3
  Beta = generate_parameters(Tt=Tt, range_ymeans=stats::qlogis(c(0.01, 0.05)))

  correct_dimensions <- c(N, (length(Beta) - 1)*(Tt+1) + 2) # Tt+1 periods of variables all variables expect U0, plus U0 and uid

  df_normal <- generate_data(N=N, Tt=Tt, Beta=Beta)
  df_binomial <- generate_data(N=N, Tt=Tt, Beta=Beta, ylink='rbinom_logit', binomial_n = round(runif(N, min=20, max=100)))

  expect_equal(dim(df_normal), correct_dimensions)
  expect_equal(dim(df_binomial), correct_dimensions + c(0, 1)) #has additional binomial_n column, but same # rows
})

test_that('data dimensions are correct (long)', {
  N=100
  Tt=3
  Beta = generate_parameters(Tt=Tt, range_ymeans = stats::qlogis(c(0.01, 0.05)))
  df_normal <- generate_data(N=N, Tt=Tt, Beta=Beta, long=TRUE)
  expect_equal(dim(df_normal), c(N*(Tt+1), length(Beta) + 2))
})
