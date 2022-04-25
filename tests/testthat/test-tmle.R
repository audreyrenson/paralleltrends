test_that("tmle wide and long return same result (currently only checking dim)", {
  Tt = 3
  Beta = generate_parameters(Tt=Tt, range_ymeans=stats::qlogis(c(0.01, 0.05)))
  N = 100

  df_wide = generate_data(N, Tt, Beta)
  df_wide_interv = df_wide %>% dplyr::mutate(dplyr::across(dplyr::starts_with('A'), ~ 0))
  truth = generate_data(N*100, Tt, Beta, TRUE) %>% dplyr::select(dplyr::starts_with('Y')) %>%
    colMeans() %>%
    diff()

  df_long = pivot_longer_gendata(df_wide) %>%
    dplyr::mutate(W22=W2^2) %>%
    dplyr::group_by(.data[['uid']]) %>%
    append_lags(3, c('A','W1','W2','W22', 'Y'), default=NA) %>%
    dplyr::ungroup() %>%
    make_time_dummies('t')


  df_long_interv = df_long %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with('A'), ~0))

  tmle_wide = tmle(df_obs=df_wide,
                   df_interv=df_wide_interv,
                   den_formula = 'A{t}~A{t-1}*(W1{t}+W2{t}+I(W2{t}^2))',
                   inside_formula_t = '~A{t}*(W1{t}+W2{t}+I(W2{t}^2))',
                   inside_formula_tmin1 = '~A{t}*(W1{t-1}+W2{t-1}+I(W2{t-1}^2))',
                   outside_formula = '~A{k}*(W1{k}+W2{k}+I(W2{k}^2))',
                   Tt=Tt,
                   models=FALSE,
                   suppress_rank_warnings=TRUE,
                   long=FALSE)

  tmle_long = tmle(df_obs=df_long,
                   df_interv=df_long_interv,
                   den_formula = 'A ~ ({tvars})*(W1+W2+W22)*A_lag1',
                   inside_formula_t = 'Y~-1 + {tvars}*(W1+W2+W22)*A - A*(W1+W2+W22)',
                   inside_formula_tmin1 ='Y_lag1~ -1 + {tvars}*(W1_lag1+W2_lag1+W22_lag1)*A - (W1_lag1+W2_lag1+W22_lag1)*A',
                   outside_formula = 'preds ~ -1 + {tvars}*(W1_lag{n}+W2_lag{n}+W22_lag{n})*A_lag{n} -(W1_lag{n}+W2_lag{n}+W22_lag{n})*A_lag{n} - t{n}:A_lag{n} - t{n}:(W1_lag{n}+W2_lag{n}+W22_lag{n}):A_lag{n} ',
                   Tt=Tt,
                   t_col=df_long$t,
                   id = 'uid',
                   models=FALSE,
                   suppress_rank_warnings=TRUE,
                   long=TRUE)


  expect_equal(dim(tmle_wide), dim(tmle_long))
})
