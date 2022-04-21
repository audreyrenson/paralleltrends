safe_predict = function(...) {
  withCallingHandlers({
    predict(...)
  }, warning = function(w) {
    #the non-integer successes warning happens anytime you weight a binomial or other discrete glm, and isn't actually a problem
    if (startsWith(conditionMessage(w), "prediction from a rank-deficient fit may be misleading"))
      invokeRestart("muffleWarning")
  })
}
