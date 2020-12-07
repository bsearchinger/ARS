# Tests for ars functions
library(testthat)

################################################################################
# Tests for approxD: Derivative Function

# Symbolic Differentiation Checks
dtest_log <- Deriv::Deriv(f = log, nderiv = 1)
dtest_log2 <- Deriv::Deriv(f = log, nderiv = 2)
dtest_exp <- Deriv::Deriv(exp)
dtest_norm <- Deriv::Deriv(dnorm)
dtest_norm2 <- Deriv::Deriv(dnorm, nderiv = 2)

test_that(
  "Derivative function evaluates correctly",{
    expect_equivalent(
      approxD(
        f = log,
        x = 2,
        f_params = NULL,
        n = 1),
      dtest_log(2)[1])
    expect_equivalent(
      approxD(
        f = log,
        x = 2,
        f_params = NULL,
        n = 2),
      dtest_log2(2)[1])
    expect_equivalent(
      approxD(
        f = exp,
        x = 2,
        f_params = NULL,
        n = 1),
      dtest_exp(2)[1])
    expect_equivalent(
      approxD(
        f = dnorm,
        x = 0.5,
        f_params = list(
          mean = 0,
          sd = 1
        ),
        n = 1),
      dtest_norm(0.5)[1])
    expect_true(
      sign(approxD(
        f = dnorm,
        x = 0.5,
        f_params = list(
          mean = 0,
          sd = 1
        ),
        n = 2)) == sign(
          dtest_norm2(0.5)[1])
    )
  }
)

test_that(
  "Derivative Function throws an error", {
    expect_error(
      approxD(
        f = dnorm,
        f_params = list(
          mean = 0,
          lambda = 1),
        x = 0.5,
        n = 1))
  }
)
