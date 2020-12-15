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

################################################################################
# Tests for tanIntersect: Computing intersections of tangent lines
ztest <- c(-0.5, -0.1, 0, 0.1, 0.5)
xtest <- seq(0, 1, length.out = 20)
test_that("output of tanIntersect is symmetric for symmetric distribution",{
  expect_true(
    all(
      tanIntersect(ztest, dnorm)$z[1:2] == -1*rev(tanIntersect(ztest, dnorm)$z[3:4])
      ) == TRUE
    )
})

test_that("output of tanIntersect is corrent dimensions", {
  expect_true(
    length(
      tanIntersect(
        x_abs = xtest,
        f = dexp,
        f_params = list(rate = 2))$z
      ) == 19
  )
  expect_true(
    length(
      tanIntersect(
        x_abs = xtest,
        f = dexp,
        f_params = list(rate = 2)
        )) == 4
  )
})

################################################################################
# Tests for checkThat: argument checking function

# check that takes a quosure
dexp_quo <- rlang::quo(dexp)
dnorm_quo <- rlang::quo(dnorm)
dgamma_quo <- rlang::quo(dgamma)

test_that("checkThat runs silently for correct imputs", {
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, sd = 1),
      starting_values = c(-0.5, -0.1, 0, 0.1, 0.5),
      sample_size = 30
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params = NULL,
      starting_values = c(-0.5, -0.1, 0, 0.1, 0.5),
      sample_size = 30
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params =  list(mean = 0, sd = 1),
      starting_values = NULL,
      sample_size = 30
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params =  NULL,
      starting_values = NULL,
      sample_size = 30
    )
  )
})

test_that("checkThat catches missing/incorrect density arguments", {
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, df = 1),
      starting_values = NULL,
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dgamma_quo,
      f_params = list(rate = 1),
      starting_values = NULL,
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dgamma_quo,
      f_params = list(scale = 1, lambda = 1),
      starting_values = NULL,
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dgamma_quo,
      f_params = NULL,
      starting_values = NULL,
      sample_size = 30
    )
  )
})

test_that("checkThat catches invalid starting values", {
  expect_error(
    checkThat(
      f = dexp_quo,
      f_params = NULL,
      starting_values = c(-5, -1, 0, 1, 5),
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dexp_quo,
      f_params = list(rate = 1),
      starting_values = c(5:10),
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 5),
      starting_values = c(0:4),
      sample_size = 30
    )
  )
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 5),
      starting_values = c(1),
      sample_size = 30
    )
  )
})

test_that("checkThat catches invalid sample_size", {
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, sd = 1),
      starting_values = c(-3:3),
      sample_size = 10.5
    )
  )
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, sd = 1),
      starting_values = c(-3:3),
      sample_size = "1"
    )
  )
  expect_error(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, sd = 1),
      starting_values = c(-3:3),
      sample_size = c(10, 20)
    )
  )
})



