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
# Tests for tanIntersect: Computing Intersections of Tangent Lines
ztest <- c(-0.5, -0.1, 0, 0.1, 0.5)
xtest <- seq(0, 1, length.out = 20)
test_that("output of tanIntersect is symmetric for symmetric distribution",{
  expect_true(
    all(
      tanIntersect(ztest, dnorm)$z[1:2] == -1*rev(tanIntersect(ztest, dnorm)$z[3:4])
    ) == TRUE
  )
})

test_that("output of tanIntersect is correct dimensions", {
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
# Tests for checkThat: Argument Checking Function

# Note: checkThat takes a quosure as an argument
dexp_quo <- rlang::quo(dexp)
dnorm_quo <- rlang::quo(dnorm)
dgamma_quo <- rlang::quo(dgamma)

test_that("checkThat runs silently for correct inputs", {
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params = list(mean = 0, sd = 1),
      starting_values = c(-0.5, -0.1, 0, 0.1, 0.5),
      sample_size = 30,
      supp = c(-Inf, Inf)
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params = NULL,
      starting_values = c(-0.5, -0.1, 0, 0.1, 0.5),
      sample_size = 30,
      supp = c(-Inf, Inf)
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params =  list(mean = 0, sd = 1),
      starting_values = NULL,
      sample_size = 30,
      supp = c(-Inf, Inf)
    )
  )
  expect_silent(
    checkThat(
      f = dnorm_quo,
      f_params =  NULL,
      starting_values = NULL,
      sample_size = 30,
      supp = c(-Inf, Inf)
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

################################################################################
# Tests for cavitySearch: Log-concavity Check Function
test_that("cavitySearch detects non-log-concavity", {
  expect_false(
    cavitySearch(
      f = dt, f_params = list(df = 5, ncp = 0), x = -4)
  )
  expect_false(
    all(
      cavitySearch(
        f = dt, f_params = list(df = 5, ncp = 0), x = c(-4, -1, 0, 1))) == TRUE
  )
})


################################################################################
# Tests for upperHull: Upper Hull Function

ztest <- c(-0.5, -0.1, 0, 0.1, 0.5)
test_that("output of upperHull is symmetric for symmetric distribution",{
  expect_true(
    all(
      upperHull(-5:5, ztest, dnorm) == rev(upperHull(-5:5, ztest, dnorm))
    ) == TRUE
  )
})

test_that("upperHull bounds the log-distribution from above", {
  expect_true(
    all(
      upperHull(-10:10/10, ztest, dnorm) >= dnorm(-10:10/10, log = T)
    )
  )
  expect_true(
    all(
      upperHull(0:50/10,
                1:10/3,
                dgamma,
                f_params = list(shape = 4, rate = 2)) >= dgamma(0:50/10, 4, 2, log = T)
    )
  )
})


################################################################################
# Tests for sampleEnv: Sampling from Piecewise Envelope

test_that("sampleEnv returns the correct sample size", {
  expect_length(
    sampleEnv(10, ztest, dnorm),
    10
  )
  expect_length(
    sampleEnv(10, 1:10/3, dgamma, f_params = list(shape = 4, rate = 2)),
    10
  )
})

test_that("sampleEnv works with several known distributions", {
  expect_silent(
    sampleEnv(10, ztest, dnorm)
  )
  expect_silent(
    sampleEnv(10, 1:10/3, dchisq, f_params = list(df = 4))
  )
  expect_silent(
    sampleEnv(10, 1:10/3, dgamma, f_params = list(shape = 4, rate = 2))
  )
})

test_that("sampleEnv returns finite values", {
  expect_true(
    all(
      is.finite(sampleEnv(10, ztest, dnorm))
    )
  )
})

################################################################################
# Tests for lowerHull: Lower Hull Function

xtest <- c(-0.25, -0.1, 0.1, 0.25)
test_that("output of lowerHull is -Inf for out-of-bounds x", {
  expect_true(
    lowerHull(-5, ztest, dnorm) == -Inf
  )
})

test_that("lowerHull is less than or equal to upperHull", {
  expect_true(
    all(
      signif(lowerHull(xtest, ztest, dnorm), digits = 7) <=
        signif(upperHull(xtest, ztest, dnorm), digits = 7)
    )
  )
  expect_true(
    all(
      signif(lowerHull(0:50/10,
                       1:10/3,
                       dgamma,
                       f_params = list(shape = 4, rate = 2)), digits = 7) <=
        signif(upperHull(0:50/10,
                         1:10/3,
                         dgamma,
                         f_params = list(shape = 4, rate = 2)), digits = 7)
    )
  )
})

test_that("lowerHull bounds the log-distribution from below", {
  expect_true(
    all(
      signif(lowerHull(-10:10/10, ztest, dnorm), digits = 7) <=
        signif(dnorm(-10:10/10, log = T), digits = 7)
    )
  )
  expect_true(
    all(
      signif(lowerHull(0:50/10,
                       1:10/3,
                       dgamma,
                       f_params = list(shape = 4, rate = 2)), digits = 7) <=
        signif(dgamma(0:50/10, 4, 2, log = T), digits = 7)
    )
  )
})


################################################################################
# General Tests for ars: Main Wrapper Function

expected_norm_0_1 <- c(-0.8002421, 0.4852841, 0.7472756, -0.2869187, 1.5645132)
expected_gamma_2_3 <- c(0.2668276, 0.7871653, 0.9318842, 0.8362395, 0.4302574)
expected_gamma_5_3 <- c(1.067018, 1.954670, 2.175522, 2.030368, 1.370136)
expected_trunc_norm <- c(-0.32586650, -0.24809180, -0.49334810, -0.06272284, -1.27120300)

# Tests for reproducible results
set.seed(1)
norm_abs_0_1 <- runif(10, min = -3, max = 3)
test_that("ars returns expected results for normal(0,1)", {
  expect_equal(
    round(ars(
      n = 5,
      x_abs = norm_abs_0_1,
      f = dnorm,
      f_params = list(
        mean = 0,
        sd = 1),
      supp = c(-Inf, Inf)
    )$vals, digits = 7),
    expected_norm_0_1
  )
})

set.seed(1)
gamma_abs_2_3 <- runif(10, min = 0, max = 1.5)
test_that("ars returns expected results for gamma(2,3)", {
  expect_equal(
    round(ars(
      n = 5,
      x_abs = gamma_abs_2_3,
      f = dgamma,
      f_params = list(
        shape = 2,
        rate = 3),
      supp = c(0, Inf)
    )$vals, digits = 7),
    expected_gamma_2_3
  )
})

set.seed(1)
gamma_abs_5_3 <- runif(10, min = 0.5, max = 3)
test_that("ars returns expected results for gamma(5,3)", {
  expect_equal(
    signif(ars(
      n = 5,
      x_abs = gamma_abs_5_3,
      f = dgamma,
      f_params = list(
        shape = 5,
        rate = 3),
      supp = c(0.5, Inf)
    )$vals, digits = 7),
    expected_gamma_5_3
  )
})

set.seed(1)
norm_abs_trunc <- rnorm(100)
norm_abs_trunc <- norm_abs_trunc[norm_abs_trunc < 0]
norm_abs_trunc <- sample(norm_abs_trunc, 10)
norm_abs_trunc <- c(norm_abs_trunc, 0)
test_that("ars returns expected results for truncated normal(0,1)", {
  expect_equal(
    signif(ars(
      n = 5,
      x_abs = norm_abs_trunc,
      f = dnorm,
      f_params = list(
        mean = 0,
        sd = 1),
      supp = c(-Inf, 0)
    )$vals, digits = 7),
    expected_trunc_norm
  )
})

# Test for detecting non-log-concavity
test_that("Non-log-concave distirbutions are caught", {
  expect_error(
    ars(
      n = 5,
      x_abs = c(-10, -9, -5, 0, 5, 9, 10),
      f = dt,
      f_params = list(
        df = 10,
        ncp = 0),
      supp = c(-Inf, Inf)
    )
  )
  expect_error(
    ars(
      n = 5,
      x_abs = c(0, 0.5, 0.75, 1, 1.25),
      f = dgamma,
      f_params = list(
        shape = 0.75,
        rate = 2),
      supp = c(0, Inf)
    )
  )
  expect_error(
    ars(
      n = 5,
      x_abs = c(0, 0.5, 0.75, 1, 1.25),
      f = dchisq,
      f_params = list(
        df = 1),
      supp = c(0, Inf)
    )
  )
})

# Test to catch mis-matched initial x_abs and support arguments
test_that("Mismatching initial values and support are detected", {
  expect_error(
    ars(
      n = 5,
      x_abs = c(-5, -4, -3, 0, 1, 2),
      f = dnorm,
      f_params = list(
        mean = 0,
        sd = 1),
      supp = c(-3.5, 1.5)
    )
  )
})

set.seed(1)
# Test that fucntion behaves normally
test_that("ars() runs silently", {
  # supplying all fields validly
  expect_silent(
    ars(
      n = 25,
      x_abs = c(-3, -0.95, -0.2, 0.1, 0.25, 1.3, 1.5),
      f = dnorm,
      f_params = list(
        mean = 0,
        sd = 1),
      supp = c(-Inf, Inf)
    )
  )
  # Removing f_params and support
  expect_silent(
    ars(
      n = 25,
      x_abs = c(-3, -0.95, -0.2, 0.1, 0.25, 1.3, 1.5),
      f = dnorm
    )
  )
  # Removing starting values
  expect_silent(
    ars(
      n = 25,
      f = dnorm
    )
  )
  # generating values based on constrained support
  expect_silent(
    ars(
      n = 25,
      f = dnorm,
      supp = c(-2, 2)
    )
  )
})

