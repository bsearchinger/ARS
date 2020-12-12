# Modular Functions for Adaptive Rejection Sampling Algorithm

# Check Conditions Function ####################################################
#' @name checkThat
#'
#' @title Check arguments of the call to ars.
#'
#' @description \code{checkThat} performs tests and checks to ensure
#'  functionality of the adaptive rejection sampling algorithm.  For example, if
#'  no starting values are supplied, \code{checkThat} will determine viable
#'  starting values for the initialization step.
#'
#' @param f A function quosure representing the target sampling distribution.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @param starting_values The initial values for the k abscissae.  If NULL,
#'  values will be determined using the bounds of the log-distribution where the
#'  first derivative is positive.
#'
#' @param sample_size The desired number of samples in the final output.
#'
#' @return A warning or error if certain conditions are detected
#'
#' @import assertthat
#' @importFrom rlang exec quo_get_expr
#' @importFrom stringr str_replace
#' @importFrom stats runif
#'
#' @keywords internal
checkThat <- function(f, f_params, starting_values, sample_size){

  assert_that(
    all(is.numeric(sample_size),
        length(sample_size) == 1,
        sample_size > 0,
        sample_size %% 1 == 0),
    msg = 'Sample size must be a single positive integer value.'
  )

  # Check function call and function parameters
  fxp <- rlang::quo_get_expr(f)

  # TODO: put this at beginning of final wrapper function
  #assert_that(
  #  is.function(f),
  #  msg = "f must be a density function."
  #)
  #if (!is.null(f_params)){
  #  assert_that(
  #    names(f_params) == names(formals(f)),
  #    msg = "Names in f_params do not match arguments for f."
  #  )
  #}

  # Check that starting values are valid
  if (!is.null(starting_values)){
    assert_that(
      length(starting_values) >= 2,
      msg = "Must have at least 2 initial abscissae."
    )
    # Check that abscissae have both positive and negative derivative values
    f_args <- list(starting_values)

    if(!is.null(f_params)){
      f_args <- append(
        f_args,
        values = f_params
      )
    }
    # # Evaluate density at the sample at starting values
    gx <- rlang::exec(
      fxp, !!!f_args
    )
    # Calculate derivatives of the density
    dgdx <- approxD(
      f = f,
      f_params = f_params,
      x = starting_values
    )

    # Deriv of log
    dhdx <- dgdx/gx

    positive_vals <- starting_values[dhdx > 0]
    negative_vals <- starting_values[dhdx < 0]

    assert_that(
      length(positive_vals) > 0,
      msg = paste(
      'No starting values have positive log-derivative.
      Consider expanding lower bound')
      )
    assert_that(
      length(negative_vals) > 0,
      msg = paste(
      'No starting values have negative log-derivative.
      Consider expanding upper bound')
    )

    initial_abs <- starting_values
  }
  # Generate starting values if none are provided
  else{
    # Change density function to sample function
    fun_expr <- rlang::as_label(f)
    rfun <- stringr::str_replace(
      fun_expr, 'd', 'r'
      )
    rfun_args = list(100)
    if(!is.null(f_params)){
      rfun_args <- append(
        rfun_args,
        values = f_params
      )
    }
    # Sample 100 values from provided density
    rsample <- rlang::exec(
      rfun, !!!rfun_args
      )

    # Evaluate density at the sample
    rfun_args[[1]] <- rsample
    gx <- rlang::exec(
      fxp, !!!rfun_args
    )

    # Calculate derivatives from the sample
    dgdx <- approxD(
      f = f,
      f_params = f_params,
      x = rsample
      )

    # Deriv of log
    dhdx <- dgdx/gx

    # Subset between positive and negative derivatives
    positive_vals <- rsample[dhdx > 0]
    negative_vals <- rsample[dhdx < 0]
    pos_len <- length(positive_vals)
    neg_len <- length(negative_vals)

    # TODO:clean this up a bit. Check for zero length
    stopifnot(
      exprs = {
        pos_len > 0
        neg_len > 0
        }
      )

    # Use the mean value on either side to determine the boundary
    lower_bound <- mean(positive_vals)
    upper_bound <- mean(negative_vals)

    # Generate 10 samples from a uniform in this interval
    initial_abs <- runif(
      10, lower_bound, upper_bound
    )
  }
  return(initial_abs)
}


# Derivative Approximation #################################################
#' @name approxD
#'
#' @title Approximate Function Derivatives
#'
#' @description \code{approxD} approximates the derivative of a specified
#' function at a specified value using a central finite difference approach.
#'
#' @param f The function to evaluate, as a quosure.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @param x The value at which to evaluate the derivative of \code{f}.
#'
#' @param n The order of derivative to calculate, 1 or 2.  Default is 1.
#'
#' @param h The finite difference to use in the approximation.  The default
#' value is \code{sqrt(.Machine$double.eps)} by which the absolute value of
#' \code{x} will be multiplied unless \code{x} is 0.
#'
#' @return The approximated derivative of the function at the specified value.
#'
#' @importFrom rlang exec quo_get_expr
#'
#' @keywords internal
approxD <- function(f,
                   f_params = NULL,
                   x,
                   n = 1,
                   h = sqrt(.Machine$double.eps)
                   ) {
  # Check that proper arguments are supplied
  #arg_names <- names(formals(f))
  #f_names <- names(f_params)

  #assertthat::assert_that(
  #  all(f_names %in% arg_names) == TRUE,
  #  msg = "Incorrect arguments supplied to density function."
  #)
  f <- rlang::quo_get_expr(f)

  dx <- abs(x) * h

  if (any(dx == 0)){
    dx[which(dx == 0)] <- h
  }
  fplus_args <- list(x + dx)
  #names(fplus_args) <- arg_names[1]
  fplus_args <- append(
    fplus_args,
    values = f_params
  )

  fminus_args <- list(x - dx)
  #names(fminus_args) <- arg_names[1]
  fminus_args <- append(
    fminus_args,
    values = f_params
  )

  fxp <- rlang::exec(f, !!!fplus_args)
  fxm <- rlang::exec(f, !!!fminus_args)


  if (n == 1) {
    fxdx <- (fxp - fxm)/(2*dx)

    return(fxdx)
  }
  else{
    fx_args <- list(x)
    #names(fx_args) <- arg_names[1]
    fx_args <- append(
      fx_args,
      values = f_params
    )

    fx <- rlang::exec(f, !!!fx_args)

    dfdx2 <- (fxp - 2*fx + fxm)/(dx^2)

    return(dfdx2)

  }
}

# Tangent Intersection Function ################################################
#' @name tanIntersect
#'
#' @title Compute x-values of Tangent Intersections
#'
#' @description \code{tanIntersect} calculates the x-values corresponding to the
#' intersections of subsequent tangent lines of a function.
#'
#' @param x_abs A \code{numeric} vector of length \code{k} for \code{k > 1}.
#'
#' @param f A function representing the sampling distribution, as a quosure.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @return A named \code{list} containing the original abscissae values,
#'  the intersection points of the tangent lines, the h values
#'  of the original abscissae, and the derivatives to be used in further functions.
#'
#' @import assertthat
#' @importFrom rlang exec quo_get_expr
#'
#' @keywords internal
tanIntersect <- function(x_abs, f, f_params = NULL) {
  assertthat::assert_that(
    length(x_abs) > 1,
    msg = 'Must have more than 1 abscissae.'
  )

  fxp <- rlang::quo_get_expr(f)

  j_plus_1 <- length(x_abs)
  j <- j_plus_1 - 1

  xl <- list(x_abs)

  if (!is.null(f_params)){
    xl <- append(
      xl,
      values = f_params)
    }

  gx <- rlang::exec(
    fxp, !!!xl)

  hx <- log(gx)

  dhdx <- (1/gx) * approxD(f = f, x = x_abs)

  # Gilks Equation (1)
  z <- (( hx[2:j_plus_1] - hx[1:j] ) - ( x_abs[2:j_plus_1] * dhdx[2:j_plus_1] ) + (
    x_abs[1:j] * dhdx[1:j] )) / ( dhdx[1:j] - dhdx[2:j_plus_1] )

  out <- list(
    x_abs = x_abs,
    z = z,
    hx = hx,
    dhdx = dhdx
  )

  return(out)
}

# Upper Hull Function #######################################################
#' @name upperHull
#'
#' @title Compute the Upper Hull
#'
#' @description \code{upperHull} calculates the y-values corresponding to the
#' upper hull formed from the tangents to a function
#'
#' @param x A \code{numeric} vector representing x values
#'
#' @param x_abs A \code{numeric} vector of length \code{k} for \code{k > 1},
#' representing the abscissae.
#'
#' @param f A function representing the target sampling distribution.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @param supp The support of \code{f}, as a two-membered \code{numeric} vector.
#' Default is \code{c(-Inf, Inf)}
#'
#' @return A named \code{list} containing the y-values of the upper hull.
#'
#' @import assertthat
#' @importFrom rlang exec
#'
#' @keywords internal
upperHull <- function(x,
                      x_abs,
                      f,
                      f_params = NULL,
                      supp = c(-Inf,Inf)
                      ) {
  # Find tangent intersections
	tans <- tanIntersect(x_abs, f, f_params)

	# Compute the upper hull
	z <- tans$z
	z <- z[z >= supp[1] & z <= supp[2]]
	idx <- findInterval(x, c(supp[1], z, supp[2])) # indices for each piece
	u_out <- tans$hx[idx] + (x - x_abs[idx])*tans$dhdx[idx]

	return(u_out)
}


# Sample From Rejection Envelope ############################################
#' @name sampleEnv
#'
#' @title Sample from the Rejection Envelope
#'
#' @description \code{sampleEnv} samples from the rejection envelope formed
#' by the piecewise upper hull for the given function
#'
#' @param n The number of samples
#'
#' @param x_abs A \code{numeric} vector of length \code{k} for \code{k > 1},
#' representing the abscissae.
#'
#' @param f A function representing the target sampling distribution.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @param supp The support of \code{f}, as a two-membered \code{numeric} vector.
#' Default is \code{c(-Inf, Inf)}
#'
#' @return A \code{numeric} vector containing the sampled values.
#'
#' @import assertthat
#' @importFrom rlang exec
#' @importFrom stats integrate runif
#'
#' @keywords internal
sampleEnv <- function(n,
                      x_abs,
                      f,
                      f_params = NULL,
                      supp = c(-Inf, Inf)) {
  # Store the values for the tangent intersections
  tans <- tanIntersect(x_abs, f, f_params)
  z <- tans$z
  z <- z[z >= supp[1] & z <= supp[2]]

  # Find the normalization constant
  norm <- integrate(function(a) exp(upperHull(a, x_abs, f, f_params, supp)),
                    lower = supp[1],
                    upper = supp[2])$value

  # Sample from this CDF using the inverse-CDF method
  # Done by taking CDF at each intersection and then splitting into pieces
  cdf <- function(x) integrate(function(a) 1/norm * exp(upperHull(a,x_abs,f,f_params, supp)),
                               lower = supp[1],
                               upper = x)$value
  z_cdf <- sapply(z, cdf)
  zb <- c(0,z_cdf,1)
  lb <- c(supp[1], z)
  unif_samp <- runif(n)
  idx <- findInterval(unif_samp, zb)

  samp <- (unif_samp - zb[idx]) * tans$dhdx[idx] * norm/exp(tans$hx[idx]) + exp((lb[idx] - x_abs[idx])*tans$dhdx[idx])
  samp <- x_abs[idx] + log(samp)/tans$dhdx[idx]
  samp <- samp[!is.na(samp)]

  while (length(samp)<n) {
    unif_samp <- runif(n)
    samp_1 <- (unif_samp - zb[idx]) * tans$dhdx[idx] * norm/exp(tans$hx[idx]) + exp((lb[idx] - x_abs[idx])*tans$dhdx[idx])
    samp_1 <- x_abs[idx] + log(samp)/tans$dhdx[idx]
    samp_1 <- samp_1[!is.na(samp_1)]
    samp <- c(samp, samp_1)
  }

  samp <- samp[1:n]

  return(samp)
}
