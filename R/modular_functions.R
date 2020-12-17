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
#' @return if all checks pass, a named \code{list} is returned containing an
#'  expression for the density function, and the initial abscissae values.
#'
#' @import assertthat
#' @importFrom rlang exec quo_get_expr as_label
#' @importFrom stringr str_replace
#' @importFrom stats runif
#'
#' @keywords internal
checkThat <- function(f, f_params, starting_values, sample_size){

  # Check that f_params is a named list
  if (!is.null(f_params)){
    assert_that(
      !is.null(names(f_params)),
      msg = 'Density parameters must be in a named list.'
    )
  }

  # Check Sample Size Argument
  assert_that(
    all(is.numeric(sample_size),
        length(sample_size) == 1,
        sample_size > 0,
        sample_size %% 1 == 0),
    msg = 'Sample size must be a single positive integer value.'
  )

  # Check for necessary default arguments to sampling density
  d <- rlang::as_label(f)
  d_args <- formals(d)
  defaults <- sapply(d_args, is.symbol)
  default_names <- names(defaults[defaults == TRUE])
  if (sum(defaults) > 1){
    default_names <- default_names[-1]
    if (!is.null(f_params)){
      assert_that(
        all(default_names %in% names(f_params)),
        msg = 'Some required arguments of denisty function are missing'
      )
    }
    else{
      stop(
        paste0(
          "Must provide non-default density arguments: ",
          paste(default_names, collapse = " ")
        )
      )
    }
  }
  else{
    assert_that(
      default_names[1] == 'x',
      msg = "Some required arguments of denisty function are missing"
    )
  }


  # Check function call and function parameters
  fxp <- rlang::quo_get_expr(f)

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
    # # Evaluate density at starting values
    gx <- rlang::exec(
      fxp, !!!f_args
    )
    # Calculate derivatives of the density
    dgdx <- approxD(
      f = fxp,
      f_params = f_params,
      x = starting_values
    )

    # Deriv of log
    dhdx <- dgdx/gx
    assert_that(
      any(is.nan(dhdx)) == FALSE,
      msg = "Some starting values are not in the range of the density."
    )

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
      f = fxp,
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

  out <- list(
    f_expr = fxp,
    initial_abs = initial_abs
    )

  return(out)
}


# Derivative Approximation #################################################
#' @name approxD
#'
#' @title Approximate Function Derivatives
#'
#' @description \code{approxD} approximates the derivative of a specified
#' function at a specified value using a central finite difference approach.
#'
#' @param f The function to evaluate, as an expression.
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

  dx <- abs(x) * h

  if (any(dx == 0)){
    dx[which(dx == 0)] <- h
  }
  fplus_args <- list(x + dx)
  fplus_args <- append(
    fplus_args,
    values = f_params
  )

  fminus_args <- list(x - dx)
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
#' @param f A function representing the sampling distribution, as an expression.
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

  x_abs <- sort(x_abs)

  # Set up index
  j_plus_1 <- length(x_abs)
  j <- j_plus_1 - 1

  # Collect density function arguments
  xl <- list(x_abs)

  if (!is.null(f_params)){
    xl <- append(
      xl,
      values = f_params)
  }
  # Evaluate Density
  gx <- rlang::exec(
    f, !!!xl)

  # evaluate Derivatives
  hx <- log(gx)

  dhdx <- (1/gx) * approxD(f = f, f_params = f_params, x = x_abs)

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
#' @param z The tangent intersection points at given abscissae.
#'
#' @param hx The log of the density \code{f} at point x
#'
#' @param dhdx The derivative of the log of the density \code{f} at point x
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
                      supp = c(-Inf,Inf),
                      z = NULL,
                      hx = NULL,
                      dhdx = NULL) {
  # Find tangent intersections if not provided
  if(is.null(z) | is.null(hx) | is.null(dhdx)){
    tans <- tanIntersect(x_abs, f, f_params)
    z <- tans$z
    hx <- tans$hx
    dhdx <- tans$dhdx
  }

	# Compute the upper hull
	z <- z[z >= supp[1] & z <= supp[2]]
	idx <- findInterval(x, c(supp[1], z, supp[2])) # indices for each piece
	u_out <- hx[idx] + (x - x_abs[idx]) * dhdx[idx]

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
#' @param z The tangent intersection points at given abscissae.
#'
#' @param hx The log of the density \code{f} at point x
#'
#' @param dhdx The derivative of the log of the density \code{f} at point x
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
                      supp = c(-Inf, Inf),
                      z = NULL,
                      hx = NULL,
                      dhdx = NULL) {
  # Store the values for the tangent intersections if not provided
  if(is.null(z) | is.null(hx) | is.null(dhdx)){
    tans <- tanIntersect(x_abs, f, f_params)
    z <- tans$z
    hx <- tans$hx
    dhdx <- tans$dhdx
  }

  # Only select z on the valid interval
  z <- z[z >= supp[1] & z <= supp[2]]

  # Find the normalization constant
  norm <- integrate(function(a) exp(upperHull(a, x_abs, f, f_params, supp)),
                    lower = supp[1],
                    upper = supp[2])$value

  # Sample from this CDF using the inverse-CDF method
  # Done by taking CDF at each intersection and then splitting into pieces
  cdf <- function(x) integrate(function(a) 1/norm *
                               exp(upperHull(a, x_abs, f, f_params, supp)),
                               lower = supp[1], upper = x)$value
  z_cdf <- sapply(z, cdf)
  zb <- c(0,z_cdf,1)
  lb <- c(supp[1], z)
  unif_samp <- runif(n)
  idx <- findInterval(unif_samp, zb)

  samp <- (unif_samp - zb[idx]) * dhdx[idx] * norm/exp(hx[idx]) +
           exp((lb[idx] - x_abs[idx])*dhdx[idx])
  samp <- x_abs[idx] + log(samp)/dhdx[idx]
  samp <- samp[!is.na(samp)]

  while (length(samp)<n) {
    unif_samp <- runif(n)
    samp_1 <- (unif_samp - zb[idx]) * dhdx[idx] * norm/exp(hx[idx]) +
               exp((lb[idx] - x_abs[idx])*dhdx[idx])
    samp_1 <- x_abs[idx] + log(samp)/dhdx[idx]
    samp_1 <- samp_1[!is.na(samp_1)]
    samp <- c(samp, samp_1)
  }

  samp <- samp[1:n]

  return(samp)
}

# Lower Hull Function ##########################################################
#' @name lowerHull
#'
#' @title Compute the Lower Hull
#'
#' @description \code{lowerHull} calculates the y-values corresponding to the
#' lower hull formed from the tangents to a function.
#'
#' @param x A \code{numeric} vector representing x values.
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
#' @param hx The log of the density \code{f} at point x
#'
#' @return A named \code{list} containing the y-values of the lower hull.
#'
#' @importFrom rlang exec
#'
#' @keywords internal

lowerHull <- function(x,
                      x_abs,
                      f,
                      f_params = NULL,
                      supp = c(-Inf, Inf),
                      hx = NULL){
  # Find hx for x_abs if not provided
  if(is.null(hx)){
    tans <- tanIntersect(x_abs, f, f_params)

    hx <- tans$hx
  }

  # anonymous function used to vectorize calculations of lk for multiple x
  lambda <- function(x){
    # return -Inf if x is outside of the range of x_abs
    if(x < min(x_abs) || x > max(x_abs)){
      lk <- -Inf
    }

    # find lk for x in the interval [x[j], x[j+1]]
    else{
      j <- findInterval(x, x_abs, rightmost.closed = T)

      lk <- ((x_abs[j+1] - x) * hx[j] + (x - x_abs[j]) * hx[j+1]) /
        (x_abs[j+1] - x_abs[j])
    }
    return(lk)
  }
  sapply(x, lambda)
}

# Adaptive Rejection Sampling Function #########################################
#' @name ars
#'
#' @title Perform Adaptive Rejection Sampling
#'
#' @description \code{ars} calculates a given number of \code{n} random samples
#' from a given density \code{f}, initialized at the given abscissae
#' \code{x_abs}.
#'
#' @param n The number of samples to return
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

ars <- function(n, x_abs, f, f_params = NULL, supp = c(-Inf, Inf)){

  # Create quosure
  fquo <- rlang::enquo(f)
  if (is.character(f)){
    fsym <- rlang::sym(f)
    fquo <- rlang::enquo(fsym)
  }

  # Perform argument checks
  checkList <- checkThat(
    f = fquo,
    f_params = f_params,
    starting_values = x_abs,
    sample_size = n
  )

  # Extract initial abscissae and function expression after checks
  x_abs <- checkList$initial_abs
  f <- checkList$f_expr

  # initialize values
  vals <- c()

  k <- i <- m <- 1

  tans <- tanIntersect(x_abs, f, f_params)

  z <- tans$z

  hx <- tans$hx

  dhdx <- tans$dhdx

  while(length(vals) < n){
    # find hx for each iteration. Only needs to recalculate if x_abs changes

    # sample x* from s(x)
    x_star <- sampleEnv(1, x_abs, f, f_params, supp, z, hx, dhdx)

    # sample w from unif(0,1)
    w <- runif(1)

    # find lower and upper hull values for x*
    lx_star <- lowerHull(x_star, x_abs, f, f_params, supp, hx)

    ux_star <- upperHull(x_star, x_abs, f, f_params, supp, z, hx, dhdx)

    # initial test for acceptance
    if(w <= exp(lx_star - ux_star)){
      vals <- append(vals, x_star)
    }
    # calculate h(x*) and h'(x*) if x* isn't at first accepted
    ## we're currently not using the updated h'(x*)
    else{
      xl_star <- append(x_star, f_params)

      gx <- rlang::exec(f, !!!xl_star)

      hx_star <- log(gx)

      dhdx_star <- (1/gx) * approxD(f = f, f_params = f_params, x = x_abs)

      # update x_abs
      x_abs <- sort(append(x_abs, x_star))

      # recalculate z, hx, dhdx with updated x_abs
      tans <- tanIntersect(x_abs, f, f_params)

      z <- tans$z

      hx <- tans$hx

      dhdx <- tans$dhdx

      # find index of x* for second rejection test
      id <- which(x_abs == x_star)

      # iterate k, the count of abscissae added
      k <- k + 1

      # perform second rejection test
      if(w <= exp(hx[id] - ux_star)){
        vals <- append(vals, x_star)
      }

      else{
        # iterate m, the count of rejected samples
        m <- m + 1
      }
    }
    # iterate i, the count of total iterations of the ars function
    i <- i + 1
  }
  return(c(vals, "k" = k,"i" = i, "m" = m))
}
