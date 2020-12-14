# Modular Functions for Adaptive Rejection Sampling Algorithm

# Derivative Approximation #################################################
#' @name approxD
#'
#' @title Approximate Function Derivatives
#'
#' @description \code{approxD} approximates the derivative of a specified
#' function at a specified value using a central finite difference approach.
#'
#' @param f The function to evaluate
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
#' @importFrom rlang exec
#'
#' @keywords internal
approxD <- function(f,
                   f_params = NULL,
                   x,
                   n = 1,
                   h = sqrt(.Machine$double.eps)
                   ) {
  # Check that proper arguments are supplied
  arg_names <- names(formals(f))
  f_names <- names(f_params)

  assertthat::assert_that(
    all(f_names %in% arg_names) == TRUE,
    msg = "Incorrect arguments supplied to density function."
  )

  dx <- abs(x) * h

  if (any(dx == 0)){
    dx[which(dx == 0)] <- h
  }
  fplus_args <- list(x + dx)
  names(fplus_args) <- arg_names[1]
  fplus_args <- append(
    fplus_args,
    values = f_params
  )

  fminus_args <- list(x - dx)
  names(fminus_args) <- arg_names[1]
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
    names(fx_args) <- arg_names[1]
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
#' @param f A function representing the target sampling distribution.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @return A named \code{list} containing the original abscissae values,
#'  the intersection points of the tangent lines, the h values
#'  of the original abscissae, and the derivatives to be used in further functions.
#'
#' @import assertthat
#' @importFrom rlang exec
#'
#' @keywords internal
tanIntersect <- function(x_abs, f, f_params = NULL) {
  assertthat::assert_that(
    length(x_abs) > 1,
    msg = 'Must have more than 1 abscissae.'
  )
  j_plus_1 <- length(x_abs)
  j <- j_plus_1 - 1

  xl <- list(x_abs)

  if (!is.null(f_params)){
    xl <- append(
      xl,
      values = f_params)}

  gx <- rlang::exec(
    f, !!!xl)

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


lowerHull <- function(x, x_abs, f, f_params = NULL, supp = c(-Inf, Inf)){
  # need to sort x_abs to find the correct x[j] and x[j+1]
  x_abs <- sort(x_abs)

  # find hx for sorted x_abs
  tans <- tanIntersect(x_abs, f, f_params)

  hx <- tans$hx

  # return -Inf if x is outside of the range of x_abs
  if(x < min(x_abs) || x > max(x_abs)){
    lk <- -Inf
  }
  # find lk for x in the interval [x[j], x[j+1]]
  else{
    j <- which(x == sort(append(x_abs, x))) - 1

    lk <- ((x_abs[j+1] - x) * hx[j] + (x - x_abs[j]) * hx[j+1]) /
      (x_abs[j+1] - x_abs[j])
  }
  return(lk)
}

ars <- function(f, f_params = NULL, supp = c(-Inf, Inf), x_abs, n){
  # initialize values
  vals <- list()

  k <- 1

  while(length(vals) < n){
    # find hx for each iteration. Only needs to recalculate if x_abs changes
    x_abs <- sort(x_abs)

    tans <- tanIntersect(x_abs, f, f_params)

    hx <- tans$hx

    # sample x* from s(x)
    x_star <- sampleEnv(1, x_abs, f, f_params, supp)

    # sample w from unif(0,1)
    w <- runif(1)

    # find lower and upper hull values for x*
    lx_star <- lowerHull(x_star, x_abs, f, f_params, supp)

    ux_star <- upperHull(x_star, x_abs, f, f_params, supp)

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

      dhdx_star <- (1/gx) * approxD(f = f, x = x_abs)

      # update x_abs
      x_abs <- sort(append(x_abs, x_star))

      # iterate k
      k <- k + 1

      # perform second rejection test
      if(w <= exp(hx_star - ux_star)){
        vals <- append(vals, x_star)
      }
    }
  }
  return(vals)
}
