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
#' @param x A \code{numeric} vector of length \code{k} for \code{k > 1}.
#'
#' @param f A function representing the target sampling distribution.
#'
#' @param f_params A \code{list} of any associated parameters of \code{f}.
#'
#' @return A named \code{list} containing the intersection points, the h values
#'  of the original x's, and the derivatives to be used in further functions.
#'
#' @import assertthat
#' @importFrom rlang exec
#'
#' @keywords internal
tanIntersect <- function(x, f, f_params = NULL) {
  assertthat::assert_that(
    length(x) > 1,
    msg = 'Must have more than 1 abscissae.'
  )
  j_plus_1 <- length(x)
  j <- j_plus_1 - 1

  xl <- list(x)

  if (!is.null(f_params)){
    xl <- append(
      xl,
      values = f_params)}

  gx <- rlang::exec(
    f, !!!xl)

  hx <- log(gx)

  dhdx <- (1/gx) * approxD(f = f, x = x)

  z <- (( hx[2:j_plus_1] - hx[1:j] ) - ( x[2:j_plus_1] * dhdx[2:j_plus_1] ) + (
    x[1:j] * dhdx[1:j] )) / ( dhdx[1:j] - dhdx[2:j_plus_1] )

  out <- list(
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
#' @return A \code{numeric} vector containing the y-values of the upper hull.
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
	idx <- findInterval(x, c(supp[1],tans$z,supp[2])) # indices for each piece
	u_out <- tans$hx[idx] + (x - x_abs[idx])*tans$dhdx[idx]
	
	return(u_out)
}