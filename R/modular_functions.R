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
#' @param f_params A \code{list} of any associated parameters of \code{fun}.
#'
#' @param x The value at which to evaluate the derivative of \code{fun}.
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
    msg = "Incorrect arguments supplid to density function."
  )

  if (x == 0){
    dx <- h
  }
  else{
    dx <- abs(x) * h
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
#' @title Compute Tangent Lines and Intersections
#'
#' @description \code{tanIntersect} calculates the intersections of subsequent
#' tangent lines of a function.
#'
#' @param x A \code{numeric} vector of length \code{k} for \code{k > 1}.
#'
#' @param f A function representing the target sampling distribution.
#'
#' @return A \code{numeric} vector of \code{k - 1} intersection points.
#'
#' @import assertthat
#'
#' @keywords internal
tanIntersect <- function(x, f) {
  assertthat::assert_that(
    length(x) > 1,
    msg = 'Must have more than 1 abscissae.'
  )
  k <- 1:length(x)
}
