#' @title expit
#' @description calculate the expit 
#' @param x a number from -inf to inf
#' @importFrom a numeric number
#' @return a positive numeric number
#' @export
#'
#' @examples \dontrun{
#' expit(2)
#' }

expit = function(x){exp(x)/(1+exp(x))}

