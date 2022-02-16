#' @export
RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}
