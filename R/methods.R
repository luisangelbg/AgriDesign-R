# Standard methods, so that fitted designs work with the usual R functions.

#' Standard methods for fitted designs
#'
#' `residuals()`, `fitted()`, `nobs()` and `anova()` work on the result of
#' [design_anova()] as on a linear model. `anova()` returns the analysis of
#' variance as a data frame of class `"anova"`, with the error stratum
#' against which each term is tested. `as.data.frame()` on the result of
#' [compare_means()] returns its table of means with their letters.
#'
#' @param object,x A fitted design from [design_anova()] (or, for
#'   `as.data.frame()`, a comparison from [compare_means()]).
#' @param type For `residuals()`: `"response"` (observed minus fitted) or
#'   `"studentized"` (internally studentized, as [stats::rstandard()]).
#' @param row.names,optional,... Passed on or not used.
#' @return A numeric vector (residuals, fitted values), a count (`nobs`), or a
#'   data frame.
#' @examples
#' fit <- design_anova(maize_nitrogen, "Yield_t_ha", "Nitrogen_kg_ha", block = "Block")
#' anova(fit)
#' head(residuals(fit, type = "studentized"))
#' nobs(fit)
#' as.data.frame(compare_means(fit, "Nitrogen_kg_ha"))
#' @name ad_methods
#' @importFrom stats residuals fitted nobs anova
NULL

#' @rdname ad_methods
#' @export
residuals.ad_anova <- function(object, type = c("response", "studentized"), ...) {
  type <- match.arg(type)
  e <- object$qr_fit$residuals
  if (type == "studentized") e <- e / (sqrt(object$residual[["ms"]]) * sqrt(pmax(1e-12, 1 - .hat_values(object$qr_fit))))
  stats::setNames(as.numeric(e), rownames(object$data))
}

#' @rdname ad_methods
#' @export
fitted.ad_anova <- function(object, ...) {
  stats::setNames(as.numeric(object$y - object$qr_fit$residuals), rownames(object$data))
}

#' @rdname ad_methods
#' @export
nobs.ad_anova <- function(object, ...) object$n

#' @rdname ad_methods
#' @export
anova.ad_anova <- function(object, ...) {
  t <- object$table
  out <- data.frame(Df = c(t$df, object$residual[["df"]]), `Sum Sq` = c(t$ss, object$residual[["ss"]]),
                    `Mean Sq` = c(t$ms, object$residual[["ms"]]), `F value` = c(t$F, NA), `Pr(>F)` = c(t$p, NA),
                    `Tested against` = c(ifelse(t$is_error, "", t$error), ""),
                    check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- c(t$source, "Residuals")
  structure(out, heading = paste0("Analysis of variance: ", object$design_name, "\nResponse: ", object$response, "\n"),
            class = c("anova", "data.frame"))
}

#' @rdname ad_methods
#' @export
as.data.frame.ad_comparison <- function(x, row.names = NULL, optional = FALSE, ...) {
  m <- x$means
  rownames(m) <- row.names
  m
}
