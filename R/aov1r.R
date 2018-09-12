
#' @name aov1r
#' @rdname aov1r
#' @title One-way random effect ANOVA
#' @description Fits a one-way random effect ANOVA model.
#'
#' @param formula a formula of the form \code{y~group}
#' @param data optional dataframe
#' @param x output of \code{summary}
#' @param ... ignored
#'
#' @return \code{aov1r} returns an object of class \code{aov1r};
#' @import data.table
#' @importFrom lazyeval f_eval_lhs f_eval_rhs f_lhs f_rhs
#'
#' @examples
#' dat <- simAV1R(I=2, J=3, mu=10, sigmab=1, sigmaw=1)
#' fit <- aov1r(y ~ group, data=dat)
#' summary(fit)
NULL

#' @rdname aov1r
#' @export
aov1r <- function(formula, data=NULL){
  DT <- data.table(y = lazyeval::f_eval_lhs(formula, data=data),
                   group = lazyeval::f_eval_rhs(formula, data=data))
  DT[, means := mean(y), by="group"]
  ssw <- with(DT, crossprod(y-means)[1L,1L])
  DT2 <- DT[, list(means = means[1L], Ji = .N), by="group"]
  DT2[, Mean:=mean(means)]
  I <- nrow(DT2)
  Jh <- I/sum(1/DT2[["Ji"]])
  ssb <- Jh*with(DT2, crossprod(Mean-means)[1L,1L])
  terms <- c(y = as.character(lazyeval::f_lhs(formula)),
             group = as.character(lazyeval::f_rhs(formula)))
  N <- nrow(data)
  out <- list(
    "Sums of squares" = c(ssw=ssw, ssb=ssb),
    "Variance components" = c(sigma2w = ssw/(N-I), sigma2b = (ssb/(I-1)-ssw/(N-I))/Jh),
    "Design" = c(I=I, Jh=Jh, N=N),
    "grandmean" = DT2$Mean[1L],
    "groupmeans" =
      setNames(
        as.data.frame(DT2[, .SD, .SDcols=c("group", "means")]),
        c("group", "mean")),
    "data" = data[, terms],
    "terms" = terms
  )
  class(out) <- "aov1r"
  return(out)
}

#' @rdname aov1r
#' @export
summary.aov1r <- function(object, ...){
  out <- list()
  class(out) <- "summary.aov1r"
  out[["Response"]] <- object$terms[["y"]]
  out[["Factor"]] <- object$terms[["group"]]
  return(out)
}

#' @rdname aov1r
#' @export
print.summary.aov1r <- function(x, ...){
  for(foo in names(x)){
    cat(foo, ": ", x[[foo]], "\n", sep="")
  }
}

#' @title Prediction interval
#' @description Prediction interval for the one-way random effect ANOVA.
#'
#' @param object an output of \code{\link{aov1r}}
#' @param level confidence level
#' @param ... ignored
#'
#' @return A vector of length two, the bounds of the prediction interval.
#' @export
#' @importFrom stats qt
#'
#' @examples
#' dat <- simAV1R(I=2, J=3, mu=10, sigmab=1, sigmaw=1)
#' fit <- aov1r(y ~ group, data=dat)
#' predict(fit)
predict.aov1r <- function(object, level=0.95, ...){
  I <- object[["Design"]][["I"]]
  J <- object[["Design"]][["Jh"]]
  N <- object[["Design"]][["N"]]
  SSb <- object[["Sums of squares"]][["ssb"]]
  SSw <- object[["Sums of squares"]][["ssw"]]
  a <- (1/J*(1+1/I))/(I-1)
  b <- (J-1)/J/(N-I)
  v <- a*SSb+b*SSw # estimates the variance of (Ynew-Ybar)
  nu <- v^2/((a*SSb)^2/(I-1)+(b*SSw)^2/(N-I)) # Satterthwaite degrees of freedom
  alpha.over.two <- (1-level)/2
  bounds <-  object[["grandmean"]] +
    c(-1,1)*sqrt(v)*qt(1-alpha.over.two, nu)
  names(bounds) <- paste0(100*c(alpha.over.two, 1-alpha.over.two), "%")
  attr(bounds, "std.error") <- sqrt(v)
  attr(bounds, "df") <- nu
  return(bounds)
}
