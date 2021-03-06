
#' @keywords internal
#' @importFrom stats terms xtabs


.multiCA.test <- function(x, scores, outcomes){
  K <- nrow(x)
  full <- length(outcomes) == K  #full test
  
  nidot <- apply(x, 2, sum)
  n <- sum(nidot)
  
  cbar <- sum(nidot * scores)/n
  
  s2 <- sum(nidot * (scores - cbar)^2)
  pdot <- prop.table(rowSums(x))[outcomes]
  nonz <- (pdot > 0)
  
  if (!any(nonz)) return(1)
  
  X <- x[outcomes, ,drop=FALSE] %*% (scores - cbar)

  #individual tests
  CAT <- X[nonz]^2 / (pdot[nonz] * (1-pdot[nonz])) / s2 
  CAT.p.value <- pchisq(CAT, df=1, lower.tail=FALSE)
  
  #overall test
  if (full || sum(pdot) >= 1){
    Tt <- ( sum(X[nonz]^2 / pdot[nonz])) / s2
  } else {
    Tt <- (sum(X)^2 / (1-sum(pdot)) + sum(X[nonz]^2 / pdot[nonz])) / s2
  }

  df <- length(outcomes) - full
  p.value <- pchisq(Tt, df=df, lower.tail=FALSE)

  res <- list(statistic = Tt, parameter = df, p.value = p.value, 
              indiv.statistics = CAT, indiv.p.value = CAT.p.value)
  return(res)
}

#'@rdname multiCA.test
#'@method multiCA.test default
#'@param scores non-decreaseing numeric vector of the same length as the number of ordered groups. Defaults to linearly increasing values
#'@param outcomes integer or character vector defining the set of outcomes (by row index or row name) over which the trend should be tested. Defaults to all outcomes.
#'@param p.adjust.method character string defining the correction method for individual outcome p-values. Defaults to "closed.set" when \code{length(outcomes)<=3}, and "Holm-Shaffer" otherwise.
#'@export
#' @importFrom utils str

multiCA.test.default <- function(x, scores=1:ncol(x), outcomes=1:nrow(x),
  p.adjust.method=c("none","closed.set","Holm-Shaffer"),...){
  if (!is.matrix(x)) {
    cat(str(x))
    stop("x should be a two-dimensional matrix")
}
  if (length(scores) != ncol(x)) stop("The length of the score vector should equal the number of columns of x")

  testres <- .multiCA.test(x=x, scores=scores, outcomes=outcomes)
 
  Tt <- c(W = testres$statistic)
  df <- c(df = testres$parameter)

  p.value <- testres$p.value
  null.value <- 0
  names(null.value) <- sprintf("slope for outcomes %s", deparse(substitute(outcomes)))

  res <- list(statistic = Tt, parameter = df, p.value = p.value, 
              method="Multinomial Cochran-Armitage trend test",
              alternative="two.sided",
              null.value=null.value,
              data.name = deparse(substitute(x)))
  class(res) <- "htest"

  
    if (missing(p.adjust.method)){
      if (length(outcomes)<=3) p.adjust.method <- "closed.set"
      else p.adjust.method <- "Holm-Shaffer"
    } else {
      p.adjust.method <- match.arg(p.adjust.method)
    }

    full.set <- (length(outcomes) == nrow(x)) 
    if (p.adjust.method=="none") {
      indiv.res <- testres$indiv.p.value
    } else if (p.adjust.method=="closed.set") {
      
        mytest <- function(hypotheses){
          .multiCA.test(x, scores, hypotheses)$p.value
        }
        indiv.res <- .p.adjust.closed(mytest, outcomes, remove=full.set)  
      
    } else if (p.adjust.method=="Holm-Shaffer") {
      
          s <- seq_along(testres$indiv.p.value)
          if (full.set) s[2] <- 3
          o <- order(testres$indiv.p.value)
          ro <- order(o)
          indiv.res <- pmin(1, cummax((length(outcomes) - s + 1L) * testres$indiv.p.value[o]))[ro]
      
    } 
    attr(indiv.res, "method") <- p.adjust.method
  

  return(list(overall = res, individual = indiv.res))  
}

#'@rdname multiCA.test
#'@method multiCA.test formula
#'@param formula a formula of the form \code{outcome ~ group} where \code{outcome} is a factor representing the cateogrical outcome and \code{group} is the grouping variable over which the trend is tested.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param na.action      a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#'@param weights an integer-valued variable representing the number of times each \code{outcome} - \code{group} combination was observed.
#'@export
#' @importFrom stats terms xtabs

multiCA.test.formula <- function(formula, data, subset, na.action,  weights, ...){
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    responsevar <- attr(attr(mf, "terms"), "response")
    response <- mf[[responsevar]]
    weightvar <- which(names(mf)=="(weights)")
    w <- if(length(weightvar) > 0)  mf[[weightvar]] else rep(1L, nrow(mf))
    g <- factor(mf[,-c(responsevar, weightvar)])

    tab <- xtabs(w ~ response + g)
    multiCA.test(tab, ...)
}

#'Internal functions
#'
#' These internal functions perform the closed set p-value adjustment calculation
#' for the multivariate Cochran-Armitage trend test. The logical constraint
#' on the possible number of true null hypotheses is incorporated.
#'
#'
#' @name internal
#' @importFrom bitops bitAnd
#' @keywords internal
.bit2boolean <- function (x, N) 
{
  base <- 2^(1:N - 1)
  bitAnd(x, base) != 0
}

#' @param test function that performs the local test. The function should accept a subvector of the hypotheses argument as input, and return a p-value.
#' @param  hypotheses identifiers of the collection of elementary hypotheses. 
#' @param remove logical indicator of whether hypotheses of length N-1 should be removed
#' @param  ...  additional parameters to the 'test' function
#' @return  numeric vector of adjusted p-values for each hypothesis
#' @keywords internal
#' @name internal
.p.adjust.closed <- function (test, hypotheses, remove=FALSE, ...) 
{
  N <- length(hypotheses)
  Nmax <- log2(.Machine$integer.max + 1)
  if (N > Nmax) 
    stop("no more than ", Nmax, " hypotheses supported in full closed testing.\n Use a shortcut-based test.")
  closure <- 1:(2^N - 1)
  base <- 2^(1:N - 1)
  offspring <- function(x) {
    res <- bitAnd(x, closure)
    res[res != 0]
  }
  lengths <- rowSums(sapply(base, function(bs) bitAnd(closure, bs) != 0))
  
  idx <- sort.list(lengths, decreasing = TRUE)
  closure <- closure[idx]
  lengths <- lengths[idx]
  if (remove)  closure <- closure[lengths != (N-1)]
  
  adjusted <- numeric(2^N - 1)
  for (i in closure) {
    if (adjusted[i] < 1) {
      localtest <- test(hypotheses[.bit2boolean(i,N)], ...)
      if (localtest > adjusted[i]) {
        offs <- offspring(i)
        adjusted[offs] <- pmax(adjusted[offs], localtest)
      }
    }
  }
  
  out <- adjusted[base]
  names(out) <- hypotheses
  return(out)
}

#' Non-centrality parameter for chi-square distribution
#'
#' Calculates the non-centrality parameter for a chi-square distribution for a given 
#' quantile. This is often needed for sample size calculation for chi-square based tests.
#'
#'@details The function is modeled after the SAS function CNONCT. If \code{p} is larger 
#' than the cumulative probability of the central chi-square distribution at \code{x}, then
#' there is no solution and NA is returned.
#'
#'@param x a numeric value at which the distribution was evaluated
#'@param p a numeric value giving the cumulative probability at \code{x}
#'@param df an integer giving the degrees of freedom of the chi-square variable
#'@examples
#' (ncp <- cnonct(qchisq(0.95, df=10), 0.8, df=10))
#' ## check
#' pchisq(qchisq(0.95, df=10), df=10, ncp=ncp)  ## 0.8
#'@export
#'@importFrom stats pchisq uniroot

cnonct <- function(x, p, df){
 
  if (pchisq(x, df=df) < p) return(NA)

  f <- function(ncp){pchisq(x, df=df, ncp=pmax(0,ncp)) - p}

  res <- uniroot(f, interval=c(0, 100), extendInt="downX", tol=.Machine$double.eps^0.5)
  res$root
}

#' Power calculations for the multinomial Cochran-Armitage trend test
#'
#' Given the probabilities of outcomes, compute the power of the overall multinomial 
#' Cochran-Armitage trend test or determine the sample size to obtain a target power. 
#'
#'@details 
#' The distribution of the outcomes can be specified in two ways: either the full matrix of 
#' outcome probabilities \code{pmatrix} can be specified, or exactly two of the parameters 
#' \code{p.ave}, \code{slopes}, \code{p.start}, and \code{p.end} must be specified, while 
#' 
#' @param N integer, the total sample size of the study. If \code{NULL} then \code{power} needs to be specified.
#' @param power target power. If \code{NULL} then \code{N} needs to be specified.
#' @param pmatrix numeric matrix of hypothesized outcome probabilities in each group,  with #' the outcomes as rows and ordered groups as columns. The columns should add up to 1. 
#' @param p.ave numeric vector of average probability of each outcome over the groups  
#' weighted by \code{n.prop}.
#' @param p.start,p.end numeric vectors of the probability of each outcome for the  
#' first / last ordered group
#' @param slopes numeric vector of the hypothesized slope of each outcome when regressed  
#' against the column \code{scores} wiht weights \code{n.prop}
#' @param scores non-decreasing numeric vector of the same length as the number of ordered groups  
#' giving the trend test scores. Defaults to linearly increasing values.
#' @param n.prop numeric vector describing relative sample sizes of the ordered groups.  
#' Will be normalized to sum to 1. Defaults to equal sample sizes.
#' @param G integer, number of ordered groups
#' @param sig.level significance level
#' @return object of class "power.htest"
#'
#' @examples
#' power.multiCA.test(power=0.8, p.start=c(0.1,0.2,0.3,0.4), p.end=c(0.4, 0.3, 0.2, 0.1), 
#'                      G=5, n.prop=c(3,2,1,2,3))
#'
#' ## Power of stroke study with 100 subjects per year and observed trends
#' data(stroke)
#' strk.mat <- xtabs(Freq ~ Type + Year, data=stroke)
#' power.multiCA.test(N=900, pmatrix=prop.table(strk.mat, margin=2))
#' @export
#' @importFrom stats pchisq qchisq weighted.mean


power.multiCA.test <- function(N=NULL, power=NULL, pmatrix=NULL, p.ave=NULL, p.start=NULL, 
                               p.end=NULL, slopes=NULL, scores=1:G, n.prop=rep(1, G),
                               G=length(p.ave), sig.level=0.05){
  if (sum(sapply(list(N, power), is.null)) != 1) 
        stop("exactly one of 'N',  and 'power' must be NULL")
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")


  if (!is.null(pmatrix)){
    K <- nrow(pmatrix)
    G <- ncol(pmatrix)
    if (!isTRUE(all.equal(colSums(pmatrix), rep(1, G), 
                          check.attributes=FALSE, use.names=FALSE))) 
      stop("pmatrix should have column sums of 1.")
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- as.vector(pmatrix %*% (n.prop * (scores-cbar))) / s2
    p.ave <- as.vector(pmatrix %*% n.prop)
  }
  else {
   if (sum(sapply(list(p.ave, slopes, p.start, p.end), is.null)) != 2) 
        stop("Either  pmatrix, or exactly two of 'p.ave', 'slopes', 'p.start', and 'p.end' must be specified (ie not NULL)")

  if (!is.null(p.ave) & !is.null(slopes)){
    if (length(p.ave) != length(slopes)) 
      stop("p.ave and slopes should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
  }
  else if (!is.null(p.ave) & !is.null(p.start)){
    if (length(p.ave) != length(p.start)) 
      stop("p.ave and p.start should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.start - p.ave) / (scores[1] - cbar)
  }
  else if (!is.null(p.ave) & !is.null(p.end)){
    if (length(p.ave) != length(p.end)) 
      stop("p.ave and p.end should have the same length")
    K <- length(p.ave)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.end - p.ave) / (scores[G] - cbar)
  }
  else if (!is.null(p.start) & !is.null(p.end)){
    if (length(p.start) != length(p.end)) 
      stop("p.start and p.end should have the same length")
    K <- length(p.start)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    slopes <- (p.end - p.start) / (scores[G] - scores[1])
    p.ave <- p.start - slopes * (scores[1] - cbar)
  }
  else if (!is.null(p.start) & !is.null(slopes)){
    if (length(p.start) != length(slopes)) 
      stop("p.start and slopes should have the same length")
    K <- length(p.start)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    p.ave <- p.start - slopes * (scores[1] - cbar)
  }
  else if (!is.null(p.end) & !is.null(slopes)){
    if (length(p.end) != length(slopes)) 
      stop("p.end and slopes should have the same length")
    K <- length(p.end)
    
        if (missing(G)){
          if (!missing(scores)) G <- length(scores)
          else if (!missing(n.prop)) G <- length(n.prop)
          else stop("The number of groups G needs to be specified explicitly or implicitly through the dimensions of pmatrix, the scores, or the n.prop vector.")
        }
        if (sum(n.prop) != 1) n.prop <- n.prop/sum(n.prop)
        cbar <- weighted.mean(scores, w=n.prop)
        s2 <- sum(n.prop * (scores-cbar)^2)
    
    p.ave <- p.end - slopes * (scores[G] - cbar)
  }
    
  
    if (!isTRUE(all.equal(sum(slopes), 0, check.attributes=FALSE, use.names=FALSE))) 
        stop("Implied or specified values of slopes should sum to 0.")
    if (!isTRUE(all.equal(sum(p.ave), 1, check.attributes=FALSE, use.names=FALSE))) 
        stop("Implied or specified values of p.ave should sum to 1.")
    check <- outer(1:K, 1:G, function(j,i)p.ave[j] + slopes[j]*(scores[i]-cbar))
    if (!all(check >= 0) || !(all(check <=1)))
      stop("The parameters do not define a valid probability matrix")  
  
}

  
  df <- K - 1
  crit <- qchisq(sig.level, df=df, lower.tail=FALSE)
  ncp0 <- sum(slopes^2 / p.ave) * s2 
  if (missing(power)){
    ncp <- ncp0 * N
    power <- pchisq(crit, df=df, ncp=ncp, lower.tail=FALSE)
   } 
   else {
     ncp <- cnonct(crit, p=1-power, df=df)
     N <-  ncp / ncp0
   }

   res <- structure(list(n = N, n.prop = n.prop, p.ave=p.ave, slopes = slopes, G = G, 
                         sig.level = sig.level, power = power,  
                         method = "Multinomial Cochran-Armitage trend test"), 
                     class = "power.htest")
   res
   }
