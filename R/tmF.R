# - formula is the regression equation as you'd specify in e.g. the lm function
# - GR gives the name of your treatment variable (here TR)
# - Lowest GR value is by default assumed to be the comparator group
# - Default trimming fraction trF is adaptive, another trimming fraction can be specified,
#   but nothing smaller than the largest observed dropout proportion.
# - If no trimming fraction is specified, and there is no dropout, then the default is 50% trimming
# - Trimming side is user specified (LOWer value trimming or HIGHer value trimming)
# - n_perm is the number of permutation to obtain p value/ 95% CI. Default value is 1000
# - The adjusted estimate can only be computed for 50% trimming. Default is no calculation of the adjusted
# - estimate (adj_est=FALSE).
# details. The trimmed means estimate is subject to two assumptions: the strong MNAR assumption
# requires that all dropouts (unobserved outcome values) are located in the fraction of the distribution
# that is trimmed away; the location shift assumption requires the group variances of the full sample
# to be equal. The adjusted trimmed means estimator relaxes the latter, but assumes normally distributed
# outcomes.

#' @name tm
#' @title Fitting Trimmed Mean Linear Models:
#'
#' @description \code{tm} performs a trimmed means analysis for data with a continuous outcome/response and a binary
#' treatment/exposure variable. Outcomes are sorted and trimmed per treatment group, and a linear
#' regression is fitted using [`lm`].
#'
#' @param formula an object of class [`formula`], specifying the model, of the form
#' \code{outcome ~ terms}, where \code{terms} must include the binary treatment variable, with additional
#' variables optional.
#' @param GR a string denoting the name of the binary treatment variable. This function assumes the
#' lowest value to be the comparator/reference group
#' @param trF a number between 0 and 1, specifying the trimming fraction: the proportion of the data that is trimmed away
#' for each treatment group. \code{trF} should be equal to or greater than the largest observed
#' dropout proportion. If left unspecified, a default trimming fraction of 0.5 is assumed.
#' @param side specifies if higher value trimming (`"HIGH"`) or lower value trimming (`"LOW"`) should be performed.
#' @param n_perm the number of permutations performed to obtain the p-value and 95% confidence intervals
#' for the estimates. Default is 1000.
#' @param adj_est logical. If \code{TRUE} the adjusted trimmed means estimate is computed. The default is `FALSE`.
#' @param data a data frame containing the variables in the model. \code{data} should contain at least the following:
#' a numeric outcome variable and a binary treatment variable (numeric, character or factor).
#'
#' @section Details: The trimmed means estimate is subject to two assumptions: the strong MNAR assumption requires
#' that all dropouts (unobserved outcome values) are located in the fraction of the distribution
#' that is trimmed away; the location shift assumption requires the group variances of the full sample
#' to be equal. The adjusted trimmed means estimator relaxes the latter, but assumes normally
#' distributed outcomes. The adjustment is performed on the group with the smallest dropout proportion.
#'
#' The p-value and 95% confidence intervals for the trimmed means estimate and the adjusted trimmed means
#' estimate are obtained in a permutation approach.
#'
#' @return \code{tm} returns an object of class `tm`.
#' The function \code{summary} is used to obtain a summary of the results. The generic accessor function
#' \code{coefficients} extracts the regression coefficients with corresponding p-values and 95% confidence intervals.
#'
#' An object of class "\code{tm}" is a list containing the following components:
#' \item{call}{the matched call}
#' \item{n}{the number of observations per treatment group}
#' \item{dropout}{the proportion of dropout per treatment group}
#' \item{trimfrac}{the proportion of data that was trimmed away per treatment group}
#' \item{trimside}{specifies if lower or higher value trimming was performed}
#' \item{n_after_trimming}{the number of observations per treatment group after trimming}
#' \item{coefficients}{an array of coefficients with corresponding p-values and 95% confidence intervals}
#' \item{Analysis_details}{reiterates trimming fraction and side, and, for adjest=TRUE specifies if the adjustment was performed on the comparator or treatment group.}
#' \item{SD_outcome}{an array of the standard deviation per treatment group, for the observed outcomes and for the trimmed outcomes}
#'
#' @examples
#' set.seed(123456)
#' test_dat <- as.data.frame(cbind(c(rep(0, 500), rep(1, 500)),
#'                           c(sort(rnorm(500, 0, 1)), sort(rnorm(500, 1, 1.5))),
#'                           rbinom(1000, 2, 0.4), rnorm(1000, 0, 1)))
#' colnames(test_dat) <- c("TR", "Y", "U", "U2")
#' test_dat$Y[1:200] <- NA
#' # Note that we usually recommend setting n_perm to a larger value, e.g., 1000
#' tm_obj <- tm(formula= Y ~ TR + U + U2,
#'              GR = "TR", trF = 0.5, side = "LOW",
#'              n_perm = 100, adj_est = TRUE, data = test_dat)
#' print(tm_obj)
#' summary(tm_obj)
#' @export
tm <- function(formula, GR, trF=NULL, side=c("LOW","HIGH"), n_perm=1000, adj_est=FALSE, data){

  cl <- match.call()

  vn <- all.vars(formula)

  if(!(GR %in% vn)){stop("TR variable not in data")}
  if (is.numeric(data[,vn[1]])==FALSE){ stop("Y non-numeric")}
  if (length(unique(data[,GR]))>2){ stop("TR non-binary")}

  TR <- as.factor(data[,GR])
  CG <- sort(levels(TR))[1]
  TG <- sort(levels(TR))[2]

  data.CG <- data[which(data[,GR]==CG),]
  data.TG <- data[which(data[,GR]==TG),]

  CG.drop <- sum(is.na(data.CG[,vn[1]]))/nrow(data.CG)
  TG.drop <- sum(is.na(data.TG[,vn[1]]))/nrow(data.TG)

  drop <- max(CG.drop,TG.drop)

  if(is.null(trF)){
    trF <- drop
    if(drop==0){
      trF=0.5
    }
  } else {
    if (drop>trF){ stop("Trimming fraction smaller than largest dropout proportion")}
  }

  rems.TG <- ceiling(nrow(data.TG)*trF)
  rems.CG <- ceiling(nrow(data.CG)*trF)

  if(side=="LOW"){
    data.TG[is.na(data.TG[,vn[1]]),vn[1]] <- -Inf
    data.CG[is.na(data.CG[,vn[1]]),vn[1]] <- -Inf
    data.CG <- data.CG[order(data.CG[,vn[1]]),]
    data.TG <- data.TG[order(data.TG[,vn[1]]),]
    data.TGtrim <- data.TG[-(1:rems.TG),]
    data.CGtrim <- data.CG[-(1:rems.CG),]
  }


  if(side=="HIGH"){
    data.TG[is.na(data.TG[,vn[1]]),vn[1]] <- Inf
    data.CG[is.na(data.CG[,vn[1]]),vn[1]] <- Inf
    data.CG <- data.CG[order(data.CG[,vn[1]]),]
    data.TG <- data.TG[order(data.TG[,vn[1]]),]
    data.TGtrim <- data.TG[-((nrow(data.TG)-rems.TG):nrow(data.TG)),]
    data.CGtrim <- data.CG[-((nrow(data.CG)-rems.CG):nrow(data.CG)),]
  }

  data.trim <- rbind(data.TGtrim,data.CGtrim)
  data.trim$TR <- as.factor(data.trim$TR)

  perm.func <- function(data.trim, var, n_perm){

    perm.f <- function(data.trim, var){
      new.tr <- sample(data.trim[,var], replace=FALSE)
      data.trim.perm <- data.trim
      data.trim.perm[,var] <- new.tr
      if(var==GR){
        var1 <- paste(var, TG, sep="")
      } else {var1 <- var}
      perm.est <- summary(stats::lm(formula, data=data.trim.perm))$coefficients[var1,1]
      return(perm.est)}

    perm.testing <- replicate(n_perm, perm.f(data.trim, var))
    lm.obj <- stats::lm(formula,data.trim)
    if(var==GR){
      var1 <- paste(var, TG, sep="")
    } else {var1 <- var}
    beta_t <- summary(lm.obj)$coefficients[var1,1]
    Pval <- (length(perm.testing)-sum(beta_t>perm.testing))/length(perm.testing)
    sd.perm <- stats::sd(perm.testing)
    conf.int <- c(beta_t-sd.perm*1.96, beta_t+sd.perm*1.96)
    out <- c(beta_t, Pval, conf.int)
    names(out) <- c("Estimate", "Pval", "95% CI LO", "95% CI HI")
    return(out)
  }

  perm.out <- list()
  for (i in seq_along(vn[-1])){
    perm.out[[i]] <- perm.func(data.trim, vn[1+i], n_perm)
  }

  perm.out <- do.call("rbind", perm.out)
  rownames(perm.out) <- vn[2:length(vn)]

  # specifying adjusted estimator on rescaling group with least dropout
  if(adj_est==TRUE){
    if(trF != 0.5){stop("Adjusted estimate can only be computed for 50% trimming")}

    if(CG.drop <= TG.drop){
      dat.resc <- data.CGtrim
      dat.oth <- data.TGtrim
      proc <- paste("Comparator (", ") group rescaled", sep=CG)
    } else {
      dat.resc <- data.TGtrim
      dat.oth <- data.CGtrim
      proc <- paste("Treatment (", ") group rescaled", sep=TG)
    }

    SD.oth <- sqrt((stats::sd(dat.oth[, vn[1]])^2)/(1-2/pi))

    x1 <- dat.resc[, vn[1]]

    if (side=="LOW"){
      x2 <- min(x1) - (abs((min(x1)-x1)))}
    if (side=="HIGH"){
      x2 <- min(x1) + (abs((min(x1)+x1)))}

    x3a <- c(x1,x2)
    x3 <- x3a - mean(x3a)
    x4 <- x3/(stats::sd(x3)/SD.oth)
    x4 <- x4 - mean(x4)
    x5 <- x4 + mean(x3a)
    x6 <- x5[seq_along(x1)]
    dat.resc[,vn[1]] <- x6
    dat.trim.resc <- rbind(dat.resc,dat.oth)
    dat.trim.resc$TR <- as.factor(dat.trim.resc$TR)
    perm.out.TR.adj <- t(data.frame(perm.func(dat.trim.resc, GR, n_perm)))
    rownames(perm.out.TR.adj) <- paste(GR, "adj", sep="")

  } else {
    proc <- NA
  }

  stats::sd(c(rep(NA,10), stats::rnorm(100,0,1)))
  sdY.trim <- c(stats::sd(data.TGtrim[,vn[1]]),stats::sd(data.CGtrim[,vn[1]]))
  sdY.obs <- c(stats::sd(data[which(data[,GR]==TG),vn[1]], na.rm=TRUE),stats::sd(data[which(data[,GR]==CG),vn[1]], na.rm=TRUE))
  sd_tab <- t(matrix(c(sdY.obs,sdY.trim), nrow=2, ncol=2))
  rownames(sd_tab) <- c(paste("Observed outcomes", vn[1], sep=" "), paste("Trimmed outcomes", vn[1], sep=" "))
  colnames(sd_tab) <- c(TG,CG)

  props <- t(matrix(
    c(TG.drop, CG.drop)))
  rownames(props) <- c("p")
  colnames(props) <- c(TG,CG)

  names(trF) <- "trimming fraction"

  n_tot <- t(matrix(c(sum(TR==TG),sum(TR==CG))))
  rownames(n_tot) <- c("N")
  colnames(n_tot) <- c(TG,CG)

  n_post_trim <- t(matrix(c(nrow(data.TGtrim), nrow(data.CGtrim))))
  rownames(n_post_trim) <- c("N")
  colnames(n_post_trim) <- c(TG,CG)

  if(side=="LOW"){
    trimside <- "lower value trimming"
  }
  if(side=="HIGH"){
    trimside <- "higher value trimming"
  }
  names(trimside) <- "trimming side"

  tp <- paste(paste(round(trF*100,1), "%", sep=""), trimside, sep=" ")
  tp1 <- paste(as.character(proc), "for adjusted TM estimate", sep=" ")
  an.det <- (matrix(c(tp, tp1), nrow=2,ncol=1))
  rownames(an.det) <- c("", "")
  colnames(an.det) <- c("Analysis details")

  final.out <- list()
  final.out$call <- cl
  final.out$n <- n_tot
  final.out$dropout <- props
  final.out$trimfrac <- trF
  final.out$trimside <- trimside
  final.out$`n_after_trimming` <- n_post_trim
  if(adj_est==TRUE){
    final.out$coefficients <- rbind(perm.out, perm.out.TR.adj)
  } else {
    final.out$coefficients <- perm.out
  }
  final.out$`Analysis_details` <- an.det
  final.out$`SD_outcome` <- sd_tab

  class(final.out) <- "tm"

  return(final.out)}


#' @export
print.tm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if (length(stats::coef(x))) {
    cat("Coefficients:\n")
    coefs <- x$coefficients
    stats::printCoefmat(coefs, digits = digits, na.print = "NA", ...)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}


#' @name summary.tm
#' @title Summarizing Trimmed Means Linear Model fits:
#'
#' @description \code{summary} method for class "\code{tm}".
#'
#' @param object an object of class "\code{tm}"
#' @param ... user specified arguments
#'
#' @return \code{summary.tm} returns an list of summary statistics of the fitted trimmed means linear
#' model in \code{object}, with components
#' \item{call}{the matched call}
#' \item{n}{the number of observations per treatment group}
#' \item{dropout}{the proportion of dropout per treatment group}
#' \item{trimfrac}{the proportion of data that was trimmed away per treatment group}
#' \item{trimside}{specifies if lower or higher value trimming was performed}
#' \item{n_after_trimming}{the number of observations per treatment group after trimming}
#' \item{coefficients}{an array of coefficients with corresponding p-values and 95% confidence intervals}
#' \item{Analysis_details}{reiterates trimming fraction and side, and, for `adjest=TRUE` specifies if the adjustment was performed on the comparator or treatment group.}
#' \item{SD_outcome}{an array of the standard deviation per treatment group, for the observed outcomes and for the trimmed outcomes}
#'
#' @seealso [`tm`]. The function [`coef`]
#' extracts the array of regression coefficients with corresponding p-values and 95% confidence intervals.
#'
#' @examples
#' set.seed(123456)
#' test_dat <- as.data.frame(cbind(c(rep(0,500),rep(1,500)),
#'                           c(sort(rnorm(500,0,1)),sort(rnorm(500,1,1.5))),
#'                           rbinom(1000,2,0.4), rnorm(1000,0,1)))
#' colnames(test_dat) <- c("TR", "Y", "U", "U2")
#' test_dat$Y[1:200] <- NA
#' tm_obj <- tm(formula= Y ~ TR + U + U2, GR = "TR", trF = 0.5,
#'              side = "LOW", n_perm = 1000, adj_est = TRUE, data = test_dat)
#' summary(tm_obj)
#' coef(tm_obj)
#' @export
summary.tm <- function (object, ...)
{
  ans <- object
  class(ans) <- "summary.tm"
  ans
}

#' @export
print.summary.tm <- function (x,
                              digits = max(3L, getOption("digits") - 3L),
                              ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("\nAnalysis details:\n")
  cat(x$`Analysis_details`[2], "\n\n", sep = "")

  if (length(stats::coef(x))) {
    cat("Coefficients:\n")
    coefs <- x$coefficients
    stats::printCoefmat(coefs, digits = digits, na.print = "NA", ...)
  } else {
    cat("No coefficients\n")
  }

  cat("\n\nDropout:\n")
  cat(format(x$dropout[[2]] * 100, digits = digits), "%", sep="")

  cat("\n\nSample size after trimming:\n")
  cat(format(x$n_after_trimming[[1]], digits = digits))

  cat("\n\nTrimming fraction: \n",
      format(x$trimfrac*100, digits = digits),
      "% ", x$trimside, "\n\n", sep = "")

  invisible(x)
}
