# The TM_bias function calculates, for a given dataset with binary treatment variable GR,
# the strong MNAR and location shift assumption biases, under assumption of normally distributed
# outcomes. Trimming fraction (trF) and side of trimming (LOWer or HIGHer value trimming) are user
# specified. The biases can be computed for any dropout spread, ranging, for each treatment group,
# from the observed dropout proportion (e.g., 0.25), to dropout across the entirety of the distribution
# (spread_TG/spread_CG=1). spread_TG indicates the dropout spread of the treatment group, spread_CG
# of the comparator group. By default, the lowest  value in the treatment variable GR is defined
# as the comparator group.
# If a dropout spread is not specified, the worst-case scenario is assumed, in which dropout is located on
# the side of the  distribution opposite from the one that is being trimmed. For example, for lower value
# trimming, the worst-case scenario would involve lower value dropout in the TG group and higher value dropout
# in the CG group, and vice versa. For both scenarios, the bias components and total bias are calculated.
# The TM_bias function returns the bias components, the total bias, the TM estimate, the adjusted TM estimate,
# and the observed and inferred treatment group SDs. Full sample SDs are calculated from the observed SDs,
# dropout proportions and specified dropout spread, under the assumption of normality.

#' @name tm_bias
#' @title Calculating Bias For Trimmed Mean Linear Models:
#'
#' @description \code{tm_bias} calculates the bias and the bias-adjusted estimate for a trimmed means analysis ([`tm`]) of a given
#' dataset, for a user-specified trimming fraction and dropout spread. \code{tm_bias} calculates, under assumption
#' of normally distributed outcomes, the bias components resulting from violation of the
#' location shift assumption and violation of the strong MNAR assumption.
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
#' @param spread_TG a number between 0 and 1, specifying the dropout spread for the treatment group.
#' \code{spread_TG} should be equal to or greater than the observed dropout proportion. If left unspecified,
#' the worst-case scenario is assumed, in which dropout is located on the side of the distribution opposite from the one
#' that is being trimmed (\code{spread_TG="max_bias"}).
#' @param spread_CG a number between 0 and 1, specifying the dropout spread for the comparator group.
#' \code{spread_CG} should be equal to or greater than the observed dropout proportion. If left unspecified,
#' the worst-case scenario is assumed, in which dropout is located on the side of the distribution opposite from the one
#' that is being trimmed (\code{spread_CG="max_bias"}).
#' @param data a data frame containing the variables in the model. \code{data} should contain at least the following:
#' a numeric outcome variable and a binary treatment variable (numeric, character or factor).
#'
#' @section Details: The trimmed means estimate is subject to two assumptions: the strong MNAR assumption requires
#' that all dropouts (unobserved outcome values) are located in the fraction of the distribution
#' that is trimmed away; the location shift assumption requires the group variances of the full sample
#' to be equal. The bias resulting from the violation of either assumption can be calculated under assumption
#' of normally distributed outcomes.
#'
#' Obtaining the strong MNAR assumption bias requires an additional assumption about
#' the distribution of the dropouts: it is assumed that the dropouts are spread homogeneously across the specified
#' dropout spread. For example, under lower value trimming (\code{side="LOW"}), and a treatment group dropout
#' spread of 0.6 (\code{spread_TG=0.6}), any value in the bottom 60% of the treatment group outcome distribution
#' is equally likely to be missing.
#'
#' The specified dropout spread for a given treatment group has implications for the unobserved full sample
#' variance that is inferred from the observed data. For example, for an observed dropout of 0.4 and an
#' assumed dropout spread of 0.5, the inferred full sample variance will be larger than for an assumed
#' dropout spread of e.g., 0.8.
#'
#' In addition to calculating the bias for a user-specified dropout spread, \code{tm_bias} also calculates
#' the maximal bias. For example, for lower value trimming (\code{side="LOW"}), the worst-case scenario would
#' involve lower value dropout in the treatment group (TG) and higher value dropout in the comparator group (CG),
#' and vice versa. Bias components are calculated for both scenarios. If the dropout spread
#' (\code{spread_TG}, \code{spread_CG}) is left unspecified for either treatment group, the function will
#' return only these quantities.
#'
#' @return \code{tm_bias} returns an object of class `tm_bias`.
#'
#' An object of class "\code{tm_bias}" is a list containing the following components:
#' \item{call}{the matched call}
#' \item{bias_components}{an array of bias components, including location shift assumption bias (LS),
#' Strong MNAR bias in the treatment group (TG) and the comparator group (CG)}
#' \item{total_bias}{the sum of all bias components}
#' \item{TM_estimate}{the trimmed means estimate of the treatment effect}
#' \item{bias_adj_TM_estimate}{the bias adjusted trimmed means estimate }
#' \item{analysis_details}{the user-specified trimming fraction, trimming side, and dropout spread in the
#' treatment (TG) and comparator groups (CG)}
#' \item{observed_TG_SD}{observed standard deviation of the treatment group (TG) outcome}
#' \item{observed_CG_SD}{observed standard deviation of the comparator group (CG) outcome}
#' \item{inferred_TG_SD}{inferred full sample standard deviation of the treatment group (TG) outcome}
#' \item{max_bias_CG}{an array of bias components, total bias, the bias adjusted estimate, and inferred full sample
#' group standard deviations, calculated under the assumption of worst-case scenario dropout, with dropout in the comparator group (CG) on the opposite
#' side of the distribution from the one that is being trimmed}
#' \item{max_bias_TG}{an array of bias components, total bias, the bias adjusted estimate, and inferred full sample
#' group standard deviations, calculated under the assumption of worst-case scenario dropout, with dropout in the treatment group (TG) on the opposite
#' side of the distribution from the one that is being trimmed}
#'
#' @examples
#' test_dat <- as.data.frame(cbind(c(rep(0, 500), rep(1, 500)),
#'   c(sort(rnorm(500, 0, 1)), sort(rnorm(500, 1, 1.5)))))
#'
#' colnames(test_dat) <- c("TR", "Y")
#'
#' test_dat$Y[which(test_dat$TR == 0)[1:150]] <- NA
#' test_dat$Y[which(test_dat$TR == 1)[sample(seq(1, 400),
#'             200, replace = FALSE)]] <- NA
#'
#' tm_bias_obj <- tm_bias(formula = Y ~ TR, "TR", trF = 0.5,
#'                        side = "LOW", spread_TG = 0.4,
#'                        spread_CG = 0.6, data = test_dat)
#'
#' print(tm_bias_obj)
#'
#' @export
#' @importFrom stats var
tm_bias <- function(formula, GR, trF, side=c("LOW", "HIGH"), spread_TG="max_bias", spread_CG="max_bias",data){

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

  if (spread_TG==1){
    spread_TG <- 0.999
  }
  if (spread_CG==1){
    spread_CG <- 0.999
  }
  if (spread_TG < TG.drop ){stop("Treatment Gr spread smaller than dropout proportion")}
  if (spread_CG < CG.drop ){stop("Comparator Gr spread smaller than dropout proportion")}


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

  TG_var <- var(data[which(data[,GR]==TG),vn[1]], na.rm=TRUE)
  CG_var <- var(data[which(data[,GR]==CG),vn[1]], na.rm=TRUE)

  lm.obj <- lm(formula,data.trim)
  beta_t <- summary(lm.obj)$coefficients[paste(GR,TG,sep=""),1]
  beta_t <- matrix(beta_t)
  colnames(beta_t) <- ""; rownames(beta_t) <- ""

  SD.func.extr <- function(spread, dr, obs.var){

    c <- spread
    fr1 <- (c-dr)/(1-dr)
    fr2 <- (1-c)/(1-dr)

    A <- (1- ((qnorm(c)*dnorm(qnorm(c)))/(c)) -
            ((dnorm(qnorm(c))/(c)))^2)
    B <- (1- ((-qnorm(c)*dnorm(qnorm(c)))/(1-c))-
            ((-dnorm(qnorm(c))/(1-c)))^2)
    C <- -(dnorm(qnorm(c)))/(c)
    D <- (dnorm(qnorm(c)))/(1-c)

    fsd <- function(sigf) {((fr1*A*sigf^2 + fr2*B*sigf^2 +
                               fr1*fr2*(C*sigf-D*sigf)^2))-
        obs.var}

    fullSD <- uniroot(fsd, lower=0.1, upper=30)$root
    return(fullSD)}


  LSA.bias <- function(TRSD, PLSD, trF, side){
    if (side=="HIGH"){
      bias <- -(TRSD - PLSD) * (dnorm(qnorm(1-trF))-dnorm(qnorm(0)))/
        (1-trF)}
    if(side=="LOW"){
      bias <- (TRSD - PLSD) * (dnorm(qnorm(1-trF))-dnorm(qnorm(0)))/
        (1-trF)
    }
    return(bias)
  }

  bias.drop <- function(SD=1,DSP=0.9, plD=0.2, trimf=0.55,
                        group="TG", side){

    DS <- 1-DSP
    calc.frac <- ((1-DS) * ((1-DS)-trimf) / ((1-DS)-plD))

    if (DSP <= trimf){
      calc.frac <- 0
    }
    calc.frac

    b.prop <- 1-trimf
    c.prop <- DS

    bst.prop <- DS+calc.frac
    a <- qnorm(DS+calc.frac,0,1)
    b <- qnorm(DS,0,1)
    a1 <- qnorm(1-trimf,0,1)

    F1 <- (((1-trimf)-DS)/(1-trimf))

    if(DSP<=trimf){
      a <- qnorm(1-trimf,0,1)
      b <- qnorm(0,0,1)
      a1 <- qnorm(1-trimf,0,1)
      F1 <- 1
    }

    elPL <- ((((dnorm(a1)-dnorm(b))/
                 (pnorm(a1)-pnorm(b)))*-SD)-
               ((dnorm(a)-dnorm(b))/
                  (pnorm(a)-pnorm(b))*-SD)) *F1
    elTR <- (((dnorm(a)-dnorm(b))/
                (pnorm(a)-pnorm(b))*-SD) -
               ((((dnorm(a1)-dnorm(b))/
                    (pnorm(a1)-pnorm(b)))*-SD))) *F1

    if(group=="TG"){
      el <- elTR
    } else { if (group=="CG"){
      el <- elPL}
    }

    if(group!="CG" & group!="TG"){
      el <- NA}

    if(side=="HIGH"){
      el <- el
    }
    if(side=="LOW"){
      el <- -1*el
    }

    return(el)}



  if(is.numeric(spread_CG) & is.numeric(spread_TG)){

    TG_SD_inf <- SD.func.extr(spread=spread_TG, dr=TG.drop,
                              obs.var=TG_var)
    CG_SD_inf <- SD.func.extr(spread=spread_CG, dr=CG.drop,
                              obs.var=CG_var)

    LSA.bias.out <- LSA.bias(TG_SD_inf, CG_SD_inf, trF, side=side)

    bias.drop.out.TG <- bias.drop(TG_SD_inf, spread_TG, TG.drop, trF, group="TG", side=side)
    bias.drop.out.CG <- bias.drop(CG_SD_inf, spread_CG, CG.drop, trF, group="CG", side=side)

    tot_bias <- LSA.bias.out + bias.drop.out.TG + bias.drop.out.CG

    bias_adj_est <- beta_t - tot_bias

    bias.tab <- matrix(c("LS", paste("Strong MNAR TG group (", paste(TG, ")", sep=""), sep=""),
                         paste("Strong MNAR CG group (", paste(CG, ")", sep=""), sep=""),
                         LSA.bias.out,bias.drop.out.TG,bias.drop.out.CG), ncol=2, nrow=3)
    colnames(bias.tab) <- c("Bias type", "Bias")
    rownames(bias.tab) <- c("","","")

    tot_bias <- matrix(tot_bias)
    colnames(tot_bias) <- ""; rownames(tot_bias) <- ""
    bias_adj_est <- matrix(bias_adj_est)
    colnames(bias_adj_est) <- ""; rownames(bias_adj_est) <- ""

  } else {
    bias.tab <- NULL
    TG_SD_inf <- matrix(NA) ; rownames(TG_SD_inf) <- colnames(TG_SD_inf) <- ""
    CG_SD_inf <- matrix(NA) ; rownames(CG_SD_inf) <- colnames(CG_SD_inf) <- ""
    tot_bias <- matrix(NA)
    colnames(tot_bias) <- ""; rownames(tot_bias) <- ""
    bias_adj_est <- matrix(NA)
    colnames(bias_adj_est) <- ""; rownames(bias_adj_est) <- ""
  }


  TG_SD_max <- SD.func.extr(spread=TG.drop, dr=TG.drop,
                            obs.var=TG_var)
  CG_SD_max <- SD.func.extr(spread=CG.drop, dr=CG.drop,
                            obs.var=CG_var)

  bias.max <- function(SD, dr, trF, viol.group, side){

    if(side=="LOW" & viol.group=="CG"){
      bias <- SD/(1-dr)* (dnorm(qnorm(trF))--dnorm(qnorm(1-dr))-dnorm(qnorm(1-dr-(1-trF))))}
    if(side=="LOW" & viol.group=="TG"){
      bias <- -SD/(1-dr)* (dnorm(qnorm(trF))--dnorm(qnorm(1-dr))-dnorm(qnorm(1-dr-(1-trF))))}
    if(side=="HIGH" & viol.group=="CG"){
      bias <- -SD/(1-dr)* (dnorm(qnorm(trF))--dnorm(qnorm(1-dr))-dnorm(qnorm(1-dr-(1-trF))))}
    if(side=="HIGH" & viol.group=="TG"){
      bias <- SD/(1-dr)* (dnorm(qnorm(trF))--dnorm(qnorm(1-dr))-dnorm(qnorm(1-dr-(1-trF))))}

    return(bias)
  }

  bias.max.viol.CG <- bias.max(CG_SD_max,CG.drop,trF, "CG", side)
  bias.max.viol.TG <- bias.max(TG_SD_max,TG.drop,trF, "TG", side)
  bias.max.LS.comp <- LSA.bias(TG_SD_max, CG_SD_max, trF, side)


  maximum_bias_CG <- matrix(c(bias.max.viol.CG,
                              bias.max.LS.comp,
                              (bias.max.viol.CG+bias.max.LS.comp),
                              beta_t-(bias.max.viol.CG+bias.max.LS.comp),
                              TG_SD_max, CG_SD_max), ncol=1)
  colnames(maximum_bias_CG) <- ""
  rownames(maximum_bias_CG) <- c("Strong MNAR bias", "LS bias", "Total bias", "Bias adjusted estimate",
                                 paste("Inferred SD TG group (", paste(TG, ")", sep=""),  sep=""),
                                 paste("Inferred SD CG group (", paste(CG, ")", sep=""),  sep=""))

  maximum_bias_TG <- matrix(c(bias.max.viol.TG,
                              bias.max.LS.comp,
                              (bias.max.viol.TG+bias.max.LS.comp),
                              beta_t-(bias.max.viol.TG+bias.max.LS.comp),
                              TG_SD_max, CG_SD_max), ncol=1)
  colnames(maximum_bias_TG) <- ""
  rownames(maximum_bias_TG) <- c("Strong MNAR bias", "LS bias", "Total bias", "Bias adjusted estimate",
                                 paste("Inferred SD TG group (", paste(TG, ")", sep=""),  sep=""),
                                 paste("Inferred SD CG group (", paste(CG, ")", sep=""),  sep=""))




  if(side=="LOW"){
    trimside <- "lower"

  }
  if(side=="HIGH"){
    trimside <- "higher"
  }
  names(trimside) <- "trimming side"

  if(spread_CG==0.999){
    spread_CG <- 1
  }

  if(spread_TG==0.999){
    spread_TG <- 1
  }

  if(is.numeric(spread_CG) & is.numeric(spread_TG)){
    analysis_details <- matrix(c(paste(paste(round(trF*100,1), "%", sep=""), "trimming,", sep=" "),
                                 paste(paste("under assumption of", trimside, sep=" "), "value dropout,", sep=" "),
                                 paste(paste("with", paste(round(spread_TG*100,1), "% dropout spread in the TG group (", sep=""), sep=" "),
                                       paste(TG, "),", sep=""), sep=""),
                                 paste(paste("and", paste(round(spread_CG*100,1), "% dropout spread in the CG group (", sep=""), sep=" "),
                                       paste(CG, ")", sep=""), sep="")), ncol=1)
    colnames(analysis_details) <- ""
    rownames(analysis_details) <- c("","","","")
  } else {
    analysis_details <- matrix(c(paste(paste(paste(round(trF*100,1), "%", sep=""), "trimming,", sep=" "),
                                       paste(paste("under assumption of", trimside, sep=" "), "value dropout", sep=" "),sep=" ")))
    colnames(analysis_details) <- ""
    rownames(analysis_details) <- ""
  }


  if(is.numeric(spread_CG) & is.numeric(spread_TG)){
    infSDTG <- matrix(c(paste(paste(paste(paste("Inferred TG group (", paste(paste(TG, ")", sep=""),"full sample SD for", sep=" "), sep=" "),
                                          paste(round(TG.drop*100,1), "%", sep=""), sep=" "),"dropout in the", sep=" "),
                              paste(paste(trimside, paste(round(spread_TG*100,1),"%", sep=""), sep=" "), "of the distribution", sep=" "), sep=" "),
                        TG_SD_inf), ncol=1)
    colnames(infSDTG) <- ""; rownames(infSDTG) <- c("","")

    infSDCG <- matrix(c(paste(paste(paste(paste("Inferred CG group (", paste(paste(CG, ")", sep=""),"full sample SD for", sep=" "), sep=" "),
                                          paste(round(CG.drop*100,1), "%", sep=""), sep=" "),"dropout in the", sep=" "),
                              paste(paste(trimside, paste(round(spread_CG*100,1),"%", sep=""), sep=" "), "of the distribution", sep=" "), sep=" "),
                        CG_SD_inf), ncol=1)
    colnames(infSDCG) <- ""; rownames(infSDCG) <- c("","")
  } else {
    infSDTG <- infSDCG <- NA
  }

  obsSDTG <- matrix(c(paste(paste(paste("Observed TG group (", paste(paste(TG, ")", sep=""),"SD for", sep=" "), sep=" "),
                                  paste(round(TG.drop*100,1), "%", sep=""), sep=" "),"dropout", sep=" "),
                      sqrt(TG_var)), ncol=1)
  colnames(obsSDTG) <- ""; rownames(obsSDTG) <- c("","")

  obsSDCG <- matrix(c(paste(paste(paste("Observed CG group (", paste(paste(CG, ")", sep=""),"SD for", sep=" "), sep=" "),
                                  paste(round(CG.drop*100,1), "%", sep=""), sep=" "),"dropout", sep=" "),
                      sqrt(CG_var)), ncol=1)
  colnames(obsSDCG) <- ""; rownames(obsSDCG) <- c("","")

  out_fin <- list()
  out_fin$call <- cl
  out_fin$bias_components <- bias.tab
  out_fin$total_bias <- tot_bias
  out_fin$TM_estimate <- beta_t
  out_fin$bias_adj_TM_estimate <- bias_adj_est
  out_fin$analysis_details <- analysis_details
  out_fin$observed_TG_SD <- obsSDTG
  out_fin$observed_CG_SD <- obsSDCG
  out_fin$inferred_TG_SD <- infSDTG
  out_fin$inferred_CG_SD <- infSDCG
  out_fin$max_bias_CG <- maximum_bias_CG
  out_fin$max_bias_TG <- maximum_bias_TG


  class(out_fin) <- "tm_bias"

  return(out_fin)

}

#' @export
print.tm_bias <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Analysis details:\n")
  print.default(x$analysis_details, quote=FALSE)

  cat("\n")

  if (length(x$bias_components)) {
    cat("Bias components:\n")
    x$bias_components[,2] <- format(as.numeric(x$bias_components[,2]), digits=digits)
    print.default(x$bias_components, print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No bias components: dropout spread has not been specified and only the maximum bias has been computed \n")
  cat("\n")

  cat("Total bias:\n")
  print.default(format(x$total_bias, digits=digits), quote=FALSE)

  cat("\n")

  cat("TM estimate:\n")
  print.default(format(x$TM_estimate, digits=digits), quote=FALSE)

  cat("\n")

  cat("Bias adjusted TM estimate:\n")
  print.default(format(x$bias_adj_TM_estimate, digits=digits), quote=FALSE)

  cat("\n")

  cat("Observed SDs:\n")

  x$observed_TG_SD[2,] <- format(as.numeric(x$observed_TG_SD[2,]), digits=digits)
  print.default(x$observed_TG_SD, quote=FALSE)
  x$observed_CG_SD[2,] <- format(as.numeric(x$observed_CG_SD[2,]), digits=digits)
  print.default(format(x$observed_CG_SD, digits=digits), quote=FALSE)

  cat("\n")

  cat("Inferred SDs:\n")
  if(!(is.null(dim(x$inferred_TG_SD)))){
    x$inferred_TG_SD[2,] <- format(as.numeric(x$inferred_TG_SD[2,]), digits=digits)
    print.default(format(x$inferred_TG_SD, digits=digits), quote=FALSE)
  } else {
    print.default(x$inferred_TG_SD, quote=FALSE)
  }
  if(!(is.null(dim(x$inferred_CG_SD)))){
    x$inferred_CG_SD[2,] <- format(as.numeric(x$inferred_CG_SD[2,]), digits=digits)
    print.default(format(x$inferred_CG_SD, digits=digits), quote=FALSE)
  } else {
    print.default(x$inferred_CG_SD, quote=FALSE)
  }

  cat("\n")

  max.bias.CG.title <- paste("Bias under maximal violation of the strong MNAR assumption in the",
                             substr(x$analysis_details[4,], (nchar(x$analysis_details[4,])-1)-10,
                                    (nchar(x$analysis_details[4,]))), sep=" ")

  cat(max.bias.CG.title)
  print.default(format(x$max_bias_CG, digits=digits), quote=FALSE)

  cat("\n")

  max.bias.TG.title <- paste("Bias under maximal violation of the strong MNAR assumption in the",
                             substr(x$analysis_details[3,], (nchar(x$analysis_details[3,])-1)-11,
                                    (nchar(x$analysis_details[3,])-1)), sep=" ")
  cat(max.bias.TG.title)
  print.default(format(x$max_bias_TG, digits=digits), quote=FALSE)

  invisible(x)
}


