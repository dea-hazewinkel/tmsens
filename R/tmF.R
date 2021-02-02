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

#'@import "stats"


tmF <- function(formula, GR, trF=NULL, side=c("LOW","HIGH"), n_perm=1000, adj_est=FALSE, data){

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
      perm.est <- summary(lm(formula, data=data.trim.perm))$coefficients[var1,1]
      return(perm.est)}

    perm.testing <- replicate(n_perm, perm.f(data.trim, var))
    lm.obj <- lm(formula,data.trim)
    if(var==GR){
      var1 <- paste(var, TG, sep="")
    } else {var1 <- var}
    beta_t <- summary(lm.obj)$coefficients[var1,1]
    Pval <- (length(perm.testing)-sum(beta_t>perm.testing))/length(perm.testing)
    sd.perm <- sd(perm.testing)
    conf.int <- c(beta_t-sd.perm*1.96, beta_t+sd.perm*1.96)
    out <- c(beta_t, Pval, conf.int)
    names(out) <- c("Estimate", "Pval", "95% CI LO", "95% CI HI")
    return(out)
  }

  perm.out <- list()
  for (i in 1:length(vn[-1])){
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

    SD.oth <- sqrt((sd(dat.oth[, vn[1]])^2)/(1-2/pi))

    x1 <- dat.resc[, vn[1]]

    if (side=="LOW"){
      x2 <- min(x1) - (abs((min(x1)-x1)))}
    if (side=="HIGH"){
      x2 <- min(x1) + (abs((min(x1)+x1)))}

    x3 <- c(x1,x2)
    x3 <- x3 - mean(x3)
    x4 <- x3/(sd(x3)/SD.oth)
    x4 <- x4 - mean(x4)
    x5 <- x4 + mean(x3)
    x6 <- x5[1:length(x1)]
    dat.resc[,vn[1]] <- x6
    dat.trim.resc <- rbind(dat.resc,dat.oth)
    dat.trim.resc$TR <- as.factor(dat.trim.resc$TR)
    perm.out.TR.adj <- t(data.frame(perm.func(dat.trim.resc, GR, n_perm)))
    rownames(perm.out.TR.adj) <- paste(GR, "adj", sep="")

  } else {
    proc <- NA
  }

  sd(c(rep(NA,10), rnorm(100,0,1)))
  sdY.trim <- c(sd(data.TGtrim[,vn[1]]),sd(data.CGtrim[,vn[1]]))
  sdY.obs <- c(sd(data[which(data[,GR]==TG),vn[1]], na.rm=TRUE),sd(data[which(data[,GR]==CG),vn[1]], na.rm=TRUE))
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
  final.out$`SD outcome` <- sd_tab

  class(final.out) <- "tmF"

  return(final.out)}


# print function (modelled after lm.print)
print.tmF <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x)[,1], digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")

  cat(paste(paste(round(x$trimfrac*100,1), "%", sep=""), x$trimside, sep=" "))

  invisible(x)
}

# summary function (not sure?)
summary.tmF <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  ans <- x
  class(ans) <- "summary.tmF"
  return(ans)
}


