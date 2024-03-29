% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TM_bias.R
\name{tm_bias}
\alias{tm_bias}
\title{Calculating Bias For Trimmed Mean Linear Models:}
\usage{
tm_bias(
  formula,
  GR,
  trF,
  side = c("LOW", "HIGH"),
  spread_TG = "max_bias",
  spread_CG = "max_bias",
  data
)
}
\arguments{
\item{formula}{an object of class \code{\link{formula}}, specifying the model, of the form
\code{outcome ~ terms}, where \code{terms} must include the binary treatment variable, with additional
variables optional.}

\item{GR}{a string denoting the name of the binary treatment variable. This function assumes the
lowest value to be the comparator/reference group}

\item{trF}{a number between 0 and 1, specifying the trimming fraction: the proportion of the data that is trimmed away
for each treatment group. \code{trF} should be equal to or greater than the largest observed
dropout proportion. If left unspecified, a default trimming fraction of 0.5 is assumed.}

\item{side}{specifies if higher value trimming (\code{"HIGH"}) or lower value trimming (\code{"LOW"}) should be performed.}

\item{spread_TG}{a number between 0 and 1, specifying the dropout spread for the treatment group.
\code{spread_TG} should be equal to or greater than the observed dropout proportion. If left unspecified,
the worst-case scenario is assumed, in which dropout is located on the side of the distribution opposite from the one
that is being trimmed (\code{spread_TG="max_bias"}).}

\item{spread_CG}{a number between 0 and 1, specifying the dropout spread for the comparator group.
\code{spread_CG} should be equal to or greater than the observed dropout proportion. If left unspecified,
the worst-case scenario is assumed, in which dropout is located on the side of the distribution opposite from the one
that is being trimmed (\code{spread_CG="max_bias"}).}

\item{data}{a data frame containing the variables in the model. \code{data} should contain at least the following:
a numeric outcome variable and a binary treatment variable (numeric, character or factor).}
}
\value{
\code{tm_bias} returns an object of class \code{tm_bias}.

An object of class "\code{tm_bias}" is a list containing the following components:
\item{call}{the matched call}
\item{bias_components}{an array of bias components, including location shift assumption bias (LS),
Strong MNAR bias in the treatment group (TG) and the comparator group (CG)}
\item{total_bias}{the sum of all bias components}
\item{TM_estimate}{the trimmed means estimate of the treatment effect}
\item{bias_adj_TM_estimate}{the bias adjusted trimmed means estimate }
\item{analysis_details}{the user-specified trimming fraction, trimming side, and dropout spread in the
treatment (TG) and comparator groups (CG)}
\item{observed_TG_SD}{observed standard deviation of the treatment group (TG) outcome}
\item{observed_CG_SD}{observed standard deviation of the comparator group (CG) outcome}
\item{inferred_TG_SD}{inferred full sample standard deviation of the treatment group (TG) outcome}
\item{max_bias_CG}{an array of bias components, total bias, the bias adjusted estimate, and inferred full sample
group standard deviations, calculated under the assumption of worst-case scenario dropout, with dropout in the comparator group (CG) on the opposite
side of the distribution from the one that is being trimmed}
\item{max_bias_TG}{an array of bias components, total bias, the bias adjusted estimate, and inferred full sample
group standard deviations, calculated under the assumption of worst-case scenario dropout, with dropout in the treatment group (TG) on the opposite
side of the distribution from the one that is being trimmed}
}
\description{
\code{tm_bias} calculates the bias and the bias-adjusted estimate for a trimmed means analysis (\code{\link{tm}}) of a given
dataset, for a user-specified trimming fraction and dropout spread. \code{tm_bias} calculates, under assumption
of normally distributed outcomes, the bias components resulting from violation of the
location shift assumption and violation of the strong MNAR assumption.
}
\section{Details}{
 The trimmed means estimate is subject to two assumptions: the strong MNAR assumption requires
that all dropouts (unobserved outcome values) are located in the fraction of the distribution
that is trimmed away; the location shift assumption requires the group variances of the full sample
to be equal. The bias resulting from the violation of either assumption can be calculated under assumption
of normally distributed outcomes.

Obtaining the strong MNAR assumption bias requires an additional assumption about
the distribution of the dropouts: it is assumed that the dropouts are spread homogeneously across the specified
dropout spread. For example, under lower value trimming (\code{side="LOW"}), and a treatment group dropout
spread of 0.6 (\code{spread_TG=0.6}), any value in the bottom 60\% of the treatment group outcome distribution
is equally likely to be missing.

The specified dropout spread for a given treatment group has implications for the unobserved full sample
variance that is inferred from the observed data. For example, for an observed dropout of 0.4 and an
assumed dropout spread of 0.5, the inferred full sample variance will be larger than for an assumed
dropout spread of e.g., 0.8.

In addition to calculating the bias for a user-specified dropout spread, \code{tm_bias} also calculates
the maximal bias. For example, for lower value trimming (\code{side="LOW"}), the worst-case scenario would
involve lower value dropout in the treatment group (TG) and higher value dropout in the comparator group (CG),
and vice versa. Bias components are calculated for both scenarios. If the dropout spread
(\code{spread_TG}, \code{spread_CG}) is left unspecified for either treatment group, the function will
return only these quantities.
}

\examples{
test_dat <- as.data.frame(cbind(c(rep(0, 500), rep(1, 500)),
  c(sort(rnorm(500, 0, 1)), sort(rnorm(500, 1, 1.5)))))
colnames(test_dat) <- c("TR", "Y")
test_dat$Y[which(test_dat$TR == 0)[1:150]] <- NA
test_dat$Y[which(test_dat$TR == 1)[sample(seq(1, 400), 200, replace = FALSE)]] <- NA
tm_bias_obj <- tm_bias(formula = Y ~ TR, "TR", trF = 0.5,
                       side = "LOW", spread_TG = 0.4,
                       spread_CG = 0.6, data = test_dat)
print(tm_bias_obj)
}
