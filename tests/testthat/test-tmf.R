context("tmf")


set.seed(407774)
test_dat <- as.data.frame(cbind(c(rep(0,500),rep(1,500)),
                                c(sort(rnorm(500,0,1)),sort(rnorm(500,1,1.5))),
                                rbinom(1000,2,0.4), rnorm(1000,0,1)))
colnames(test_dat) <- c("TR", "Y", "U", "U2")
colMeans(test_dat)

test_dat0 <- test_dat
test_dat$Y[1:200] <- NA
test_dat2 <- test_dat; test_dat2$Y[1:10] <- "Oops"
test_dat3 <- test_dat; test_dat3$TR[1:10] <- 3


# checking TM estimate and adjusted TM estimate
expect_equal(round(as.numeric(tm(Y ~ TR + U + U2, GR="TR", trF=0.5, side="LOW",
                      n_perm=1000, adj_est=TRUE, data=test_dat)$coefficients[c(1,4),1]),4),
             round(c(1.482352,1.032575 ),4))

# checking default 0.5 trimming when no dropout
expect_equal(as.numeric(tm(Y ~ TR + U + U2, GR="TR", side="LOW",
                                 n_perm=1000, adj_est=TRUE, data=test_dat0)$trimfrac),
             0.5)

# checking default adaptive trimming under dropout
expect_equal(as.numeric(tm(Y ~ TR + U + U2, GR="TR", side="LOW",
                           n_perm=1000, adj_est=FALSE, data=test_dat)$trimfrac),
             sum(is.na(test_dat$Y[test_dat$TR==0]))/length(which(test_dat$TR==0)))


# checking error messages
expect_error(tm(Y ~ TR + U + U2, GR="Trt", trF=0.5, side="LOW", n_perm=1000, adj_est=TRUE, data=test_dat),
             "TR variable not in data")

expect_error(tm(Y ~ TR + U + U2, GR="TR", trF=0.4, side="LOW", n_perm=1000, adj_est=TRUE, data=test_dat),
             "Adjusted estimate can only be computed for 50% trimming")

expect_error(tm(Y ~ TR + U + U2, GR="TR", trF=0.5, side="LOW", n_perm=1000, adj_est=TRUE, data=test_dat2),
             "Y non-numeric")

expect_error(tm(Y ~ TR + U + U2, GR="TR", trF=0.5, side="LOW", n_perm=1000, adj_est=TRUE, data=test_dat3),
             "TR non-binary")

expect_error(tm(Y ~ TR + U + U2, GR="TR", trF=0.3, side="LOW", n_perm=1000, adj_est=FALSE, data=test_dat),
             "Trimming fraction smaller than largest dropout proportion")


