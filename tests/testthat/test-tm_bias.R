context("tm_bias")

skip_on_cran()

set.seed(4077)
B_test_dat <- as.data.frame(cbind(c(rep(0,500),rep(1,500)),
                                c(sort(rnorm(500,0,1)),sort(rnorm(500,1,1.5)))))
colnames(B_test_dat) <- c("TR", "Y")

B_test_dat$Y[which(B_test_dat$TR==0)[1:150]] <- NA
B_test_dat$Y[which(B_test_dat$TR==1)[sample(seq(1,400),
                                        200, replace=FALSE)]] <- NA

B_test_dat2 <- B_test_dat; B_test_dat2$Y[1:10] <- "Oops"
B_test_dat3 <- B_test_dat; B_test_dat3$TR[1:10] <- 3



# checking total bias
expect_equal(round(as.numeric(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                         side="LOW", spread_TG=0.4,
                         spread_CG=0.8, data=B_test_dat)$total_bias),4),
             round(1.634035,4))

# checking TM estimate
expect_equal(round(as.numeric(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                                      side="LOW", spread_TG=0.4,
                                      spread_CG=0.8, data=B_test_dat)$TM_estimate),4),
             round(1.204756,4))

# checking total bias under maximal violation of strong MNAR assumption in CG group
expect_equal(round(as.numeric(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                                      side="LOW", spread_TG=0.4,
                                      spread_CG=0.8, data=B_test_dat)$max_bias_CG[4]),4),
             round(-0.7331188,4))


# checking errors
expect_error(tm_bias(formula= Y ~ TR, "Trt", trF=0.5,
                                      side="LOW", spread_TG=0.4,
                                      spread_CG=0.8, data=B_test_dat),
             "TR variable not in data")

expect_error(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                     side="LOW", spread_TG=0.4,
                     spread_CG=0.8, data=B_test_dat2),
             "Y non-numeric")

expect_error(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                     side="LOW", spread_TG=0.4,
                     spread_CG=0.8, data=B_test_dat3),
             "TR non-binary")

expect_error(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                     side="LOW", spread_TG=0.4,
                     spread_CG=0.2, data=B_test_dat),
             "Comparator Gr spread smaller than dropout proportion")

expect_error(tm_bias(formula= Y ~ TR, "TR", trF=0.5,
                     side="LOW", spread_TG=0.1,
                     spread_CG=0.8, data=B_test_dat),
             "Treatment Gr spread smaller than dropout proportion")

expect_error(tm_bias(formula= Y ~ TR, "TR", trF=0.35,
                     side="LOW", spread_TG=0.4,
                     spread_CG=0.3, data=B_test_dat),
             "Trimming fraction smaller than largest dropout proportion")
