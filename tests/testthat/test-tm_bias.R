context("tm_bias")

#test_dat <- as.data.frame(cbind(c(rep(0,500),rep(1,500)),
#c(sort(rnorm(500,0,1)),sort(rnorm(500,1,1.5))),
#rbinom(1000,2,0.4), rnorm(1000,0,1)))
#colnames(test_dat) <- c("TR", "Y", "U", "U2")
#test_dat$Y[which(test_dat$TR==0)[1:150]] <- NA
#test_dat$Y[which(test_dat$TR==1)[sample(seq(1,400),200, replace=FALSE)]] <- NA
#out1 <- tm(Y ~ TR + U + U2, GR="TR", trF=0.3, side="LOW", n_perm=1000, adj_est=FALSE, data=test_dat)
#out1$trimfrac
#sum(is.na(test_dat$Y[test_dat$TR==0]))/length(which(test_dat$TR==0))
