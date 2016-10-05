setwd(normalizePath(dirname(R.utils::commandArgs(asValues=TRUE)$"f")))
source("../../../scripts/h2o-r-test-setup.R")
##
# Test out the cor() functionality
##


test.cor <- function() {
  #1 column case
  print("Running one column case with no NA")
  n <- 20
  g <- runif(n)
  h <- runif(n)
  run.cor.tests(g,h,one_row=TRUE)

#  print("Running one column case with an NA")
#  g[4] <- NA
#  run.cor.tests(g,h,one_row=TRUE, has_nas=TRUE)

  #Matrices
  print("Running with no NA rows")
  g <- matrix(runif(n),nrow=4)
  h <- matrix(runif(n),nrow=4)
  run.cor.tests(g,h)

  #Run with NAs. R's native cor() does not have a na.rm arg, so we must test it this way
#  print("Running with 1 NA Row")
#  g[2,3] <- NA
#  run.cor.tests(g,h,has_nas=TRUE)

  print("Running with 2 NA Row")
  g[2,3] <- NA
  g[1,4] <- NA
  run.cor.tests(g,h,has_nas=TRUE)

  print("Running with 3 NA Row")
  g[2,3] <- NA
  g[1,4] <- NA
  h[2,3] <- NA
  run.cor.tests(g,h,has_nas=TRUE)

  print("Running with 4 NA Row")
  g[2,3] <- NA
  g[1,4] <- NA
  h[2,3] <- NA
  h[1,4] <- NA
  run.cor.tests(g,h,has_nas=TRUE)

}

run.cor.tests <- function (g,h,one_row=FALSE,has_nas=FALSE) {
  h2o_g <- as.h2o(g)
  h2o_h <- as.h2o(h)
  h2o_g = h2o.rbind(h2o_g,h2o_g,h2o_g) #Ensure across chunks
  h2o_h = h2o.rbind(h2o_h,h2o_h,h2o_h) #Ensure across chunks
  uses <- c("everything","all.obs","complete.obs") #Add complete.obs method
  if (has_nas) uses <- uses[-2]
  for (na.rm in c(FALSE, TRUE)) {
    for (use in uses) {
      if (one_row) {
        h2o_cor <- h2o.cor(as.h2o(g),as.h2o(h), na.rm = na.rm,use = use)
      } else {
        h2o_cor <- h2o.cor(h2o_g, h2o_h, na.rm = na.rm, use = use)
      }

      R_cor <- stats::cor(g, h, use = use)
      print("H2O cor():")
      print(paste0("Use: ",use))
      print(h2o_cor)

      print("R cor():")
      print(paste0("Use: ",use))
      print(R_cor)
      h2o_and_R_equal(h2o_cor, R_cor)
    }
  }
}


doTest("Test out the cor() functionality", test.cor)
