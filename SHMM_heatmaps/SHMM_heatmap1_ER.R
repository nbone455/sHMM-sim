# true model threshold 

## Mk vs ThreshMK


# load packages and functions
library(ape)
library(treeplyr)
library(textshape)
library(vctrs)
library(stringr)
library("phytools")
require(dplyr)
require(corHMM)
library(geiger)
library(ggplot2)
library(gplots)
library(phytools)

nsim = 50

setwd("~/Documents/hiThresh power/R")


#old matrix maker functions from corhmm
rate.mat.maker <- function (rate.cat, hrm = TRUE, ntraits = NULL, nstates = NULL, model = c("ER", "SYM", "ARD")){
  if (hrm == TRUE) {
    k = 2
    mat1 <- matrix(NA, k * rate.cat, k * rate.cat)
    mat2 <- matrix(NA, k * rate.cat, k * rate.cat)
    vec.tmp1 <- rep(c(0, 1), rate.cat)
    vec.tmp2 <- rep(1:rate.cat, rep(2, rate.cat)) - 1
    for (i in 1:(k * rate.cat)) {
      mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
      mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
    }
    matFINAL <- mat1 + mat2
    rate.mat.index <- matrix(NA, k * rate.cat, k * rate.cat)
    np <- k + (rate.cat - 1) * 6
    index <- matFINAL == 1
    rate.mat.index[index] <- 1:np
    if (rate.cat == 1) {
      rownames(rate.mat.index) <- c("(0)", "(1)")
      colnames(rate.mat.index) <- c("(0)", "(1)")
    }
    if (rate.cat == 2) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)")
    }
    if (rate.cat == 3) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
    }
    if (rate.cat == 4) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)")
    }
    if (rate.cat == 5) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)")
    }
    if (rate.cat == 6) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
    }
    if (rate.cat == 7) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)",
                                    "(0,R7)", "(1,R7)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)",
                                    "(0,R7)", "(1,R7)")
    }
  }
  if (hrm == FALSE) {
    k = ntraits
    nl = 2
    if (ntraits == 1) {
      k <- 1
      nl <- nstates
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl, nl)
        tmp2 <- cbind(1:(nl^k), 1:(nl^k))
        index <- matrix(TRUE, nl^k, nl^k)
        diag(index) <- FALSE
        if (model == "ER") {
          np <- 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- nl * (nl - 1)/2
          sel <- col(rate.mat.index) < row(rate.mat.index)
          rate.mat.index <- t(rate.mat.index)
          rate.mat.index[sel] <- 1:np
          rate.mat.index[upper.tri(rate.mat.index)] = t(rate.mat.index)[upper.tri(rate.mat.index)]
        }
        if (model == "ARD") {
          np <- nl * (nl - 1)
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 2) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 0, 1, 1)
      vec.tmp2 <- c(0, 1, 0, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
      }
      matFINAL <- mat1 + mat2
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 4
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 4, 6)] <- rate.mat.index[index][c(3,
                                                                          5, 7, 8)] <- 1:np
        }
        if (model == "ARD") {
          np <- 8
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 3) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      mat3 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 1, 0, 0, 1, 1, 0, 1)
      vec.tmp2 <- c(0, 0, 1, 0, 1, 0, 1, 1)
      vec.tmp3 <- c(0, 0, 0, 1, 0, 1, 1, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
        mat3[i, ] <- abs(vec.tmp3 - vec.tmp3[i])
      }
      matFINAL <- mat1 + mat2 + mat3
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 12
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 3, 5, 6, 8, 9,
                                  11, 12, 15, 18, 21)] <- rate.mat.index[index][c(4,
                                                                                  7, 10, 13, 16, 14, 19, 17, 20, 22, 23, 24)] <- 1:np
        }
        if (model == "ARD") {
          np <- 24
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 1) {
      rownames(rate.mat.index) <- as.character(1:nl)
      colnames(rate.mat.index) <- as.character(1:nl)
    }
    if (ntraits == 2) {
      rownames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)",
                                    "(1,1)")
      colnames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)",
                                    "(1,1)")
    }
    if (ntraits == 3) {
      rownames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)",
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)",
                                    "(1,1,1)")
      colnames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)",
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)",
                                    "(1,1,1)")
    }
  }
  return(rate.mat.index)
}

rate.par.drop <- function(rate.mat.index=NULL,drop.par=NULL){
  if(is.null(rate.mat.index)){
    stop("Rate matrix needed.  See mat.maker to create one.\n")
  }
  if(is.null(drop.par)){
    cat("No parameters indicated to drop.  Original matrix returned.\n")
    return(rate.mat.index)
  }
  if(max(rate.mat.index,na.rm=TRUE) < max(drop.par,na.rm=TRUE)){
    cat("Some parameters selected for dropping were not in the original matrix.\n")
  }
  drop.par <- unique(drop.par) # in case parameters listed more than once in drop vector
  drop.par <- drop.par[order(drop.par)]
  max <- max(rate.mat.index,na.rm=TRUE)
  for(drop.which in 1:length(drop.par)){
    drop.locs <- which(rate.mat.index == drop.par[drop.which],arr.ind=TRUE)
    rate.mat.index[drop.locs] <- NA
  }
  max <- max - length(drop.par)
  exclude <- which(is.na(rate.mat.index))
  rate.mat.index[-exclude] <- 1:max
  
  return(rate.mat.index)
}

#root function
root.obs <- function(prop, bins, cutoff=3.1){
  qq <- qnorm(prop, 0, 1)
  x <- pnorm(seq(-1*cutoff, cutoff, length.out=bins-1), qq,1)
  x <- c(x[2:(length(x)-1)],1) - c(0,x[2:(length(x)-1)])
  P <- c(0, x,0)
  P <- P/sum(P) #Make sure it sums to 1
  return(P)
}

#proportion function
prop <- function(data) {
  x = table(unlist(data[,2]))
  prop_1 <- x[1] / (x[1]+x[2])
  return(prop_1)
}

#BIC function
bic <- function(loglik, k, n) {
  ##where k is number of states 
  ##and n is number of taxa (up for debate lol)
  stopifnot(n > 0, k >= 0)
  bic_v <- -2*loglik + k*log(n)
  return(bic_v) 
} 

# making all matrices

## lists to contain the 200 sims for each Markov model and the threshold model
matrix_listARD <- replicate(n=8, expr=list())
matrix_listER <- replicate(n=8, expr=list())
matrix_listSYM <- replicate(n=8, expr=list())
matrix_listThreshER <- replicate(n=8, expr=list())
matrix_listThreshARD <- replicate(n=8, expr=list())

## number of traits (1-7, might need to expand this later but I will need to fix how I store the matrices)
for (i in  1:8){
  #ARD
  m1 <- rate.mat.maker(ntraits = 1, nstates = i, model = "ARD", hrm = F)
  matrix_listARD[[i]]<- m1
  
  #ER
  m2 <- rate.mat.maker(ntraits = 1, nstates=i, model="ER", hrm = F)
  matrix_listER[[i]] <- m2
  
  #SYM
  m3 <- rate.mat.maker(ntraits = 1, nstates=i, model ="SYM", hrm=F)
  matrix_listSYM[[i]] <- m3
  
  #need to modify this for a threshold model
  m4 <- rate.mat.maker(ntraits = 1, nstates=i, model = "SYM", hrm=F)
  matrix_listThreshER[[i]] <- m4
}


## two doesnt work for a thresh model (where would the absorbing states be :^) ) 

two <- matrix_listThreshER[[2]]

## three

three <- matrix( c(NA,1,NA,1,NA,1,NA,1,NA), nrow = 3, ncol= 3)
rownames(three) <- c("1","2","3")
colnames(three) <- c("1","2","3")
matrix_listThreshER[[3]] <- three
## four

four <- matrix(c(NA,1,NA,NA,NA,NA,1,NA,NA,1,NA,NA,NA,NA,1,NA), nrow=4, ncol=4)
rownames(four) <- c("1","2","3","4")
colnames(four) <- c("1","2","3","4")
matrix_listThreshER[[4]] <- four

## five -- don't really need this maybe 
five <- matrix(c(NA,1,NA,NA,NA,1,NA,1,NA,NA,NA,1,NA,1,NA,NA,NA,1,NA,1,NA,NA,NA,1,NA), nrow =5, ncol =5)
rownames(five) <- c("1","2","3","4","5")
colnames(five) <- c("1","2","3","4","5")
matrix_listThreshER[[5]] <- five

## six

six <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,NA,NA,NA,NA,NA,1,NA),
              nrow=6, ncol=6)
rownames(six) <- c("1","2","3","4","5","6")
colnames(six) <- c("1","2","3","4","5","6")
matrix_listThreshER[[6]] <- six

## seven
seven <- matrix(c(NA,1,NA,NA,NA,NA,NA,
                  1,NA,1,NA,NA,NA,NA
                  ,NA,1,NA,1,NA,NA,NA
                  ,NA,NA,1,NA,1,NA,NA,
                  NA,NA,NA,1,NA,1,NA,
                  NA,NA,NA,NA,1,NA,1,
                  NA,NA,NA,NA,NA,1,NA),
                ncol=7, nrow=7)
rownames(seven) <- c("1","2","3","4","5","6","7")
colnames(seven) <- c("1","2","3","4","5", "6","7")
matrix_listThreshER[[7]] <- seven


## eight
eight <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,
                  NA,NA,1,NA,NA,NA,NA, NA
                  ,NA,1,NA,1,NA,NA,NA,NA
                  ,NA,NA,1,NA,1,NA,NA,NA,
                  NA,NA,NA,1,NA,1,NA,NA,
                  NA,NA,NA,NA,1,NA,1,NA,
                  NA,NA,NA,NA,NA,1,NA,NA
                  ,NA,NA,NA,NA,NA,NA,1,NA
),
ncol=8, nrow=8)
rownames(eight) <- c("1","2","3","4","5","6","7","8")
colnames(eight) <- c("1","2","3","4","5", "6","7","8")
matrix_listThreshER[[8]] <- eight

#fresh ARD thresh

## two -- doesnt work for a thresh model

two <- matrix_listThreshARD[[2]]

## three

three_temp <- matrix( c(NA,1,NA,2,NA,3,NA,4,NA), nrow = 3, ncol= 3)
rownames(three_temp) <- c("1","2","3")
colnames(three_temp) <- c("1","2","3")
matrix_listThreshARD[[3]] <- three_temp
## four

four_temp <- matrix(c(NA,1,NA,NA,NA,NA,3,NA,NA,4,NA,NA,NA,NA,6,NA), nrow=4, ncol=4)
rownames(four_temp) <- c("1","2","3","4")
colnames(four_temp) <- c("1","2","3","4")
matrix_listThreshARD[[4]] <- four_temp

## five 
five_temp <- matrix(c(NA,1,NA,NA,NA,2,NA,3,NA,NA,NA,4,NA,5,NA,NA,NA,6,NA,7,NA,NA,NA,8,NA), nrow =5, ncol =5)
rownames(five_temp) <- c("1","2","3","4","5")
colnames(five_temp) <- c("1","2","3","4","5")
matrix_listThreshARD[[5]] <- five_temp

## six

six_temp <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,3,NA,NA,NA,NA,4,NA,5,NA,NA,NA,NA,6,NA,7,NA,NA,NA,NA,8,NA,NA,NA,NA,NA,NA,10,NA),
                   nrow=6, ncol=6)
rownames(six_temp) <- c("1","2","3","4","5","6")
colnames(six_temp) <- c("1","2","3","4","5","6")
matrix_listThreshARD[[6]] <- six_temp

## seven
seven_temp <- matrix(c(NA,1,NA,NA,NA,NA,NA,
                       2,NA,3,NA,NA,NA,NA
                       ,NA,4,NA,5,NA,NA,NA
                       ,NA,NA,6,NA,7,NA,NA,
                       NA,NA,NA,8,NA,9,NA,
                       NA,NA,NA,NA,10,NA,11,
                       NA,NA,NA,NA,NA,12,NA),
                     ncol=7, nrow=7)
rownames(seven_temp) <- c("1","2","3","4","5","6","7")
colnames(seven_temp) <- c("1","2","3","4","5", "6","7")
matrix_listThreshARD[[7]] <- seven_temp

## eight
eight_temp <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,
                       NA,NA,3,NA,NA,NA,NA, NA
                       ,NA,4,NA,5,NA,NA,NA,NA
                       ,NA,NA,6,NA,7,NA,NA,NA,
                       NA,NA,NA,8,NA,9,NA,NA,
                       NA,NA,NA,NA,10,NA,11,NA,
                       NA,NA,NA,NA,NA,12,NA,NA,
                       NA,NA,NA,NA,NA,NA,14,NA
),
ncol=8, nrow=8)
rownames(eight_temp) <- c("1","2","3","4","5","6","7","8")
colnames(eight_temp) <- c("1","2","3","4","5", "6","7","8")
matrix_listThreshARD[[8]] <- eight_temp

#hidden thresh models

thresh_four <- matrix_listThreshER[[4]]
thresh_six <- matrix_listThreshER[[6]]
thresh_eight <- matrix_listThreshER[[8]]
# with absorbing state
#thresh_ten <-  matrix(ncol = 10, nrow = 10,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                                        10,NA,11,NA,NA,NA,NA,NA,NA,NA,
#                                                        NA,20,NA,21,NA,NA,NA,NA,NA,NA,
#                                                        NA,NA,30,NA,31,NA,NA,NA,NA,NA,
#                                                        NA,NA,NA,40,NA,42,NA,NA,NA,NA,
#                                                        NA,NA,NA,NA,51,NA,52,NA,NA,NA,
#                                                        NA,NA,NA,NA,NA,61,NA,62,NA,NA,
#                                                        NA,NA,NA,NA,NA,NA,71,NA,72,NA,
#                                                        NA,NA,NA,NA,NA,NA,NA,81,NA,82,
#                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
#ER
thresh_ten <-  matrix(ncol = 10, nrow = 10,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

thresh_12 <- matrix(ncol = 12, nrow = 12,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

thresh_14 <- matrix(ncol = 14, nrow = 14,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))



thresh_16 <-  matrix(ncol = 16, nrow = 16,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))


mk_four <- matrix_listER[[4]]
mk_six <- matrix_listER[[6]]
mk_eight <- matrix_listER[[8]]

#mk_ten <- matrix(ncol = 10, nrow = 10,byrow = T, c(NA,1,2,3,4,5,6,7,8,9,
#                                                   10,NA,11,12,13,14,15,16,17,18,
#                                                  19,20,NA,21,22,23,24,25,26,27,
#                                                   28,29,30,NA,31,32,33,34,35,36,
#                                                   37,38,39,40,NA,42,43,44,45,46,
#                                                   47,48,49,50,51,NA,52,53,54,55,
#                                                   56,57,58,59,60,61,NA,62,63,64,
#                                                   65,66,67,68,69,70,71,NA,72,73,
#                                                   74,75,76,77,78,79,80,81,NA,82,
#                                                   83,84,85,86,87,88,89,90,91,NA))

mk_ten <- rate.mat.maker(ntraits =1, nstates =10, model = "ER", hrm=F)
mk_twel <- rate.mat.maker(ntraits = 1, nstates = 12, model = "ER", hrm = F)
mk_14 <- rate.mat.maker(ntraits = 1, nstates = 14, model = "ER", hrm = F)
mk_16 <- rate.mat.maker(ntraits=1, nstates = 16, model = "ER", hrm = F)
mk_20 <- rate.mat.maker(ntraits=1, nstates = 20, model = "ER", hrm = F)


sim_aic <- list()
sim_bic <- list()

for (i in 1:nsim) {
  count <- as.character(i)
  #trees
  tree50 <- sim.bdtree(b=1,d = 0,n = 50)
  tree100 <- sim.bdtree(b=1, d=0, n=100)
  tree200 <- sim.bdtree(b=1,d=0, n=200)
  #data
  a0 <- runif(1,-1,1)
  thresh50 <- fastBM(tree50, a=a0, sig2 = 1)
  th50_df <- data.frame(thresh50)
  th50_four <- case_when(
    th50_df[,1] > 0 ~ "1&2",
    th50_df[,1] < 0 ~ "3&4")
  
  while(length(unique(th50_four)) == 1) {
    thresh50 <- fastBM(phy = tree50, a = a0, sig2 =1)
    th50_df <- data.frame(thresh50)
    th50_four <- case_when(
      th50_df[,1] > 0 ~ "1&2",
      th50_df[,1] < 0 ~ "3&4")
  }
  
  thresh100 <- fastBM(tree100, a=a0, sig2 = 1)
  thresh200 <- fastBM(tree200, a=a0, sig2 = 1)
  
  #splitting
  #50
  th50_df <- data.frame(thresh50)
  th50_four <- case_when(
    th50_df[,1] > 0 ~ "1&2",
    th50_df[,1] < 0 ~ "3&4")
  th50_six <- case_when(
    th50_df[,1] > 0 ~ "1&2&3",
    th50_df[,1] < 0 ~ "4&5&6"
  )
  th50_eight <- case_when(
    th50_df[,1] > 0 ~ "1&2&3&4",
    th50_df[,1] < 0 ~ "5&6&7&8"
  )
  th50_ten <- case_when(
    th50_df[,1] > 0 ~ "1&2&3&4&5",
    th50_df[,1] < 0 ~ "6&7&8&9&10"
  )
  th50_12 <- case_when(
    th50_df[,1] > 0 ~ "1&2&3&4&5&6",
    th50_df[,1] < 0 ~ "7&8&9&10&11&12"
  )
  th50_14 <- case_when(
    th50_df[,1] > 0 ~ "1&2&3&4&5&6&7",
    th50_df[,1] < 0 ~ "8&9&10&11&12&13&14"
  )
  th50_16 <- case_when(
    th50_df[,1] >0 ~ "1&2&3&4&5&6&7&8",
    th50_df[,1] <0 ~ "9&10&11&12&13&14&15&16"
  )
  th50_20 <- case_when(
    th50_df[,1] >0 ~ "1&2&3&4&5&6&7&8&9&10",
    th50_df[,1] <0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  
  
  th50_four <- data.frame(tree50$tip.label, th50_four)
  th50_six <- data.frame(tree50$tip.label, th50_six)
  th50_eight <- data.frame(tree50$tip.label, th50_eight)
  th50_ten <- data.frame(tree50$tip.label, th50_ten)
  th50_twel <- data.frame(tree50$tip.label, th50_12)
  th50_14 <- data.frame(tree50$tip.label, th50_14)
  th50_16 <- data.frame(tree50$tip.label, th50_16)
  th50_20 <- data.frame(tree50$tip.label, th50_20)
  
  #100 species
  th100_df <- data.frame(thresh100)
  th100_four <- case_when(
    th100_df[,1] > 0 ~ "1&2",
    th100_df[,1] < 0 ~ "3&4")
  th100_six <- case_when(
    th100_df[,1] > 0 ~ "1&2&3",
    th100_df[,1] < 0 ~ "4&5&6"
  )
  th100_eight <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4",
    th100_df[,1] < 0 ~ "5&6&7&8"
  )
  th100_ten <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5",
    th100_df[,1] < 0 ~ "6&7&8&9&10"
  )
  th100_12 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6",
    th100_df[,1] < 0 ~ "7&8&9&10&11&12"
  )
  th100_14 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6&7",
    th100_df[,1] < 0 ~ "8&9&10&11&12&13&14"
  )
  th100_16 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6&7&8",
    th100_df[,1] < 0 ~ "9&10&11&12&13&14&15&16"
  )
  th100_12 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6",
    th100_df[,1] < 0 ~ "7&8&9&10&11&12"
  )
  th100_14 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6&7",
    th100_df[,1] < 0 ~ "8&9&10&11&12&13&14"
  )
  th100_20 <- case_when(
    th100_df[,1] > 0 ~ "1&2&3&4&5&6&7&8&9&10",
    th100_df[,1] < 0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  th100_four <- data.frame(tree100$tip.label, th100_four)
  th100_six <- data.frame(tree100$tip.label, th100_six)
  th100_eight <- data.frame(tree100$tip.label, th100_eight)
  th100_ten <- data.frame(tree100$tip.label, th100_ten)
  th100_twel <- data.frame(tree100$tip.label, th100_12)
  th100_14 <- data.frame(tree100$tip.label, th100_14)
  th100_16 <- data.frame(tree100$tip.label, th100_16)
  th100_20 <- data.frame(tree100$tip.label, th100_20)
  
  #200 species
  th200_df <- data.frame(thresh200)
  th200_four <- case_when(
    th200_df[,1] > 0 ~ "1&2",
    th200_df[,1] < 0 ~ "3&4")
  th200_six <- case_when(
    th200_df[,1] > 0 ~ "1&2&3",
    th200_df[,1] < 0 ~ "4&5&6"
  )
  th200_eight <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4",
    th200_df[,1] < 0 ~ "5&6&7&8"
  )
  th200_ten <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4&5",
    th200_df[,1] < 0 ~ "6&7&8&9&10"
  )
  th200_12 <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4&5&6",
    th200_df[,1] < 0 ~ "7&8&9&10&11&12"
  )
  th200_14 <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4&5&6&7",
    th200_df[,1] < 0 ~ "8&9&10&11&12&13&14"
  )
  th200_16 <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4&5&6&7&8",
    th200_df[,1] < 0 ~ "9&10&11&12&13&14&15&16"
  )
  
  th200_20 <- case_when(
    th200_df[,1] > 0 ~ "1&2&3&4&5&6&7&8&9&10",
    th200_df[,1] < 0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  th200_four <- data.frame(tree200$tip.label, th200_four)
  th200_six <- data.frame(tree200$tip.label, th200_six)
  th200_eight <- data.frame(tree200$tip.label, th200_eight)
  th200_ten <- data.frame(tree200$tip.label, th200_ten)
  th200_twel <- data.frame(tree200$tip.label, th200_12)
  th200_14 <- data.frame(tree200$tip.label, th200_14)
  th200_16 <- data.frame(tree200$tip.label, th200_16)
  th200_20 <- data.frame(tree200$tip.label, th200_20)
  
 
# MODEL RUNS 
  
  # 4 thresh/Mk
  
  ## 50
  four50_Mk <- rayDISC(phy = tree50, data = th50_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th50_four), 4), model = "ER", node.states = "none")
  
  four50_thresh <- rayDISC(phy= tree50, data = th50_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th50_four), 4), model = "ER", node.states = "none")
  ## 100
  four100_Mk <- rayDISC(phy= tree100, data = th100_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th100_four), 4), model = "ER", node.states = "none")
  
  four100_thresh <- rayDISC(phy= tree100, data = th100_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th100_four), 4), model = "ER", node.states = "none")
  
  ## 200
  four200_Mk <- rayDISC(phy = tree200, data = th200_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th200_four), 4), model = "ER", node.states = "none")
  
  four200_thresh <- rayDISC(phy = tree200, data = th200_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th200_four), 4), model = "ER", node.states = "none")
  
  # 6 thresh/Mk
  
  ## 50
  six50_Mk <- rayDISC(phy = tree50, data = th50_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th50_six), 6), node.states = "none", model = "ER")
  six50_thresh <- rayDISC(phy= tree50, data = th50_six, rate.mat = thresh_six, root.p =  root.obs(prop = prop(th50_six), 6), node.states = "none", model = "ER")
  
  ## 100
  six100_Mk <- rayDISC(phy= tree100, data = th100_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th100_six), 6), node.states = "none", model = "ER")
  six100_thresh <- rayDISC(phy= tree100, data = th100_six, rate.mat = thresh_six, root.p = root.obs(prop = prop(th100_six), 6), node.states = "none", model = "ER")
  
  ## 200
  six200_Mk <- rayDISC(phy = tree200, data = th200_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th200_six), 6), node.states = "none", model = "ER")
  six200_thresh <- rayDISC( phy = tree200, data = th200_six, rate.mat = thresh_six, root.p = root.obs(prop = prop(th200_six), 6), model = "ER")
  
  params <- sapply(na.omit(c(six200_thresh$index.mat)), 
                   function(x) na.omit(c(six200_thresh$solution))[x])
  ancRECON(rate.cat = 1, phy = tree200, data = th200_six, p = params, root.p = root.obs(prop = prop(th200_six), 6), model = "ER", method = "joint")
  # 8 thresh/mk

  ## 50
  eight50_Mk <- rayDISC(phy = tree50, data = th50_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th50_eight), 8), model = "ER", node.states = "none")
  eight50_thresh <- rayDISC(phy= tree50, data = th50_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th50_eight), 8), model = "ER", node.states = "none")
  
  ## 100
  eight100_Mk <- rayDISC(phy= tree100, data = th100_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th100_eight), 8), model = "ER", node.states = "none")
  eight100_thresh <- rayDISC(phy= tree100, data = th100_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th100_eight), 8), model = "ER", node.states = "none")
  
  ## 200
  eight200_Mk <- rayDISC(phy = tree200, data = th200_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th200_eight), 8), model = "ER", node.states = "none")
  eight200_thresh <- rayDISC(phy = tree200, data = th200_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th200_eight), 8), model = "ER", node.states = "none")
  
  
  ## 50
  ten50_Mk <- rayDISC(phy = tree50, data = th50_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th50_ten), 10), model = "ER", node.states = "none")
  ten50_thresh <- rayDISC(phy= tree50, data = th50_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th50_ten), 10), model = "ER", node.states = "none")
  
  ## 100
  ten100_Mk <- rayDISC(phy= tree100, data = th100_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th100_ten), 10), model = "ER", node.states = "none")
  ten100_thresh <- rayDISC(phy= tree100, data = th100_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th100_ten), 10), model = "ER", node.states = "none")
  
  ## 200
  ten200_Mk <- rayDISC(phy = tree200, data = th200_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th200_ten), 10), model = "ER", node.states =  "none")
  ten200_thresh <- rayDISC(phy = tree200, data = th200_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th200_ten), 10), model = "ER", node.states =  "none")
  
  ## 50
  twel50_Mk <- rayDISC(phy = tree50, data = th50_twel, rate.mat = mk_twel, root.p = root.obs(prop = prop(th50_twel), 12), model = "ER", node.states = "none")
  twel50_thresh <- rayDISC(phy= tree50, data = th50_twel, rate.mat = thresh_12, root.p = root.obs(prop = prop(th50_twel), 12), model = "ER", node.states = "none")
  
  ## 100
  twel100_Mk <- rayDISC(phy= tree100, data = th100_twel, rate.mat = mk_twel, root.p = root.obs(prop = prop(th100_twel), 12), model = "ER", node.states = "none")
  twel100_thresh <- rayDISC(phy= tree100, data = th100_twel, rate.mat = thresh_12, root.p = root.obs(prop = prop(th100_twel), 12), model = "ER", node.states = "none")
  
  ## 200
  twel200_Mk <- rayDISC(phy = tree200, data = th200_twel, rate.mat = mk_twel, root.p = root.obs(prop = prop(th200_twel), 12), model = "ER", node.states =  "none")
  twel200_thresh <- rayDISC(phy = tree200, data = th200_twel, rate.mat = thresh_12, root.p = root.obs(prop = prop(th200_twel), 12), model = "ER", node.states =  "none")
  
  ## 50
  fourteen50_Mk <- rayDISC(phy = tree50, data = th50_14, rate.mat = mk_14, root.p = root.obs(prop = prop(th50_14), 14), model = "ER", node.states = "none")
  fourteen50_thresh <- rayDISC(phy= tree50, data = th50_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th50_14), 14), model = "ER", node.states = "none")
  
  ## 100
  fourteen100_Mk <- rayDISC(phy= tree100, data = th100_14, rate.mat = mk_14, root.p = root.obs(prop = prop(th100_14), 14), model = "ER", node.states = "none")
  fourteen100_thresh <- rayDISC(phy= tree100, data = th100_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th100_14), 14), model = "ER", node.states = "none")
  
  ## 200
  fourteen200_Mk <- rayDISC(phy = tree200, data = th200_14, rate.mat = mk_14, root.p = root.obs(prop = prop(th200_14), 14), model = "ER", node.states =  "none")
  fourteen200_thresh <- rayDISC(phy = tree200, data = th200_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th200_14), 14), model = "ER", node.states =  "none")
  
  ## 50
  sixteen50_Mk <- rayDISC(phy = tree50, data = th50_16, rate.mat = mk_16, root.p = root.obs(prop = prop(th50_16), 16), model = "ER", node.states =  "none")
  sixteen50_thresh <- rayDISC(phy= tree50, data = th50_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th50_16), 16), model = "ER", node.states =  "none")
  
  ## 100
  sixteen100_Mk <- rayDISC(phy= tree100, data = th100_16, rate.mat = mk_16, root.p = root.obs(prop = prop(th100_16), 16), model = "ER", node.states =  "none")
  sixteen100_thresh <- rayDISC(phy= tree100, data = th100_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th100_16), 16), model = "ER", node.states =  "none")
  
  ## 200
  sixteen200_Mk <- rayDISC(phy = tree200, data = th200_16, rate.mat = mk_16, root.p=root.obs(prop = prop(th200_16), 16), model = "ER", node.states =  "none")
  sixteen200_thresh <- rayDISC(phy = tree200, data = th200_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th200_16),16), model = "ER", node.states =  "none")
  
  #delta AIC and BIC
  all_aic <- matrix(nrow = 3, ncol = 7)
  
  #four AIC
  all_aic[1,1] <- delta50_four <- four50_thresh$AIC - four50_Mk$AIC
  all_aic[2,1] <- delta100_four <- four100_thresh$AIC - four100_Mk$AIC
  all_aic[3,1] <- delta200_four <- four200_thresh$AIC - four200_Mk$AIC
  
  #six AIC
  all_aic[1,2] <- delta50_six <- six50_thresh$AIC - six50_Mk$AIC
  all_aic[2,2] <- delta100_six <- six100_thresh$AIC - six100_Mk$AIC
  all_aic[3,2] <- delta200_six <- six200_thresh$AIC - six200_Mk$AIC
  
  #eight AIC
  all_aic[1,3] <- delta50_eight <- eight50_thresh$AIC - eight50_Mk$AIC
  all_aic[2,3] <- delta100_eight <- eight100_thresh$AIC - eight100_Mk$AIC
  all_aic[3,3] <- delta200_eight <- eight200_thresh$AIC - eight200_Mk$AIC

  #ten AIC
  all_aic[1,4] <- delta50_ten <- ten50_thresh$AIC - ten50_Mk$AIC
  all_aic[2,4] <- delta100_ten <- ten100_thresh$AIC - ten100_Mk$AIC
  all_aic[3,4] <- delta200_ten <- ten200_thresh$AIC - ten200_Mk$AIC
  
  #12 AIC
  all_aic[1,5] <- delta50_12 <- twel50_thresh$AIC - twel50_Mk$AIC
  all_aic[2,5] <- delta100_12 <- twel100_thresh$AIC - twel100_Mk$AIC
  all_aic[3,5] <- delta200_12 <- twel200_thresh$AIC - twel200_Mk$AIC
  
  #14 AIC
  all_aic[1,6] <- delta50_14 <- fourteen50_thresh$AIC - fourteen50_Mk$AIC
  all_aic[2,6] <- delta100_14 <- fourteen100_thresh$AIC - fourteen100_Mk$AIC
  all_aic[3,6] <- delta200_14 <- fourteen200_thresh$AIC - fourteen200_Mk$AIC
  
  #sixteen AIC
  all_aic[1,7] <- delta50_sixteen <- sixteen50_thresh$AIC - sixteen50_Mk$AIC
  all_aic[2,7] <- delta100_sixteen <- sixteen100_thresh$AIC - sixteen100_Mk$AIC
  all_aic[3,7] <- delta200_sixteen <- sixteen200_thresh$AIC - sixteen200_Mk$AIC
  
  allaic_df <- data.frame(all_aic)

  #heatmap.2(x=all_bic)
  print(i)
  sim_aic[[i]] <- all_aic
}

## AIC 
all_sims <- sim_aic
saveRDS(all_sims, file = "./output/heatmap1_aic_results.rds")



avg_sims <- do.call(cbind, all_sims)
avg_sims <- array(avg_sims, dim=c(dim(all_sims[[1]]), length(all_sims)))
avg_sims <- apply(avg_sims, c(1, 2), mean, na.rm = TRUE)
median_sims <- apply(avg_sims, c(1,2), median, na.rm=  TRUE)

h1 <- heatmap.2(x = avg_sims, labRow = c(50,100,500), labCol = c(4,6,8,10,12,14,16), xlab = "Hidden states", main = "Heatmap 1: Threshold data, deltaAIC = threshAIC- MkAIC",
                   ylab = "sample size", dendrogram = 'none',Rowv=FALSE, Colv=FALSE,trace='none')

saveRDS(h1, file = "heatmap1")
